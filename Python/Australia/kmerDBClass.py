'''
Keeps track of countfiles, distances and their
'''

import os
import sqlite3
import distances

class DistanceDataBase:

	@staticmethod
	def canonical_key(word1, word2):
		'''
		Creates a unique, well defined identifier based on two unique identifiers

		Args:
		word1 - Unique keyword (str)
		word2 - Unique keyword (str)
		'''
		w1 = os.path.split(word1)[1]
		w2 = os.path.split(word2)[1]
		return '$'.join(sorted([w1, w2]))

	# DB/Table structure
	count_file_table = [
		'count_file TEXT PRIMARY KEY',			
		'package_id TEXT NOT NULL'
	]

	distance_table = [
		'file_pair TEXT PRIMARY KEY',
		'count_file1 TEXT NOT NULL',
		'count_file2 TEXT NOT NULL',
		'cosine REAL',
		'jensenshannon REAL',
		'braycurtis REAL',
		'angular REAL'
	]

	# Not in use
	resource_table = [
		'resource_id TEXT PRIMARY KEY',
		'package_id TEXT NOT NULL',
		'downloaded INTEGER NOT NULL',
		'trimmed INTEGER NOT NULL',
		'counted INTEGER NOT NULL'
	]

	tables = {
		'count_files':count_file_table,
		'distance_table':distance_table,
		'resource_table':resource_table
	}

	def __init__(self, name, test=False):
		self.name = name
		self.test = test
		if not os.path.isfile(self.name):
			print('No database found. Initiating new database...')
			self.initiate_database()

	def connect(self):
		''' Connects to database and returns connection.

		Keyword arguments:
		test -- If true, opens connection in memory (default False)
		'''
		if self.test:
			con = sqlite3.connect(':memory:')
		else:
			con = sqlite3.connect(self.name)
		con.row_factory = sqlite3.Row
		return con

	def initiate_database(self):
		'''
		Sets up a database for keeping track of kmer counts and distances between
		samples
		'''	
		connection = self.connect()
		c = connection.cursor()
		for table in self.tables:
			command = f'CREATE TABLE {table} ('
			for column in self.tables[table]:
				command += f'{column}, '
			command = command.strip()
			command = command[:-1] + ')'
			c.execute(command)
		connection.close()
		print('Database initiated.')

	def add_file(self, count_file, package_id):
		'''
		Adds a row to the database table "count_table"

		Args:
		count_file - (str) Unique identifier and name of count file
		package_id - (str) id of resource package the count file was created from
		'''
		connection = self.connect()
		command = '''
			INSERT INTO count_files (count_file, package_id) VALUES (?, ?)
		'''
		values = (count_file, package_id)

		with connection:
			c = connection.cursor()
			c.execute(command, values)
		connection.close()

	def get_file(self, count_file=None, package_id=None):
		'''
		Returns a dictionary containing kmer count files with a given package id
		and or the package id of a given kmer count file

		Args:
		count_file - (str, default=None) Name of kmer count file
		package_id - (str, default=None) Search database for kmer count files
			with given package id
		'''
		connection = self.connect()
		def execute(connection, command, values):
			with connection:
				cursor = connection.cursor()
				cursor.execute(command, values)
				rows = cursor.fetchall()
				return rows

		result = {'count_files':[], 'package_id':package_id}
		if count_file is not None:
			command = '''
				SELECT
					package_id
				FROM
					count_files
				WHERE
					count_file = ?
			'''
			values = (count_file,)
			rows = execute(connection, command, values)
			result['package_id'] = rows[0][0]
		if package_id is not None:
			command = '''
				SELECT
					count_file
				FROM
					count_files
				WHERE
					package_id = ?
			'''
			values = (package_id,)
			rows = execute(connection, command, values)
			for row in rows:
				for field in row:
					result['count_files'].append(field)
		return result 

	def get_file_pairs(self):
		'''Returns a list of all pairs of files entred in the database'''

		connection = self.connect()
		command = '''
			SELECT
				count_file
			FROM
				count_files
		'''
		with connection:
			cursor = connection.cursor()
			cursor.execute(command)
			rows = cursor.fetchall()
		connection.close()
		file_entries = [row[0] for row in rows]
		file_pairs = []
		for i in range(len(file_entries) - 1):
			for j in range(i + 1 , len(file_entries)):
				file_pairs.append(
					self.canonical_key(file_entries[i], file_entries[j]))
		return file_pairs

	def add_distances(self):
		'''
		Adds a rows to the database table "distance_table" for each pair of
		files in the "count_files" table
		'''
		# Get all file pairs in DB
		connection = self.connect()
		with connection:
			cursor = connection.cursor()
			cursor.execute('''SELECT file_pair FROM distance_table''')
			db_pairs = [row[0] for row in cursor.fetchall()]
		# Get all possible pairs of files in DB
		total_pairs = self.get_file_pairs()
		print('Total file pairs: ', len(total_pairs))
		print('Pairs in data base: ', len(db_pairs))
		to_update = [pair for pair in total_pairs if pair not in db_pairs]
		print('Pairs to add: ', len(to_update))
		if len(to_update) > 0:
			command = '''
			INSERT INTO
				distance_table(
					file_pair,
					count_file1,
					count_file2
					)
			VALUES(?,?,?)
			'''
			values = []
			for pair in to_update:
				# Split canonical key into file names
				files = pair.split('$')
				values.append((pair, files[0], files[1]))

			with connection:
				c = connection.cursor()
				c.executemany(command, values)

	def get_unknown_distances(self, invert=False):
		''' Finds all distance entries in database that have not yet
		been calculated

		Args:
		invert - (bol, default False) If true, returns all known distances
		'''
		switch = ''
		if invert:
			switch = 'NOT'
		connection = self.connect()
		command = f'''
			SELECT 
				file_pair, count_file1, count_file2
			FROM
				distance_table
			WHERE
				cosine IS {switch} NULL
		'''
		with connection:
			c = connection.cursor()
			c.execute(command)
			rows = c.fetchall()

		connection.close()

		result = []
		for row in rows:
			result.append(
				{column[0]:row[i] for i, column in enumerate(c.description)})

		return result

	def get_count_files(self):
		''' Returns a list of dictionaries with one
		entry for each kmer count file entered in database
		'''
	
		connection = self.connect()
		command = f'''
			SELECT 
				count_file
			FROM
				count_files
		'''
		with connection:
			c = connection.cursor()
			c.execute(command)
			rows = c.fetchall()

		connection.close()

		result = []
		for row in rows:
			result.append(
				{column[0]:row[i] for i, column in enumerate(c.description)})

		return result		

	def set_distances(self, distance_dict, multi=False):
		''' Updates entry in distance table of the database

		Args:
		distance_dict - (dic, lst) dictionary of distances or list of such
		multi - (bol, default=False) - If true, function will assume a list of dicts
			as the first argument
		'''
		command = '''
			UPDATE 
				distance_table
			SET
				cosine = :cosine,
				jensenshannon = :jensenshannon,
				braycurtis = :braycurtis,
				angular = :angular
			WHERE
				file_pair = :file_pair
		'''
		connection = self.connect()
		with connection:
			cursor = connection.cursor()
			if multi:
				cursor.executemany(command, distance_dict)
			else:
				cursor.execute(command, distance_dict)

	def get_distance(self, file1, file2):
		'''
		Search database for calculated distances between specified files.

		Args:
		file1 - Database identifier for file (str)
		file1 - Database identifier for file (str)
		'''

		file_pair = self.canonical_key(file1, file2)

		command = '''
		SELECT
			cosine, jensenshannon, braycurtis, angular
		FROM
			distance_table
		WHERE
			file_pair = ?
			'''
		connection = self.connect()
		c = connection.cursor()
		c.execute(command, (file_pair,))
		row = c.fetchone()
		connection.close()
		if row is None:
			return None
		else:
			result = {column[0]:row[i] for i, column in enumerate(c.description)}
			return result

	def testtest(self):
		print(self.name)

def main():
	db = DistanceDataBase('kmerdatabase.db', test=False)
	db.testtest()

if __name__ == '__main__':
	main()
	