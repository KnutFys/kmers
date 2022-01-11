# Standard library imports
import os
import gzip
import concurrent.futures
import random
import sys
from time import time
# Local imports
import readcounter
import kmer_DB as DB
import settings

def file_name_resolver(filename, suffix=None, exclude=[]):
	def name_acceptable(name):
		if name in exclude:
			return False
		if os.path.isfile(os.path.join(sample_path, name)):
			return False
		con = DB.connect()
		c = con.cursor()
		command = 'SELECT path FROM sample_files WHERE path = ?'
		c.execute(command, (name,))
		row = c.fetchall()
		con.close()
		if len(row) > 0:
			return False
		return True
	sample_path = settings.get_setting('sample_path')
	base_name = filename.split('.')[0]
	if suffix is None:
		suffix = filename.split('.')[1]
	new_name = base_name + f'.{suffix}'
	i = 1
	acceptable = name_acceptable(new_name)
	while not acceptable and i < 1000:
		new_name = base_name + f'_{i}' + f'.{suffix}'
		acceptable = name_acceptable(new_name)
		i += 1
	return new_name

def extract_lots(files, sample_sizes):
	'''Extracts one sample for each sample size from each file

	Arguments:
	files			-- (lst) List of files to extract from
	sample_sizes	-- (lst) List of sample sizes to draw
	'''

	all_database_entries = []
	sample_names = []
	# Get read counts
	con = DB.connect()
	c = con.cursor()
	command = 'SELECT path, read_count FROM raw_files WHERE path IN ('
	for f in files:
		command += '?,'
	command = command[:-1] + ')'
	values = tuple([os.path.split(f)[1] for f in files])
	c.execute(command, values)
	rows = c.fetchall()
	read_counts = {r['path']:r['read_count'] for r in rows}
	con.close()
	# For each file create a dictionary for each
	# samplesize that can be inserted into the database
	for f in files:
		database_entries = []
		filename = os.path.split(f)[1]
		count = read_counts[filename]
		sample = filename.split('-')[-1].split('_')
		sample_id = sample[0]
		strand = sample[2]
		sample = f'{sample_id}_{strand}'
		if 'trimmed' in filename:
			trimmed = True
			file_path = os.path.join(settings.get_setting('trimmed_data'),f)
		else:
			trimmed = False
			sample += '_raw'
			file_path = os.path.join(settings.get_setting('raw_data'),f)		
		for size in sample_sizes:		
			out_file = sample + f'_{size}'
			out_file = file_name_resolver(
				out_file, exclude=sample_names, suffix='sample')
			entry = {
				'table':'sample_files',
				'path':out_file,
				'sample_id':sample_id,
				'raw_file':filename,
				'sample_size':size,
				'strand':strand
				}
			database_entries.append(entry)
			sample_names.append(out_file)
		all_database_entries.append((database_entries, file_path, count))

	with concurrent.futures.ProcessPoolExecutor() as executor:
		for r in executor.map(extract_multiple_samples, all_database_entries):
			pass
	
	print('Extraction completed successfully')

def extract_multiple_samples(samples):
	'''Extract multiple random samples from a zipped fastq
	file. Adds entries to database.

	Args:
	samples		-- (list)  Structured as follows
				[list of database entries, file to extract from, read count]
	'''
	sample_sizes = []
	out_files = []
	raw_files = []
	for s in samples[0]:
		raw_files.append(s['raw_file'])
		out_files.append(s['path'])
		sample_sizes.append(s['sample_size'])
	in_file = samples[1]
	reads_in_file = samples[2]

	output_directory = settings.get_setting('sample_path')
	if not os.path.isdir(output_directory):
		os.mkdir(output_directory)
	# Create random index dictionary
	read_lists = {}
	# Create dictionary to store extracted reads
	reads = {}
	# No need to read past the last read we are looking for
	last_sample = 0
	# Draw random index set for each sample size
	for i, sample_size in enumerate(sample_sizes):
		read_list = random.sample(range(reads_in_file), sample_size)
		read_list.sort()
		if read_list[-1] > last_sample:
			last_sample = read_list[-1]
		
		read_lists[out_files[i]] = read_list
		reads[out_files[i]] = []

	###################################
	# Keeps track of the what read is being processed. Each fourth
	# line in a fastq file is a new read
	read_number = 0
	sample_name = os.path.split(in_file)[1]
	with gzip.open(in_file, "rt") as in_handle:
		# Loop through file 4 lines at a time and check if read is to be
		# extracted. Selected reads exported to fasta format  with
		# read decription as >new read number - old readnumber
		for i in range(last_sample):
			line1 = next(in_handle)
			line2 = next(in_handle)
			line3 = next(in_handle)
			line4 = next(in_handle)
			for key in read_lists:
				if read_lists[key]:
					if read_number == read_lists[key][0]:
						r = (f'>{sample_name}_{read_number}', f'{line2.strip()}')
						read_lists[key].pop(0)
						reads[key].append(r)
			read_number += 1
	for key in reads:
		out_file = os.path.join(output_directory, key)
		with open(out_file, 'w') as fh:
			for read in reads[key]:
				fh.write(f'{read[0]}\n')
				fh.write(f'{read[1]}\n')
	con = DB.connect()
	print('Accessing database.')
	with con:
		for s in samples[0]:
			DB.insert_dictionary(s, con)
	print('Closing database')
	con.close()
	
def main():
	settings.set_defaults()
	test()

def test():
	'''
	

	'''
	t0 = time()
	f1 = '1-SV001_S68_R1_001_trimmed.fq.gz'
	f2 = '1-SV001_S68_R2_001_trimmed.fq.gz'
	two_files = [f1, f2]
	sample_sizes = [
		1600000, 1600000, 1600000, 1600000, 1600000,
		1600000, 1600000, 1600000, 1600000, 1600000,
		1600000
	]
	extract_lots(two_files,sample_sizes)
	print('Running time: ', time()-t0)

if __name__ == '__main__':
	main()