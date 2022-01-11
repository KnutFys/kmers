# Standard library imports
import sys
import os
import json
from time import time
import subprocess
import unittest
# Imports
import pandas as pd # Not needed?
import numpy as np
# Local imports
import kmerDBClass # # Stores inter sample distances
import dbSearch # Connects to AMB database
import resultViewer # Examines AMB search results
import processgenome # Handles genome processing via trimgalore and jellyfish
import distances	# Calculates distances

class MetaDataFrame:
	'''
	
	'''
	field_size = 20
	cutoff = field_size - 5

	def __init__(self, samples=[], fields=[]):
		self.fields = fields
		self.samples = samples
		self.only_valid = True

	def iterate_samples(self, fields=None):
		''' Iterates over samples and yields samples with valid entries for
		specified fields
		'''
		if fields is None:
			fields = self.fields

		for sample in self.samples:
			try:
				for field in fields:
					value = sample.fields.get(field, None)
					if value is None or value == '':
						raise ValueError('Invalid field entry')
				yield sample
			except ValueError as e:
				continue


	def to_csv(self, filename, sep=',', fields=None, pole_file=None):
		''' Saves data in csv format 
		
		Args:
		filename - (str) Basename for ou put files
		sep - (str, default ",") Separator to use when writing file
		fields - (lst, None) List of fields to include, If None uses all fields
		add_pole - (str, None) If specified, attempts to load pole distances and 
					add these to the meta data.
		'''

		filename = filename + '_meta.csv'
		if fields is None:
			fields = self.fields

		if pole_file is not None:
			with open(pole_file, 'r') as fh:
				pole_distances = json.load(fh)
			measures = ['angular', 'jensenshannon', 'braycurtis', 'kmer_count']

		def line_iterator(samples):
			for sample in samples:
				try:
					words = [sample.name]
					for field in fields:
						value = sample.fields.get(field, None)
						if value is None or value == '':
							raise ValueError('Invalid field entry')
						words.append(value)
					if pole_file is not None:
						for measure in measures:
							value = pole_distances[sample.name][measure]
							words.append(value)
					yield words
				except ValueError as e:
					continue

		if not os.path.isfile(filename):
			header = ['id'] + self.fields
			if pole_file is not None:
				header = header + measures
			header = sep.join(header) + '\n'
			with open(filename, 'w') as fh:
				fh.write(header)
				for words in line_iterator(self.iterate_samples(fields)):
					words = [str(word) for word in words]
					line = sep.join(words) + '\n'
					fh.write(line)
		else:
			print('Filename in use.')

	def __str__(self):
		''' Prints out meta data in table form'''

		def style_formatter(to_display):
			''' Formats cell structure '''

			if not isinstance(to_display, str):
				to_display = f'{type(to_display)}'
			if to_display == '':
				to_display = 'No entry'
			return f'{to_display:{self.field_size}}|'

		separator = '_'*((len(self.fields)+1)*(self.field_size + 1))

		header = style_formatter('Id')
		for f in self.fields:
			header += style_formatter(f)
		print('')
		print(separator)
		print(header)
		print(separator)
		
		samples = self.samples
		if self.only_valid:
			samples = self.iterate_samples()
		for s in samples:
			if len(s.name) > self.cutoff:
				line = style_formatter(s.name[:self.cutoff]+'...')
			else:
				line = style_formatter(s.name)
			for f in self.fields:
				line += style_formatter(s.fields.get(f, ' --- '))
			print(line)
		print(separator)
		print(header)
		return separator

class DistanceMatrix:
	def __init__(self, samples, database, remove_missing=True):
		self.samples = samples
		self.database = database
		self.create_arrays()

	def create_arrays(self):
		self.distances = {
			'angular':np.zeros((len(self.samples), len(self.samples))),
			'cosine':np.zeros((len(self.samples), len(self.samples))),
			'braycurtis':np.zeros((len(self.samples), len(self.samples))),
			'jensenshannon':np.zeros((len(self.samples), len(self.samples)))
		}
		# Collect distances from database
		for i in range(len(self.samples)-1):
			for j in range(i+1, len(self.samples)):
				file1 = self.samples[i].name + '.count'
				file2 = self.samples[j].name + '.count'
				distance = self.database.get_distance(file1, file2)
				if distance is None:
					distance = {
						'angular':None,
						'jensenshannon':None,
						'braycurtis':None,
						'cosine':None
						}
				for distance_measure, matrix in self.distances.items():
					matrix[i,j] = distance[distance_measure]
					matrix[j,i] = distance[distance_measure]

	def to_csv(self, filename, sep=','):
		for distance_measure in self.distances:
			out_file = filename + f'_{distance_measure}.csv'
			with open(out_file, 'w') as fh:
				matrix = self.distances[distance_measure]
				number_of_samples = matrix.shape[0]
				header = ['id']
				for sample in self.samples:
					header.append(sample.name)
				header = sep.join(header) + '\n'
				fh.write(header)
				for i in range(number_of_samples):
					line = [self.samples[i].name]
					for j in range(number_of_samples):
						line.append(f'{matrix[i,j]}')
					line = sep.join(line) + '\n'
					fh.write(line)

class Sample:

	def __init__(self, name, package_id, R1=None, R2=None):
		self.name = name
		self.package = package_id
		self.R1 = R1
		self.R2 = R2
		self.fields = {}
		self.resource = None
	
	def add_field(self, field, value):
		self.fields[field] = value

class Project:
	'''
	
	'''

	def __init__(self, name, path=None):
		self.name = name
		if path is None:
			self.path = os.path.join(
				os.getcwd(),'Projects', name, name+'.prj')
		else:
			self.path = path
		# Load existing project or create new if no project exists
		try:
			self.project = self.load_project()
		except FileNotFoundError as e:
			print('Initiating new project')
			self.new_project()

		self.packages = {}
		self.samples = {}
		self.search_results = {}
		self.load_search_results()
		self.fields = set()
		for field in self.project.get('fields', []):
			self.add_field(field)
		self.database_file = os.path.join(
			self.project['directories']['base_dir'], self.project['dbname'])
		self.database = kmerDBClass.DistanceDataBase(self.database_file)

	def set_api(self, api):
		''' Stores API key for accessing AMB database. Needed for downloads and
		searches

		Args:
		api - (str) Api. Contact
		'''	
		self.project['api'] = api
		self.save_project()

	def new_project(self):
		'''
		Create directories structure and store paths

		Args:
		base_dir - (str) Path to parent directory where project directories will
			be created
		'''
		
		base_dir = os.path.join(os.getcwd(), 'Projects', self.name)
		directories = {
			'base_dir':base_dir,
			'count_dir':os.path.join(base_dir, 'Counts'),
			'download_dir':os.path.join(base_dir, 'Downloads'),
			'trimmed_dir':os.path.join(base_dir, 'Trimmed'),
			'search_dir':os.path.join(base_dir, 'Searchresults')
		}
		# Check if project folder exists and exit if it does
		if os.path.isdir(base_dir):
			sys.exit(f'Conflict with existing project at {base_dir}')
		for key in directories:
			if not os.path.isdir(directories[key]):
				os.makedirs(directories[key])

		self.project = {
			'name':self.name,
			'dbname':self.name + '.db',
			'projectfile':self.name + '.prj',
			'api':None,
			'directories':directories,
			'searches':[]
		}
		self.save_project()

	def view_project(self):
		'''Prints info on project settings and files'''

		def purtify(word):
			print('-'*20)
			print(word)
			print('-'*20)

		purtify('Project settings: ')
		for key, item in self.project.items():
			if key not in ['searches', 'directories', 'fields']:
				print(key, ' : ', item)
		purtify('Project Directories: ')
		for key, item in self.project['directories'].items():
			print(key, ' : ', item)
		purtify('Mandatory fields')
		for field in self.fields:
			print(field)
		purtify(f'Searches: {len(self.project["searches"])}')
		for s in self.project['searches']:
			r = resultViewer.SearchResult(s['file'])
			r.view_search()
		purtify(f'Samples: {len(self.samples)}')
		print(f'Packages: {len(self.packages)}')
		if len(self.fields) > 0:
			valid_samples = [s for s in self.iterate_samples()]
			valid_packages = {s.package:None for s in valid_samples}

			print(f'Packages with valid fields: {len(valid_packages)}')
			print(f'Samples with valid fields: {len(valid_samples)}')
		samples_in_db = self.database.get_count_files()
		print('Samples in database: ', len(samples_in_db))
		
		purtify('Count status')
		distances_in_db = self.database.get_unknown_distances(invert=True)		
		sample_pairs = int((
			len(samples_in_db)*len(samples_in_db) - len(samples_in_db))/2)
		print('Possible distances: ', sample_pairs)
		print('Distances in database: ', len(distances_in_db))

		purtify('Download status')
		sample_types = {
			'Pelagic':0, 'Coastal water':0, 'Sediment':0, 'Unknown':0, 'Soil':0}
		unk_sample_types = {
			'Pelagic':0, 'Coastal water':0, 'Sediment':0, 'Unknown':0, 'Soil':0}

		for name, sample in self.samples.items():
			sample_type = sample.fields.get('sample_type', 'Unknown')
			if sample_type is None or sample_type == '':
				sample_type = 'Unknown'
			if self.check_status(sample.name.replace('_RX', '_R1'))['counted']:
				sample_types[sample_type] += 1
			else:
				unk_sample_types[sample_type] += 1
		for key in sample_types:
			print(key, ' Counted:', sample_types[key], 'Not counted:', unk_sample_types[key])


	def view_searches(self):
		searches = [
			resultViewer.SearchResult(
				s['file']) for s in self.project['searches']]
		for s in searches:
			print('-'*25)
			print('Query: ', s.query)
			for f in s.fields:
				s.add_field(f)
		

	def create_sample_metadata_csv(self, fields=None):
		file_name = os.path.join(
			self.project['directories']['base_dir'],'sample_meta_data.csv')

		if fields is None:
			fields = ['geo_loc_name', 'sample_type']
		
		header = ['Sample', 'Package'] + fields + ['Downloaded']
		header = ','.join(header)

		with open(file_name, 'w') as fh:
			fh.write(header + '\n')
			for name, sample in self.samples.items():
				line = [name, sample.package]
				for field in fields:
					value = self.packages[sample.package].get(field, 'Unknown')
					line.append(value)
				
				line.append(str(self.check_status(sample.R1)['counted']))
				line = ','.join(line)
				fh.write(line + '\n')

	def find_complete_samples(self, fields=None):
		''' Adds "complete" attribute to the project pointing to a dictionary
		of samples where all samples have valid field entries and distances
		between all samples are known.

		Args:
		fields - (lst, default None) List of fields to include. 
			If None, use all project fields.
		'''

		if fields is None:
			fields = self.fields

		print('Project samples: ', len(self.samples))
		# Find samples with valid field entries
		t0 = time()
		self.complete = {
			s.name:s for s in self.iterate_samples(fields=fields, only_valid=True)}
		t1 = time()
		print(t1-t0)
		print('Samples with valid fields: ', len(self.complete))
		# Find processed samples
		processed_samples = [
			s['count_file'].split('.')[0] for s in self.database.get_count_files()]
		to_remove = [s for s in self.complete if s not in processed_samples]
		for s in to_remove:
			removed = self.complete.pop(s)
		print('Processed samples: ', len(self.complete))
		# Find samples that lack distances
		unknown_pairs = [u['file_pair'] for u in self.database.get_unknown_distances(
			invert=False)]
		unknown_pairs = set(unknown_pairs)

		# Score by number of missing distances
		def prune_samples(unknown_pairs, samples):
			''' Recursively removes samples with fewest known distances in order
			to get a set of samples with a complete set of inter sample distances

			Args:
			unknown_pairs - (lst) list of canonical file pairs with unknown distance
			samples - (dic) Dictionary of samples to be pruned 
			'''
			relevant_pairs = []
			for sample1 in samples:
				for sample2 in samples:
					if sample1 == sample2:
						continue
					relevant_pairs.append(self.database.canonical_key(
						sample1 + '.count', sample2 + '.count'))
			relevant_pairs = set(relevant_pairs)
			relevant_unknown = unknown_pairs.intersection(relevant_pairs)

			fewest_known = [None, 0]
			for name, s in samples.items():
				score = 0
				for pair in relevant_unknown:
					if name in pair:
						score+=1
				if score > fewest_known[1]:
					fewest_known = [name, score]
			if fewest_known[1] == 0:
				return samples
			else:
				removed = samples.pop(fewest_known[0])
				samples = prune_samples(unknown_pairs, samples)
				return samples

		self.complete = prune_samples(unknown_pairs, self.complete)
		print('Samples with valid fields and complete distance set: ',
			len(self.complete))

	def iterate_samples(self, fields=None, only_valid=True):
		''' Iterates over samples and yields samples with valid entries for
		specified fields
		'''
		if fields is None:
			fields = self.fields

		for name, sample in self.samples.items():
			try:
				for field in fields:
					value = sample.fields.get(field, None)
					if (value is None or value == '') and only_valid:
						raise ValueError('Invalid field entry')
				yield sample
			except ValueError as e:
				continue

	def add_field(self, field):
		''' Adds a metadata field to the project and 
		unpacks relevant metadata for each sample in project.
		'''
		self.fields.add(field)
		for name, smpl in self.samples.items():
			smpl.fields[field] = self.packages[smpl.package].get(field, None)
		fields = list(self.fields)
		if fields is None:
			fields = []
		self.project['fields'] = list(self.fields)
		
	def remove_field(self, field):
		''' Remove mandatory field from project'''
		try:
			self.fields.remove(field)
		except KeyError:
			pass
		try:
			self.project['fields'] = list(self.fields)
		except TypeError:
			# If set is empty, list conversion fails
			self.project['fields'] = []

	def create_data_frame(self, fields=None, only_valid=True):
		''' Create pandas data frame with meta data entries

		Args:
		fields - (lst, default None) List of field names, if None use all fields.
		only_valid - (bol, default True) Only include samples with valid entries for
			all fields.
		'''
		
		if fields is None:
			fields = [field_name for field_name in self.fields]	
		samples = self.iterate_samples(fields, only_valid=only_valid)

		# Create records/dictionaries for selected samples
		data = []
		for sample in samples:		
			record = {'ID':sample.name}
			for field in fields:
				record[field] = sample.fields[field]
			data.append(record)		
		self.data_frame = pd.DataFrame.from_records(data)

	def create_data_table(self, fields=None, only_valid=True, samples=None):
		''' Create a table containing selected meta data fields
		for project samples

		Args:
		fields - (lst) List of meta data fields to include
		only_valid - (bol) If true, only include samples with valid entries
			for all fields.
		'''
		if fields is None:
			fields = [name for name in self.fields]
		if samples is None:
			samples = [sample for sample in self.iterate_samples(
				fields=fields, only_valid=only_valid)]
		self.data_table = MetaDataFrame(samples=samples, fields=fields)

	def create_distance_matrices(self, fields=None, only_valid=True, samples=None):
		if fields is None:
			fields = [name for name in self.fields]
		if samples is None:
			samples = [sample for sample in self.iterate_samples(
				fields=fields, only_valid=only_valid)]
		self.dm = DistanceMatrix(samples, self.database)

	def to_csv(self, filename, fields=None):
		''' Creates meta data and distance distance csv files.

		Args:
		filename - (str) base for output file names
		fields - (lst) list of fields to add to metadata
		'''
		print('Finding complete samples...')
		self.find_complete_samples(fields=fields)
		samples = [s for name, s in self.complete.items()]
		print('Creating data tables...')
		self.create_data_table(samples=samples, fields=fields)
		self.create_distance_matrices(samples=samples, fields=fields)
		out_csv_file = os.path.join(
			self.project['directories']['base_dir'], filename)
		pole_file = os.path.join(
			self.project['directories']['base_dir'], 'pole_dist.json')
		self.dm.to_csv(out_csv_file)
		self.data_table.to_csv(
			out_csv_file, fields=fields, pole_file=pole_file)

	def plot_R(self, filename, col, shape, measure):
		''' Calls an R script that calculates and plots an 
		NMDS representation of exported distance matrix.

		Args:
		filename - (str) Base filename, metadata and distance matrix files
			will be created based on this name
		col - (str) Field to use when coloing samples in plot
		shape - (str) Field to use when determining shape representation for
			samples in plot.
		measure - (str) What distance measure to plot
		'''

		rscript = os.path.join(os.getcwd(), 'rforpoject.R')
		base_dir = self.project['directories']['base_dir']
		R_command = ['--vanilla', rscript, base_dir, filename, measure, col, shape]
		
		try:
			subprocess.run(R_command, check=True , shell=False,
			executable='/usr/bin/Rscript')
		except subprocess.CalledProcessError as e:
			print('Error running R script.')
			print(e)
			return 1
			
	def load_project(self):
		''' Load settings for project 
	
		Args:
		path -- (str )path to project file
		'''

		try:
			with open(self.path, 'r') as fh:
				project = json.load(fh)
			return project
		except FileNotFoundError as e:
			print('Could not find project file')
			raise FileNotFoundError

	def save_project(self):
		''' Save settings for project '''
		try:
			with open(self.path, 'w') as fh:
				json.dump(self.project, fh)
		except Exception as e:
			print('Error saving project settings')
			print(e)

	def load_search_results(self):
		''' Loads all search results stored in project or found in 
		search folder '''

		# Find searches	performed from this project	
		searches = self.project['searches']
		# Look for seaches added to search folder
		search_path = self.project['directories']['search_dir']
		local_files = os.listdir(search_path)
		local_files = [
			os.path.join(search_path,f) for f in local_files]
		project_search_files = [s['file'] for s in searches]
		external_searches = [
			f for f in local_files if f not in project_search_files]
		for s in external_searches:
			external_search = {'query':'Unknown', 'file':s}
			searches.append(external_search)
		
		self.searches = searches
		# Load all searches and unpack
		for search in self.searches:
			with open(search['file'], 'r') as fh:
				search_result = json.load(fh)
			for package in search_result['results']:
				self.packages[package['id']] = package
		self.unpack_search()

	def unpack_search(self):
		''' Unpacks resources and packages and links resources
		belonging to same sample '''

		# Pair resources
		for package_name, package in self.packages.items():

			try:
				# Skip controls
				if 'control' in package_name:
					continue

				paired_resources = {}
				for resource in package['resources']:
					# Skip checksums
					if not resource['name'].endswith('.fastq.gz'):
						continue
					base_name = resource['name'].split('.')[0]
					if '_R1_' in base_name:
						base_name = base_name.replace('_R1_', '_RX_')
						strand = 'R1'
					elif '_R2_' in base_name:
						base_name = base_name.replace('_R2_', '_RX_')
						strand = 'R2'
					else:
						#print(resource['name'])
						raise ValueError('Unknown naming scheme')

					pair = paired_resources.get(base_name, {})
					pair[strand] = resource
					pair['name'] = base_name
					pair['package'] = package_name
					paired_resources[base_name] = pair

				for name, pair in paired_resources.items():
					smpl = Sample(name, pair['package'],
						R1=pair['R1'], R2=pair['R2'] )
					self.samples[smpl.name] = smpl
			except ValueError as e:
				#print(e)
				continue

	def search(self, query, rows=1000):
		''' Searches AMB database for metagenomes relevant
		to query.
		'''
		# Find unused filename
		search_file_number = 1
		free_name_found = False
		while not free_name_found:
			search_file = os.path.join(
				self.project['directories']['search_dir'],
				f'search{search_file_number}.json'	
			)
			if os.path.isfile(search_file):
				search_file_number += 1
				continue
			free_name_found = True

		s = {'query':query, 'file':search_file}
		dbSearch._apikey = self.project['api']
		result = dbSearch.searchDB(query=query, rows=rows)
		with open(search_file, 'w') as fh:
			json.dump(result, fh)

		self.project['searches'].append(s)
		self.save_project()

	def update_db(self):
		'''
		Add count files and empty distance records to database
		'''

		# Find local count files
		count_path = self.project['directories']['count_dir']
		count_files = [f for f in os.listdir(count_path) if f.endswith('.count')]
		# Add count files to database
		total = len(count_files)
		known = 0
		added = 0
		for f in count_files:
			try:
				sample_name = f.split('.count')[0]
				self.database.add_file(f, self.samples[sample_name].package)
				added += 1
			except kmerDBClass.sqlite3.IntegrityError:
				known += 1
		self.database.add_distances()
		print('Total count files: ', total)
		print('Count files in data base: ', known)
		print('Files added to data base: ', added)
		# Add distance pairs to database

	def calculate(self, batchsize=1, workers=1):

		# Get uncalculated distances in database
		unknown = self.database.get_unknown_distances()
		if len(unknown) < 1:
			print('No distances to calculate. Try updating database.')
			return
		# Set number of workers/cpus for multiprocessing
		self.workers = workers
		print('Uncalculated distances: ', len(unknown))
		# Set full paths for count files
		for record in unknown:
			record['f1'] = os.path.join(
				self.project['directories']['count_dir'], record['count_file1'])
			record['f2'] = os.path.join(
				self.project['directories']['count_dir'], record['count_file2'])
		
		batch = [[record['f1'], record['f2'], record] for record in unknown]
		results = distances.multi_batch(
			batch, batchsize=batchsize, max_workers=self.workers)
		distance_dicts = []
		for r in results:
			distance_dict = r[0][2]
			for distance_measure in r[1]:
				distance_dict[distance_measure] = r[1][distance_measure]
			distance_dicts.append(distance_dict)
		self.database.set_distances(distance_dicts, multi=True)

	def count(self, from_process=False):
		''' Counts kmers for trimmed files assosciated with the project. IE
		files in the 'trimmed' folder.
		'''

		count_path = self.project['directories']['count_dir']
		trimmed_path = self.project['directories']['trimmed_dir']
		trimmed_files = [f for f in os.listdir(trimmed_path) if '_trimmed' in f]
		counted_pairs = 0
		for trimmed_pair in self.pairs(trimmed_files):
			# Skip counted files
			if self.check_status(trimmed_pair[0])['counted']:
				print(f'{trimmed_pair[0]}\nand\n{trimmed_pair[1]}')
				print('allready counted')
				continue
			print('Counting kmers for file pair:')
			print(f'{trimmed_pair[0]}\n{trimmed_pair[1]}')

			input_files = [os.path.join(trimmed_path, f) for f in trimmed_pair]
			# Format final output file name 
			count_file = self.get_file_names(trimmed_pair[0])['count']
			output_file = os.path.join(count_path, count_file) 

			processgenome.jellyfish(input_files, output_file)
			counted_pairs += 1
			if from_process:
				return counted_pairs
			return True
		return False

	def clean_up_downloads(self):
		''' Remove processed files '''
	
		# Remove trimmed and/or  counted files from
		# download folder
		path = self.project['directories']['download_dir']
		files = [f for f in os.listdir(path) if f.endswith('q.gz')]
		for f in files:
			if self.check_status(f)['trimmed']:
				print(f'{f} has been trimmed')
				file_path = os.path.join(path, f)
				os.remove(file_path)
				print(f'{file_path} deleted...')
			else:
				print(f'{f} not yet trimmed')

	def clean_up_trimmed(self):
		''' Remove files from trimmed folder if kmers have been counted'''
	
		path = self.project['directories']['trimmed_dir']
		files = [f for f in os.listdir(path) if '_trimmed.fq.gz' in f]
		for f in files:
			if self.check_status(f)['counted']:
				print(f'Kmers for {f} have been counted')
				file_path = os.path.join(path, f)
				os.remove(file_path)
				print(f'{file_path} deleted...')
			else:
				print(f'Kmers for {f} not yet counted')

	def check_status(self, fileorresource):
		''' Checks if metagenome has been downloaded and or counted. Assumes
		that if a file has been trimmed, it has been downloaded, if it has been
		counted it has been trimmed.

		Args:
		fileorresource - (dic OR str) AMB resource dictionary or filename
		'''

		checks = {
			'downloaded':True,
			'counted':True,
			'trimmed':True
		}
		# Get standardized filenames for count files and trimmed files
		file_names = self.get_file_names(fileorresource)

		if os.path.isfile(os.path.join(
			self.project['directories']['count_dir'], file_names['count'])):
			return checks

		checks['counted'] = False
		if os.path.isfile(os.path.join(
			self.project['directories']['trimmed_dir'], file_names['trimmed'])):
			return checks

		checks['trimmed'] = False
		if os.path.isfile(os.path.join(
			self.project['directories']['download_dir'], file_names['download'])):
			return checks
		checks['downloaded'] = False
		return checks

	def check_file_availability(self, missing=False, available=False):
		'''	Checks whether or not kmer count files entered in database
		can be found in the project count folder and prints result.

		Args:
		missing - (bol, False) If True, prints all missing files
		available - (bol, False) If True, prints all available files
		'''
		database_files = self.database.get_count_files()
		local_file_available = []
		local_file_unavailable = []
		for f in database_files:
			path = os.path.join(
				self.project['directories']['count_dir'], f['count_file'])
			if os.path.isfile(path):
				local_file_available.append(f['count_file'])
			else:
				local_file_unavailable.append(f['count_file'])

		if missing:
			print('Missing files')
			print('-'*20)
			for f in local_file_unavailable:
				print(f)
			print('-'*20)
		if available:
			print('Missing files')
			print('-'*20)
			for f in local_file_unavailable:
				print(f)
			print('-'*20)
		print('Total kmer count files in database: ', len(database_files))
		print('Local kmer count files: ', len(local_file_available))
		print('Missing kmer count files: ', len(local_file_unavailable))

	def download(self, batchsize=1, from_process=False, forced_fields=None):
		'''
		Downloads metagenomes in search results

		Args:
		args - (namespace) Parsed command line options
		from_process - (bol) Switch determining if the function is being called
			from the 'process' function
		forcefield - (dict) Dictionary of fields with mandatory values 
			IE {field:mandatory_value}
		'''

		# Download up to 'args.batchsize' files
		downloaded = 0
		for package_id, package in self.packages.items():
			# Ignore packages lacking specified field entries
			try:
				for field in self.project['fields']:
					field_value = package.get(field, None)
					if field_value is None or field_value == '':
						raise ValueError(
							f'{package_id} lacking entry for field {field}')
				if forced_fields is not None:
					for forced_field, value in forced_fields.items():
						field_value = package.get(forced_field, None)
						if field_value != value:
							raise ValueError(
								f'Wrong entry for forced field {forced_field}')

			except ValueError as e:
				print(e)
				continue

			for resource in package['resources']:
				# Ignore checksum files
				if not resource['name'].endswith('.fastq.gz'):
					#print('Skipping checksum file...')
					continue
				# Ignore counted, trimmed or downloaded files
				if not self.check_status(resource)['downloaded']:
					# Download file
					download_path = self.project['directories']['download_dir']
					download_file = resource['name']
					dbSearch.downloader(
						resource, os.path.join(download_path, download_file))
					downloaded += 1
				# Check if anything has been downloaded
				if downloaded == 0:
					continue
				# Check if called from process function
				if from_process:
					return downloaded
				if downloaded >= batchsize:
					break
			if downloaded >= batchsize:
				break
		return 0

	def trim(self, batchsize=1, from_process=False):
		'''
		Trim downloaded metagenome

		'''
		# Find existing files
		counted_files = os.listdir(self.project['directories']['count_dir'])
		trimmed_files = os.listdir(self.project['directories']['trimmed_dir'])
		downloaded_files = os.listdir(self.project['directories']['download_dir'])
		downloaded_files = [f for f in downloaded_files if f.endswith('fastq.gz')]

		trimmed = 0
		for file_pair in self.pairs(downloaded_files):
			queue = []
			# Add untrimmed files to queue
			for file_name in file_pair:
				if not self.check_status(file_name)['trimmed']:
					queue.append(file_name)
			if len(queue) < 1:
				continue
			# Set full path
			queue = [os.path.join(
				self.project['directories']['download_dir'], f) for f in queue]

			# Call trim_galore
			for file_path in queue:
				processgenome.trim_galore(file_path,
					trimmed_path=self.project['directories']['trimmed_dir'])
				trimmed += 1

			# If called from process function, return
			if from_process:
				return trimmed

			if trimmed >= batchsize:
				break
		if trimmed < batchsize:
			print(f'{trimmed} files were available for trimming')
			return trimmed

	def pole_distance(self, batchsize=2, export=True, add_fields=None):
		''' Calculate distances to pole for each sample 
		
		Args:
		batchsize - (int) Number of calculations to do
		export - (bol) If true skips calculations and creates a csv file of
			allready calculated distances
		add_fields - (lst) If supplied, adds specified fields to the csv file
		'''
		add_fields = ['longitude', 'latitude']
		pole_path = os.path.join(
			self.project['directories']['base_dir'], 'pole_dist.json')
		try:
			with open(pole_path, 'r') as fh:
				pole_distances = json.load(fh)
		except FileNotFoundError as e:
			print('No file found...')
			print('Creating pole distance file...')
			pole_distances = {}
		
		if export:
			#Exports csv file
			sample_names = [name for name in pole_distances]
			fields = [field_name for field_name in pole_distances[sample_names[0]]]
			header = ['Sample']
			header += fields
			if add_fields is not None:
				header += add_fields
			with open(pole_path.replace('.json', '.csv'), 'w') as fh:
				fh.write(f'{",".join(header)}\n')
				for sample, record in pole_distances.items():
					line = [sample]
					for field in fields:
						line.append(str(record.get(field, 'Unknown')))
					if add_fields is not None:
						for field in add_fields:
							package_id = self.samples[sample].package
							word = self.packages[package_id].get(field, 'Unknown')
							line.append(word)
					fh.write(f'{",".join(line)}\n')
			# Exit if export was set to true
			return

		sample_names = [key for key in self.samples]
		count_files = [f['count_file'] for f in self.database.get_count_files()]
		count_files = [f for f in count_files if f.split('.')[0] not in pole_distances]
		if batchsize is None:
			batchsize = len(count_files)
		batchsize = min([len(count_files), batchsize])
		for i, f in enumerate(count_files[:batchsize]):
			print(f'Calculating distance {i+1}/{batchsize}', end='\r')

			sample_name = f.split('.')[0]
			path = os.path.join(
				self.project['directories']['count_dir'], f)
			kmer_counts = distances.load_jellyfish_counts(path)
			# Calculate pole distances
			dist = distances.calculate_distance(kmer_counts, 'pole11')

			sample = self.samples[sample_name]
			record = {measure:value for measure, value in dist.items()}
			record['package'] = sample.package

			# Count kmers
			kmer_count = 0
			for key, value in kmer_counts.items():
				kmer_count += value
			record['kmer_count'] = kmer_count

			# Add fields and 
			for field, value in sample.fields.items():
				record[field] = value
			pole_distances[sample_name] = record

		with open(pole_path, 'w') as fh:
			json.dump(pole_distances, fh)
		print(f'Total pole distances in file: {len(pole_distances)}')

	@staticmethod
	def get_file_names(fileorresource):
		''' Returns standardized file names for
		kmer count files, downloaded fastq files and 
		trimmed fastq files.

		Args:
		fileorresource - (dic OR str) AMB resource dictionary
			or Name of count file, trimmed fasta file or
			fastq file with standardized name.
		'''

		if isinstance(fileorresource, dict):
			file_name = fileorresource['name']
		else:
			file_name = fileorresource
		# Find common root name by removing trimmed tag and replacing
		# strand tag by RX
		# Check if assumed naming scheme holds
		if file_name.count('_R1') + file_name.count('_R2') != 1:
			raise Exception(f'{file_name} has an unresolvable naming scheme')
		base = file_name.replace('_trimmed', '').split('.')[0]
		base = base.replace('_R1', '_RX').replace('_R2', '_RX')
		# Find paired file name
		if '_R1' in file_name:
			paired_file = file_name.replace('_R1', '_R2')
		if '_R2' in file_name:
			paired_file = file_name.replace('_R2', '_R1')
		count = base + '.count'
		download = file_name
		trimmed = file_name.split('.')[0] + '_trimmed.fq.gz'
	
		return {'count':count, 'download':download,
			'trimmed':trimmed, 'base':base, 'paired':paired_file}

	@staticmethod
	def pairs(file_list):
		'''
		Iterator that pairs paired end read files. AMB naming
		scheme is assumed. IE paired files have identical names 
		except for strand indication _R1_ or _R2_. Alternative
		naming scheme uses _R1 and _R2

		Args:
		file_list - (lst) List of file names
		'''
		for f in file_list:
			# Check for first naming scheme
			if '_R1_' in f:
				paired_file = f.replace('_R1_', '_R2_')
			# Check for second naming scheme
			elif '_R1.' in f:
				paired_file = f.replace('_R1.', '_R2.')
			# Skip non matching
			else:
				continue

			yield [f, paired_file]

def main():
	pass

if __name__ == '__main__':
	main()