##########################
#       IMPORTS
##########################

# Standard library
import os
import sys
import json
import concurrent.futures
import gzip
# Local imports
import zipsampler
import readcounter
import calculations
import settings
import kmer_DB as DB
from time import time
#######################
#  K-mer manipulation
#######################
def add_nuc(kmers,alfabet):
	''' Takes a list of kmers and returns a list of all k+1-mers extended
	 by alphabet

	Args:
	kmers - (lst) List of kmers
	alfabet - List of letters/nucleotides
	'''
	new_kmers = []
	for s in kmers:
		for n in alfabet:
			new_kmers.append(s+n)
	return new_kmers

def canonical(kmer):
	''' Returns canonic nucleotide sequence, that is the first of
	the lexicographically ordered pair of the sequence and its reverse
	compliment

	Args:
	kmer - (str) A string using only letters AGCT
 	'''
	translation = 'string'.maketrans('ACGT','TGCA')
	reverse = kmer[::-1]
	canonic = sorted([kmer,reverse.translate(translation)])[0]
	return canonic

def count_kmers(seq, offset, k=None, **kwargs):
	'''Counts K-mers in a single sequence and returns a dict {kmer:count}
	for each offset.
	
	Arguments:
	offset			-- (lst) List of offsets to use
	k				-- (int) Kmer length

	Kwargs:	
	k_dictionaries	-- (lst) List of dictionaries containing kmer counts
					for each offset. If supplied the function will fill these
					dictionaries rather than start from emty ones.
	'''

	if k is None:
		k = settings.get_setting('k')
	if 'k_dictionaries' in kwargs:
		k_dictionaries = kwargs['k_dictionaries']
	else:
		k_dictionaries = [{} for o in offset]
	for ofs, k_dictionary in zip(offset, k_dictionaries):
		for i in range(len(seq)-k-ofs):
			kmer = seq[i+ofs:i+k+ofs]
			kmer = canonical(kmer)
			if kmer in k_dictionary:
				k_dictionary[kmer] += 1
			else:
				k_dictionary[kmer] = 1
	return k_dictionaries

def convert_jf_file(filename):
	''' Creates a count file from a jellyfish dump file and store it
	in the database.

	Arguments:
	filename		-- (str) filename of jellyfish dump file
	'''
	if not filename.endswith('.dump'):
		print('Invalid filetype')
		return
	count_path = settings.get_setting('count_path')
	out_file = os.path.split(filename)[1].replace('.dump', '_ofs_0.count')
	sample_id = out_file.split('_')[0]
	k = int(out_file.split('_k')[1].split('_ofs')[0])
	sample_file = f'full_{out_file.split("_")[1]}'
	database_entry = {
		'table':'count_files',
		'path':out_file,			
		'sample_id':sample_id,
		'sample_file':sample_file,
		'offset':0,
		'kmer':k,
		'strand':'both'
	}
	out_file = os.path.join(count_path, out_file)
	path = database_entry["path"]
	condition = f'path="{path}"'
	rows = DB.search(table='count_files', condition=condition)

	if len(rows) < 1:
		print('Not in database')
		in_database = False
	else:
		print('In database')
		in_database = True
	if os.path.isfile(out_file):
		print('Alleady converted')
	else:
		print('Converting...')
		k_dictionary = {}
		with open(filename, 'r') as fh:
			for line in fh:
				count = int(line.strip().replace('>', ''))
				kmer = next(fh).strip()
				k_dictionary[kmer] = count
		with open(out_file, 'w') as fh:
			json.dump(k_dictionary, fh)
		if not in_database:
			con = DB.connect()
			with con:
				DB.insert_dictionary(database_entry, con)
			con.close()
			print('Database entry added')

def count_kmers_in_file(filename, offset=None, k=None):
	'''Counts kmers in fasta files. Used to count kmers in subsamples.

	Args:
	filename - (str)  Path to fasta file 
	'''
	
	if k is None:
		k = settings.get_setting('k')
	if offset is None:
		offset = settings.get_setting('offset')
	if not isinstance(offset, list):
		offset = [offset]
	# Standardize output file name
	name_base = os.path.split(filename)[1]
	filename_out = os.path.join(settings.get_setting('count_path'), name_base)
	files_out = [
		f'{filename_out.split(".")[0]}_k{k}_ofs_{ofs}.count' for ofs in offset
	]
	# Check if file has been counted for each offset
	not_counted = [] 
	for f, ofs in zip(files_out, offset):
		if os.path.isfile(f):
			print(f'File {filename} allready counted with offset {ofs}')
		else:
			not_counted.append((f, ofs))
	# Create database entries
	database_entries = []
	sample_id = name_base.split('_')[0]
	strand = name_base.split('_')[1]
	for f_o in not_counted:
		entry = {
			'table':'count_files',
			'path':os.path.split(f_o[0])[1],
			'sample_id':sample_id,
			'sample_file':name_base,
			'offset':f_o[1],
			'kmer':k,
			'strand':strand
			}
		database_entries.append(entry)

	print(f'Starting counting of {k}-mers in files:')
	for e in database_entries:
		print(f'{e["path"]}')
	k_dictionaries = [{} for _ in not_counted]
	offset = [f_o[1] for f_o in not_counted]
	filename = os.path.join(
		settings.get_setting('sample_path'), filename)
	with open(filename,'rt') as file_handle:
		for line in file_handle:
			# Skip fasta headers
			if line.startswith('>'):
				pass
			else:
				seq = line.strip()
				count_kmers(
					seq, offset=offset, k=k, k_dictionaries=k_dictionaries)
	for f_o, d in zip(not_counted, k_dictionaries):
		with open(f_o[0], 'w') as file_handle:
			json.dump(d, file_handle)

	return database_entries
	
def multi_count_kmers(file_list):
	''' Multiprocessing wrapper for function count_kmers_file(filename)

	Args:
	file_list - (lst) List of paths to fasta files
	'''
	database_entries = []
	with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
		for r in executor.map(count_kmers_in_file, file_list):
			for entry in r:
				database_entries.append(entry)
		
	con = DB.connect()
	with con:
		for e in database_entries:
			DB.insert_dictionary(e, con)
	con.close()
	print('Counting finished')

def load_list(filename):
	try:
		if os.path.isfile:
			with open(filename, 'r') as fh:
				loaded_list = [f.strip() for f in fh]
		else:
			raise Exception(f'{filename} not found.')
	except Exception as e:
		print(e)
		sys.exit(2)
	return loaded_list
###################################################################
# Runtime
###################################################################
def compute_distance_pair(count_file_1, count_file_2, cutoff=None):
	''' Checks database for known distances between files and
	calculates unknown distances before updating database.

	Args:
	count_file_1 - (str) Path to kmer count file
	count_file_2 - (str) Path to kmer count file
	cutoff - (int, None) Ignore kmers with counts < cutoff
	'''

	count_path = settings.get_setting('count_path')
	count_file_1 = os.path.split(count_file_1)[1]
	count_file_2 = os.path.split(count_file_2)[1]
	distances = DB.check_distance(count_file_1, count_file_2)
	if distances is None:
		print('Not in database')
		in_database = False
		distances = {
		'cosine':None,
		'jensenshannon':None,
		'braycurtis':None
		}
		all_known = False
	else:
		print('In database')
		in_database = True
		all_known = True
		for d in distances:
			if distances[d] is None:
				all_known = False
	if not all_known:
		con = DB.connect()
		f1 = os.path.join(count_path, count_file_1)
		f2 = os.path.join(count_path, count_file_2)
		distance_to_update = calculations.calculate_distance(
			f1, f2, cutoff=cutoff, dist=distances)
		if in_database:
			with con:
				DB.update_distance(f1, f2, con, **distance_to_update)
		else:
			distance_to_insert = DB.create_distance_dic(
				file1=f1, file2=f2, **distance_to_update)
			with con:
				DB.insert_dictionary(distance_to_insert, con)
		con.close()

def compute(**kwargs):
	# Check max number of runs
	max_runs = kwargs.get('max', None)

	samples = kwargs.get('samples', None)
	offset = kwargs.get('offset', 0)
	sample_size = kwargs.get('sample_size', 100000)
	kmer = settings.get_setting('k')
	count_path = settings.get_setting('count_path')

	sample_ids = [r['sample_id'] for r in DB.get_sample_ids()]
	print(sample_ids)
	files = DB.get_count_files(offset=offset, samples=samples, kmer=kmer)
	files = [f for f in files if f'_{sample_size}_' in f]
	#Filter
	filtered = []
	for s in sample_ids:
		filtered += [f for f in files if s in f][2:4]
	files = [os.path.join(count_path, f) for f in filtered]

	file_pairs = []
	cutoff = settings.get_setting('cutoff')
	known = 0
	unknown_filepair = 0
	db_col_names = {
			'cs':'cosine',
			'js':'jensenshannon',
			'bc':'braycurtis'
		}
	for i in range(len(files)-1):
		for j in range(i+1, len(files)):
			#Ignore pairs from same sample		
			if files[i].split('_')[0] != files[j].split('_')[0]:
				f1 = os.path.join(count_path, files[i])
				f2 = os.path.join(count_path, files[j])
				distances = DB.check_distance(f1, f2)
				if distances is None:
					unknown_filepair += 1
					distances = {
						'cosine':None,
						'jensenshannon':None,
						'braycurtis':None
					}
					file_pairs.append(
						(f1, f2, distances))
				else:
					known += 1
	total = unknown_filepair + known
	print('Total distance pairs to compute: ', total)
	print('Distances allready known: ', known)
	print('Distances to calculate: ', unknown_filepair)
	if max_runs is not None:
		print('Max runs set to ', max_runs)
		if max_runs < len(file_pairs):
			total = max_runs
			file_pairs = file_pairs[:max_runs]

	con = DB.connect()
	for i, f in enumerate(file_pairs):
		print(f'Calculating distance {i+1} of {len(file_pairs)}')
		distance_to_update = calculations.calculate_distance(
			f[0], f[1], cutoff=cutoff, dist=f[2])
	
		distance_to_insert = DB.create_distance_dic(
				file1=f[0], file2=f[1], **distance_to_update)
		with con:
			DB.insert_dictionary(distance_to_insert, con)
		print('Database updated')
	
	con.close()
	print('Computation done')

def extract(**kwargs):
	####################################
	# Not implemented
	####################################
	''' Extracts reads from zipped fastq files.

	kwargs:
	sample_size	--	(int) Number of reads to extract
	filename	--	(str) Path to zipped fastq file
	'''
	try:
		pass
	except Exception as e:
		print('Error resolving arguments')
		print(e)
		sys.exit(2)	
	zipsampler.extract_lots(files)

def count_all_uncounted(k=11, sample=None):
	''' Counts kmers for all uncounted fasta files in database. For counting
	kmers in subsamples.

	Args:
	k - (int) Kmer size to use
	sample - (str) Only count sub samples starting with specified string
	'''
	unc = DB.uncounted(k=k)	
	settings.set_setting('k', k)
	if sample is not None:
		unc = [u for u in unc if u.startswith(sample)]
	multi_count_kmers(unc)

def combine_counts(file1, file2, outfile):
	'''Combines counts from two samples'''

	with open(file1, 'r') as fh:
		counts1 = json.load(fh)
	with open(file2, 'r') as fh:
		counts2 = json.load(fh)
	kmers1 = set(counts1.keys())
	kmers2 = set(counts2.keys())
	kmers = kmers1.union(kmers2)
	total_counts = {}
	for k in kmers:
		total_counts[k] = counts1.get(k, 0) + counts2.get(k, 0)
	with open(outfile, 'w') as fh:
		json.dump(total_counts, fh)
	
def test(**kwargs):
	'''
	Convergence of alpha diversity?
	'''

	k=11
	count_path = settings.get_setting('count_path')
	files = os.listdir(count_path)
	pole = [os.path.join(
		count_path,f) for f in files if f'k{k}' in f and 'pole' in f][0]
	outfiles = [os.path.join(
		count_path, f'Combined{i+1}Samples.count') for i in range(25)]

def argument_handler():
	command = sys.argv.pop(1)
	arguments = {}

	def single(arg):
		if arg in sys.argv:
			pos = sys.argv.index(arg)
			value = sys.argv.pop(pos+1)
			_ = sys.argv.pop(pos)
		return value
	def multi(arg):
		value = []
		pos = sys.argv.index(arg)
		end_pos = len(sys.argv) - 1
		for word in sys.argv[pos+1:]:
			if word.startswith('-'):
				end_pos = sys.argv.index(word)-1
				break
		for i in range(end_pos-pos):
			value.append(sys.argv.pop(pos+1))
		_ = sys.argv.pop(pos)
		return value

	try:
		if '-k' in sys.argv:
			settings.set_setting('k', int(single('-k')))
		if '-offset' in sys.argv:
			offset = int(single('-offset'))
			arguments['offset'] = offset
			settings.set_setting('offset', offset)
		if '-cutoff' in sys.argv:
			settings.set_setting('cutoff', int(single('-cutoff')))
		if '-p' in sys.argv:
			arguments['path'] = single('-p')
		if '-f' in sys.argv:
			arguments['file'] = single('-f')
		if '-samples' in sys.argv:
			arguments['samples'] = multi('-samples')
		if '-ft' in sys.argv:
			arguments['filetype'] = single('-ft')
		if '-dm' in sys.argv:
			settings.set_setting('distance_measure', single('-dm'))
		if '-size' in sys.argv:
			arguments['sample_size'] = int(single('-size'))
		if '-sizes' in sys.argv:
			arguments['sample_sizes'] = multi('-sizes')
		if '-c' in sys.argv:
			arguments['contains'] = single('-c')
		if '-e' in sys.argv:
			arguments['end'] = single('-e')
		if '-s' in sys.argv:
			arguments['start'] = single('-s')
		if '-metric' in sys.argv:
			arguments['metric'] = single('-metric')
		if '-max' in sys.argv:
			arguments['max'] = int(single('-max'))

	except Exception as e:
		print(e)
		print('Failed to resolve arguments')
		sys.exit(2)

	commands = {
		'count':count,
		'compute':compute,
		'extract':extract,
		'test':test,
	}
	try:
		commands[command](**arguments)
	except Exception as e:
		print(f'Failed to run command: {command}')
		print(e)
		sys.exit(2)

def main():
	# Set defaults	
	settings.set_defaults()
	# Commands available from command line
	argument_handler()	

if __name__ == "__main__":
	main()

