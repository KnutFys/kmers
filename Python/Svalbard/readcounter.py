# Standard imports
import os
import gzip
import concurrent.futures

# Path to file containing filenames and their read counts
_count_file = 'readcounts.kmer'
# Linux Path to drive where samples are stored 
_drive = '/media/knut/Storage'

def read_counter(filename):
	# Counts reads in zippe fasta file
	print(f'Opening file {filename}')
	counter = 0
	with gzip.open(filename,"rt") as infile:
		for line in infile:
			counter += 1
	return int(counter/4)

def multi_read(file_list):
	# Multiprocessing for function read_counter(filename)
	results = []
	with concurrent.futures.ProcessPoolExecutor() as executor:
		for r in executor.map(read_counter, file_list):
			results.append(r)
	print('Counting finished')
	results = list(zip(file_list,results))
	return results

def count(file_list, read_count_file=_count_file):
	known_counts = load_read_counts(read_count_file)
	known_files = set([c[0] for c in known_counts])
	new_files = set(file_list)
	for f in known_files.intersection(new_files):
		print(f'{f} allready counted. Skipped')
	new_files = list(new_files.difference(known_files))
	for f in new_files:
		print(f'{f} will be counted')
	results = multi_read(new_files)
	if len(results) > 0:
		known_counts += results
	save_read_counts(known_counts, read_count_file)
	return known_counts

def load_read_counts(filename=_count_file):
	read_counts = []
	with open(filename, 'r') as fh:
		for line in fh:
			read_counts.append(tuple(line.strip().split()))
	read_counts = [(convert_path(x[0]), x[1]) for x in read_counts]
	return read_counts

def save_read_counts(counts, filename):
	with open(filename, 'w') as fh:
		for c in counts:
			fh.write(f'{c[0]} {c[1]}\n')

def convert_path(path):
	''' Converts windows like paths to Linux. Reads were originally
	counted on a windows system. This is an ad-hoc solution
	to use the counts on a linux system.
	'''
	if os.name == 'posix':
		path = path.split(':')
		path = path[-1]
		path = path.split('\\')
		path = os.path.join(_drive, *path)
	elif os.name == 'nt':
		pass
	return path

def main():
	pass

if __name__ == '__main__':
	main()
