##########################
#       IMPORTS
##########################

# Standard library
import os
import sys
import json
# Third party imports
import concurrent.futures
import numpy as np
from scipy.spatial import distance
from time import time
import kmer_DB as DB
import settings
#############################################################################
###
#############################################################################
def load_jellyfish_counts(filename):
	'''Reads Jellyfish outputted fasta files and returns
	a dictionary of kmer counts'''
	kd = {}
	test_counter = 0
	with open(filename, 'r') as fh:
		for line in fh:
			test_counter += 1
			if line.startswith('>'):
				value = line.strip().replace('>', '')
				key = next(fh).strip()
				kd[key] = int(value)
	if len(kd) == test_counter:
		return kd
	else:
		raise Exception(f'Something went wrong reading {filename}')

def load_count_matrix(file1, file2, cutoff=None):
	''' Create a matrix of kmer counts for computation
	
	Arguments:
	file1, file2		-- (str) files of kmer counts
	cutoff				-- (int) removes kmers with kmer count < cutoff
	'''
	if file1.endswith('.dump'):
		kd1 = load_jellyfish_counts(file1)
	else:
		with open(file1, 'r') as fh:		
			kd1 = json.load(fh)
	if file2.endswith('.dump'):
		kd2 =  load_jellyfish_counts(file2)
	else:
		with open(file2, 'r') as fh:
			kd2 = json.load(fh)
	if cutoff is not None and cutoff != 0:
		kd1 = filter_low(kd1, cutoff)
		kd2 = filter_low(kd2, cutoff)
	#print('Performing set operations....')
	kmers = set(kd1.keys()).union(set(kd2.keys()))
	#print('Removing failed base calls (N)....')
	invalid_kmers = set([kmer for kmer in kmers if 'N' in kmer])
	kmers = kmers.difference(invalid_kmers)
	#print(f'Removed {len(invalid_kmers)} kmers containing failed base calls.')
	kmer_matrix = np.zeros((len(kmers), 2), dtype=np.uint32)
	#print('Creating matrix....')
	for i, kmer in enumerate(kmers):
		value_one = kd1.get(kmer)
		value_two = kd2.get(kmer)
		if value_one is not None:
			kmer_matrix[i, 0] = np.uint32(value_one)
		if value_two is not None:
			kmer_matrix[i, 1] = np.uint32(value_two)
	return kmer_matrix

#############################################################################
### Differnce calculations
#############################################################################
def denormalize_kmer_spectrum(k_dictionary):
	normalization_factor = k_dictionary.pop('normalization_factor',1)
	for key in k_dictionary:
		k_dictionary[key] = k_dictionary[key]*normalization_factor

def normalize_kmer_spectrum(k_dictionary):
	total = sum(k_dictionary.values())
	for key in k_dictionary:
		k_dictionary[key] = k_dictionary[key]/total
	k_dictionary['normalization_factor'] = total

def spectrum_difference(k_dictionaries):
	'''Return the difference in kmer count for 2 kmer dictionaries'''
	kmers = set()
	for d in k_dictionaries:
		kmers = kmers.union(set(d))
	for kmer in kmers:
		for d in k_dictionaries:
			if kmer not in d:
				d[kmer] = 0
	if 'normalization_factor' in kmers:
		kmers.remove('normalization_factor')
	difference = {}
	for kmer in kmers:
		difference[kmer] = abs(k_dictionaries[0][kmer]-k_dictionaries[1][kmer])
	return difference

def calculate_pairwise_L1_differences(kmer_spectra):
	# L1 calculation
	n = len(kmer_spectra)
	d = []
	matrix = np.zeros((n,n))
	for i in range(n-1):
		for j in range(i+1,n):
			k_dictionaries = [kmer_spectra[i],kmer_spectra[j]]
			d.append(spectrum_difference(k_dictionaries))
	def generator_list(gen_lst):
		for element in gen_lst:
			yield element
	gen_lst = generator_list(d)
	for i in range(n-1):
		for j in range(i+1,n):
			matrix[i,j] = sum(next(gen_lst).values())
	return d, matrix
##########################################################################
### Cosine calculations
########################################################################
def cosine_similarity(v1, v2):

	similarity = v1.dot(v2)
	return similarity

###################################################################
### Entropy calculations
####################################################################
def leibler_kullback(p, q):
	discontinuity_count = 0
	entropy = np.zeros(len(p))
	for i in range(len(p)):
		# Check continuity
		if q[i] == 0:
			if not p[i] == 0:
				discontinuity_count += 1
		else:
			entropy[i] = p[i]*np.log(p[i]/q[i])
	if discontinuity_count > 0:
		print('Leibler-Kullback not defined.')
		print(f'Discontinuity in {discontinuity_count} kmers')
	return np.sum(entropy)

def jensen_shannon(file1, file2, normalize=True, cutoff=None):
	# Not used. Using scipy
	print('Calculating Jensen-Shannon...')
	kmer_matrix = load_count_matrix(file1, file2, cutoff=cutoff)
	if normalize:
		print('Normalizing....')
		kmer_matrix = kmer_matrix / np.linalg.norm(kmer_matrix, axis=0)
	print('Calculating....')
	p = kmer_matrix[:, 0]
	q = kmer_matrix[:,1]
	return distance.jensenshannon(p,q)
#################################################################
# Filters
##############################################################
def filter_low(kmer_dictionary, cutoff):
	'''Removes kmers with counts below cutoff from dictionary
	and returns it.

	Args:
	kmer_dictionary		-- (dic) Dictionary {kmer:count}
	cutoff				-- (int) Remove kmers with lower count than cutoff
	'''
	filtered_dictionary = {}
	removed = 0
	for key, value in kmer_dictionary.items():
		if value < cutoff:
			filtered_dictionary[key] = 0
			removed += 1
		else:
			filtered_dictionary[key] = value
	if removed > 0:
		print(f'Removed {removed} kmers with counts below {cutoff}')		
	return filtered_dictionary
####################################################################
# Runtime
###################################################################
def calculate_distance(
	file1, file2, cutoff, measure=None, dist=None):
	'''Calculates distances between two count files. 

	Arguments:
	file1		-- (str) Full path to count file
	file2		-- (str) Full path to count file
	cutoff		-- (int) Filter out kmers with counts below cutoff
	measure		-- (str, default: None) Measure to calculate
	dist		-- (dict, default: None) distance dictionary
	'''
	measures = {
		'cs':cosine_similarity,
		'js':distance.jensenshannon,
		'cd':distance.cosine,
		'bc':distance.braycurtis
	}
	kmer_matrix = load_count_matrix(file1, file2, cutoff=cutoff)
	if dist is None:
		if measure == 'cs':
			kmer_matrix = kmer_matrix / np.linalg.norm(kmer_matrix, axis=0)
		if measure == 'bc':
			kmer_matrix = kmer_matrix / np.sum(kmer_matrix, axis=0)
		return measures[measure](kmer_matrix[:,0], kmer_matrix[:,1])
	else:
		calculations = {
			'cosine':cosine_similarity,
			'jensenshannon':distance.jensenshannon,
			'braycurtis':distance.braycurtis
		}
		new_values = {}
		for c in calculations:
			if dist[c] is None:
				if c == 'cosine':
					print('Calculating cosine similarity')
					norm_matrix = kmer_matrix / np.linalg.norm(
						kmer_matrix, axis=0)
					new_values[c] = calculations[c](
						norm_matrix[:,0], norm_matrix[:,1])
				if c == 'braycurtis':
					print('Calculating braycurtis dissimilarity')
					norm_matrix = kmer_matrix / np.sum(kmer_matrix, axis=0)
					new_values[c] = calculations[c](
						norm_matrix[:,0], norm_matrix[:,1])
				if c == 'jensenshannon':
					print('Calculating JS dissimilarity')
					new_values[c] = calculations[c](
						kmer_matrix[:,0], kmer_matrix[:,1])
		return new_values

def multi(fffp):
	'''	Helper function for multiprocessing in multi_batch function

	'''
	d = calculate_distance(fffp[0], fffp[1], 0, measure='bc')
	return (fffp[2], d)

def multi_batch(counts, files):
	''' Ad hoc function to recalculated erroneus bray-curtis  distances. Not
	called from other scripts.
	
	files - (lst) List of files to 
	'''
	path = settings.get_setting('count_path')
	'''
	files = DB.find_unupdated_distances()
	files = [f for f in files if 'raw' not in f['file_pair']]
	#files = [f for f in files if 'ofs_10' not in f['file_pair']]
	#files = [f for f in files if 'SV002' in f['file_pair']]
	files = [f for f in files if 'k13' not in f['file_pair']]
	#files = [f for f in files if '_trimmed_' in f['file_pair']]
	'''
	print(f'Counting {counts} of {len(files)} uncounted files...')
	if len(files) == 0:
		print('No files to count')
		return None
	counts = min([counts, len(files)])

	max_workers = 2
	tuple_list = []
	for i in range(counts):
		file1 = os.path.join(path, files[i]['file1'])
		file2 = os.path.join(path, files[i]['file2'])
		file_pair = files[i]['file_pair']
		print(file_pair)
		tuple_list.append([file1, file2, file_pair])
	results = []
	i = 0
	with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
		for r in executor.map(multi, tuple_list):
			i +=1
			sys.stdout.write(f'Finished with {i} of {counts} counts.\r')
			sys.stdout.flush()
			results.append(r)
	for r in results:
		print(r[0], ' ', r[1])
		DB.update_distance2(r[0], r[1])

def main():
	batches = 800
	t0 = time()
	counts = 100
	settings.set_defaults()

	######### FIX BRAY 800 +1600 ###############
	count_path = settings.get_setting('count_path')
	files = os.listdir(count_path)
	k = 13
	trimmed = [f for f in files if f'_k{k}_' in f]
	samples = [f for f in trimmed if '1600000' in f or '800000' in f]
	sv001 = [f for f in trimmed if 'SV001_tr' in f][0]
	sv002 = [f for f in trimmed if 'SV002_tr' in f][0]
	file_pairs = []
	counter = 1
	for sa in samples:
		file_pairs.append([sv001,sa])
		file_pairs.append([sv002,sa])
	print(len(file_pairs))
	t0 = time()
	file_pairs = [{
		'file_pair':DB.canonical_key(fp[0],fp[1]),
		'file1':fp[0],
		'file2':fp[1]
		} for fp in file_pairs]
	multi_batch(counts, file_pairs)
	###############################################
	#for b in range(batches):
	#	print(f'Starting batch {b+1} of {batches}')
	#	multi_batch(counts)
	#	t = ((time()-t0)/(b+1))*(batches - 1 -b)
	#	print(f'Estimated time remaining {t}')
	#	if t - t0 > 25000:
	#		break
	#t1 = time() - t0
	#print(f'Running time: {t1}')

if __name__ == "__main__":
	main()
