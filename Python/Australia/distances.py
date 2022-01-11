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

def add_nuc(kmers,alfabet):
    # Takes a list of kmers and returns a list of k+1-mers extended
    # by alphabet
    new_kmers = []
    for s in kmers:
        for n in alfabet:
            new_kmers.append(s+n)
    return new_kmers

def create_kmer_index(k=3):
    # Creates a complete list of k-mers based on the alphabet alf
    alf = set(['A','C','G','T'])
    translation = 'string'.maketrans('ACGT','TGCA')
    kmers = [x for x in alf]
    for i in range(k-1):
        kmers = add_nuc(kmers,alf)
    reverse = [x[::-1] for x in kmers]
    reverse = [x.translate(translation) for x in reverse]
    kmers = [sorted([x[0],x[1]])[0] for x in zip(kmers,reverse)]
    kmers = list(set(kmers))
    kmers = sorted(kmers, key=str.upper)
    return kmers
#############################################################################
### Load kmer tables
#############################################################################
def load_jellyfish_counts(filename):
	'''Reads Jellyfish outputted fasta files and returns
	a dictionary of kmer counts
	
	args:
	filename - (str) Path to jellyfish dump file
	'''
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

def load_count_matrix(file1, file2, cutoff=None, from_dict=False):
	''' Create a matrix of kmer counts for computation
	
	Arguments:
	file1, file2		-- (str) files of kmer counts
	cutoff				-- (int) removes kmers with kmer count < cutoff
	'''
	#Check if dictionaries are passed
	if isinstance(file1, dict):
		kd1 = file1
	else:
		kd1 = load_jellyfish_counts(file1)	
	if isinstance(file2, dict):
		kd2 = file2
	else:
		# Check if pole is to be used as second file
		if 'pole' in file2:
			# Creates a sample of evenly distributed k-mers
			try:
				k = file2.replace('pole','')
				k = int(k)
				kd2 = {kmer:1 for kmer in create_kmer_index(k)}
			except:
				kd2 =  load_jellyfish_counts(file2)
		else:
			kd2 =  load_jellyfish_counts(file2)

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

##########################################################################
### Cosine calculations
########################################################################
def cosine_similarity(v1, v2):
	''' Calculates cosine similarity for unit vectors IE 
	dot product '''

	similarity = v1.dot(v2)
	return similarity

def cosine_similarity_to_angular(similarity):
	'''Converts cosine similarity to angular distance'''

	if isinstance(similarity, list):
		for s in similarity:
			s['cosine'] = cosine_similarity_to_angular(s['cosine'])
		return similarity
	else:
		return (np.arccos(similarity))*(2/np.pi)
###################################################################
### Entropy calculations NOT USED
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
	# Not used
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
# Distance calculation battery
###################################################################
def calculate_distance(file1, file2, cutoff=None):
	'''Calculates distances between two count files. 

	Arguments:
	file1		-- (str) Path to jellyfish dump file
	file2		-- (str) Path to jellyfish dump file
	cutoff		-- (int) Ignore kmers appearing less than 'cutoff' times
	from_dict   -- (bol) If true, uses kmer dict instead of jellyfish files
	'''

	kmer_matrix = load_count_matrix(file1, file2, cutoff=cutoff)
	
	calculations = {
		'cosine':cosine_similarity,
		'jensenshannon':distance.jensenshannon,
		'braycurtis':distance.braycurtis
	}
	calculated_distances = {}
	for c in calculations:
		if c == 'cosine':
			# Normalize to unit vectors
			norm_matrix = kmer_matrix / np.linalg.norm(
				kmer_matrix, axis=0)
			calculated_distances[c] = calculations[c](
				norm_matrix[:,0], norm_matrix[:,1])
			# Add angular distance to dict
			calculated_distances['angular'] = cosine_similarity_to_angular(
				calculated_distances[c])
		if c == 'braycurtis':
			# Normalize to relative abundance/frequency vectors
			norm_matrix = kmer_matrix / np.sum(kmer_matrix, axis=0)
			calculated_distances[c] = calculations[c](
				norm_matrix[:,0], norm_matrix[:,1])
		if c == 'jensenshannon':
			# Normalization automatically carried out by scipy.spatial
			calculated_distances[c] = calculations[c](
				kmer_matrix[:,0], kmer_matrix[:,1])
	return calculated_distances

##########################################
# Multi core
##########################################
def multi(file_pair):
	d = calculate_distance(file_pair[0], file_pair[1])
	return (file_pair, d)

def multi_batch(file_pairs, batchsize=1000, max_workers=2):
	''' Calculates Bray-Curtis, Jensen shannon and
	cosine distances for a list of file pairs

	Args:
	file_pairs - (lst) List of pairs of paths to jellyfish dump files
	batchsize - (int, default 1000) max number of pairs to calculate
	'''

	print(f'Calculating distances for {batchsize} of {len(file_pairs)}')
	print('file pairs...')
	if len(file_pairs) == 0:
		print('No distances to calculate')
		return None
	batchsize = min([batchsize, len(file_pairs)])
	file_pairs = file_pairs[:batchsize]
	
	results = []
	i = 0
	with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
		for r in executor.map(multi, file_pairs):
			i +=1
			sys.stdout.write(f'Finished with {i} of {batchsize} calculations.\r')
			sys.stdout.flush()
			results.append(r)
	return results

def main():
	pass

if __name__ == "__main__":
	main()
