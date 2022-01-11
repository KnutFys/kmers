# Standard library imports
import sys
import os
import json
from time import time
# Third party imports
import pandas as pd
import numpy as np
import skbio
from sklearn.manifold import MDS
from matplotlib import pyplot as plt
import seaborn as sns
# Local imports
import kmer_DB as DB
import outputter
import settings
import features as feat
import calculations

###################################################
### Helper functions
##################################################
def load_kmer_dictionary(filename):
	'''Returns a kmer_dictionary read from file'''

	path = settings.get_setting('count_path')
	filename = os.path.join(path, filename)
	with open(filename, 'r') as fh:
		kmer_dict = json.load(fh)
	return kmer_dict

def cosine_similarity_to_angular(similarity):
	'''Converts cosine similarity to angular distance'''

	if isinstance(similarity, list):
		for s in similarity:
			s['cosine'] = cosine_similarity_to_angular(s['cosine'])
		return similarity
	else:
		return (np.arccos(similarity))*(2/np.pi)

def r_export(dfs,names=None):
	'''Saves pandas dataframes to csv for import in R'''
	
	def save(df, name):
		path = '/media/knut/Storage/Rfiles'
		df.to_csv(os.path.join(path, name))
	if isinstance(dfs, list):
		if len(dfs) != len(names):
			return False
		else:
			for df, name in zip(dfs, names):
				save(df, name)
	if isinstance(dfs, dict):
		for key, item in dfs.items():
			save(item, f'{key}.csv')
	else:
		save(dfs, names)

###################################################
### Distance Matrix creation
##################################################
def one_from_each():
	'''Returns three distance matrices. Angular, Jensen-Shannon and Bray-Curtis.
	based on one random subsample for each sample.
	'''
	def get_rand_dist(samples, n=1, force_same_strand=False):
		include = [
			f'{samples[0]}',
			f'{samples[1]}',
			f'\_k{k}\_',
			f'\_{sample_size}\_'
		]
		exclude = ['\_raw\_', f'\_ofs\_{10}']
		if force_same_strand:
			exclude += ['\_R2\_']
		rows = DB.filter_distances(include=include, exclude=exclude)
		return rng.choice(rows, size=n)[0]

	rng = np.random.default_rng()

	sample_size = settings.get_setting('sample_size')
	k = settings.get_setting('k')
	
	samples = [s['sample_id'] for s in DB.get_sample_ids()]
	angular = np.zeros((len(samples), len(samples)))
	bray_curtis = np.zeros((len(samples), len(samples)))
	jensen_shannon = np.zeros((len(samples), len(samples)))
	for i in range(len(samples)-1):
		for j in range(i+1, len(samples)):
			dist = get_rand_dist([samples[i], samples[j]])
			angular[i, j] = cosine_similarity_to_angular(
				dist['cosine'])
			jensen_shannon[i, j] = dist['jensenshannon']
			bray_curtis[i, j] = dist['braycurtis']
			angular[j, i] = angular[i, j]
			jensen_shannon[j, i] = jensen_shannon[i, j]
			bray_curtis[j, i] = bray_curtis[i, j]
	angular = skbio.stats.distance.DistanceMatrix(
		angular, samples)
	jensen_shannon = skbio.stats.distance.DistanceMatrix(
		jensen_shannon, samples)
	bray_curtis = skbio.stats.distance.DistanceMatrix(
		bray_curtis, samples)

	return angular, jensen_shannon, bray_curtis

def full_counts():
	'''Returns three distance matrices. Angular, Jensen-Shannon and Bray-Curtis.
	based on all reads for each sample.
	'''
	k = settings.get_setting('k')
	count_files = DB.get_count_files(kmer=k, strand='both')
	count_files = sorted(count_files, key=lambda f: f.split('_')[0])
	sample_ids = [f.split('_')[0] for f in count_files]
	angular = np.zeros((len(count_files), len(count_files)))
	bray_curtis = np.zeros((len(count_files), len(count_files)))
	jensen_shannon = np.zeros((len(count_files), len(count_files)))
	for i in range(len(count_files)-1):
		for j in range(i+1, len(count_files)):
			dist = DB.get_distance(count_files[i], count_files[j])
			angular[i, j] = cosine_similarity_to_angular(
				dist['cosine'])
			jensen_shannon[i, j] = dist['jensenshannon']
			bray_curtis[i, j] = dist['braycurtis']
			angular[j, i] = angular[i, j]
			jensen_shannon[j, i] = jensen_shannon[i, j]
			bray_curtis[j, i] = bray_curtis[i, j]
	angular = skbio.stats.distance.DistanceMatrix(
		angular, sample_ids)
	jensen_shannon = skbio.stats.distance.DistanceMatrix(
		jensen_shannon, sample_ids)
	bray_curtis = skbio.stats.distance.DistanceMatrix(
		bray_curtis, sample_ids)
	return angular, jensen_shannon, bray_curtis

def pole_distances():
	'''Returns a distance matrix based on angular distance with an ideal
	"pole" sample added.
	'''
	k = settings.get_setting('k')
	count_files = DB.get_count_files(kmer=k)
	count_files = [f for f in count_files if 'trimmed' in f or 'pole' in f]
	count_files = sorted(count_files, key=lambda f: f.split('_')[0])
	sample_ids = [f.split('_')[0] for f in count_files]
	angular = np.zeros((len(count_files), len(count_files)))
	for i in range(len(count_files)-1):
		for j in range(i+1, len(count_files)):
			dist = DB.get_distance(count_files[i], count_files[j])
			angular[i, j] = cosine_similarity_to_angular(
				dist['cosine'])
			angular[j, i] = angular[i, j]
	angular = skbio.stats.distance.DistanceMatrix(
		angular, sample_ids)
	return angular

def feature_matrix(df, features, standardize=False, group_loc=False):
	'''Returns distance matrix based on selected features'''

	#Standardize
	if standardize:
		df[features] = df[features].apply(lambda x: (x - x.mean())/x.std())
	#Create distance matrix
	size = df.shape[0]
	distance = np.zeros((size, size))
	for i in range(size-1):
		for j in range(i+1, size):
			x = [df.iloc[[i]][f] for f in features]
			y =	[df.iloc[[j]][f] for f in features]
			x = np.array(x)
			y = np.array(y)
			distance[i, j] = np.linalg.norm(x-y)
			distance[j, i] = distance[i, j]
	sample_ids = df['sample_id'].tolist()
	distance = skbio.stats.distance.DistanceMatrix(
		distance, sample_ids)
	return distance

def all_matrices(standardize=False):
	'''Returns a dictionary with all distance matrices.'''

	default_features = [
		'toc',
		'doc',
		'tn',
		'dn',
		'totP',
		'dp',
		'altitude',
		'dnavolume',
		'rnavolume',
		'temperature',
		'cond',
		'pH',
		'pigmentvolume'
	]
	names = ['Angular', 'Jensen-Shannon', 'Bray-Curtis']
	matrices = {name:m for name, m in zip(names, list(full_counts()))}
	matrices['Spatial'] = spatial_distance()
	for f in default_features:
		matrices[f] = feature_matrix([f], standardize=standardize)
	return matrices
###################################################
### Distance Matrix calculations
##################################################
def multidimensional_scaling(dm, samples, n=2):
	'''Tries to find a nD representation of a distance matrix
	DEPRECATED

	args:
	matrix		-- (obj) skbio DistanceMatrix
	samples		-- (lst) List of sample names for annotations
	n			-- (int, default 2) Dimensions used for embedding
	'''
	#Fetch metadata
	con = DB.connect()
	df = pd.read_sql_query('SELECT * FROM sample_table', con)
	con.close()

	#Perform scaling
	embedding = MDS(
		n_components=n,
		dissimilarity='precomputed',
		metric=False,
		n_init=100
		)
	transform = embedding.fit_transform(dm.data)
	print('Stress: ', embedding.stress_)
	#Plot
	fig = plt.figure(1)
	ax = plt.axes([0., 0., 1., 1.])	
	X = transform[:, 0]
	Y = transform[:, 1]
	sns.scatterplot(x=X, y=Y, data=df, hue='loacality')
	for x, y, s in zip(X, Y, samples):
		point = f'{s}'
		plt.annotate(point, (x,y))
	
	plt.show()

def compare_mantel(features=None, distance_measure='Bray-Curtis'):
	'''DEPRECATED'''

	default_features = [
		'toc',
		'doc',
		'tn',
		'dn',
		'totP',
		'dp',
		'altitude',
		'dnavolume',
		'rnavolume',
		'temperature',
		'cond',
		'pH',
		'pigmentvolume'
	]
	if features is None:
		features = default_features
	exclude_samples = set()
	matrices = all_matrices()
	matrices['Many features'] = feature_matrix(
		features=['temperature','pH'])
	features += ['Many features']
	tests = [[distance_measure, 'Spatial']]
	tests += [[distance_measure, f] for f in features]
	tests += [['Angular', f] for f in features]
	tests += [['Jensen-Shannon', f] for f in features]
	tests += [['Spatial', f] for f in features]
	tests += [['Angular', 'Spatial']]
	tests += [['Jensen-Shannon', 'Spatial']]
		
	summary = {
		'Matrix A':[],
		'Matrix B':[],
		'Correlation':[],
		'p-value':[],
		'N':[]
	}
	
	for t in tests:
		ids_1 = set(matrices[t[0]].ids)
		ids_2 = set(matrices[t[1]].ids)
		ids = list(ids_1.intersection(ids_2).difference(exclude_samples))
		d1 = matrices[t[0]].filter(ids, strict=True)
		d2 = matrices[t[1]].filter(ids, strict=True)
		val = (skbio.stats.distance.mantel(d1, d2))
		summary['Matrix A'].append(t[0])
		summary['Matrix B'].append(t[1])
		summary['Correlation'].append(val[0])
		summary['p-value'].append(val[1])
		summary['N'].append(val[2])
	df = pd.DataFrame.from_dict(summary)
	c = ['r','g','b','y']
	fig, ax = plt.subplots()
	for i, f in enumerate(['pH','temperature', 'tn','dn']):
		linear(matrices['Angular'], matrices[f], ax,c=c[i])
	plt.show()
	print(df.sort_values(['Matrix A', 'p-value']))
	
###################################################
### Plots
##################################################
def plot_all_offsets_and_ks():
	''' Create plot for thesis, Discriminatory capacity'''

	ks = [5,7,9,11,13]
	offset = 0
	strands = ['R1', 'R2', None]
	hue = ['Ang','JS','BC','Ang2','JS2','BC2']
	for strand in strands:
		dataframes = []
		for k in ks:
			dataframes.append([all_vs_full_DF(
				k=k, offset=offset, show=False, raw=False,strand=strand),
			f'Kmer size {k}'])
		if strand is None:
			strand = 'both'
		title = f'Samples from {strand} strands'
		filename = f'All_K_strand_{strand}.png'
		outputter.multiplotter(
			dataframes, hue, mode='line', rows=3, title=title, filename=filename)
		outputter.multiplotter(
			dataframes, hue, mode='box', rows=3, title=title, filename=filename)

def prepare_plot(ds, offset=None, save=False):
	filename = None
	if save:
		if offset is None:
			filename = f'PlotADJSBC_SV001SV002_k{k}.png'
		else: 
			filename = f'PlotADJSBC_SV001SV002_k{k}_{offset}.png'
	hue = ['Ang','JS','BC','Ang2','JS2','BC2']
	outputter.plotter(ds, hue, filename=filename, mode=mode)

##################################################
### kmer distribution
##################################################
def kmer_bar():
	files = DB.get_count_files(offset=0, kmer=11)
	files = [f for f in files if '_trimmed_' in f]
	f = files[0]
	print(f)
	kmer_dict = load_kmer_dictionary(f)	
	kmer_counts = list(kmer_dict.values())	
	sns.histplot(data=kmer_counts, bins=20, log_scale=True)
	plt.show()

def alphas():
	files = DB.get_count_files(offset=0, kmer=11)
	p = [f for f in files if 'pole' in f][0]
	files = [f for f in files if '_trimmed_' in f]
	distances = DB.get_distances_to_file(p)
	distances = [d for d in distances if 'trimmed' in d['file_pair']]
	shannon = []
	js = []
	cos = []
	bray = []
	simpsons = []
	all_measures = []
	for f in files:
		kmer_dict = load_kmer_dictionary(f)	
		kmer_counts = [kmer_dict[key] for key in kmer_dict if 'N' not in key]
		sh = skbio.diversity.alpha.shannon(np.array(kmer_counts))
		shannon.append(sh)
		sm = skbio.diversity.alpha.simpson_e(np.array(kmer_counts))
		simpsons.append(sm)
		d = [d for d in distances if f in d['file_pair']][0]
		cos.append(d['cosine'])
		bray.append(d['braycurtis'])
		ang = 1 - cosine_similarity_to_angular(d['cosine'])
		js.append(d['jensenshannon'])
		print(d['cosine'], ' ', sh, ' ', sm)
		sample_id = f.split('_')[0]
		all_measures.append([d['cosine'], d['jensenshannon'],
		d['braycurtis'], sm, sh, ang, sample_id])
	sort_by_cos = sorted(all_measures, key=lambda x: x[0])
	#print(sort_by_cos)
	df = pd.DataFrame(sort_by_cos, columns=[
		'Cosine', 'Js', 'BrayCurtis','Simpson', 'Shannon', 'Angular', 'Sample'])
	df.to_csv('alphadiversity.csv')
	sns.lineplot(x='Cosine', y='Simpson', data = df)
	plt.show()

def plotalpha(feature_to_plot=None):
	if feature_to_plot is None:
		feature_to_plot = 'dn'
	df = pd.read_csv('alphadiversity.csv')
	df.set_index('Sample')
	rows = DB.get_sample_ids()
	samples = pd.Series({s['sample_id']:s[feature_to_plot] for s in rows})
	df = df.set_index('Sample')
	df = pd.concat([df, pd.DataFrame(samples, columns=[feature_to_plot])], axis=1)
	melted = pd.melt(df, id_vars=feature_to_plot, value_vars=[
		'Simpson', 'Cosine', 'Angular'])
	sns.lineplot(data=melted,x=feature_to_plot, y='value', hue='variable')
	plt.show()
###################################################
### Spreadsheets
##################################################
def convert_to_sheets(dfs, k, offset, sheets, trimmed_or_raw):
	dist_within_sample = dfs[:3]
	dist_between_samples = dfs[3:]
	names = [
		f'Angular distance, k{k}, ofs{offset}',
		f'Jensen-Shannon divergence, k{k}, ofs{offset}',
		f'Bray-Curtis distance, k{k}, ofs{offset}'
		]
	sheet_name = f'{trimmed_or_raw} Distances k{k}, SV001-SV001'
	sheets[sheet_name] = []
	for d, n in zip(dist_within_sample, names):
			sheets[sheet_name].append((d, d.describe(), n))
	sheet_name = f'{trimmed_or_raw} Distances k{k}, SV001-SV002'
	sheets[sheet_name] = []
	for d, n in zip(dist_between_samples, names):
		sheets[sheet_name].append((d, d.describe(), n))
	return sheets

def save_sheets_all_offsets_and_ks():
	ks = [5,7,9,11,13]
	offsets = [0]
	trimmed_or_raw = ['Trimmed']
	sheets = {}
	for trim in trimmed_or_raw:
		for offset in offsets:
			dataframes = []
			for k in ks:
				if trim == 'Raw':
					dfs = all_vs_full_DF(k=k, offset=offset, show=False, raw=True)
				else:
					dfs =all_vs_full_DF(k=k, offset=offset, show=False, raw=False)
				sheets = convert_to_sheets(dfs, k, offset, sheets, trim)

	outputter.excel_writer(sheets, 'AllSheetsSV1SV2.xlsx')
	
def all_vs_full_DF(k=11, offset=0, show=True, raw=False, strand=None):
	'''Returns a list of dataframes with distances for sample SV001 and SV002

	Arguments:
	k			-- (int, default 11) kmer size
	offset		-- (int, default 10) only use counts with given offset 
	show		-- (bol, default True) If true, prints dataframes to stout
	'''

	sample_sizes = [
		10000, 20000, 50000, 100000,
		200000, 400000, 800000, 1600000]
	sv1_cs = {}
	sv1_js = {}
	sv1_bc = {}
	sv2_cs = {}
	sv2_js = {}
	sv2_bc = {}

	exclude = []
	if offset is not None:
		include = [f'0___\_k{k}\_ofs\_{offset}']
	else:
		include = [f'\_k{k}']
	if raw:
		include += ['\_raw\_']
	else:
		exclude += ['\_raw\_']
	if strand is not None:
		include += [f'\_{strand}\_']
	for ss in sample_sizes:
		include1 = include + [f'SV001\_trimmed', f'\_{ss}\_'] 
		include2 = include + [f'SV002\_trimmed', f'\_{ss}\_']
		exclude1 = exclude
		exclude2 = exclude
		d1 = DB.filter_distances(include=include1, exclude=exclude1)
		sv1_cs[ss] = [cosine_similarity_to_angular(d['cosine']) for d in d1]
		sv1_js[ss] = [d['jensenshannon'] for d in d1]
		sv1_bc[ss] = [d['braycurtis'] for d in d1]
		d2 = DB.filter_distances(include=include2, exclude=exclude2)
		sv2_cs[ss] = [cosine_similarity_to_angular(d['cosine']) for d in d2]
		sv2_js[ss] = [d['jensenshannon'] for d in d2]
		sv2_bc[ss] = [d['braycurtis'] for d in d2]
		print(len(d1),' : ', len(d2))
	# Descriptions for writing to screen
	descript = [
		'Angular distance between SV001 and SV001',
		'Jensen-Shannon distance between SV001 and SV001',
		'Bray-Curtis distance between SV001 and SV001',
		'Angular distance between SV001 and SV002',
		'Jensen-Shannon distance between SV001 and SV002',
		'Bray-Curtis distance between SV001 and SV002'
	]
	# Convert dictionaries to dataframes
	ds = [sv1_cs, sv1_js, sv1_bc, sv2_cs, sv2_js, sv2_bc]
	ds = [pd.DataFrame.from_dict(d) for d in ds]

	if show:
		for d, des in zip(ds, descript):
			print(des, '\n', '-'*80)
			#print(d, '\n', '-'*80)
			print(d.describe(),'\n', '-'*80)
	return ds

def anosim(dm, df, col):
	res = skbio.stats.distance.anosim(dm, df, column = col)
	return res
###################################################
### Runtime
##################################################
def export_all_dm_R():
	'''Creates distannce matrices from sub samples and saves csv'''

	path = './dms'
	ks = [5, 7, 9, 11, 13]
	sample_sizes = [100000]
	for k in ks:
		for ss in sample_sizes:
			settings.set_setting('k', k)
			settings.set_setting('sample_size', ss)
			ang, js, bc = one_from_each()
			ang.to_data_frame().to_csv(os.path.join(path, f'{k}_{ss}ang.csv'))
			js.to_data_frame().to_csv(os.path.join(path, f'{k}_{ss}_js.csv'))
			bc.to_data_frame().to_csv(os.path.join(path, f'{k}_{ss}_bc.csv'))

def evenness_as_richness():
	'''Plot evenness of kmers as a function of number of combined
	communities'''

	#Select distance measure and kmer length
	measure = 'cosine'
	k = 11
	#Get pole count file
	count_path = settings.get_setting('count_path')
	files = os.listdir(count_path)
	pole = [f for f in files if f'k{k}' in f and 'pole' in f][0]
	#Get count files for SV001, SV001+SV002, SV001+SV002+SV003 + ...
	outfiles = [f'Combined{i+1}Samples.count' for i in range(25)]
	number_of_samples = range(26)[1:]
	#Get distances to pole
	distances = DB.filter_distances(include=['pole'])
	distances = [d for d in distances if 'Combined' in d['file_pair']]
	angular = [(1-cosine_similarity_to_angular(d['cosine'])) for d in distances]
	bc = [d['braycurtis'] for d in distances]
	js = [d['jensenshannon'] for d in distances]
	plt.plot(number_of_samples, angular)
	plt.plot(number_of_samples, js)
	plt.plot(number_of_samples, bc)
	plt.show()

def test():
	def subplotter(n):
		''' Returns double columned plot axes'''
		rows = int(n//2)
		if rows < n/2:
			rows += 1
		fig, ax = plt.subplots(rows, 2)
		return fig, ax

	default_features = [
		'toc',
		'doc',
		'tn',
		'dn',
		'totP',
		'dp',
		'altitude',
		'dnavolume',
		'rnavolume',
		'temperature',
		'cond',
		'pH',
		'pigmentvolume'
	]
	#####################################################
	#alphas()
	#path = '/media/knut/Storage/'
	#angular = pole_distances()
	#df = angular.to_data_frame()
	#df.to_csv('AngularPole11.csv')
	#print(angular)
	#plotalpha(feature_to_plot='glacial')
	#settings.set_setting('k', 13)
	#an, js ,bc = full_counts()
	#an.to_data_frame().to_csv('AngularFull13.csv')
	#js.to_data_frame().to_csv('JensenShannon13.csv')
	#bc.to_data_frame().to_csv('BrayCurtisFull13.csv')
	#plot_all_offsets_and_ks()
	#evenness_as_richness()
	alphas()
	plotalpha(feature_to_plot='temperature')

def argument_handler():
	if '-k' in sys.argv:
		try:
			k = int(sys.argv[sys.argv.index('-k')+1])
			settings.set_setting('k', k)
		except IndexError as e:
			print(e)
	if '-o' in sys.argv:
		try:
			offset = int(sys.argv[sys.argv.index('-o')+1])
			settings.set_setting('offset', offset)
		except IndexError as e:
			print(e)

if __name__ == '__main__':
	settings.set_defaults()
	try:
		command_line_functions = {
			'test':test,
			'mantel':compare_mantel
		}
		choose_function = sys.argv.pop(1)
	except IndexError as e:
		print(e)
		print('No arguments supplied. Test...')
	argument_handler()
	command_line_functions[choose_function]()

