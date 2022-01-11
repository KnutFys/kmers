###########################
# Deprecated
###########################


# Standard library imports
import sys
import os
from math import radians
# Third party imports
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.metrics.pairwise import haversine_distances
import matplotlib.pyplot as plt
import seaborn as sns
import skbio
# Local imports
import kmer_DB as DB
import settings

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

def r_export(dfs,names):
	'''Saves pandas dataframes to csv for import in R'''

	def save(df, name):
		path = '/media/knut/Storage/Rfiles'
		dfs.to_csv(os.path.join(path, names))
	if isinstance(dfs, list):
		if len(dfs) != len(names):
			return False
		else:
			for df, name in zip(dfs, names):
				save(df, name)
	else:
		save(dfs, names)

def spatial_distance(samples=None, normalize=True):
	'''Returns a matrix of the physical distance between sample locations.

	Arguments:
	samples		-- (lst) List of samples. Assumes atleast 2 samples
	'''

	if samples is None:
		samples = [r['sample_id'] for r in DB.get_sample_ids()]
	condition = '('
	for s in samples:
		condition += '?,'
	condition = condition[:-1] + ')'
	command = f'''
	SELECT 
		north, east
	FROM
		sample_table
	WHERE
		sample_id IN {condition}
	'''
	values = tuple(samples)
	con = DB.connect()
	with con:
		c = con.cursor()
		c.execute(command, values)
		rows = c.fetchall()
	con.close()
	latitude_longditude = [[row['north'], row['east']] for row in rows]
	earth_radius = 6371
	latitude_longditude_in_rad = [[radians(lat_lon[0]), radians(
		lat_lon[1])] for lat_lon in latitude_longditude]
	distances = haversine_distances(latitude_longditude_in_rad)
	distances = distances * earth_radius
	if normalize:
		distances = distances / np.max(distances)
	
	distance_matrix = skbio.stats.distance.DistanceMatrix(
			distances, samples)
	return distance_matrix

def standatrdized_all_features(transform='log'):
	''' Standardized features'''
	#Feature types
	#Continuous -> standardize or log
	std_features = ['temperature', 'pH', 'north', 'east', 'cond']
	log_features = ['tn', 'dn', 'totP', 'dp', 'toc','doc','cond']
	scale_features = ['temperature', 'altitude'] + log_features
	#No difference -> drop
	equal_features = ['gasvolume', 'gasheadspace']
	#Categorical like?
	class_feaures = []
	#Exclude
	excluded_features = [
	'cond',
	'pH',
	'name',
	'loacality',
	'north',
	'east',
	'dnavolume',
	'rnavolume',
	'pigmentvolume',
	'comments',
	'date'
	]
	#Fetch metadata
	con = DB.connect()
	df = pd.read_sql_query('SELECT * FROM sample_table', con)
	con.close()

	df.at[14, 'tn'] = 0.0001
	if transform == 'log':
		df[log_features] = df[log_features].apply(lambda x: (np.log(x)))
		#df[scale_features] = df[scale_features].apply(lambda x: (x/np.max(x)))
	if transform == 'std':
		df[std_features] = df[std_features].apply(lambda x: (x - x.mean())/x.std())

	df = df.drop(equal_features, axis=1)
	df = df.drop(excluded_features, axis=1)
	#Change index to sample id
	df.set_index('sample_id', inplace=True, drop=True)

	#Drop samples with NaN
	#df2 = df.dropna(subset=c_features)
	df.to_csv(f'/media/knut/Storage/Rfiles/metadata_{transform}.csv')

def create_feature_matrix(features=None, standardize=False):
	'''Reaturns distance matrix based on selected features and a list
	of sample ids
	'''
	
	if features is None:
		features = default_features

	samples = get_samples_with_features(features=features)

	#Standardize features
	min_maxed_features = {}
	for f in features:
		min_maxed_features[f] = min_max([s[f] for s in samples])
	
	#Create matrix
	distance = np.zeros((len(samples),len(samples)))
	for i in range(len(samples)-1):
		for j in range(i+1, len(samples)):
			x = [min_maxed_features[f][i] for f in features]
			y =	[min_maxed_features[f][j] for f in features]
			x = np.array(x)
			y = np.array(y)
			distance[i, j] = np.linalg.norm(x-y)
			distance[j, i] = distance[i, j]
	sample_ids = [s['sample_id'] for s in samples]

	return distance, sample_ids

def create_feature_dataframe(features=None):
	'''Returns a dataframe with all samples with valid
	feature values
	'''
	if features is None:
		features = default_features
	samples = get_samples_with_features(features=features)
	data = {'Id':[s['sample_id'] for s in samples]}
	for feature in features:
		data[feature] = [s[feature] for s in samples]
	df = pd.DataFrame(data)
	df = df.set_index('Id')
	
	return df

def run_PCA(features=None):
	n_components = 2
	df = create_feature_dataframe()

	y = df.index.values.tolist()
	ssc = StandardScaler()
	ssc.fit(df)
	X = ssc.transform(df)
	pca = PCA(n_components=n_components)
	pca.fit(X)
	X_transformed = pca.fit_transform(X)

	biplot(
		X_transformed[:,0:2],
		np.transpose(pca.components_[0:2,:]),
		y
	)
	X_transformed = pca.fit_transform(X)
	#print(pca.explained_variance_ratio_)
	#print(pca.components_)

	#fig, axes = plt.subplots(1,2)
	#axes[0].scatter(X[:,0], X[:,1])
	#axes[0].set_xlabel('x1')
	#axes[0].set_ylabel('x2')
	#axes[0].set_title('Before PCA')
	#axes[1].scatter(X_transformed[:,0], X_transformed[:,1])
	#b = [X, X_transformed]
	#for a , p in zip(axes, b):
	#	for i, l in enumerate(y):
	#		a.annotate(l, (p[i,0], p[i,1]))
	#axes[1].set_xlabel('PC1')
	#axes[1].set_ylabel('PC2')
	#axes[1].set_title('After PCA')
	plt.show()

def biplot(score, coeff , y):
	'''
	Original Author: Serafeim Loukas, serafeim.loukas@epfl.ch
	Inputs:
	   score: the projected data
	   coeff: the eigenvectors (PCs)
	   y: the class labels
	'''    
	xs = score[:,0] # projection on PC1
	ys = score[:,1] # projection on PC2
	n = coeff.shape[0] # number of variables
	plt.figure(figsize=(10,8), dpi=100)
	classes = np.unique(y)

	for s,l in enumerate(classes):
		plt.scatter(xs[y==l],ys[y==l])
	for i in range(n):
		#plot as arrows the variable scores (each variable has a score for PC1 and one for PC2)
		plt.arrow(
			0, 0, coeff[i,0],
			coeff[i,1], color='k',
			alpha = 0.9,linestyle='-',
			linewidth = 1.5, overhang=0.2
			)
		plt.text(
			coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1),
			color = 'k', ha = 'center', va = 'center',fontsize=10)

	plt.xlabel("PC{}".format(1), size=14)
	plt.ylabel("PC{}".format(2), size=14)
	limx= int(xs.max()) + 1
	limy= int(ys.max()) + 1
	plt.xlim([-limx,limx])
	plt.ylim([-limy,limy])
	plt.grid()
	plt.tick_params(axis='both', which='both', labelsize=14)

def create_h_cluster(dm, clusters=2):
	'''Clusters samples by distance matrix. Resturn cluster labels

	'''
	if not isinstance(dm, list):
		dm = [dm]
	
	cluster = AgglomerativeClustering(
			affinity='precomputed', n_clusters=clusters, linkage='complete')
	
	results = []
	for i, d in enumerate(dm):
		cluster.fit(d.data)
		results.append(cluster.labels_)
	return results

def spatial_plot(include_samples, plotname, clusters, ax=None):
	'''Plots clusters on map'''

	columns = DB.get_table_cols('sample_table')
	sample_df = pd.DataFrame(DB.get_sample_ids(), columns=columns)
	#Drop samples not specified
	sample_df = sample_df[sample_df.sample_id.isin(include_samples)]
	sample_df['hue'] = clusters
	coords = pd.melt(sample_df, id_vars=['east', 'north'], value_vars=['hue'])
	sns.scatterplot(x='east', y='north', data=coords, hue='value', ax=ax)
	if ax is None:
		plt.show()
	else:
		ax.set_title(plotname)
		
def create_s_cluster(dm, clusters=2, delta=1):
	'''Clusters samples by distance matrix. Resturn cluster labels

	'''
	if not isinstance(dm, list):
		dm = [dm]
	
	cluster = SpectralClustering(
			affinity='precomputed', n_clusters=clusters, assign_labels='discretize')
	
	results = []
	for i, d in enumerate(dm):
		similarity = np.exp(- d.data ** 2 / (2. * delta ** 2))
		cluster.fit(similarity)
		results.append(cluster.labels_)
	return results

def create_dendrogram(dms):
	'''Plots dendrogram based on distance matrix'''

	def dendrify(dm,ax):
		X = dm.condensed_form()
		linked = linkage(X, 'single', optimal_ordering=True)
		labelList = list(dm.ids)
		dendrogram(linked,
			orientation='top',
			labels=labelList,
			distance_sort='descending',
			show_leaf_counts=True,
			ax=ax)

	if not isinstance(dms, list):
		fig, ax = plt.subplots(1)
		dendrify(dms, ax)
	else:	
		fig, ax = plt.subplots(len(dms))
		for i, dm in enumerate(dms):
			dendrify(dm, ax[i])
	plt.show()

if __name__ == '__main__':
	#print(create_feature_dataframe())
	samples = DB.get_sample_ids()
	samples = [s['sample_id'] for s in samples]
	spatial_plot()