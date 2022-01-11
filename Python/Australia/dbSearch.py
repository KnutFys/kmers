# Create connection to Bioplatforms Australia, based on example at 
# https://usersupport.bioplatforms.com/programmatic_access.html

# Standard library imports
import os
import json
import time
import hashlib
# Third party imports
import ckanapi
import requests
# Local imports

# Settings
_apikey = None # Needed to search database, but is personal
_downloads = os.path.join(os.getcwd(), 'Downloads')
_counts = os.path.join(_downloads, 'Counts')


def check_integrity(resource, filename):
	'''
	Checks file integrity using md5 checksum
	'''

	def md5(fname):
		''' 
		function written by quantumSoup
		https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
		'''
		hash_md5 = hashlib.md5()
		with open(fname, "rb") as f:
			for chunk in iter(lambda: f.read(4096), b""):
				hash_md5.update(chunk)
		return hash_md5.hexdigest()

	if md5(filename) == resource['md5']:
		return True
	return False

def view_package(package, fields=None):
	'''Displays selected information from search result package'''

	id_fields = [
		'title',
		'id',
		'sample_type',
		'data_type',
		'sample_database_file',
		'notes',
		'description',
		'env_medium',
		'information',
	]
	url_fields = [
		'base_url',
		'url'
	]
	geo_fields = [
		'geo_loc_name',
		'latitude',
		'longitude'
	]
	msc_fields = [
		'sample_metadata_ingest_file',
		'number_of_raw_reads'
	]

	if fields is None:
		fields = [
			id_fields,
			geo_fields,
			url_fields,
			msc_fields
		]

	for f in fields:
		print('---------------------------')
		for key in f:
			print(f'{key} : {package.get(key, "-----")}')

def view_resource(resource, fields=None):

	fields = [[
		'id',
		'package_id',
		'size',
		'format',
		'read',
		'name',
		'md5'
	]]

	for f in fields:
		print('---------------------------')
		for key in f:
			print(f'{key} : {resource.get(key, "-----")}')

def searchDB(query, rows=1000):
	'''
	Query microbiom data base based on examples at 
	https://usersupport.bioplatforms.com/programmatic_access.html
	returns 

	'''

	remote = ckanapi.RemoteCKAN(
		'https://data.bioplatforms.com', apikey=_apikey
	)

	result = remote.action.package_search(
		q=query,
		rows = rows,
		include_private = True
	)
	result['Query'] = query
	result['Date'] = time_stamp()
	return result

def downloader(resource, filename, stream=False):
	'''
	Dowloads metagenomes from australian microbiom

	Args:
	resource - search result resource from search result package from AMB
	destination - filename to save as, full path 
	'''
	print(f'Attempting to download {resource["name"]}')

	remote = ckanapi.RemoteCKAN(
		'https://data.bioplatforms.com', apikey=_apikey
	)

	try:
		response = requests.get(
			resource['url'],
			headers={'Authorization': remote.apikey},
			stream=stream
			)
		print(resource['url'])
	except Exception as e:
		print('Error downloading file', e)

	try:
		with open(filename, 'wb') as fh:
			if stream:
				# Code snippet copied from 
				# https://usersupport.bioplatforms.com/programmatic_access.html
				for chunk in response.iter_content(chunk_size=1024):
					if chunk:
						fh.write(chunk)
			else:
				fh.write(response.content)
	except Exception as e:
		print('Error writing to file', e)

	if check_integrity(resource, filename):
		print('File downloaded')
	else:
		print('Download failed file integrity check')
		new_filename = filename.split('.')[0] + '.corrupt'
		os.rename(filename, new_filename)
		print('Resource information on failed download:')
		view_resource(resource)


def main():
	test_search()
	test_load()
	test_download()

def examine_sample_types(result_package):
	'''
	Views sample types present in search results

	Args:
	result_package - A dictionary of search results returned from AMB
	'''
	types = {}
	for r in result_package['results']:
		types[r.get('sample_type', 'No entry')] = types.get(
			r.get('sample_type', 'No entry'), 0) + 1
	for t in types:
		print(t, f': {types[t]}')

def test_load():
	with open('testSearch.json', 'r') as fh:
		results = json.load(fh)
	examine_sample_types(results)
	show_structure(results)

def show_structure(result_package):
	'''
	Displayes structure of result package dictionary

	Args:
	result_package - A dictionary of search results returned from AMB
	'''
	for key in result_package:
		print('Key: ', key)
		if key == 'count':
			print('Metagenomes found: ', result_package[key])
		else:
			print('Sub categories: ', len(result_package[key]))
	print('Search query: ', result_package['Query'])
	print('Search date: ', result_package['Date'])

def test_search():
	''' Test of search terms 
	data_type -- MGE = metagenomes ??
	library_construction_protocol
	omics -- metagenomics
	ph -- 
	sample_type -- Coastal water
	'''
	terms = [
		'data_type:MGE',
		'sample_type:"Coastal water"',
		'omics:metagenomics'
	]

	results = searchDB(terms[2])
	with open('testSearch.json', 'w') as fh:
		json.dump(results, fh)

	examine_sample_types(results)

def time_stamp():
	return time.ctime()

def test_download():
	pass

if __name__ == '__main__':
	main()