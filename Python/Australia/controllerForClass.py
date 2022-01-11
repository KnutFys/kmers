'''
Main script for cleaning fasta files and counting kmers

'''
# Standard library imports
import os
import sys
import argparse
# Local imports
import prjclass # Organizes project and exports csv files to R

def link_packages_to_files(packages):
	links = {}
	for package in packages:
		package_id = package.get('id', 'Unknown')
		for resource in package['resources']:
			if not resource['name'].endswith('.fastq.gz'):
				continue
			file_names = get_file_names(resource=resource)
			links[file_names['count']] = package
	return links

def find_paired(resource, package):
	'''
	Searches package for other read. IE for R1 read looks for R2 read

	Args:
	resource - (dic) Resource to be paired
	package - (dic) Package the resource belongs to
	'''
	names = get_file_names(resource=resource)
	paired_resource = None
	for res in package['resources']:
		if res['name'] == names['paired']:
			return res
	
def download_old(args, from_process=False):
	'''
	Downloads metagenomes in search results

	Args:
	args - (namespace) Parsed command line options
	from_process - (bol) Switch determining if the function is being called
		from the 'process' function
	'''
	# Set mandatory fields. Samples that do not have valid entries
	# in all specified fields will be ignored.
	mandatory_fields = ['sample_type']
	# TEMP Need better way to selec search ################
	project = load_project(args.name)
	
	packages = result['results']

	# Download up to 'args.batchsize' files
	downloaded = 0
	for package_id, package in project.packages.items():
		# Ignore packages lacking specified field entries
		try:
			for field in mandatory_fields:
				field_value = package.get(field, None)
				if field_value is None or field_value == '':
					raise ValueError(f'Lacking entry for field {field}')

		except ValueError as e:
			print(e)
			continue

		for resource in package['resources']:
			# Ignore checksum files
			if not resource['name'].endswith('.fastq.gz'):
				#print('Skipping checksum file...')
				continue
			# Find paired resource
			paired_resource = find_paired(resource, package)
			# Ignore counted, trimmed or downloaded files
			if not check_status(project, resource)['downloaded']:
				# Download file
				download_path = project['directories']['download_dir']
				download_file = resource['name']
				dbSearch.downloader(
					resource, os.path.join(download_path, download_file))
				downloaded += 1
			if not check_status(project, paired_resource)['downloaded']:
				# Download file
				download_path = project['directories']['download_dir']
				download_file = paired_resource['name']
				dbSearch.downloader(
					paired_resource, os.path.join(download_path, download_file))
				downloaded += 1
			# Check if anything has been downloaded
			if downloaded == 0:
				continue
			# Check if called from process function
			if from_process:
				return downloaded
			if downloaded >= args.batchsize:
				break
		if downloaded >= args.batchsize:
			break

def download(args, from_process=False):
	# Load project info
	project = prjclass.Project(args.name)
	if args.forcefield is None:
		project.download(batchsize=args.batchsize)
	else:
		fields = args.forcefield.split(' AND ')
		forced_fields = {f.split('=')[0]:f.split('=')[1] for f in fields}
		project.download(batchsize=args.batchsize, forced_fields=forced_fields)

def trim(args, from_process=False):
	project = prjclass.Project(args.name)
	project.trim(batchsize=args.batchsize)

def count(args, from_process=False):
	
	# Load project info
	project = prjclass.Project(args.name)
	if from_process:
		if project.count():
			return 1
		else:
			return 0
	counted_pairs = 0
	for i in range(args.batchsize):
		if not project.count():
			print('No more uncounted files detected...')
			break
		counted_pairs += 1
	print(f'Counted kmers for {counted_pairs} of {args.batchsize} pairs')
	
def clean_up(args):
	''' Remove processed files '''
	project = prjclass.Project(args.name)

	# Remove trimmed and/or  counted files from
	# download folder
	if args.clean_downloads or args.clean_all:
		project.clean_up_downloads()

	# Remove counted files from trimmed folder
	if args.clean_trimmed or args.clean_all:
		project.clean_up_trimmed()

def process(args):
	batchsize = args.batchsize
	downloaded = 0
	trimmed = 0
	counted = 0
	project = prjclass.Project(args.name)
	for i in range(batchsize):
		if args.forcefield is None:
			downloaded += project.download(
				batchsize=batchsize, from_process=True)
		else:
			fields = args.forcefield.split(' AND ')
			forced_fields = {f.split('=')[0]:f.split('=')[1] for f in fields}
			downloaded += project.download(
				batchsize=batchsize, from_process=True, forced_fields=forced_fields)			
		trimmed += project.trim(
			batchsize=batchsize, from_process=True)
		counted += project.count(from_process=True)
		print('Downloaded file: ', downloaded)
		print('Trimmed files: ', trimmed)
		print('Counted files: ', counted)
		project.clean_up_downloads()
		project.clean_up_trimmed()
	print('Done')

def calculate(args):
	'''
	Calculates distances between counted metagenomes and stores
	result in project database
	'''
	workers = 1
	if args.multi:
		workers = 3

	project = prjclass.Project(args.name)
	if args.update:
		project.update_db()

	project.calculate(batchsize=args.batchsize, workers=workers)

def search(args):
	''' Search for metagenomes and populate database '''

	project = prjclass.Project(args.name)
	project.search(args.query)

def view_project(args):
	'''
	Display information about selected project
	'''
	project = prjclass.Project(args.name)
	if args.detailed_searches:
		project.view_searches()
	elif args.detailed_fields:
		project.create_sample_metadata_csv()
	else:
		project.view_project()
		project.check_file_availability()

def new_project(args):
	project = prjclass.Project(args.name)
	print(f'Project {args.name} initiated...')

def export(args):
	fields = None
	samples = None
	if args.out_file is None:
		filename = args.name
	else:
		filename = args.out_file

	project = prjclass.Project(args.name)
	project.add_field('sample_type')
	project.add_field('geo_loc_name')
	project.add_field('longitude')
	project.add_field('latitude')
	project.add_field('collection_date')
	project.add_field('collection_date')
	project.to_csv(filename, fields=None)

def plot(args):
	if args.out_file is None:
		filename = args.name
	else:
		filename = args.out_file
	project = prjclass.Project(args.name)
	project.plot_R(filename, args.col, args.shape, args.measure)

def set_api(args):
	project = prjclass.Project(args.name)
	project.set_api(args.api)
	print(f'Api for project {args.name} has been set')

def fields(args):
	''' Adds or removes mandatory fields from project. Samples without
	valid entry for any mandatory field will be ignored'''

	project = prjclass.Project(args.name)
	if args.add:
		for f in args.fields:
			project.add_field(f)
			print(f'Field "{f}" for project "{args.name}" has been set to mandatory')
	if not args.add:
		for f in args.fields:
			print(f)
			project.remove_field(f)
			print(f'Field "{f}" for project "{args.name}" is no longer mandatory')		
	project.save_project()

def update_db(args):
	project = prjclass.Project(args.name)
	project.update_db()

def pole_distances(args):
	''' Calculates distances from each sample to idealized pole sample and
	stores results in a json file along with meta data from project fields.
	'''
	project = prjclass.Project(args.name)
	project.pole_distance(batchsize=args.batchsize, export=args.export)

def test_function(args):
	'''
	Purely for developmental purposes.
	'''
	project = prjclass.Project(args.name)
	project.pole_distance(batchsize=args.batchsize)

def main():
	help_text = {
		'name':'''Project name. Example: <myproject> . Related files
		will be named accordingly. I.E. the database file will
		be named <myproject.db>, project settings <myproject.prj>''',
		'rows':'Maximum number of rows returned in search',
		'search':'''Search AMB for metagenomes''',
		'query':'''Search term formatted as specified by the AMB database. Syntax
		searchfield:value
		If query contains white spaces, enclose entire query in quotation marks.
		Logical operators can be used:
		"sample_type:marine AND sample_type:soil"
		''',
		'update':'''Check if new files have been downloaded, if so add them
		to the data base'''
	}
	#################################
	# Main parser
	#################################
	parser = argparse.ArgumentParser(
	description='''
	Downloads metagnomes and calculates the similarity using kmer counts
	'''
	)
	parser.add_argument(
		'-t', '--test', help='Test run', action='store_true', default=False
	)
	parser.add_argument(
		'name', help=help_text['name']
	)
	subparsers = parser.add_subparsers(help='Choose command')
	#################################
	# Subparser for 'search' command
	#################################
	parser_search = subparsers.add_parser(
		'search', help=help_text['search']
	)
	parser_search.add_argument(
		'query', type=str, help=help_text['query']
	)
	parser_search.add_argument(
		'-r', '--rows', type=int, help='Maximum number of rows returned in search',
		dest='rows'
	)
	parser_search.set_defaults(command=search, rows=1000)
	###################################
	# Subparser for 'process' command
	###################################
	parser_process = subparsers.add_parser(
		'process', help='Process metagenome. Downloads, trims and counts kmers'
	)
	parser_process.add_argument(
		'-n', type=int, help='Number of metagenomes to process', dest='batchsize'
	)
	parser_process.add_argument(
		'-ff', '--forcefield',
		help='Only download samples with given field value', dest='forcefield')
	parser_process.set_defaults(command=process, batchsize=1, forcefield=None)
	###################################
	# Subparser for 'download' command
	###################################
	parser_download = subparsers.add_parser(
		'download', help='Download metagenome'
	)
	parser_download.add_argument(
		'-n', type=int, help='Number of metagenomes to download',
		dest='batchsize'
	)
	parser_download.add_argument(
		'-ff', '--forcefield',
		help='Only download samples with given field value', dest='forcefield')
	parser_download.set_defaults(command=download, batchsize=1, forcefield=None)
	###################################
	# Subparser for 'trim' command
	###################################
	parser_trim = subparsers.add_parser(
		'trim', help='Trim downloaded metagenomes'
	)
	parser_trim.add_argument(
		'-n', type=int, help='Number of metagenomes to process', dest='batchsize'
	)
	parser_trim.set_defaults(command=trim, batchsize=1)
	###################################
	# Subparser for 'count' command
	###################################
	parser_count = subparsers.add_parser(
		'count', help='Count kmers in trimmed metagenomes'
	)
	parser_count.add_argument(
		'-n', type=int, help='Number of metagenomes to process', dest='batchsize'
	)
	parser_count.add_argument(
		'-c', '--clean', help='Remove download and trimmed files', dest='clean',
		action='store_true'
	)
	parser_count.set_defaults(command=count, batchsize=1, clean=False)
	###################################
	# Subparser for 'clean' command
	###################################
	parser_clean = subparsers.add_parser(
		'clean', help='Remove processed metagenomes from drive'
	)
	parser_clean.add_argument('-d', '--downloads', action='store_true',
		dest='clean_downloads', help='Only remove download files')
	parser_clean.add_argument('-t', '--trimmed', action='store_true',
		dest='clean_trimmed', help='Only remove trimmed files')
	parser_clean.add_argument('-a', '--all', action='store_true',
		dest='clean_all', 
		help='Remove download files and trimmed files')
	parser_clean.set_defaults(command=clean_up, clean_downloads=False,
		clean_trimmed=False, clean_all=False)
	###################################
	# Subparser for 'new' command
	###################################
	parser_new = subparsers.add_parser(
		'new', help='Initiates a new project'
	)
	parser_new.set_defaults(command=new_project)
	###################################
	# Subparser for 'pole' command
	###################################
	parser_pole = subparsers.add_parser(
		'pole', help='Calculates distances to idealized pole sample'
	)
	parser_pole.add_argument('-n', '--batchsize', type=int,
		dest='batchsize', help='Maximum number of calculations')
	parser_pole.add_argument('-e', '--export', action='store_true', dest='export',
		help='Export pole distances as csv file')
	parser_pole.set_defaults(command=pole_distances, batchsize=None, export=False)
	###################################
	# Subparser for 'calculate' command
	###################################
	parser_calculate = subparsers.add_parser(
		'calculate', help='Calculate distances between samples'
	)
	parser_calculate.add_argument('-n', '--batchsize', type=int,
		dest='batchsize', help='Maximum number of calculations')
	parser_calculate.add_argument('-m', '--multi', 
		dest='multi', action='store_true', help='Use multiple cores')
	parser_calculate.add_argument('-u', '--update', dest='update',
		action='store_true', help=help_text['update'])
	parser_calculate.set_defaults(command=calculate, batchsize=1)
	###################################
	# Subparser for 'view' command
	###################################
	parser_view = subparsers.add_parser(
		'view', help='View project details'
	)
	parser_view.add_argument(
		'-s','--searches', help='View detailed search information',
		action='store_true', dest='detailed_searches'
	)
	parser_view.add_argument(
		'-f','--fields', help='View detailed field information',
		action='store_true', dest='detailed_fields'
	)
	parser_view.set_defaults(
		command=view_project, detailed_searches=False, detailed_fields=False
		)
	###################################
	# Subparser for 'field' command
	###################################
	parser_field = subparsers.add_parser(
		'fields', help='Manage mandatory fields'
	)
	parser_field.add_argument(
		'-a', '--add', help='Add mandatory fields', action='store_true',
		dest='add'
	)
	parser_field.add_argument(
		'-r', '--remove', help='Remove mandatory fields', action='store_false',
		dest='add'
	)
	parser_field.add_argument(
		'fields', help='Mandatory fields to be added or removed',
		nargs='+'
	)
	parser_field.set_defaults(command=fields, add=True, fields=[])
	###################################
	# Subparser for 'api' command
	###################################
	parser_api = subparsers.add_parser(
		'setapi', help='Set api to access AMB database'
	)
	parser_api.add_argument(
		'api', help='Personal api', nargs='?'
	)
	parser_api.set_defaults(command=set_api, api=None)
	###################################
	# Subparser for 'updatedb' command
	###################################
	parser_updatedb = subparsers.add_parser(
		'updatedb', help='Look for new files and update database'
	)
	parser_updatedb.set_defaults(command=update_db)
	###################################
	# Subparser for 'updatedb' command
	###################################
	parser_dev = subparsers.add_parser(
		'dev', help='Test function. Runs function currently being developed.'
	)
	parser_dev.set_defaults(command=test_function)
	###################################
	# Subparser for 'export' command
	###################################
	parser_export = subparsers.add_parser(
		'export', help='Create csv files of meta data and distance matrices.'
	)
	parser_export.add_argument(
		'-f', '--fields', dest='fields',
		help='''Mandatory metadata fields. Samples lacking data are excluded'''
	)
	parser_export.add_argument(
		'-o', '--out', dest='out_file',
		help='''Basename for output files. By default project name is used.'''
	)
	parser_export.set_defaults(command=export, out_file=None)
	###################################
	# Subparser for 'plot' command
	###################################
	parser_plot = subparsers.add_parser(
		'plot', help='Create R plot.'
	)
	parser_plot.add_argument(
		'-m', '--measure', dest='measure',
		choices=['jensenshannon', 'braycurtis', 'angular'],
		help='''Distance measure to use for plot'''
	)
	parser_plot.add_argument(
		'-c', '--colour', dest='col',
		help='''Field to use as basis for colour scheme of points in the plot.'''
	)
	parser_plot.add_argument(
		'-s', '--shape', dest='shape',
		help='''Field to use as basis for shape of points in the plot.'''
	)
	parser_plot.add_argument(
		'-o', '--out', dest='out_file',
		help='''Basename for output files. By default project name is used.'''
	)
	parser_plot.set_defaults(
		command=plot, measure='jensenshannon',
		col='geo_loc_name', shape='sample_type',
		out_file=None
	)
	##########################################
	args = parser.parse_args()

	# Set path to project file
	args.path = os.path.join(
		os.getcwd(), 'Projects', args.name, args.name + '.prj')

	args.command(args)

if __name__ == '__main__':
	main()

