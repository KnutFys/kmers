import json

class SearchResult:

	def __init__(self, filename):
		self.filename = filename
		self.load(filename)
		self.fields = [
			'sample_type', 'type', 'am_environment',
			'library_construction_protocol', 'sequence_data_type',
			'geo_loc_name']

	def load(self, filename):
		with open(filename, 'r') as fh:
			result = json.load(fh)
		
		self.query = result['Query']
		self.packages = result['results']
		self.result = result

	def add_field(self, field):
		total = len(self.packages)

		values = {}
		for package in self.packages:
			value = package.get(field, 'No entry')
			if value == '':
				value = 'No entry'
			if values == 'No entry':
				continue
			current_count = values.get(value, 0)
			values[value] = current_count + 1

		print('Field: ', field,)
		for key, item in values.items():
			print(' '*5, key, ':', item)
	
	def view_all_package_fields(self):
		field_types = {
			'string':{}, 'list':{}, 'nested':{}, 'unknown':{}
			}
		for package in self.packages:
			for field_name, field_entry in package.items():
				if isinstance(field_entry, str):
					field_types['string'][field_name] = 'Poop'
				elif isinstance(field_entry, list):
					field_types['list'][field_name] = '...'
				elif isinstance(field_entry, dict):
					field_types['nested'][field_name] = '...'
				else:
					field_types['unknown'][field_name] = '...'
		for field_type in field_types:
			print(field_type)
			print('-'*25)
			for field in field_types[field_type]:
				print(' '*5, field)

						  
	def view_search(self):
		def purtify(word):
			print(' '*5, word)
		print(f'Query: {self.result["Query"]}')
		purtify(f'Date: {self.result["Date"]}')
		purtify(f'Number of metagenomes found: {self.result["count"]}')
		purtify(f'Filename: {self.filename}')

	@staticmethod
	def view_nested(nested, offset, field=None):
		if isinstance(nested, list):
			if field is None:
				print(' '*offset, 'List of length ', len(nested))
			else:
				for item in nested:
					print(' '*offset,f'{field}: ', item.get(
						field,'Unknown '), 'Entries: ', len(item))
		elif isinstance(nested, dict):
			for key, item in nested.items():
				print(' '*offset, key)
				if isinstance(nested, dict):
					SearchResult.view_nested(item, offset+10)
		else:
			print(' '*offset, nested)



			
		

