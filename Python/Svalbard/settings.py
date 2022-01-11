# Standard imports
import os

#################################
### Settings
#################################
_settings = {}

def set_defaults():
	global _settings
	# Kmer settings
	# k - kmer length
	# offset - position to start kmer count in a sequence
	# cutoff - Ingore kmer counts below cutoff
	k = 11
	offset = 0
	cutoff = 0
	sample_size = 100000
	# File paths
	meta_path = '/media/knut/Storage/Metagenomes'
	meta_file = 'sampling_loc_nutrient_conc_knut.csv'
	sample_path = '/media/knut/Storage/Metagenomes/Trim/multisamples'
	count_path = '/media/knut/Storage/Metagenomes/Trim/multisamples/counts'
	raw_data = '/media/knut/Storage/Metagenomes'
	trimmed_data = '/media/knut/Storage/Metagenomes/Trim'
	kraken_data = '/media/knut/Storage/Kraken'

	_settings['meta_data'] = os.path.join(meta_path, meta_file)
	_settings['sample_path'] = sample_path
	_settings['count_path'] = count_path
	_settings['raw_data'] = raw_data
	_settings['trimmed_data'] = trimmed_data
	_settings['kraken_data'] = kraken_data
	_settings['k'] = k
	_settings['offset'] = offset
	_settings['cutoff'] = cutoff
	_settings['sample_size'] = sample_size

def get_setting(setting):
	global _settings
	return(_settings[setting])

def set_setting(setting, value):
	global _settings
	_settings[setting] = value

if __name__ == '__main__':
	pass







