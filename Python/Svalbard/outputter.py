# Standard library imports

# Third party imports
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
# Local imports
import kmer_DB as DB

def excel_writer(sheets, out_name):
	'''Writes to spreadsheet xlsx file


	Arguments:
	sheets      -- (dic) Dictionary with sheetnames as keys and 2 
				data frames pr key
	out_name    -- (str) Filename to write to. Will overwrite
	'''
	with pd.ExcelWriter(out_name) as writer:
		for key in sheets.keys():
			startrow = 1
			startcol = 1
			for d in sheets[key]:
				d[0].to_excel(
					writer, startrow=startrow, startcol=startcol,
					index=False, sheet_name=key)
				d[1].to_excel(
					writer, startrow=(len(d[0])+ 2), startcol=startcol-1,
					index=True, sheet_name=key)
				startcol += 8

def multiplotter(x_dicts, hue, filename=None, mode='box',rows=1, columns=2, title=''):
	''' Creates multiplots of sample sizes vs distance measure

	'''
	figure, ax = plt.subplots(
		rows, columns, figsize=(8,10), sharey='row', sharex='col')
	figure.suptitle(title, fontsize=16)
	row = 0
	col = 1
	for i, x_dict in enumerate(x_dicts):
		current_axis = ax[row, col]
		#if mode == 'line':
		#	current_axis.set_xscale('log')
		plotter(x_dict[0], hue, mode=mode, return_ax=current_axis)
		handles, labels = current_axis.get_legend_handles_labels()
		current_axis.legend([],[])
		current_axis.set_title(x_dict[1])
		if col == 0:
			col += 1
		else:
			current_axis.tick_params(left=False, labelleft=False)
			current_axis.tick_params(right=True, labelright=True)
			col = 0
			row += 1

	current_axis = ax[0, 0]
	labels = [
		'Angular distance to SV001',
		'Jensen-Shannon distance to SV001',
		'Bray-Curtis dissimilarity to SV001',
		'Angular distance to SV002',
		'Jensen-Shannon distance to SV002',
		'Bray-Curtis dissimilarity to SV002'
		]
	current_axis.get_yaxis().set_visible(False)
	current_axis.get_xaxis().set_visible(False)
	for key in current_axis.spines.keys():
		current_axis.spines[key].set_visible(False)
	current_axis.legend(handles, labels, loc='center',
		bbox_to_anchor=(0.45, 0.45, 0.1, 0.1), frameon=False)
	if filename is None:
		plt.show()
	else:
		filename = filename.replace('.png','') + f'{mode}.png'
		plt.savefig(filename)

def plotter(x_dict, hue, filename=None, mode='box', return_ax=None):
	'''Plots samples size vs distance measure

	'''
	if return_ax is None:
		figure, ax = plt.subplots(figsize=(7,5))
	else:
		ax = return_ax
	if not isinstance(x_dict, list):
		x_dict = [x_dict]
	if not isinstance(hue, list):
		hue = [hue]
	data = pd.melt(x_dict[0], var_name='Sample size', value_name='Distance',)
	data['measure'] = hue[0]
	for x, h in zip(x_dict[1:], hue[1:]):
		d = pd.melt(x, var_name='Sample size', value_name='Distance',)
		d['measure'] = h
		data = data.append(d)
	if mode == 'box':
		sns.boxplot(
			x='Sample size', y='Distance', data=data, hue='measure', ax=ax)
	if mode == 'line':
		sns.lineplot(
			x='Sample size', y='Distance', ci="sd", data=data,ax=ax, hue='measure')
	if return_ax is not None:
		return ax
	else:
		if filename is None:
			plt.show()
		else:
			filename = filename.replace('.png','') + f'{mode}.png'
			plt.savefig(filename)