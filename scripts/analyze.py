import pandas as pd
import numpy as np
from pathlib import Path
import phylotreelib as pt
from scripts.utils import create_symlink, symlink_directories
import matplotlib.pyplot as plt
from tqdm import tqdm
import itertools


def load_pyseer_hits(structure, params):
    pyseer_hits_table = structure['processed']['pyseer_hits_final_table_filtered']
    return pd.read_csv(pyseer_hits_table, sep='\t')


class AnalyzeResults:
	def __init__(self, structure, params):

		self.params = params
		self.structure = structure

		# params
		self.phenotypes = params['gwas']['phenotypes']
		self.collections_colors = params['colors_and_shapes']['collections_colors']
		self.ecod_colors = params['colors_and_shapes']['ecod_colors']
		self.recombinant_depos_colors = params['colors_and_shapes']['recombinant_depos_colors']
		self.min_precision = params['colors_and_shapes']['min_precision']
		self.min_recall = params['colors_and_shapes']['min_recall']
		self.high_precision = params['colors_and_shapes']['high_precision']

		# paths
		self.bacteria_table = structure['input_processed']['bacteria_tsv']
		self.tree = structure['input']['bacteria_tree']
		self.recombinant_depos_table = structure['input_processed']['recombinant_depos_table']
		self.prediction_depos_table = structure['input']['predictions_table']
		self.phenotypes_matrix = structure['gwas']['phenotypes_matrix'] # alwasy the same
		self.pyseer_hits_table = structure['processed']['pyseer_hits_final_table_filtered']
		self.mmseqs_dir = structure['mmseqs']['mmseqs_dir'].parent

		self.per_locus_dir = structure['analyze']['per_locus_dir']
		self.aggregated_data_dir = structure['analyze']['aggregated_data_dir']
		self.predictions_and_enzymes_table = structure['analyze']['predictions_and_enzymes_table']



	def aggregate_data(self, run=True):

		# checkpoint
		if run: pass
		else: return


		# paths
		src_dir = self.structure['input_processed']['processed_dir']
		dst_dir = self.aggregated_data_dir

		# symlink directories with input tables
		symlink_directories(src_dir, dst_dir, retain_structure=True)


		# symlink pyseer files
		src1 = self.structure['processed']['pyseer_hits_final_table_all']
		src2 = self.structure['processed']['pyseer_hits_final_table_filtered']

		for src in [src1, src2]:

			# retain file names
			fname = Path(src).name
			dst = Path(dst_dir, fname)

			create_symlink(src, dst)


	def per_locus_analysis(self, run=True):

		# checkpoint
		if run: pass
		else: return

		# load
		pyseer_hits = pd.read_csv(self.pyseer_hits_table, sep='\t')


		#################
		#### TESTING ####
		#################

		# params
		# testing_loci = ['KL3', 'KL7', 'KL64']
		# printing_testing_loci = ' '.join(testing_loci)
		# tqdm.write(f'Select capsules for testing: {printing_testing_loci}')

		# filt_test = (pyseer_hits['locus'].isin(testing_loci))
		# pyseer_hits = pyseer_hits.loc[filt_test].reset_index(drop=True)

		#################
		#### TESTING ####
		#################


		## per K locus alignments and phandango

		# groupby
		groups = pyseer_hits.groupby(['version', 'mode', 'locus'])

		# progress bar
		print('Alignments and phandango... ')

		# iterate
		for (version, mode, locus), group in groups:
			
			# symlink alignments
			self._link_alignments(version, mode, locus, group)

			# PCs table for phandango
			variants_df = self._get_variants_table(version, mode, locus, group)

			# tree for phandango
			self._get_subtree(version, mode, locus, variants_df)			


		### recall vs precision
		
		# aggregated results per K locus
		print('Precision and recall plots per capsule... ')
		groups = pyseer_hits.groupby(['mode', 'locus'])
		for (mode, locus), group in groups:
			output_pdf = Path(self.per_locus_dir, locus, mode) / f'{locus}_{mode}.pdf'
			print(output_pdf)
			self._plot_recall_precision_with_ci_per_kl(mode_locus_df=group, output_pdf=output_pdf)
			break


	def _get_subtree(self, version, mode, locus, variants_df, max_leaves=300, outgroup='398KBV'):
	
		# paths
		output_dir = Path(self.per_locus_dir) / f'{locus}/{mode}/phandango/{version}'
		subtree_file = Path(output_dir) / '_subtree.nwk'

		# checkpoint
		if Path(subtree_file).exists(): return

		# Remove color columns
		cols = [col for col in variants_df.columns if ':colour' not in col]
		variants_df = variants_df[cols]

		# Calculate the abundance of variants in each genome
		variants_abundance = variants_df.sum(axis=1)

		### Define leaves to retain

		# Genomes with at least one variant
		leaves_with_variant = variants_abundance.loc[variants_abundance > 0].index

		# Sample genomes without variants to keep the total number of leaves within max_leaves
		n_leaves_to_sample = max(0, (max_leaves - len(leaves_with_variant)))
		# if not n_leaves_to_sample:
			# print('Associated variants are in too many genomes. No genomes to downsample.')
		leaves_without_variant_sampled = variants_abundance.loc[variants_abundance == 0].sample(n=n_leaves_to_sample).index

		# Combine leaves to retain
		leaves2retain = set(leaves_with_variant).union(set(leaves_without_variant_sampled))
		leaves2retain = list(leaves2retain) + [outgroup]

		### Create subtree

		# Load the full phylogenetic tree
		with pt.Newicktreefile(self.tree) as treefile:
			subtree = treefile.readtree()

		# Identify leaves to remove
		leaves_all = subtree.leaflist()
		leaves2remove = set(leaves_all).difference(set(leaves2retain))

		# Remove unwanted leaves
		subtree.remove_leaves(leaves2remove)

		# Root the tree
		subtree.rootout(outgroup)

		# Save the subtree
		with open(subtree_file, 'w') as f:
			f.write(subtree.newick())
	

	def _get_variants_table(self, version, mode, locus, group):
		
		# paths
		variants_table = Path(self.mmseqs_dir) / f'{version}/3_binary_matrix.tsv'
		output_dir = Path(self.per_locus_dir) / f'{locus}/{mode}/phandango/{version}'
		phandango_table = Path(output_dir) / '_variants.csv'

		if phandango_table.exists():
			return pd.read_csv(phandango_table, sep='\t')

		# create		
		Path(output_dir).mkdir(exist_ok=True, parents=True)

		# parameters
		map_colours = {'phenotype': '#273746', 'PC': '#5DADE2', 'other': '#FDFEFE'}

		# load files
		# pyseer_hits = pd.read_csv(self.pyseer_hits_table, sep='\t')
		phenotypes_matrix_df = pd.read_csv(self.phenotypes_matrix, sep='\t')
		variants_df = pd.read_csv(variants_table, sep='\t')

		# convert phenotype matrix
		new_cols_dict = {col: col.split('_0')[0] for col in list(phenotypes_matrix_df.columns) if '_0' in col}
		phenotypes_matrix_df = phenotypes_matrix_df.rename(new_cols_dict, axis=1)

		# convert variants matrix: set row names from first column, drop first column, transpose, etc.
		variants_df.index = variants_df.iloc[:, 0]
		variants_df = variants_df.drop(variants_df.columns[0], axis=1)
		variants_df.index.name = ''
		variants_df = variants_df.T.reset_index().rename({'index': 'genomeID'}, axis=1)

		def get_sorted_variants(df):
			# sort and extract variants (PCs)
			order_reported_topology = self.ecod_colors.keys()
			df['reported_topology'] = pd.Categorical(df['reported_topology'], categories=order_reported_topology, ordered=True)
			df_sorted = df.sort_values(by=['reported_topology', 'precision'], ascending=[True, True])
			return list(df_sorted['PC'].unique())

		
		# select variants
		variants = get_sorted_variants(group)
		variants_cols = ['genomeID'] + variants
		variants_df = variants_df[variants_cols]

		# Selected phenotype column
		phenotype_cols = ['genomeID', locus]
		phenotype_df = phenotypes_matrix_df[phenotype_cols]

		# Merge phenotype and variant data
		df = phenotype_df.merge(variants_df, on='genomeID', how='left')

		# Assign phenotype colors
		phenotype_map_colours = {1: map_colours['phenotype'], 0: map_colours['other']}
		df[f'{locus}:colour'] = df[locus].map(phenotype_map_colours)

		# Create a mapping from PC to its reported_topology (assuming each PC has one unique topology in group)
		pc_topology_mapping = group.drop_duplicates(subset=['PC']).set_index('PC')['reported_topology'].to_dict()

		# assign a colors to PCs
		for col in variants_cols[1:]:
			topology = pc_topology_mapping.get(col, None)
			variant_color = self.ecod_colors.get(topology, map_colours['other'])
			df[f'{col}:colour'] = df[col].map(lambda val: variant_color if val == 1 else map_colours['other'])

		# compatible with subtree: set the index to genomeID and remove that column
		df.index = df['genomeID']
		df.index.name = 'genomeID'
		df = df.drop('genomeID', axis=1)

		# save the full phandango variants table
		df.to_csv(phandango_table, sep=',')

		# split the table into smaller fragments (e.g., for Phandango)
		def split_list(lst, chunk_size):
			"""Split a list into chunks of given size."""
			return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

		table_fragments_cols = split_list(variants, chunk_size=3)
		for i, frag_variants in enumerate(table_fragments_cols):
			outdir, name_extension = phandango_table.parent, phandango_table.name.split('.')
			name, extension = '.'.join(name_extension[:-1]), name_extension[-1]
			outfile = Path(outdir, f'FRAGMENT_{i}.{extension}')
			frag_variants_colour = [f'{col}:colour' for col in frag_variants]
			cols = [locus] + frag_variants + [f'{locus}:colour'] + frag_variants_colour
			df[cols].to_csv(outfile, sep='\t')

		return df

			
	def _link_alignments(self, version, mode, locus, group):

		# paths
		src_dir = Path(self.mmseqs_dir) / f'{version}/alignments'
		dst_dir = Path(self.per_locus_dir) / f'{locus}/{mode}/alignments/{version}'

		# create
		Path(dst_dir).mkdir(exist_ok=True, parents=True)

		# symlink
		pcs = list(group['PC'].unique())
		for pc in pcs:
			src = Path(src_dir) / f'{pc}.fasta'
			dst = Path(dst_dir) / f'{pc}.fasta'

			# skip existing
			if dst.exists(): return

			create_symlink(src, dst)


	def _plot_recall_precision_with_ci_per_kl(self, mode_locus_df, output_pdf):
		"""
		Create a PDF with 12 subplots (4 rows x 3 columns).
		"""

		import pandas as pd
		import matplotlib.pyplot as plt
		from matplotlib.backends.backend_pdf import PdfPages
		import numpy as np
		from adjustText import adjust_text

		# params
		mode = mode_locus_df['mode'].unique()[0]
		locus = mode_locus_df['locus'].unique()[0]

		# Initialize figure with 4 rows x 3 columns
		fig, axes = plt.subplots(4, 3, figsize=(10, 12), sharex=True, sharey=True)

		# Assumed identity and coverage values
		identity_list = ['00', '50', '80']
		coverage_list = ['50', '80']

		configs = [
			{"use_sc": False, "row_offset": 0, "label": "Isolates"},
			{"use_sc": True, "row_offset": 2, "label": "SCs"}
		]

		# Loop through the configurations
		for config in configs:
			use_sc = config["use_sc"]
			row_offset = config["row_offset"]

			if use_sc:
				precision_col = 'precision_SC'
				recall_col = 'recall_SC'
				precision_ci_low = 'precision_SC_CI5_LOWER'
				precision_ci_up = 'precision_SC_CI95_UPPER'
				recall_ci_low = 'recall_SC_CI5_LOWER'
				recall_ci_up = 'recall_SC_CI95_UPPER'
			else:
				precision_col = 'precision'
				recall_col = 'recall'
				precision_ci_low = 'precision_CI5_LOWER'
				precision_ci_up = 'precision_CI95_UPPER'
				recall_ci_low = 'recall_CI5_LOWER'
				recall_ci_up = 'recall_CI95_UPPER'

			for i, coverage in enumerate(coverage_list):
				for j, identity in enumerate(identity_list):
					filt = mode_locus_df['version'].str.contains(f'PCI{identity}') & mode_locus_df['version'].str.contains(f'C{coverage}$')
					subdf = mode_locus_df[filt]
					ax = axes[i + row_offset, j]
					texts = []

					for _, row in subdf.iterrows():
						rep_top = row['reported_topology']
						color = self.ecod_colors.get(rep_top, 'black')

						jitter_x = np.random.uniform(-0.02, 0.02)
						jitter_y = np.random.uniform(-0.02, 0.02)

						x = row[precision_col] - jitter_x
						y = row[recall_col] - jitter_y

						jittered_precision_ci_low = row[precision_ci_low] - jitter_x
						jittered_precision_ci_up = row[precision_ci_up] - jitter_x
						jittered_recall_ci_low = row[recall_ci_low] - jitter_y
						jittered_recall_ci_up = row[recall_ci_up] - jitter_y

						ax.errorbar(
							x, y,
							xerr=[[np.abs(x - jittered_precision_ci_low)], [np.abs(jittered_precision_ci_up - x)]],
							yerr=[[np.abs(y - jittered_recall_ci_low)], [np.abs(jittered_recall_ci_up - y)]],
							fmt='o',
							mfc=color,
							mec='black',
							markeredgewidth=1,
							markersize=6,
							ecolor=color,
							capsize=2,
							alpha=0.75
						)

						if mode == 'elastic_net':
							if row[precision_ci_low] >= self.min_precision or row[recall_ci_low] >= self.min_recall:
								txt = ax.text(x, y, row['PC'], fontsize=6, color='#FF7074')
								texts.append(txt)
						else: 
							txt = ax.text(x, y, row['PC'], fontsize=6, color='#FF7074')
							texts.append(txt)

					ax.axvline(x=self.min_precision, color='red', linestyle='--', linewidth=1, alpha=0.2)
					ax.axvline(x=self.high_precision, color='gray', linestyle='--', linewidth=1, alpha=0.2)
					ax.axhline(y=self.min_recall, color='red', linestyle='--', linewidth=1, alpha=0.2)

					ax.set_xlim(-0.05, 1.05)
					ax.set_ylim(-0.05, 1.05)

					if i + row_offset == 3:
						ax.set_xlabel('Precision')
					else:
						ax.set_xlabel('')

					if j == 0:
						ax.set_ylabel('Recall', labelpad=20)
					else:
						ax.set_ylabel('')

					adjust_text(texts, ax=ax, expand=(1.2, 1.5), arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

			row_center_y = (axes[row_offset, 0].get_position().y0 + axes[row_offset + 1, 0].get_position().y1) / 2
			fig.text(0.05, row_center_y, config["label"], ha='center', va='center', fontsize=12, fontweight='bold', rotation=90)

		# Align identity labels at the top of each column
		for j, identity in enumerate(identity_list):
			col_center_x = axes[0, j].get_position().x0 + axes[0, j].get_position().width / 2
			fig.text(col_center_x, 0.9, f"Identity: {identity}", ha='center', va='bottom', fontsize=12)

		# Align coverage labels to the right of each row
		for config in configs:
			row_offset = config["row_offset"]
			for i, coverage in enumerate(coverage_list):
				row_index = i + row_offset
				row_center_y = axes[row_index, -1].get_position().y0 + axes[row_index, -1].get_position().height / 2
				fig.text(0.93, row_center_y, f"Coverage: {coverage}", ha='left', va='center', fontsize=12, rotation=90)

		# Adjust layout and margins
		fig.subplots_adjust(left=0.15, right=0.9, top=0.88, bottom=0.12, hspace=0.3, wspace=0.2)

		fig.suptitle(f"{locus} {mode}", fontsize=16)

		with PdfPages(output_pdf) as pdf:
			pdf.savefig(fig)

		plt.close(fig)


	def combine_predictions_and_enzymes(self, run=True):

		# checkpoint
		if run: pass
		else: return

		# params
		recombinant_cols = ['proteinID', 'group', 'expression', 'specificity', 'K_locus_host', 'seq']
		prediction_cols = ['locus', 'PC', 'prediction_strength', 'seq']

		# tables
		bacteria_table = self.bacteria_table
		recombinant_depos_table = self.recombinant_depos_table
		prediction_depos_table = self.prediction_depos_table

		# load
		bacteria_df = pd.read_csv(bacteria_table, sep='\t')
		recombinant_df = pd.read_csv(recombinant_depos_table, sep='\t')
		prediction_df = pd.read_excel(prediction_depos_table, sheet_name='BEST_DEPO_HITS')

		# clean recombinant
		recombinant_df = recombinant_df[recombinant_cols].rename({'group': 'source'}, axis=1)
		recombinant_df['expression'] = recombinant_df['expression'].fillna('PRODUCED')
		recombinant_df.loc[recombinant_df['source'].str.contains('zdk'), 'source'] = 'PROPHAGE_ZDKLAB'
		recombinant_df.loc[recombinant_df['source'].str.contains('genscript'), 'source'] = 'GENSCRIPT'
		recombinant_df.loc[recombinant_df['source'].str.contains('lytic'), 'source'] = 'LITERATURE_SEARCH'

		recombinant_df['activity'] = 'ACTIVE'
		recombinant_df.loc[recombinant_df['specificity'].str.contains('NO_KL'), 'activity'] = 'INACTIVE'
		recombinant_df.loc[recombinant_df['expression'].str.contains('NOT_PRODUCED'), 'activity'] = np.nan
		recombinant_df.loc[recombinant_df['expression'].str.contains('NOT_PRODUCED'), 'specificity'] = np.nan
		recombinant_df.loc[recombinant_df['activity'] == 'INACTIVE', 'specificity'] = np.nan


		recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'] = \
		recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'].str.replace('K', 'KL')

		recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'] = \
		recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'specificity'].str.replace('KLN', 'KN')

		recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'proteinID'] = \
		recombinant_df.loc[recombinant_df['source'] == 'LITERATURE_SEARCH', 'proteinID'].str.rsplit('_', n=1).str[0]


		# clean prediction
		prediction_df = prediction_df[prediction_cols].rename({'locus': 'specificity', 'PC': 'proteinID'}, axis=1)
		prediction_df = prediction_df.dropna()
		prediction_df['source'] = 'PREDICTION'
		prediction_df['prediction_strength'] = prediction_df['prediction_strength'].replace({'good': 'yes', 'likely': 'maybe (+)'})


		# concat
		depos_df = pd.concat([recombinant_df, prediction_df])
		depos_df['assing_K_locus_host_when_specificity_missing'] = depos_df['specificity'].fillna(depos_df['K_locus_host'])

		# strains per K locus
		kaspah_ref_df = bacteria_df.query('collection == "kaspah_complete"').groupby("MGG_K_locus").size().reset_index().rename({'MGG_K_locus': 'specificity', 0: '# KASPAH-REF'}, axis=1)
		ksc_gwas_df = bacteria_df.groupby("MGG_K_locus").size().reset_index().rename({'MGG_K_locus': 'specificity', 0: '# GWAS_KSC'}, axis=1)
		strains_availability_per_k_locus = kaspah_ref_df.merge(ksc_gwas_df, on='specificity', how='outer')
		strains_availability_per_k_locus['assing_K_locus_host_when_specificity_missing'] = strains_availability_per_k_locus['specificity']
		strains_availability_per_k_locus = strains_availability_per_k_locus.drop('specificity', axis=1)

		# merge
		depos_df = depos_df.merge(strains_availability_per_k_locus, on='assing_K_locus_host_when_specificity_missing', how='left')

		# save
		final_cols = ['proteinID', 'source','expression', 'assing_K_locus_host_when_specificity_missing', 'specificity', 'K_locus_host','# KASPAH-REF','# GWAS_KSC','activity','prediction_strength','seq']
		depos_df[final_cols].to_csv(self.predictions_and_enzymes_table, sep='\t', index=False)









