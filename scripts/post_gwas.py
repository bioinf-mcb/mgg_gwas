import pandas as pd
from itertools import combinations

class PostGWAS:
	def __init__(self, structure, params):

		# params
		self.phenotypes = params['gwas']['phenotypes']
		self.collections_colors = params['figures']['collections_colors']
		self.ecod_colors = params['figures']['ecod_colors']
		self.recombinant_depos_colors = params['figures']['recombinant_depos_colors']
		self.min_precision_q25th_circle = params['figures']['min_precision_q25th_circle']
		self.min_precision_q25th_square = params['figures']['min_precision_q25th_square']


		# paths
		self.figures_basic_dir = structure['figures']['figures_basic_dir']


	def _load_results(self):
		
		# params
		bootstrap_version = 'BTSP20'

		# paths
		lasso_tables = list(self.figures_basic_dir.rglob(f'*PCI*C*{bootstrap_version}*/*lasso.tsv'))
		elastic_net_tables = list(self.figures_basic_dir.rglob(f'*PCI*C*{bootstrap_version}*/*elastic_net.tsv'))

		# load
		def _load_and_append_version(table):
			if 'lasso' in str(table) : i_version = -4
			elif 'elastic_net' in str(table) : i_version = -5
			else: 
				print('error')
				exit()

			df = pd.read_csv(table, sep='\t')
			version = str(table).split('_')[i_version]
			df['version'] = version
			return df

		lasso_dfs, elastic_net_dfs = [], []
		for lasso_table, elastic_net_table in zip(lasso_tables, elastic_net_tables):
			lasso_df = _load_and_append_version(lasso_table)
			elastic_net_df = _load_and_append_version(elastic_net_table)

			lasso_dfs.append(lasso_df)
			elastic_net_dfs.append(elastic_net_df)

		# concat
		lasso_df = pd.concat(lasso_dfs)
		elastic_net_df = pd.concat(elastic_net_dfs)

		# add pc version
		lasso_df['PC_version'] = lasso_df['PC'] + '_' + lasso_df['version']
		elastic_net_df['PC_version'] = elastic_net_df['PC'] + '_' + elastic_net_df['version']

		lasso_df['method'] = 'lasso'
		elastic_net_df['method'] = 'elastic_net'
		
		# att
		self.gwas_df = pd.concat([lasso_df, elastic_net_df])
		
		# # split
		# main_predictions_df = all_predictions_df.query('version == @main_version').reset_index(drop=True)
		# other_predictions_df = all_predictions_df.query('version != @main_version').reset_index(drop=True)

		return self.gwas_df

		
	def calculate_pc_similarity(self, data):

		# Prepare sets of members for each PC_version
		pc_sets = {row['PC_version']: set(row['members'].split(',')) for _, row in data.iterrows()}

		results = []

		# Compute pairwise Jaccard and overlap coefficient
		for (pc1, members1), (pc2, members2) in combinations(pc_sets.items(), 2):

			intersection = len(members1 & members2)
			union = len(members1 | members2)
			jaccard_index = round(intersection / union if union > 0 else 0, 2)

			overlap_coefficient_1 = round(intersection / len(members1) if len(members1) > 0 else 0, 2)
			overlap_coefficient_2 = round(intersection / len(members2) if len(members2) > 0 else 0, 2)

			results.append({
			'PC1': pc1,
			'PC2': pc2,
			'jaccard_similarity': jaccard_index,
			'overlap_coefficient_PC1': overlap_coefficient_1,
			'overlap_coefficient_PC2': overlap_coefficient_2
			})

		# create
		similarity_df = pd.DataFrame(results).sort_values('jaccard_similarity', ascending=True).reset_index(drop=True)

		# non-zero similarities
		filt_jaccard_similarity = similarity_df['jaccard_similarity'] > 0
		filt_overlap_coefficient_PC1 = similarity_df['overlap_coefficient_PC1'] > 0
		filt_overlap_coefficient_PC2 = similarity_df['overlap_coefficient_PC2'] > 0

		filt_non_zero = filt_jaccard_similarity | filt_overlap_coefficient_PC1 | filt_overlap_coefficient_PC2

		# self
		not_self = (similarity_df['PC1'] != similarity_df['PC2'])

		# non-zero and not self
		filt = filt_non_zero & not_self

		# apply
		similarity_df = similarity_df.loc[filt].reset_index(drop=True)

		return similarity_df


	def get_figure_3(df, sub_gap=0.2, gap=0.3):

		# Add missing loci
		missing_loci = [locus for locus in self.phenotypes if locus not in df['locus'].unique()]
		missing_df = pd.DataFrame({'locus': missing_loci})
		df = pd.concat([df, missing_df], ignore_index=True)


		# Convert 'recall' to numeric, handling NaNs
		df['recall'] = pd.to_numeric(df['recall'], errors='coerce')

		# Sort loci
		df['sort_locus'] = df['locus'].str.extract('(\d+)').astype(float)
		df = df.sort_values(['sort_locus', 'recall', 'precision'], ascending=[True, False, False]).reset_index(drop=True)

		# Parameters
		loci = list(df['locus'].unique())
		per_locus_n_points = df.groupby('locus').size().to_dict()

		# Plot
		fig, ax = plt.subplots(figsize=(20, 3))

		### Calculate X positions

		# Init
		starting_position = 0
		x_positions = []
		label_positions = {}  # Store label positions

		# Get
		for locus in loci:
		    # Points
		    locus_n_points = per_locus_n_points.get(locus, 1)
		    x_positions_locus = []
		    for j in range(locus_n_points):
		        x_positions_locus.append(starting_position)
		        x_positions.append(starting_position)
		        starting_position += sub_gap

		    # Calculate label position
		    if locus_n_points % 2 == 1:
		        # Odd number of points: label under the central point
		        label_position = np.median(x_positions_locus)
		    else:
		        # Even number of points: label between the two central points
		        mid_index = locus_n_points // 2
		        label_position = (x_positions_locus[mid_index - 1] + x_positions_locus[mid_index]) / 2

		    label_positions[locus] = label_position
		    starting_position += gap

		# Add x positions for each point
		df['x_position'] = x_positions

		# Drop rows with missing critical values
		df = df.dropna(subset=['recall', 'color', 'shape', 'multitopologies'])

		# Validate and map shapes
		shape_mapping = {'square': 's', 'circle': 'o'}
		df['marker'] = df['shape'].map(shape_mapping).fillna('o')

		# Plot points
		for idx, row in df.iterrows():
		    locus = row['locus']
		    x_position = row['x_position']
		    recall = row['recall']
		    color = row['color']
		    marker = row['marker']
		    multitopology = row['multitopologies']
		    edge_width = 2 if multitopology else 0.5

		    # Plot point
		    ax.scatter(
		            x_position, recall,
		            s=100,
		            c=color,
		            marker=marker,
		            edgecolor='black',
		            linewidths=edge_width
		    )

		### Add labels at calculated positions
		for locus, label_pos in label_positions.items():
		    ax.text(label_pos, -0.5, locus, ha='center', va='center', fontsize=10, rotation=45,fontweight='bold')

		# Set gridlines for all X positions
		ax.set_xticks(x_positions, minor=False)
		ax.grid(visible=True, which='major', axis='both', alpha=0.4)
		ax.tick_params(axis='x', which='both', size=20, bottom=False, top=False, labelbottom=False)

		# Polish
		ax.set_ylim(-0.1, 1.1)
		ax.set_ylabel('Recall', size=12, fontweight='bold')

		# Final
		plt.title('Strong predictors with Lasso regression (alpha = 0.8)\n')
		plt.tight_layout()
		# plt.savefig(fig3, dpi=400)
		# plt.close()