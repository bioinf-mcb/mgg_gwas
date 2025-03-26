import pandas as pd
from pathlib import Path
import phylotreelib as pt

import warnings

# Suppress all warnings
warnings.filterwarnings('ignore')

class Phandango:
    def __init__(self, structure, params):

        # params
        self.phenotypes = params['gwas']['phenotypes']

        # input
        self.gwas_variants = structure['gwas']['variants_tsv']
        self.pyseer_metrics_topologies = structure['processed']['pyseer_metrics_topologies']
        self.tree_file = structure['input']['bacteria_tree']
        self.phenotypes_matrix = structure['gwas']['phenotypes_matrix']

        # output
        self.phandango_dir = structure['processed']['phandango_dir']

        # load
        self.variants_df = self._reshape_variants_table()
        self.phenotypes_matrix_df = self._rename_phenotype_table_columns()
        self.pyseer_df = pd.read_csv(self.pyseer_metrics_topologies, sep='\t')
        

        # select variants params
        min_precision = None
        min_recall = None
        min_beta = None
        min_minus_log_pvalue = None
        function = None

        self.select_variants_vars = min_precision, min_recall, min_beta, min_minus_log_pvalue, function



    def run(self, phenotype_variants_map):
        """ generate phandango files """

        # create structure
        self.folders_dict = self._get_folders_dict()

        # get phandango tables per phenotype
        for phenotype in self.phenotypes:

            # per mode
            for mode in ['lasso']:

                # paths
                outdir = self.folders_dict[phenotype][mode]

                self.all_variants_outfile    = Path(outdir, '_variants.csv')
                self.subtree_outfile         = Path(outdir, '_subtree.nwk')
                # self.select_variants_outfile = Path(outdir, '3_SELECTED_VARIANTS.csv')

                # get all variants and sub-tree
                variants = self._get_locus_associated_variants(phenotype, mode)

                try: variants = phenotype_variants_map[phenotype]
                except KeyError: continue
                
                self._get_subtable(phenotype, variants)
                self._get_subtree()
                    

    
    def _get_subtree(self, max_leaves=300, outgroup='398KBV'):

        # load
        variants_df = pd.read_csv(self.all_variants_outfile, sep='\t', index_col=0)

        # remove colors
        cols = [col for col in variants_df.columns if ':colour' not in col]
        variants_df = variants_df[cols]

        # variants and phenotype abundance
        variants_abundance = variants_df.sum(axis=1)
        
        ### define leaves to retain

        # to retain always
        leaves_with_variant = variants_abundance.loc[variants_abundance > 0].index

        # to down sample
        n_leaves_to_sample = max(0, (max_leaves - len(leaves_with_variant)))
        if not n_leaves_to_sample: print('Associated variants are in too many genomes. No genomes to down sample.')
        leaves_without_variant_sampled = variants_abundance.loc[variants_abundance == 0].sample(n=n_leaves_to_sample).index

        # combined
        leaves2retain = set(leaves_with_variant).union(set(leaves_without_variant_sampled))
        leaves2retain = list(leaves2retain) + [outgroup]

        ### subtree

        # load tree
        with pt.Newicktreefile(self.tree_file) as treefile:
            subtree = treefile.readtree()

        # re-draw subtree
        leaves_all = subtree.leaflist()
        leaves2remove = set(leaves_all).difference(set(leaves2retain))

        # remove leaves
        subtree.remove_leaves(leaves2remove)

        # root
        subtree.rootout(outgroup)

        # save
        with open(self.subtree_outfile, 'w') as f:
            f.write(subtree.newick())



    def _get_subtable(self, phenotype, variants, map_colours={'phenotype': '#273746', 'PC': '#5DADE2', 'other': '#FDFEFE'}):
        """ select variants in table """

        # selected variants
        variants_cols = ['genomeID'] + variants
        selected_variants_df = self.variants_df[variants_cols]

        # selected phenotype
        phenotype_cols = ['genomeID', phenotype]
        phenotype_df = self.phenotypes_matrix_df[phenotype_cols]

        # subtable
        df = phenotype_df.merge(selected_variants_df, on='genomeID', how='left')
        
        # assign phenotype colours
        phenotype_map_colours = {1: map_colours['phenotype'], 0: map_colours['other']}
        df[f'{phenotype}:colour'] = df[phenotype].map(phenotype_map_colours)

        # assign PC colours
        variants_map_colours = {1: map_colours['PC'], 0: map_colours['other']}
        for col in variants_cols[1:]:
            df[f'{col}:colour'] = df[col].map(variants_map_colours)

        # select & save
        df.to_csv(self.all_variants_outfile, sep='\t', index=False)

        ### split into smaller tables

        def split_list(lst, chunk_size):
            """Split a list into chunks of given size using list comprehension."""
            return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

        table_fragments_cols = split_list(variants, chunk_size=3)

        for i, variants_cols in enumerate(table_fragments_cols):
            outdir, (name, extension) = self.all_variants_outfile.parent, self.all_variants_outfile.name.split('.')
            outfile = Path(outdir, f'FRAGMENT_{i}.{extension}')

            variants_cols_colours = [f'{col}:colour' for col in variants_cols]            
            cols = ['genomeID', phenotype] + variants_cols + [f'{phenotype}:colour'] + variants_cols_colours
            df[cols].to_csv(outfile, sep='\t', index=False)


    
    def _get_locus_associated_variants(self, locus, mode, min_precision=None, min_recall=None, min_beta=None, min_minus_log_pvalue=None, functions=None):
        """ select variants only for given phenotype """


        # filters
        filt_locus = (self.pyseer_df['locus'] == locus)
        filt_mode  = (self.pyseer_df['mode'] == mode)

        filt = filt_locus & filt_mode

        # select variants
        selected_variants_df = self.pyseer_df.loc[filt]


        # extract
        variants = list(selected_variants_df['PC'].unique())

        return variants



    def _get_folders_dict(self):
        
        # init
        folders_dict = {}
        for phenotype in self.phenotypes:

            # create dict
            folders_dict[phenotype] = {}

            # for method
            for i, mode_key in enumerate(['lasso', 'elastic_net', 'incorrect']):
                
                # folder   
                mode_folder = f'{i}_{mode_key.upper()}'

                # path
                path = Path(self.phandango_dir, phenotype, mode_folder)

                # create
                path.mkdir(exist_ok=True, parents=True)

                # append
                folders_dict[phenotype][mode_key] = path

        return folders_dict
            


    def _reshape_variants_table(self):

        # load
        variants_df = pd.read_csv(self.gwas_variants, sep='\t')

        # clean variants
        variants_df.index = variants_df.iloc[:, 0]
        variants_df = variants_df.drop(variants_df.columns[0], axis=1)
        variants_df.index.name = ''
        variants_df = variants_df.T.reset_index().rename({'index': 'genomeID'}, axis=1)

        return variants_df


    def _rename_phenotype_table_columns(self):

        # load
        phenotypes_matrix_df = pd.read_csv(self.phenotypes_matrix, sep='\t')

        # rename cols
        cols = phenotypes_matrix_df.columns[1:]
        rename_dict = {col: col.split('_')[0] for col in cols}

        phenotypes_matrix_df = phenotypes_matrix_df.rename(rename_dict, axis=1)

        return phenotypes_matrix_df

