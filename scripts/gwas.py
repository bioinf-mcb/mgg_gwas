from pathlib import Path
import pandas as pd
from functools import reduce
from scripts.utils import create_symlink
import subprocess
from tqdm import tqdm

class GWASWorkflow:
    def __init__(self, structure, params):

        # params
        self.phenotypes = params['gwas']['phenotypes']
        self.sample_size = params['gwas']['bootstrap_sample_size']
        self.nbootstrap = params['gwas']['nbootstrap']

        # input
        self.bacteria_tsv = structure['input_processed']['bacteria_tsv']
        self.mmseqs_variants = structure['mmseqs']['clusters_matrix']

        # output
        self.full_phenotype_files = []
        self.bootstrap_phenotype_files = []
        self.phenotypes_dir = structure['gwas']['phenotypes_dir']
        self.vectors_dir = structure['gwas']['phenotypes_vectors_dir']
        self.bootstrap_dir = structure['gwas']['bootstrap_phenotypes_dir']
        
        self.phenotypes_matrix = structure['gwas']['phenotypes_matrix']
        self.bootstrap_phenotypes_matrix = structure['gwas']['bootstrap_phenotypes_matrix']

        self.gwas_variants = structure['gwas']['variants_tsv']
        self.gwas_n_variants = structure['gwas']['n_variants_tsv']
        self.covariates_tsv = structure['gwas']['covariates_tsv']
        
        self.lasso_dir = structure['gwas']['lasso_dir']
        self.elastic_net_dir = structure['gwas']['elastic_net_dir']
        self.incorrect_dir = structure['gwas']['incorrect_dir']

        self.lasso_sh = structure['gwas']['lasso_sh']
        self.elastic_net_sh = structure['gwas']['elastic_net_sh']
        self.incorrect_sh = structure['gwas']['incorrect_sh']
        
        # manual alpha adjustment
        self.manual_alpha_adjustment = params['manual_alpha_adjustment']

    def get_input_files(self):

        self._get_phenotypes()
        self._get_variants()
        self._get_covariates()


    def _get_phenotypes(self): 
        """ generates phenotype files with and without bootstrapping """
        
        # prompt
        print(f'{len(self.phenotypes)} phenotypes {self.nbootstrap}x bootstrap with replacement [sample size = {self.sample_size}] [genomes without prophages excluded]')

        # generate phenotype files 
        for locus in self.phenotypes:
            print(f"{locus}.. ", end='')
            self._phenotypes_vectors(locus)
            self._phenotypes_bootstrap(locus)

        # concatenate phenotypes to matrix
        self._phenotypes_vectors_table()

        # concatenate bootstrap phenotypes to matrix
        self._phenotypes_bootstrap_table()

        print('\nDone!')


    def _phenotypes_vectors(self, locus):
        """ generate files for each locus with its presence/absence and save """

        # paths
        outfile = Path(self.vectors_dir, f'{locus}_0.tsv')

        # do not generate if exists
        if Path(outfile).exists(): pass
        else: 
            # presence/absence of a phenotype
            phenotype_df = self._get_phenotype(locus, outfile)

        # create list of all full phenotype files (pyseer)
        self.full_phenotype_files.append(outfile)



    def _phenotypes_bootstrap(self, locus):
        """ generate files for each locus with its sampled presence/absence and save"""
        
        # params
        width = len(str(self.nbootstrap))

        # output
        outdir = Path(self.bootstrap_dir, locus)

        # create
        Path(outdir).mkdir(exist_ok=True, parents=True)

        # bootstrap phenotype sampling
        for i in range(1, self.nbootstrap + 1):
            
            # path
            fname = f'{locus}_{str(i).zfill(width)}.tsv'
            outfile = Path(outdir, fname)

            # do not generate if exists
            if Path(outfile).exists(): pass
            else:
                # bootstrap
                phenotype_df = self._get_phenotype(locus, outfile, sample_size=self.sample_size)

            # create list of all bootstrap phenotype files (pyseer)
            self.bootstrap_phenotype_files.append(outfile)


    def _phenotypes_vectors_table(self):
        """ merge phenotypes presence/absence to one matrix """

        dfs = []
        for i, vector in enumerate(Path(self.vectors_dir).glob('*tsv')):
            if i == 0: vector_df = pd.read_csv(vector, sep='\t', usecols=[0,1])
            else: vector_df = pd.read_csv(vector, sep='\t', usecols=[1])
            dfs.append(vector_df)

        # matrix
        phenotypes_matrix = pd.concat(dfs, axis=1)

        # order
        cols = ['genomeID'] + self.phenotypes

        # rename
        rename_dict = {p: f'{p}_0' for p in self.phenotypes}
        phenotypes_matrix = phenotypes_matrix[cols].rename(rename_dict, axis=1)

        # save
        phenotypes_matrix.to_csv(self.phenotypes_matrix, sep='\t', index=False)



    def _phenotypes_bootstrap_table(self):
        """ merge sampeld phenotypes presence/absence to one matrix """ 

        # per bootstrap file
        dfs = []
        for file in self.bootstrap_phenotype_files:

            # params
            locus, i = Path(file).stem.split('_')
            locus_i = f'{locus}_{int(i)}'

            # load
            vector_df = pd.read_csv(file, sep='\t')
            
            # rename to locus to bootstrap
            vector_df = vector_df.rename({locus: locus_i}, axis=1)
            
            # append
            dfs.append(vector_df)

        # combine
        phenotypes_matrix = reduce(lambda left, right: pd.merge(left, right, on='genomeID', how='outer'), dfs)

        # save
        phenotypes_matrix.to_csv(self.bootstrap_phenotypes_matrix, sep='\t', index=False)



    def _get_variants(self):
        """ symbolic link to variants """
        create_symlink(self.mmseqs_variants, self.gwas_variants)


    def _get_covariates(self, i_sc_col=15):
        
        # load
        covariates_df = pd.read_csv(self.bacteria_tsv, sep='\t', usecols=[0, i_sc_col])

        # save
        covariates_df.to_csv(self.covariates_tsv, sep='\t', index=False)


    def _pyseer_command(self, variants_tsv, phenotype_tsv, covariates_tsv, correction, alpha, cor_filter, min_af, max_af, outfile, cpu=9):

        # components
        cmd = [
            'pyseer',
            '--phenotypes', f'"{phenotype_tsv}"',
            '--pres', f'"{variants_tsv}"',
            '--wg', 'enet',
            '--alpha', str(alpha),
            '--cor-filter', str(cor_filter),
            '--min-af', str(min_af),
            '--max-af', str(max_af),
            '--filter-pvalue', '1',
            '--lrt-pvalue', '1',
            '--cpu', str(cpu),
            '--uncompressed',
            '--print-samples'
        ]

        # add correction
        if correction:
            cmd.extend([
                '--lineage-clusters', f'"{covariates_tsv}"',
                '--sequence-reweighting'
            ])


        # Join the command components into a string
        cmd_str = ' '.join(cmd) + f' >> "{outfile}"'

        return cmd_str


    def get_script(self, mode, correction, alpha, cor_filter, min_af, max_af, activate_env='source ~/.zshrc; conda activate pyseer; '):
        
        # input
        variants_tsv = self.gwas_variants
        covariates_tsv = self.covariates_tsv
        phenotype_files = self.full_phenotype_files + self.bootstrap_phenotype_files

        # output
        if mode == 'lasso':
            outdir = self.lasso_dir
            pyseer_sh = self.lasso_sh
        elif mode == 'elastic_net':
            outdir = self.elastic_net_dir
            pyseer_sh = self.elastic_net_sh
        elif mode == 'incorrect': 
            outdir = self.incorrect_dir
            pyseer_sh = self.incorrect_sh
        else: 
            print('Error in get_script')
            exit()

        ######################################################################################
        # WARNIGN! For some capsules and clustering levels the alpha must be adjusted ########
        # to avoid matrix convergance error ##################################################
        ######################################################################################

        def apply_manual_alpha_adjustment(mode, variants_tsv, phenotype_tsv, default_alpha):

            # params
            PCIC = Path(variants_tsv).parent.parent.stem.split('_')[-2]
            k_locus = phenotype_tsv.stem.split('_')[0]
            keys = self.manual_alpha_adjustment.keys()

            # only for some PCICS
            if PCIC not in keys:
                return default_alpha

            # change only for lasso
            if mode == 'elastic_net':
                return default_alpha


            for key in keys:
                if key == PCIC:
                    values = self.manual_alpha_adjustment[key]
                    break
                    
            for value in values:
                if value == k_locus:
                    return 0.6

            return default_alpha

        ### prepare pyseer script
        cmd = []
        default_alpha = alpha
        for phenotype_tsv in phenotype_files:
            outfile = Path(outdir, phenotype_tsv.stem + '.tsv')

            # for some runs manually change value of alpha
            alpha = apply_manual_alpha_adjustment(mode, variants_tsv, phenotype_tsv, default_alpha=default_alpha)

            # skip existing files
            if Path(outfile).exists(): continue
            
            pyseer_vars = [variants_tsv, phenotype_tsv, covariates_tsv, correction, alpha, cor_filter, min_af, max_af, outfile]
            pyseer_cmd = self._pyseer_command(*pyseer_vars)
            cmd.append(pyseer_cmd)


        # script
        script = ' \n'.join(cmd)
        script = f'{activate_env}\n' + script
        
        # save
        with open(pyseer_sh, 'w') as f:
            f.write(script)



    def _get_phenotype(self, locus, path, sample_size=None, i_k_col=12, i_o_col=13, i_nprophages_col=43):
        """ Generate phenotype file. """

        # params
        icol = i_k_col if locus[:2] == 'KL' else i_o_col
        colname = 'MGG_K_locus' if locus[:2] == 'KL' else 'MGG_O_locus'

        # load
        bacteria_df = pd.read_csv(self.bacteria_tsv, sep='\t', usecols=[0, icol, i_nprophages_col])
        bacteria_df = bacteria_df.query("nprophages != 0").copy()
        bacteria_df = bacteria_df.drop('nprophages', axis=1)

        # binary phenotype
        bacteria_df[locus] = bacteria_df[colname].apply(lambda x: 1 if x == locus else 0)
        bacteria_df = bacteria_df.drop(colname, axis=1)
        phenotype_df = bacteria_df

        # sample
        if sample_size: 
            phenotype_df = phenotype_df.sample(n=sample_size, replace=True)
            phenotype_df = phenotype_df.drop_duplicates('genomeID')
        
        # save
        phenotype_df.to_csv(path, sep='\t', index=False)

        return phenotype_df


    def run_scripts(self):

        print('Running lasso... ', end='')
        result = subprocess.run(f'bash "{self.lasso_sh}"', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print('Done!')

        print('Running elastic net... ', end='')
        result = subprocess.run(f'bash "{self.elastic_net_sh}"', shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print('Done!')

        # print('Running elastic net scripts... ', end='')
        # run(f'bash "{self.incorrect_sh}"', shell=True)
        # print('Done!')


    def compute_n_variants(self, run=True):

        # checkpoint
        if run: pass
        else: return

        # checkpoint
        if Path(self.gwas_n_variants).exists():
            print('Number of variants per locus exist... remove existing file to run again. ')
            return None
        else: pass
        

        # load
        phenotypes_df = pd.read_csv(self.phenotypes_matrix, sep='\t')
        phenotypes_btsp_df = pd.read_csv(self.bootstrap_phenotypes_matrix, sep='\t')
        variants_df = pd.read_csv(self.gwas_variants, sep='\t').T

        # Set the first row as the header
        variants_df.columns = variants_df.iloc[0]  # Set the first row as header
        variants_df = variants_df[1:]              # Remove the first row from the DataFrame

        # Optionally, reset the index
        variants_df = variants_df.reset_index(drop=False).rename(columns={'index': 'genomeID'})

        # merge and clean
        data = phenotypes_df.merge(phenotypes_btsp_df, on='genomeID', how='outer').merge(variants_df, on='genomeID', how='outer')
        data = data.drop('genomeID', axis=1)
    
        ### compute

        # select cols
        locus_bootstrap_cols = [col for col in data.columns if not col.startswith('PC')]
        pc_cols = [col for col in data.columns if col.startswith('PC')]
        
        # init
        locus_pc_count = {}
        
        # iterate
        for locus in tqdm(locus_bootstrap_cols, desc='Calculating number of variants per run: ', ncols=50, ascii=True):
            submatrix = data.loc[~data[locus].isna(), pc_cols]
            valid_pc_cols = submatrix.loc[:, (submatrix != 0).any(axis=0) & (submatrix != 1).any(axis=0)]
            num_pcs = valid_pc_cols.shape[1]
            locus_pc_count[locus] = num_pcs

        # data frame
        n_variants_df = pd.DataFrame.from_dict(locus_pc_count, orient='index', columns=['n_PC_tested'])
        n_variants_df = n_variants_df.reset_index()
        n_variants_df = n_variants_df.rename({'index': 'locus_nbootstrap'}, axis=1)

        # save
        n_variants_df.to_csv(self.gwas_n_variants, sep='\t', index=False)
