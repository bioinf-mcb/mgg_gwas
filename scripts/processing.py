import pandas as pd
from io import StringIO
from pathlib import Path
import numpy as np
from tqdm import tqdm
import math
from scipy import stats



class Processor:
    def __init__(self, structure, params):

        # input
        self.version = params['mmseqs']['version']
        self.phenotypes = params['gwas']['phenotypes']
        self.nbootstrap = params['gwas']['nbootstrap']
        self.report_topology_map = params['report_topology_map']

        self.lasso_dir = structure['gwas']['lasso_dir']
        self.elastic_net_dir = structure['gwas']['elastic_net_dir']
        self.incorrect_dir = structure['gwas']['incorrect_dir']

        self.phenotypes_matrix = structure['gwas']['phenotypes_matrix']
        self.bootstrap_phenotypes_matrix = structure['gwas']['bootstrap_phenotypes_matrix']
        self.gwas_variants = structure['gwas']['variants_tsv']
        self.gwas_n_variants = structure['gwas']['n_variants_tsv']

        self.bacteria_table = structure['input_processed']['bacteria_tsv']
        self.clusters_mmseqs = structure['mmseqs']['clusters']
        self.clusters_functions = structure['functions']['clusters_functions']


        # output
        self.pyseer_cols = ['PC', 'pvalue', 'beta', 'locus_nbootstrap', 'mode']
        self.metrics_cols = ['PC', 'locus_nbootstrap', 'TP', 'FN', 'TN', 'FP', 'precision', 'recall', 'specificity', 'accuracy', 'F1_score', 'MCC', 'TP_SC','FN_SC','TN_SC','FP_SC','precision_SC','recall_SC','specificity_SC','accuracy_SC','F1_score_SC', 'MCC_SC', 'PC_abundance', 'locus_abundance', 'PC_freq']
        self.pyseer_bootstrap_table_cols = ['version', 'PC', 'mode', 'locus', 'nbootstrap', 'minus_log10_pvalue_corr', 'beta', 'precision', 'recall',  'precision_SC','recall_SC', 'PC_abundance', 'locus_abundance', 'PC_freq', 'TP', 'FN', 'TN', 'FP', 'specificity', 'accuracy', 'F1_score', 'MCC', 'TP_SC','FN_SC','TN_SC','FP_SC','specificity_SC','accuracy_SC','F1_score_SC', 'MCC_SC', 'pvalue', 'n_PC_tested', 'pvalue_corr', 'minus_log10_pvalue']
        self.pyseer_aggregated_cols = ['version', 'PC', 'mode', 'locus', 'nbootstrap', 'minus_log10_pvalue_corr', 'beta', 'precision', 'recall',  'precision_SC','recall_SC', 'PC_abundance', 'locus_abundance', 'PC_freq', 'TP', 'FN', 'TN', 'FP', 'specificity', 'accuracy', 'F1_score', 'MCC', 'TP_SC','FN_SC','TN_SC','FP_SC','specificity_SC','accuracy_SC','F1_score_SC', 'MCC_SC', 'pvalue', 'n_PC_tested', 'pvalue_corr', 'minus_log10_pvalue']
        
        self.pyseer_tsv = structure['processed']['raw_pyseer_table']
        self.processed_pyseer_tsv = structure['processed']['clean_pyseer_table']
        self.variants_metrics = structure['processed']['pcs_metrics_table']
        self.pyseer_bootstrap_table = structure['processed']['pyseer_bootstrap_table']

        # pyseer hits: per clustering (version)
        self.pyseer_hits_version_table_all = structure['processed']['pyseer_hits_version_table_all']
        self.pyseer_hits_version_table_filtered = structure['processed']['pyseer_hits_version_table_filtered']

        # pyseer hits: aggregated
        self.pyseer_hits_final_table_all = structure['processed']['pyseer_hits_final_table_all']
        self.pyseer_hits_final_table_filtered = structure['processed']['pyseer_hits_final_table_filtered']

        # functions: aggregated
        self.clusters_functions_best_all = structure['functions']['clusters_functions_best_all'] 


    def concatenate_pyseer(self):


        # create new file with header
        with open(self.pyseer_tsv, 'w') as f:
            f.write('\t'.join(self.pyseer_cols) + '\n') 

        self._process_and_concatenate_files(self.lasso_dir, mode='lasso')    
        self._process_and_concatenate_files(self.elastic_net_dir, mode='elastic_net')
        self._process_and_concatenate_files(self.incorrect_dir, mode='incorrect')

        # checkpoint
        with open(self.pyseer_tsv, 'r') as f:
            line1 = f.readline()
            line2 = f.readline()

            if not line2:
                print('STOP! There are no pyseer results.')
                exit()


        # add beta zero lines
        self._add_beta_zero()


    def _process_and_concatenate_files(self, pyseer_results_dir, mode):


        # per raw file (correct)
        pyseer_files = self._get_pyseer_raw_files(pyseer_results_dir)
        
        # dummy checkpoint
        self._pyseer_results_files_checkpoint(pyseer_files, mode)

        # extract from files
        for file in pyseer_files:

            # process raw table
            processed_df = self._process_raw(file, mode)

            # append to final table
            self._append_to_file(processed_df, self.pyseer_cols, self.pyseer_tsv)


    def _pyseer_results_files_checkpoint(self, pyseer_files, mode):
        
        if not len(pyseer_files): 
            print(f'WARNING! No pyseer results found for mode {mode}!')


    def _add_beta_zero(self):

        # Load data
        pyseer_df = pd.read_csv(self.pyseer_tsv, sep='\t')

        # Generate permutations for each phenotype
        permutations_dict = {
            phenotype: [f'{phenotype}_{n}' for n in range(self.nbootstrap + 1)]
            for phenotype in self.phenotypes
        }


        # Initialize dictionary to store missing data
        missing_data = {
            'PC': [],
            'locus_nbootstrap': [],
            'pvalue': [],
            'beta': [],
            'mode': []
        }
        ### identify missing bootstraps for each PC and mode
        pyseer_df['locus'] = pyseer_df['locus_nbootstrap'].str.split('_', expand=True)[0]

        for (pc, locus, mode), group in pyseer_df.groupby(['PC', 'locus', 'mode']):
            locus_bootstrap_done = set(group['locus_nbootstrap'])
            locus = next(iter(locus_bootstrap_done)).split('_')[0]
            locus_bootstrap_all = set(permutations_dict[locus])
            locus_bootstrap_missing = locus_bootstrap_all - locus_bootstrap_done

            for locus_bootstrap in locus_bootstrap_missing:
                missing_data['PC'].append(pc)
                missing_data['locus_nbootstrap'].append(locus_bootstrap)
                missing_data['pvalue'].append(np.nan)
                missing_data['beta'].append(np.nan)
                missing_data['mode'].append(mode)
        
        # clean
        pyseer_df = pyseer_df.drop(['locus'], axis=1)

        # convert
        missing_df = pd.DataFrame(missing_data)

        # combine
        combined_df = pd.concat([pyseer_df, missing_df], ignore_index=True)

        # sort
        combined_df[['locus', 'nbootstrap']] = combined_df['locus_nbootstrap'].str.split('_', expand=True)
        combined_df['nbootstrap'] = combined_df['nbootstrap'].astype(int)
        combined_df = combined_df.sort_values(['locus', 'mode', 'PC', 'nbootstrap'], ascending=[True, True, True, True])

        # clean
        combined_df = combined_df.drop(['locus', 'nbootstrap'], axis=1)

        # save
        combined_df.to_csv(self.processed_pyseer_tsv, sep='\t', index=False)


    def _process_raw(self, file, mode):

        # params
        locus, nbootstrap = Path(file).stem.split('_')
        nbootstrap = str(int(nbootstrap))

        # clean [warnings, errors]
        header_tsv = ['variant\taf\tfilter-pvalue\tlrt-pvalue\tbeta\tlineage\tk-samples\tnk-samples\tnotes\n']
        with open(file, 'r') as f:
            lines = [line for line in f.readlines() if line[:2] == 'PC']
            lines = header_tsv + lines
            raw_file_io = StringIO('\n'.join(lines))

        # load
        results_df = pd.read_csv(raw_file_io, sep='\t', usecols=[0,1,2,4])
        
        # rename cols
        colnames_map = {'variant': 'PC', 'filter-pvalue': 'pvalue'}
        results_df = results_df.rename(columns=colnames_map)

        ### compute -log10
        results_df['pvalue'] = results_df['pvalue'].astype('float')
        results_df.loc[results_df['pvalue'] == 0, 'pvalue'] = 10**-200
        results_df['pvalue'] = np.round(results_df['pvalue'], 3)
        results_df['beta'] = np.round(results_df['beta'], 3)

        # add cols
        results_df['locus_nbootstrap'] = f'{locus}_{nbootstrap}'
        results_df['mode'] = mode

        return results_df


    def _get_pyseer_raw_files(self, path):
        raw_files = list(Path(path).glob('*.tsv'))
        return raw_files


    def _append_to_file(self, df, cols, outfile):
        """Appends the dataframe to the output file in tsv format."""

        df[cols].to_csv(outfile, mode='a', sep='\t', header=False, index=False)


    def compute_metrics(self, run=True):
        
        # checkpoint
        if run: pass
        else: return

        # checkpoint
        if Path(self.variants_metrics).exists():
            print('Metrics for variants exist... remove existing files to run again. ')
            return None
        else: pass
        

        # load
        bacteria_df = pd.read_csv(self.bacteria_table, sep='\t')
        phenotypes_matrix_df = pd.read_csv(self.phenotypes_matrix, sep='\t')
        bootstrap_phenotypes_matrix_df = pd.read_csv(self.bootstrap_phenotypes_matrix, sep='\t')
        variants_df = pd.read_csv(self.gwas_variants, sep='\t')
        pyseer_df = pd.read_csv(self.pyseer_tsv, sep='\t')

        # clean
        pyseer_df[['locus', 'nbootstrap']] = pyseer_df['locus_nbootstrap'].str.split('_', expand=True)
        
        # clean & remove isolates without prophages
        filt_with_prophages = (bacteria_df['nprophages'] > 0)
        bacteria_df = bacteria_df.loc[filt_with_prophages].reset_index(drop=True)
        bacteria_df = bacteria_df[['genomeID', 'MGG_SC']].copy()

        # convert variants matrix
        variants_df.index = variants_df.iloc[:, 0]
        variants_df = variants_df.drop(variants_df.columns[0], axis=1)
        variants_df.index.name = 'PC'
        variants_df = variants_df.T.reset_index().rename({'index': 'genomeID'}, axis=1)
        
        # merge to big matrix
        phenotypes_df = bacteria_df.merge(phenotypes_matrix_df, on='genomeID', how='outer')
        phenotypes_df = phenotypes_df.merge(bootstrap_phenotypes_matrix_df, on='genomeID', how='outer')
        big_matrix_df = phenotypes_df.merge(variants_df, on='genomeID', how='outer')
        
        # maps
        variants_per_locus = pyseer_df.groupby('locus')['PC'].unique().apply(list).to_dict()
        
        bootstrap_per_locus = {}
        for locus in self.phenotypes:
            btsp_per_phenotype = []
            for btsp in list(map(str, list(range(0, self.nbootstrap+1)))):
                btsp_per_phenotype.append(f'{locus}_{btsp}')
            bootstrap_per_locus[locus] = btsp_per_phenotype


        # create new file with header
        with open(self.variants_metrics, 'w') as f:
            f.write('\t'.join(self.metrics_cols) + '\n') 


        # pairwise locus_bts vs variant
        for locus in tqdm(self.phenotypes, desc='Calculating metrics... ', ncols=50):
            
            locus_bootstraps_list = bootstrap_per_locus[locus]
            variants_list = variants_per_locus[locus]


            for locus_btsp in locus_bootstraps_list:
                for variant in variants_list:
                    cols = ['MGG_SC', locus_btsp, variant]
                    
                    subset_matrix = big_matrix_df[cols]
                    subset_matrix = subset_matrix.dropna(subset=[locus_btsp])
                    subset_sc_groups = subset_matrix.groupby('MGG_SC')

                    confusion_matrix_vars = [subset_matrix, subset_sc_groups, locus_btsp, variant]
                    confusion_matrix = self._confusion_matrix(*confusion_matrix_vars)
                    self._append_to_file(confusion_matrix, self.metrics_cols, self.variants_metrics)


    def _confusion_matrix(self, df, sc_groups, phenotype, variant):
        # Convert to boolean arrays
        phenotype_mask = (df[phenotype].values == 1)
        variant_mask = (df[variant].values == 1)
        length = len(df)

        pc_abundance = variant_mask.sum()
        phenotype_abundance = phenotype_mask.sum()
        pc_freq = pc_abundance / length if length > 0 else 0

        # Edge case if df is empty
        if not pc_abundance or not phenotype_abundance:
            data = {
                'PC': variant,
                'locus_nbootstrap': phenotype,
                'TP': np.nan, 'FN': np.nan, 'TN': np.nan, 'FP': np.nan,
                'precision': np.nan, 'recall': np.nan, 'specificity': np.nan, 'accuracy': np.nan, 'F1_score': np.nan, 'MCC': np.nan,
                'TP_SC': np.nan, 'FN_SC': np.nan, 'TN_SC': np.nan, 'FP_SC': np.nan,
                'precision_SC': np.nan, 'recall_SC': np.nan, 'specificity_SC': np.nan, 'accuracy_SC': np.nan, 'F1_score_SC': np.nan, 'MCC_SC': np.nan,
                'PC_abundance': pc_abundance,
                'locus_abundance': phenotype_abundance,
                'PC_freq': "{:.2e}".format(pc_freq)
            }
            return pd.DataFrame([data])


        # Confusion matrix masks
        TP_mask = phenotype_mask & variant_mask
        FN_mask = phenotype_mask & ~variant_mask
        TN_mask = ~phenotype_mask & ~variant_mask
        FP_mask = ~phenotype_mask & variant_mask

        # Confusion matrix counts (row-level)
        TP = TP_mask.sum()
        FN = FN_mask.sum()
        TN = TN_mask.sum()
        FP = FP_mask.sum()

        recall = _safe_round(_safe_division(TP, TP + FN), 2)
        precision = _safe_round(_safe_division(TP, TP + FP), 2)
        specificity = _safe_round(_safe_division(TN, TN + FP), 2)
        accuracy = _safe_round(_safe_division(TP + TN, TP + TN + FP + FN), 2)
        f1_score = _safe_round(_safe_division(2 * precision * recall, precision + recall), 2)
        mcc = _safe_round(_safe_division((TP * TN) - (FP * FN), math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))), 2)

        # SC-level classification
        # We'll categorize each SC into exactly one of TP, TN, FP, FN, following the logic:
        # 1. If at least one TP row: SC is TP.
        # 2. Else if no phenotype and no variant: SC is TN.
        # 3. Else if no phenotype but variant present: SC is FP.
        # 4. Else if phenotype present but no variant: SC is FN.
        # 5. Else (phenotype and variant present but no TP): SC is FN AND FP!.

        TP_SC_count = 0
        TN_SC_count = 0
        FP_SC_count = 0
        FN_SC_count = 0

        for grp_name, idx_array in sc_groups.indices.items():
            # Check conditions for this SC
            grp_TP = TP_mask[idx_array].any()
            grp_has_phenotype = phenotype_mask[idx_array].any()
            grp_has_variant = variant_mask[idx_array].any()

            if grp_TP:
                # Group is TP
                TP_SC_count += 1
            else:
                # No direct TP rows
                if not grp_has_phenotype and not grp_has_variant:
                    # TN scenario
                    TN_SC_count += 1
                elif not grp_has_phenotype and grp_has_variant:
                    # FP scenario
                    FP_SC_count += 1
                elif grp_has_phenotype and not grp_has_variant:
                    # FN scenario
                    FN_SC_count += 1
                else:
                    # grp_has_phenotype and grp_has_variant are both True, but no TP
                    # This is also FN and FP
                    FP_SC_count += 1
                    FN_SC_count += 1


        # Compute SC-level metrics
        TP_SC = TP_SC_count
        FN_SC = FN_SC_count
        TN_SC = TN_SC_count
        FP_SC = FP_SC_count

        recall_SC = _safe_round(_safe_division(TP_SC, TP_SC + FN_SC), 2)
        precision_SC = _safe_round(_safe_division(TP_SC, TP_SC + FP_SC), 2)
        specificity_SC = _safe_round(_safe_division(TN_SC, TN_SC + FP_SC), 2)
        accuracy_SC = _safe_round(_safe_division(TP_SC + TN_SC, TP_SC + TN_SC + FP_SC + FN_SC), 2)
        f1_score_SC = _safe_round(_safe_division(2 * precision_SC * recall_SC, precision_SC + recall_SC), 2)
        mcc_SC = _safe_round(_safe_division((TP_SC * TN_SC) - (FP_SC * FN_SC), math.sqrt((TP_SC + FP_SC) * (TP_SC + FN_SC) * (TN_SC + FP_SC) * (TN_SC + FN_SC))), 2)


        data = {
            'PC': variant,
            'locus_nbootstrap': phenotype,
            'TP': TP,
            'FN': FN,
            'TN': TN,
            'FP': FP,
            'precision': precision,
            'recall': recall,
            'specificity': specificity,
            'accuracy': accuracy,
            'F1_score': f1_score,
            'MCC': mcc,
            'TP_SC': TP_SC,
            'FN_SC': FN_SC,
            'TN_SC': TN_SC,
            'FP_SC': FP_SC,
            'precision_SC': precision_SC,
            'recall_SC': recall_SC,
            'specificity_SC': specificity_SC,
            'accuracy_SC': accuracy_SC,
            'F1_score_SC': f1_score_SC,
            'MCC_SC': mcc_SC,
            'PC_abundance': pc_abundance,
            'locus_abundance': phenotype_abundance,
            'PC_freq': "{:.2e}".format(pc_freq)
        }

        return pd.DataFrame([data])


    def combine_info_bootstrap(self):

        # load
        metrics_df = pd.read_csv(self.variants_metrics, sep='\t')
        btsp_pyseer_df = pd.read_csv(self.processed_pyseer_tsv, sep='\t')
        n_variants_df = pd.read_csv(self.gwas_n_variants, sep='\t')

        # version
        btsp_pyseer_df['version'] = self.version

        #################
        #### METRICS ####
        #################

        # prepare merge
        metrics_df['PC_locus_nbootstrap'] = metrics_df['PC'] + '_' + metrics_df['locus_nbootstrap']
        btsp_pyseer_df['PC_locus_nbootstrap'] = btsp_pyseer_df['PC'] + '_' + btsp_pyseer_df['locus_nbootstrap']

        # clean variants
        metrics_df = metrics_df.drop(['PC', 'locus_nbootstrap'], axis=1)

        ### merge

        # metrics
        btsp_pyseer_df = btsp_pyseer_df.merge(metrics_df, on='PC_locus_nbootstrap', how='left')

        # tested variants
        btsp_pyseer_df = btsp_pyseer_df.merge(n_variants_df, on='locus_nbootstrap', how='left')


        # clean
        btsp_pyseer_df[['PC', 'locus', 'nbootstrap']] = btsp_pyseer_df['PC_locus_nbootstrap'].str.split('_', expand=True)
        btsp_pyseer_df = btsp_pyseer_df.drop(['PC_locus_nbootstrap', 'locus_nbootstrap'], axis=1)

        # format
        btsp_pyseer_df['PC_freq'] = btsp_pyseer_df['PC_freq'].apply(lambda x: "{:.2e}".format(x))
        btsp_pyseer_df['nbootstrap'] = btsp_pyseer_df['nbootstrap'].astype(int)
        
        # sort
        btsp_pyseer_df = btsp_pyseer_df.sort_values(['locus', 'mode', 'PC', 'nbootstrap'], ascending=[True, False, True, True])

        ###############
        #### STATS ####
        ###############

        # calculate Bonferroni pvalues
        btsp_pyseer_df['pvalue_corr'] = btsp_pyseer_df['pvalue'] * btsp_pyseer_df['n_PC_tested']
        btsp_pyseer_df.loc[btsp_pyseer_df['pvalue_corr'] > 1, 'pvalue_corr'] = 1

        # convert
        btsp_pyseer_df['minus_log10_pvalue'] = -np.log10(btsp_pyseer_df['pvalue'])
        btsp_pyseer_df['minus_log10_pvalue_corr'] = -np.log10(btsp_pyseer_df['pvalue_corr'])


        ###############
        #### SAVE #####
        ###############

        btsp_pyseer_df = btsp_pyseer_df[self.pyseer_bootstrap_table_cols]
        btsp_pyseer_df.to_csv(self.pyseer_bootstrap_table, sep='\t', index=False)


    def pyseer_hits_with_CI(self):


        # load
        btsp_pyseer = pd.read_csv(self.pyseer_bootstrap_table, sep='\t')

        ###############
        #### STATS ####
        ###############

        # compute stats
        pyseer_hits_CI = self._get_confidence_intervals_per_pc(btsp_pyseer)

        #######################################
        #### SIGNIFICANT POSITIVE NON-ZERO ####
        #######################################

        # significant positive non-zero associations in btsp=0
        pyseer_hits_CI_filtered = self._significant_non_zero_positive_hits(pyseer_hits_CI)


    def _significant_non_zero_positive_hits(self, pyseer_hits):

        # process
        pyseer_hits['PC_locus_mode'] = pyseer_hits['PC'] + '_' + pyseer_hits['locus'] + '_' + pyseer_hits['mode']
        pyseer_hits = pyseer_hits.loc[pyseer_hits['mode'].isin(['lasso', 'elastic_net'])].reset_index(drop=True)

        # filter 
        filt_main_run = (pyseer_hits['nbootstrap'] == 0)        # no sampling
        filt_significance = (pyseer_hits['pvalue_corr'] < 0.05) # level of significance (alpha): 5%
        filt_beta = (pyseer_hits['beta'] > 0)                   # positive non-zero beta

        # apply
        filt = filt_main_run & filt_significance & filt_beta
        pc_locus_mode = pyseer_hits.loc[filt, 'PC_locus_mode'].to_list()
        pyseer_hits = pyseer_hits.loc[pyseer_hits['PC_locus_mode'].isin(pc_locus_mode)].reset_index(drop=True)

        # clean
        pyseer_hits = pyseer_hits.drop('PC_locus_mode', axis=1)

        # save
        pyseer_hits.to_csv(self.pyseer_hits_version_table_filtered, sep='\t', index=False)

        return pyseer_hits


    def _get_confidence_intervals_per_pc(self, pyseer_hits):


        def get_confidence_intervals(df, column, lower_quantile, upper_quantile):
        

            # not sampled, actuall run
            # in the large bootstrap this value should be close to the median
            not_sampled_run = df.loc[df['nbootstrap'] == 0]
            actual_value = not_sampled_run[column].unique()[0]

            # only bootstrap
            data = df.loc[df['nbootstrap'] != 0]

            # Drop any NaNs just in case
            data = data[column].dropna().values

            # residues
            residues = data - actual_value

            # Calculate median, lowe quantile, and upper quantile
            median = actual_value + np.median(residues)
            lower_quantile = actual_value + np.percentile(residues, lower_quantile)
            upper_quantile = actual_value + np.percentile(residues, upper_quantile)

            return round(lower_quantile, 3), round(median, 3), round(upper_quantile, 3)


        # checkpoint
        if Path(self.pyseer_hits_version_table_all).exists():
            print('Pyseer hits with statistics (CI) already computed... remove existing file to run again... loading existing file.')
            
            # load
            return pd.read_csv(self.pyseer_hits_version_table_all, sep='\t')
        else: pass


        dfs = []
        for (pc, locus, mode), group in pyseer_hits.groupby(['PC', 'locus', 'mode']):
        
            # precision & recall    
            lower_quantile, upper_quantile = 5, 95
            columns = ['precision', 'recall', 'precision_SC', 'recall_SC']
            quantile_labels = [f'CI{lower_quantile}_LOWER', 'MEDIAN', f'CI{upper_quantile}_UPPER']

            # precision | recall for main run
            not_sampled_run = group.loc[group['nbootstrap'] == 0]

            for column in columns:
                CI1, CI2, CI3 = get_confidence_intervals(group, column=column, lower_quantile=lower_quantile, upper_quantile=upper_quantile)
                for ci, label in zip([CI1, CI2, CI3], quantile_labels):
                    not_sampled_run[f'{column}_{label}'] = ci

            dfs.append(not_sampled_run)

        # concat
        hits_df = pd.concat(dfs)

        # sort
        hits_df = hits_df.sort_values(['mode', 'locus', 'precision'], ascending=[False, True, False]).reset_index(drop=True)

        # save
        hits_df.to_csv(self.pyseer_hits_version_table_all, sep='\t', index=False)

        return hits_df


    def aggregate_pyseer_hits(self):

        def glob_concat_versions_tables(versions_dir, pattern, outfile):
            
            # load
            dfs = []
            for file in versions_dir.rglob(pattern):
                df = pd.read_csv(file, sep='\t')
                dfs.append(df)

            # concatenate
            concat_dfs = pd.concat(dfs)

            # clean
            cols2drop = ['nbootstrap', 'PC_locus_mode']
            for col in cols2drop:
                try: concat_dfs = concat_dfs.drop(col, axis=1)
                except: continue

            # save
            concat_dfs.to_csv(outfile, sep='\t', index=False)


        #########################
        #### ALL PYSEER HITS ####
        #########################

        versions_dir = self.pyseer_hits_version_table_all.parent.parent
        pattern = '*pyseer_hits_all.tsv*'
        outfile = self.pyseer_hits_final_table_all
        glob_concat_versions_tables(versions_dir, pattern, outfile)


        #######################################
        #### SIGNIFICANT POSITIVE NON-ZERO ####
        #######################################

        versions_dir = self.pyseer_hits_version_table_all.parent.parent
        pattern = '*pyseer_hits_filtered.tsv*'
        outfile = self.pyseer_hits_final_table_filtered
        glob_concat_versions_tables(versions_dir, pattern, outfile)


    def aggregate_function_prediction(self):

        def glob_concat_function_tables(versions_dir, pattern, outfile):

            # load
            dfs = []
            for file in versions_dir.rglob(pattern):
                df = pd.read_csv(file, sep='\t')
                dfs.append(df)

            # concatenate
            concat_dfs = pd.concat(dfs)

            # save
            concat_dfs.to_csv(outfile, sep='\t', index=False)

        ############################
        #### ALL FUNCTIONS HITS ####
        ############################

        versions_dir = Path(self.clusters_functions).parent.parent
        pattern = '*functions.tsv*'
        outfile = self.clusters_functions_best_all
        glob_concat_function_tables(versions_dir, pattern, outfile)

        

def _safe_division(numerator, denominator):
    """
    Performs safe division and returns np.nan if the denominator is zero.

    :param numerator: The numerator in the division.
    :param denominator: The denominator in the division.
    :return: The result of the division or np.nan if the denominator is zero.
    """
    if denominator == 0:
        return np.nan
    else:
        return numerator / denominator


def _safe_round(value, digits):
    """
    Rounds the value to the specified number of digits, handling np.nan.

    :param value: The value to round.
    :param digits: Number of decimal places.
    :return: Rounded value or np.nan if the input is np.nan.
    """
    if np.isnan(value):
        return np.nan
    else:
        return round(value, digits)


def _get_report_topology(ecod_values, report_priority, default='-'):
    """
    Given a list of ECOD topology strings and a priority list,
    return the reported topology based on the first match in priority order.
    """
    
    # ecod_values is a list of 5 strings (e.g., ECOD1 ... ECOD5)
    for _, search_str, reported_value in report_priority:
        # Check if any non-filler entry contains the search string.
        # (Here we use a substring search; adjust if you need an exact match.)
        if any(search_str in ecod for ecod in ecod_values if ecod != '-'):
            return reported_value
    
    return default       








