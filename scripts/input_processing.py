import pandas as pd
from Bio import SeqIO
import re
import hashlib
import os
import json
import phylotreelib as pt
from scripts.folder_structure import FolderStructure
from pathlib import Path
from tqdm import tqdm
import numpy as np



class InputProcessor:
    def __init__(self, structure, params):
        """
        Initialize the InputProcessor with input file paths, probability threshold, and minimum protein length.
        
        :param structure: Dictionary containing paths to input files (e.g., search.tsv, PCs2proteins.tsv, etc.)
        :param prob_threshold: Probability threshold for filtering ECOD predictions (default = 70)
        :param min_protein_length: Minimum protein length for filtering proteins (default = 300)
        """

        # params
        self.min_protein_length = params['input']['min_protein_length']
        self.min_ncontigs = params['input']['min_ncontigs']
        self.retain_species = params['input']['retain_species']
        self.recombinant_depos_colors = params['colors_and_shapes']['recombinant_depos_colors']

        self.ecod_prob_threshold = params['input']['ecod_prob_threshold']
        self.phrogs_prob_threshold = params['input']['phrogs_prob_threshold']
        self.pfam_prob_threshold = params['input']['pfam_prob_threshold']
        self.alan_prob_threshold = params['input']['alan_prob_threshold']
        self.report_topology_map = params['report_topology_map']


        ### paths

        # prophage proteins (sequence)
        self.input_proteins_files = structure['input']['proteins_files']
        self.processed_proteins_file = structure['input_processed']['proteins_file']
        self.processed_proteins_length = structure['input_processed']['protein_length_tsv']
        self.protein_sequences_file = structure['input_processed']['protein_sequences_tsv']


        # protein function prediction (table)
        self.search_tsv = structure['input']['search_tsv']
        self.pcs2proteins_tsv = structure['input']['pcs2proteins_tsv']
        self.protein_length_tsv = structure['input_processed']['protein_length_tsv']
        self.pc80_map_tsv = structure['input_processed']['pc80_map_tsv']
        self.pc80_functions_tsv = structure['input_processed']['pc80_functions_tsv']
        self.pc80_functions_best_tsv = structure['input_processed']['pc80_functions_best_tsv']


        # bacteria (table)
        self.bacteria_tsv = structure['input']['bacteria_tsv']
        self.processed_bacteria_tsv = structure['input_processed']['bacteria_tsv']

        # bacteria (tree)
        self.tree_file = structure['input']['bacteria_tree']
        self.processed_tree_file = structure['input_processed']['processed_bacteria_tree_nwk']

        # prophages (table)
        self.prophages_tsv = structure['input']['prophages_tsv']
        self.processed_prophages_tsv = structure['input_processed']['prophages_tsv']


        # experimental data
        self.lysogenic_table = structure['input']['lysogenic_table']
        self.lytic_table = structure['input']['lytic_table']
        self.recombinant_depos_table = structure['input_processed']['recombinant_depos_table']


    def process_bacteria_table(self):

        # paths
        bacteria_tsv = self.bacteria_tsv
        processed_bacteria_tsv = self.processed_bacteria_tsv

        # load
        bacteria_df = pd.read_csv(bacteria_tsv, sep='\t')

        # prophage freq
        prophages_freq = self._bacteria_prophage_frequency()

        # merge
        bacteria_df = bacteria_df.merge(prophages_freq, on='genomeID', how='outer')

        # format
        bacteria_df['nprophages'] = bacteria_df['nprophages'].fillna(0)
        bacteria_df['nprophages'] = bacteria_df['nprophages'].astype(int)

        # filter
        too_many_ncontigs = (bacteria_df['ncontigs_assembly_stats'] >= self.min_ncontigs)
        filt_ksc = (bacteria_df['species_abbreviation'].isin(self.retain_species))
        
        filt = ~too_many_ncontigs & filt_ksc
        bacteria_filt_df = bacteria_df.loc[filt]

        # save
        bacteria_filt_df.to_csv(processed_bacteria_tsv, sep='\t', index=False)

        # attribuite to clean prophages
        self.remove_bacteria_genomeIDs = bacteria_df.loc[~filt, 'genomeID'].to_list()

        print('\nProcessed bacteria table saved.')
        print(f'Bacterial genomes (n={len(bacteria_df)}): {len(bacteria_filt_df)} retained (-{(~filt).sum()}, ', end='')
        print(f'where {too_many_ncontigs.sum()} have more than {self.min_ncontigs} contigs and ', end='')
        print(f'{(~filt_ksc).sum()} not in {" ".join(self.retain_species)})\n')
    

    def process_bacteria_iqtree(self, outgroup = '398KBV'):

        # paths
        bacteria_df = pd.read_csv(self.processed_bacteria_tsv, sep='\t')
        genomeIDs = set(bacteria_df['genomeID'].to_list() + [outgroup])

        # create new leaves labels
        bacteria_df['new_leaves_labels'] = bacteria_df.apply(lambda row: row['MGG_SC'] + '_' + row['MGG_K_locus'], axis=1)
        leaves_mapping = dict(zip(bacteria_df['genomeID'], bacteria_df['new_leaves_labels']))

        # load
        with pt.Newicktreefile(self.tree_file) as treefile:
            tree = treefile.readtree()

        # all leaves
        all_leaves = tree.leaflist()

        # remove some
        leaves_to_remove = set(all_leaves) - genomeIDs
        tree.remove_leaves(leaves_to_remove)

        # # new labels
        # for leaf in tree.leaflist():
        #     if leaf in leaves_mapping:
        #         leaf = leaves_mapping[leaf]

        # re-root the tree
        tree.rootout(outgroup)

        # save
        with open(self.processed_tree_file, 'w') as f:
            f.write(tree.newick())


    def process_prophage_table(self):

        # paths
        prophages_tsv = self.prophages_tsv
        processed_prophages_tsv = self.processed_prophages_tsv

        # load
        prophages_df = pd.read_csv(prophages_tsv, sep='\t')

        # filter prophages
        remove_prophages = prophages_df['genomeID'].isin(self.remove_bacteria_genomeIDs)
        prophages_filt_df = prophages_df.loc[~remove_prophages]
        
        # save
        prophages_filt_df.to_csv(processed_prophages_tsv, sep='\t', index=False)

        # attribiute to clean proteins
        self.remove_prophageIDs = prophages_df.loc[remove_prophages, 'prophageID'].to_list()

        print('Reading not filtered table of prophages!')
        print('Processed prophages table saved.')
        print(f'Prophage genomes (n={len(prophages_df)}): {len(prophages_filt_df)} retained ', end='')
        print(f'(-{remove_prophages.sum()} that were found in excluded genomes)\n')

    
    def _bacteria_prophage_frequency(self):
        """ get dataframe with frequency of prophages in bacterial genomes """

        prophages_df = pd.read_csv(self.prophages_tsv, sep='\t', usecols=[0])
        prophages_freq = prophages_df.groupby('genomeID').size().reset_index(drop=False).rename({0: 'nprophages'}, axis=1)

        return prophages_freq


    def process_prophage_proteins(self):
        """
        Process prophage proteins by:
        1. Concatenating multiple FASTA files into a single file.
        2. Filtering out proteins below the minimum length threshold.
        3. Writing protein ID and length to a TSV file.
        """

        # paths
        fasta_files = self.input_proteins_files
        processed_fasta = self.processed_proteins_file
        protein_length = self.processed_proteins_length
        protein_sequences = self.protein_sequences_file

        # filter files
        fasta_files = self._filter_prophage_protein_files(fasta_files)

        # create file with protein lengths
        n_prot_retain, n_prot_total = 0, 0
        added_proteinID = set()  # Use a set for faster lookups

        # Prepare to accumulate data before writing
        length_lines = []
        sequence_lines = []
        fasta_records = []

        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                
                # total
                n_prot_total += 1
                proteinID = record.description.split('|||')[1].strip()
                record.description = ''
                record.id = proteinID
                seq = record.seq
                seq_len = len(seq)

                # write to fasta and table
                if seq_len >= self.min_protein_length and proteinID not in added_proteinID:
                    
                    # Collect data
                    fasta_records.append(record)
                    length_lines.append(f"{proteinID}\t{seq_len}\n")
                    sequence_lines.append(f"{proteinID}\t{seq}\n")
                    
                    # Mark protein as added
                    added_proteinID.add(proteinID)
                    
                    # count
                    n_prot_retain += 1

        # Write all data at once to avoid multiple IO operations
        with open(processed_fasta, 'w') as fasta_outfile, open(protein_length, 'w') as length_file, open(protein_sequences, 'w') as sequences_file:
            length_file.write("proteinID\tlength_aa\n")
            sequences_file.write("proteinID\tseq\n")
            
            # Write the collected data to the files
            SeqIO.write(fasta_records, fasta_outfile, "fasta")
            length_file.writelines(length_lines)
            sequences_file.writelines(sequence_lines)

        print(f"Prophage proteins processed and saved. [proteins: {n_prot_retain} with length >= {self.min_protein_length}] [removed: {n_prot_total - n_prot_retain}]\n")
        print("Duplicated genes at the end of phages are removed here! A problem of mgg_annotation pipeline.")


    def _filter_prophage_protein_files(self, paths_list):
        
        clean_paths_list = []
        for path in paths_list:
            if Path(path).stem in self.remove_prophageIDs:
                continue
            else:
                clean_paths_list.append(path)

        ntotal = len(paths_list)
        nretained = len(clean_paths_list)
        nremoved = ntotal - nretained
        print("Reading already filtered prophages by confidence and completeness! ")
        print(f"Prophage protein files (n={ntotal}): retained {nretained} (-{nremoved} that were found in excluded genomes)")

        return clean_paths_list


    def process_function_predictions(self, run=True):
        """ Process HMM search/hhblits results """

        # checkpoint
        if run: pass
        else: return

        # paths
        search_file = self.search_tsv
        pcs2proteins_file = self.pcs2proteins_tsv
        protein_length = self.protein_length_tsv

        pc80_map_tsv = self.pc80_map_tsv
        pc80_functions_tsv = self.pc80_functions_tsv
        pc80_functions_best_tsv = self.pc80_functions_best_tsv

        # load
        search_df = pd.read_csv(search_file, sep='\t')
        pcs_df = pd.read_csv(pcs2proteins_file, sep='\t')
        protein_length_df = pd.read_csv(protein_length, sep='\t')
        pcs2proteins_df = pcs_df.merge(protein_length_df, on='proteinID', how='outer')

        # filter pc80 by length
        pcs2proteins_df = pcs2proteins_df.loc[~pcs2proteins_df['length_aa'].isna()].reset_index(drop=True)
        pcs2proteins_df = pcs2proteins_df.rename({'PC': 'PC80'}, axis=1)
        all_pcs80 = list(pcs2proteins_df['PC80'].unique())
        all_pcs80_df = pd.DataFrame({'PC80': all_pcs80}).dropna()

        # filter hhsearch results by probability
        filt_ecod = (search_df['db'] == "ECOD") & (search_df['prob'] >= self.ecod_prob_threshold)
        filt_phrogs = (search_df['db'] == "PHROGS") & (search_df['prob'] >= self.phrogs_prob_threshold)
        filt_pfam = (search_df['db'] == "PFAM") & (search_df['prob'] >= self.pfam_prob_threshold)
        filt_alan = (search_df['db'] == "ALANDB") & (search_df['prob'] >= self.alan_prob_threshold)
        filt_queries = search_df['query'].isin(all_pcs80)

        filt_search = filt_queries & (filt_ecod | filt_phrogs | filt_pfam | filt_alan)
        search_df = search_df.loc[filt_search].reset_index(drop=True)

        # clean
        search_df = search_df.rename({'query': 'PC80'}, axis=1)
        search_df['function'] = search_df.apply(_parse_db_annotations, axis=1)

        def _compute_reported_topology(sub_df):
            ecod_rows = sub_df[sub_df['db'] == 'ECOD']
            return _report_topology_for_group(ecod_rows, self.report_topology_map)


        reported_topologies_map = search_df.groupby('PC80').apply(_compute_reported_topology)
        search_df['reported_topology_PC80'] = search_df['PC80'].map(reported_topologies_map)

        # add missing pc80 (no hit)
        search_df = all_pcs80_df.merge(search_df, on='PC80', how='left')
        search_df[['db', 'function', 'reported_topology_PC80']] = \
        search_df[['db', 'function', 'reported_topology_PC80']].fillna('no hit')
    

        # sort
        search_df['PC-sort'] = search_df['PC80'].str.strip('PC').astype(int)
        search_df = search_df.sort_values(['PC-sort', 'db', 'function'], ascending=[True, False, False])
        search_df = search_df.drop('PC-sort', axis=1)

        # sort cols
        functions_cols = ['PC80', 'db', 'function', 'reported_topology_PC80', 'prob', 'bits']
        rest_cols = [col for col in list(search_df.columns) if col not in functions_cols]
        search_df = search_df[functions_cols + rest_cols]

        # report best hits for pc80
        best_hits_list = []
        for (pcid, db), group in search_df.groupby(['PC80', 'db']):
            if db  == 'PFAM': 
                for func, subgroup in group.groupby('function'):
                    best_hits_list.append(_get_max_bitscore_hit(subgroup, nrows=1))
                if pcid == 'PC0001':
                    print(best_hits_list)
            elif db == 'ALANDB':
                best_hits_list.append(_get_max_bitscore_hit(group, nrows=1))
            elif db == 'PHROGS':
                for func, subgroup in group.groupby('function'):
                    best_hits_list.append(_get_max_bitscore_hit(subgroup, nrows=1))
            elif db == 'ECOD':
                best_hits_list.append(_report_unique_x_levels(group))
            elif db == 'no hit':
                best_hits_list.append(group)
            else:
                raise ValueError(f'ERROR! DB NOT RECOGNIZED: {db}')

        best_hits_df = pd.concat(best_hits_list).reset_index(drop=True)

        # save
        pcs2proteins_df.to_csv(pc80_map_tsv, sep='\t', index=False)
        search_df.to_csv(pc80_functions_tsv, sep='\t', index=False)
        best_hits_df.to_csv(pc80_functions_best_tsv, sep='\t', index=False)


    def process_recombinant_depos(self, run=True):

        # checkpoint
        if run: pass
        else: return

        print('Processing recombinant depolymerases...', end='\t')
        # assign
        lysogenic_table = self.lysogenic_table
        lytic_table = self.lytic_table

        # rename
        rename_lysogenic_dict = {'protein seq': 'seq', 'expression level': 'expression', 'K_locus_specificity': 'specificity'}
        rename_lytic_dict = {'K_locus_specificity':'specificity', 'protein_seq': 'seq'}
        genscript_list = ['0367_12','0574_17','0391_11','0496_72','021_14','391_03','184_04','145_08','174_38','164_08']
        
        # cols
        cols = ['proteinID', 'report', 'group', 'expression', 'specificity', 'K_locus_host', 'seq']

        # load
        lysogenic_df = pd.read_excel(lysogenic_table, sheet_name='DEPOS_ZDK', usecols=[0,3,5,6,11]).rename(rename_lysogenic_dict, axis=1)
        lytic_df = pd.read_csv(lytic_table, sep='\t').rename(rename_lytic_dict, axis=1)

        # clean lysogenic
        lysogenic_df['specificity'] = lysogenic_df['specificity'].replace('K3/KL146', 'KL3/KL146')
        remove_na = ~((lysogenic_df['proteinID'].isna()) | (lysogenic_df['proteinID'] == 'GWAS:'))
        lysogenic_df = lysogenic_df.loc[remove_na].reset_index(drop=True)
        ugly_string = lysogenic_df['specificity'].str.contains('-')
        lysogenic_df.loc[ugly_string, 'specificity'] = 'NO_KL'

        ### lysogenic

        # clean
        clean_expression_dict = {'high': 'PRODUCED', 'medium': 'PRODUCED', 'low': 'PRODUCED', 'none': 'NOT_PRODUCED'}
        lysogenic_df['expression'] = lysogenic_df['expression'].replace(clean_expression_dict)

        # assign color
        is_produced = (lysogenic_df['expression'] == 'PRODUCED')
        is_active = ~(lysogenic_df['specificity'].str.contains('NO_KL'))
        is_genscript = lysogenic_df['proteinID'].isin(genscript_list)
        is_zdk = ~(lysogenic_df['proteinID'].isin(genscript_list))

        zdk_not_produced = is_zdk & ~(is_produced)
        zdk_produced_active = is_zdk & is_produced & is_active
        zdk_produced_inactive = is_zdk & is_produced & ~(is_active)

        genscript_not_produced = is_genscript & ~(is_produced)
        genscript_produced_active = is_genscript & is_produced & is_active
        genscript_produced_inactive = is_genscript & is_produced & ~(is_active)


        lysogenic_df.loc[zdk_not_produced, 'group'] = 'lysogenic_zdk_not_produced'
        lysogenic_df.loc[zdk_produced_active, 'group'] = 'lysogenic_zdk_produced_active'
        lysogenic_df.loc[zdk_produced_inactive, 'group'] = 'lysogenic_zdk_produced_inactive'

        lysogenic_df.loc[genscript_not_produced, 'group'] = 'lysogenic_genscript_not_produced'
        lysogenic_df.loc[genscript_produced_active, 'group'] = 'lysogenic_genscript_produced_active'
        lysogenic_df.loc[genscript_produced_inactive, 'group'] = 'lysogenic_genscript_produced_inactive'

        # report
        lysogenic_df['report'] = lysogenic_df['proteinID'] + '_' + lysogenic_df['specificity'] + '_' + lysogenic_df['expression'] + '_HOST_' + lysogenic_df['K_locus_host']

        ### clean lytic
        lytic_df['group'] = 'lytic'
        lytic_df['report'] = lytic_df['proteinID']

        # combine
        recombinant_depos_df = pd.concat([lysogenic_df, lytic_df])

        # clean
        recombinant_depos_df = recombinant_depos_df[cols]
        recombinant_depos_df.to_csv(self.recombinant_depos_table, sep='\t', index=False)
        print('Done!')



def _parse_db_annotations(row):
    
    db = row['db']

    if db == "ECOD": 
        return _parse_ecod(row)
    elif db == "PHROGS":
        return _parse_phrogs(row)
    elif db == "PFAM":
        return _parse_pfam(row)
    elif db == "ALANDB":
        return _parse_alandb(row)
    else:
        return "DB_NOT_RECOGNIZED"


def _parse_ecod(row, pattern_ecodID = r"\|\s*((?:\d+\.\d+\.\d+\.\d+)|(?:\d+\.\d+\.\d+))\s*\|"):


    # init
    name = row['name']
    ecodID, A_LEVEL, X_LEVEL, H_LEVEL, T_LEVEL, F_LEVEL, PROTEIN_LEVEL = ['EXTRACTION_ERROR'] * 7

    try:
        # initial split
        ecod_id_raw, rest = name.split(' | A: ')

        # ECOD ID
        match = re.search(pattern_ecodID, ecod_id_raw)
        ecodID = match.group(1)

        # ECOD levels
        A_LEVEL, rest = rest.split(', X: ')
        X_LEVEL, rest = rest.split(', H: ')
        H_LEVEL, rest = rest.split(', T: ')
        T_LEVEL, rest = rest.split(', F: ')
        F_LEVEL, PROTEIN_LEVEL = rest.split(' | ')
        PROTEIN_LEVEL = PROTEIN_LEVEL.strip('Protein: ')

    except: pass

    report_full = f"{ecodID} | A: {A_LEVEL} | X: {X_LEVEL} | H: {H_LEVEL} | T: {T_LEVEL} | F: {F_LEVEL} | PROTEIN: {PROTEIN_LEVEL}"
    report_some = f"{ecodID} | X: {X_LEVEL} | H: {H_LEVEL} | T: {T_LEVEL}"
    return report_some


def _parse_pfam(row):
    pfamID, func_short, func_long = row['name'].split(' ; ')
    function = f'{pfamID} | {func_long}'
    return function


def _parse_phrogs(row):
    return row['annot']


def _parse_alandb(row):
    if pd.notna(row['annot']):
        return row['annot']
    else: 
        return row['category']

     
def _report_topology_for_group(group, mapping):

    if group.empty:
        return "no hit"
    # Loop through mapping except the fallback
    for idx, search_str, reported_topo in mapping:
        # Check if any row in column 'function' contains the search string
        if group['function'].str.contains(search_str, case=False, na=False).any():
            return reported_topo
    # If no priority match was found but there are ECOD rows, return fallback
    return "other"


def _get_max_bitscore_hit(group, nrows=1):
    """Sort the group by 'bits' in descending order and return the top nrows."""
    sorted_group = group.sort_values('bits', ascending=False)
    return sorted_group.head(nrows)


def _extract_x_level(report_full):
    """
    Given a function string in the format:
    "{ecodID} | A: {A_LEVEL} | X: {X_LEVEL} | H: {H_LEVEL} | T: {T_LEVEL} | F: {F_LEVEL} | PROTEIN: {PROTEIN_LEVEL}"
    extract and return the X_LEVEL value.
    """
    parts = report_full.split('|')
    for part in parts:
        part = part.strip()
        if part.startswith("X:"):
            # Split on "X:" and strip any whitespace to get the X_LEVEL value.
            return part.split("X:")[1].strip()
    return None


def _report_unique_x_levels(group):
    """
    For a group with db 'ECOD', extract the X_LEVEL from the function column,
    group by that level, and return the row with the highest bitscore for each X_LEVEL.
    """
    # Create a copy to avoid SettingWithCopyWarning.
    group = group.copy()
    
    # Extract X_LEVEL into a new column.
    group['X_level'] = group['function'].apply(_extract_x_level)
    
    best_hits = []
    # Group by the extracted X_level and get the best hit in each subgroup.
    for x_level, subgroup in group.groupby('X_level'):
        best_hits.append(_get_max_bitscore_hit(subgroup, nrows=1))
    # Combine the best hits for each unique X_level into one DataFrame.
    return pd.concat(best_hits).drop(columns='X_level', errors='ignore')
