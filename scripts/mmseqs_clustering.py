import subprocess
import shutil
import pandas as pd
from pathlib import Path

# alignments only
import concurrent.futures
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO

class MMseqsClustering:
    def __init__(self, input_fasta, mmseqs_binary, mmseqs_dir, raw_clusters, clusters, clusters_matrix, alignments_dir, min_identity, min_coverage, max_evalue, sensitivity):
        """
        Initialize MMseqsClustering.
        """

        # attribiutres
        self.input_fasta = input_fasta
        self.mmseqs_binary = mmseqs_binary
        self.mmseqs_dir = mmseqs_dir
        self.raw_clusters = raw_clusters
        self.clusters = clusters
        self.clusters_matrix = clusters_matrix
        self.alignments_dir = alignments_dir
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.max_evalue = max_evalue
        self.sensitivity = sensitivity


    def _clean_mmseqs(self):
        """ remove the directory and all its contents """
        shutil.rmtree(self.mmseqs_dir)


    def run_mmseqs(self, tmp_dir='mmseqs_tmp_dir', run=True):
        """
        Run MMseqs2 clustering with a specified identity threshold and input fasta file.
        
        :param db_path: Temporary path for MMseqs database.
        :param min_identity: clustering identity threshold.
        :param min_coverage: clustering coverage threshold.
        """

        # checkpoint
        if run: pass
        else: return

        
        # clean previous results
        # print('Clean previous mmseqs run... ', end='')
        # self._clean_mmseqs()

        # checkpoint
        if Path(self.raw_clusters).exists():
            print('Mmseqs results exist... to run remove existing files. ')
            return None
        else: pass
        

        # paths
        tmp_dir = Path(self.mmseqs_dir, tmp_dir)


        # create
        Path(tmp_dir).mkdir(exist_ok=True, parents=True)

        # commands
        cd        = f'cd "{tmp_dir}";\n'
        createdb  = f'{self.mmseqs_binary} createdb "{self.input_fasta}" "DB";\n'
        cluster   = f'{self.mmseqs_binary} cluster --cluster-mode 2 -s {self.sensitivity} -c {self.min_coverage} --cov-mode 0 -e {self.max_evalue} --min-seq-id {self.min_identity} DB cluster tmp;\n'
        createtsv = f'{self.mmseqs_binary} createtsv DB DB cluster "{self.raw_clusters}";\n'

        # run
        cmd = f'{cd}{createdb}{cluster}{createtsv}'

        print("Running mmseqs... ", end='')
        out = subprocess.run(cmd, shell=True, capture_output=True)
        print("clusters saved.")


    def clean_clustering(self):
        """ process raw mmseqs results"""

        # load
        clusters_df = pd.read_csv(self.raw_clusters, sep='\t', header=None)
        clusters_df.columns = ['repr', 'proteinID']

        # rename PCs
        pcs = clusters_df.groupby('repr').size().sort_values(ascending=False).index
        npcs = len(pcs)
        width = len(str(npcs))
        pcs_dict = {pc: 'PC' + f'{i+1}'.zfill(width) for i, pc in enumerate(pcs)}

        # apply rename
        clusters_df['PC'] = clusters_df['repr'].map(pcs_dict)
        clusters_df[f'PC-sort'] = pd.to_numeric(clusters_df['PC'].str.strip('PC'), downcast='integer')
        clusters_df = clusters_df.sort_values(f'PC-sort', ascending=True)
        clusters_df = clusters_df.drop(f'PC-sort', axis=1)
        
        # save
        clusters_df.to_csv(self.clusters, sep='\t', index=False)

        # set att
        self.clusters_df = clusters_df
        print(f'Clean raw mmseqs protein clusters saved. [pcs: {npcs}]')


    def process_clustering(self, prophages_tsv):

        # load
        prophage2genome_df = pd.read_csv(prophages_tsv, sep='\t', usecols=[0,2])

        # prophage identifier
        self.clusters_df['prophageID'] = self.clusters_df['proteinID'].str.split('_PROTEIN', expand=True)[0]

        # map & update att
        self.clusters_df = self.clusters_df.merge(prophage2genome_df, on='prophageID', how='left')

        # save
        self.clusters_df.to_csv(self.clusters, sep='\t', index=False)
        print('Process clusters saved.')


    def compute_matrix(self, bacteria_tsv, matrix, min_pc_freq, max_pc_freq, min_pc_freq_in_SCs):
        """
        Compute presence/absence matrix for the clusters.
        Output matrix will show which sequences are present in each sample.

        Genomes with at least on prophage retained.
        PCs filtered based on their frequency.
        """

        # matrix
        matrix_df = self.clusters_df.pivot_table(index='genomeID', columns='PC', aggfunc='size', fill_value=0)
        bacteria_df = pd.read_csv(bacteria_tsv, sep='\t', usecols=[0,15,43])


        # append missing samples        
        with_prophage_genomeIDs = bacteria_df.query('nprophages != 0')['genomeID'].to_list()

        # report
        total_samples = len(bacteria_df)
        n_retained_samples = len(with_prophage_genomeIDs)
        no_prophages_samples = total_samples - n_retained_samples
        
        # filter
        matrix_df = matrix_df.reindex(with_prophage_genomeIDs, fill_value=0)

        # make binary
        matrix_df[matrix_df > 0] = 1

        # filter by SCs
        matrix_df = matrix_df.merge(bacteria_df[['genomeID', 'MGG_SC']], on='genomeID', how='left')

        # clean
        matrix_df.index = matrix_df['genomeID']
        matrix_df = matrix_df.drop('genomeID', axis=1)


        ### FILTER PCs

        # PC frequency in lineages (SCs)
        pcs = [col for col in matrix_df.columns if col.startswith('PC')]
        freq_sc_df = matrix_df.groupby('MGG_SC')[pcs].apply(lambda x: (x == 1).sum())
        freq_sc_df[freq_sc_df > 0] = 1
        freq_sc_df = freq_sc_df.sum().reset_index().rename({'index': 'PC', 0: 'n_SCs'}, axis=1)

        # clean
        matrix_df = matrix_df.drop('MGG_SC', axis=1)

        ### PC total frequency
        pc_freq = matrix_df.sum(axis=0).reset_index().rename({'index': 'PC', 0: 'pc_abundance'}, axis=1)
        pc_freq['pc_freq'] = pc_freq['pc_abundance'] / n_retained_samples
        pc_stats = pc_freq.merge(freq_sc_df, on='PC', how='outer')

        # retain variants
        filt_low_freq = (pc_stats['pc_freq'] <= min_pc_freq)
        filt_high_freq = (pc_stats['pc_freq'] >= max_pc_freq)
        filt_too_few_SCs = (pc_stats['n_SCs'] < min_pc_freq_in_SCs)
        filt_retain_variants = ~(filt_low_freq | filt_high_freq | filt_too_few_SCs)

        retain_variants = pc_stats.loc[filt_retain_variants, 'PC'].to_list()

        # report
        total_variants = matrix_df.shape[1]
        n_low_freq_variants = filt_low_freq.sum()
        n_high_freq_variants = filt_high_freq.sum()
        n_too_few_SCs = filt_too_few_SCs.sum()
        n_exclude_variants = n_low_freq_variants + n_high_freq_variants + n_too_few_SCs
        n_retained_variants = len(retain_variants)

        # apply filters
        matrix_df = matrix_df.loc[:, retain_variants]

        # transpose
        matrix_df.columns.name = 'genomeID'
        matrix_df = matrix_df.T
        
        # save
        matrix_df.to_csv(matrix, sep='\t')


        # print
        print(f'\nBinary matrix saved with {n_retained_samples} samples and {n_retained_variants} variants.')
        print(f'Samples (n={total_samples}): {n_retained_samples} retained (-{no_prophages_samples} wihout prophages)')
        print(f'Variants (n={total_variants}): {n_retained_variants} retained (-{n_exclude_variants}, where -{n_low_freq_variants} too rare (<={min_pc_freq}), -{n_high_freq_variants} too abundand (>={max_pc_freq}), -{n_too_few_SCs} in too little SCs (<{min_pc_freq_in_SCs})) \n')
        

    def get_alignments(self, max_workers=4):
        """
        Writes per-PC FASTA files in parallel.

        clusters_df     : Pandas DataFrame with columns 'PC' and 'proteinID'
        input_fasta     : Single multi-sequence FASTA file with many protein IDs
        alignments_dir  : Path to output directory
        max_workers     : Number of threads/processes to use (default=4)
        """

        # variables
        clusters_df, input_fasta, alignments_dir = self.clusters_df, self.input_fasta, self.alignments_dir

        # 1) Read the entire FASTA once into a dictionary: proteinID -> sequence
        seq_dict = {}
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_dict[record.id] = str(record.seq)

        # 2) Group the clusters by PC
        grouped = clusters_df.groupby("PC")

        # 3) For parallel execution, use ThreadPoolExecutor
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Keep track of futures
            futures = []

            for pc, group in grouped:
                futures.append(
                    executor.submit(
                        _get_alignment_for_pc,
                        pc,
                        group,
                        seq_dict,
                        alignments_dir
                    )
                )

            # 4) Optionally, track progress with tqdm
            for _ in tqdm(
                concurrent.futures.as_completed(futures), 
                total=len(futures),
                desc="?Alignments",
                ncols=50,
                ascii=True
            ):
                pass


    def map_functions(self, pc80_map_tsv, pc80_functions_best_tsv, clusters_map, clusters_functions, mmseqs_version, report_topology_map):

        # paths
        clusters = self.clusters

        pc80_map_tsv = pc80_map_tsv
        pc80_functions_best_tsv = pc80_functions_best_tsv
        clusters_map = clusters_map
        clusters_functions = clusters_functions

        # load
        mmseqs_df = pd.read_csv(clusters, sep='\t')
        pc80_map_df = pd.read_csv(pc80_map_tsv, sep='\t')
        pc80_functions_best_df = pd.read_csv(pc80_functions_best_tsv, sep='\t')

        # clean
        mmseqs_df = mmseqs_df.drop('repr', axis=1)
        pc80_map_df = pc80_map_df.drop('repr', axis=1)        
        mmseqs_df['version'] = mmseqs_version

        # map pc80 to GWAS clustering
        mapping_df = mmseqs_df.merge(pc80_map_df, on='proteinID', how='left')

        # save
        cols = ['PC', 'version', 'PC80', 'genomeID', 'prophageID', 'proteinID', 'length_aa']
        mapping_df = mapping_df[cols]
        mapping_df.to_csv(clusters_map, sep='\t', index=False)


        # retain only mapping of pc to pc80
        cols = ['PC', 'version', 'PC80']
        mapping_df = mapping_df[cols]

        mapping_df['PC-version-PC80'] = mapping_df['PC'] + '-' + mapping_df['version'] + '-' + mapping_df['PC80']
        mapping_df = mapping_df.drop_duplicates('PC-version-PC80')
        mapping_df = mapping_df.drop('PC-version-PC80', axis=1)

        # map functions
        mapping_df = mapping_df.merge(pc80_functions_best_df, on='PC80', how='left')


        # report topology for gwas clusters
        report_topology_order = [entry[2] for entry in sorted(report_topology_map, key=lambda element: element[0])]
        pc80_to_gwas_topologies_map = {}

        for pcid, group in mapping_df.groupby('PC'):
            unique_reported_topologies = group['reported_topology_PC80'].unique()
            
            # no hit
            if len(unique_reported_topologies) == 1 and unique_reported_topologies[0] == 'no hit':
                pc80_to_gwas_topologies_map[pcid] = 'no hit'
                continue
            
            # other
            if len(unique_reported_topologies) == 2 and set(unique_reported_topologies) == {'no hit', 'other'}:
                pc80_to_gwas_topologies_map[pcid] = 'other'
                continue

            # other
            if len(unique_reported_topologies) == 1 and unique_reported_topologies[0] == 'other':
                pc80_to_gwas_topologies_map[pcid] = 'other'
                continue

            
            # otherwise report topology
            for topology in report_topology_order:
                if topology in unique_reported_topologies:
                    pc80_to_gwas_topologies_map[pcid] = topology
                    break

        # apply
        mapping_df['reported_topology_PC'] = mapping_df['PC'].map(pc80_to_gwas_topologies_map)

        # save
        cols = ['PC','version','PC80','db','function','reported_topology_PC80','reported_topology_PC','prob','bits','target','pvalue','ident','qcov','tcov','qstart','qend','qlength','tstart','tend','tlength','evalue','name','color','annot','category','phrog/alan_profile']
        mapping_df = mapping_df[cols]
        mapping_df.to_csv(clusters_functions, sep='\t', index=False)
        


def _get_alignment_for_pc(pc, group, seq_dict, alignments_dir):
        """
        Create a FASTA file for the given PC and its corresponding cluster group.

        pc              : The PC identifier (string or int).
        group           : A subset of the dataframe rows for this PC.
        seq_dict        : Dict mapping proteinID -> sequence.
        alignments_dir  : Directory where PC FASTA files should be written.
        """
        fasta_file = Path(alignments_dir) / f"{pc}.fasta"
        
        # If file already exists, skip creation
        if fasta_file.exists():
            return
        
        # Build FASTA text
        fasta_parts = []
        for row in group.itertuples():
            seq = seq_dict.get(row.proteinID)
            if seq is not None:
                fasta_parts.append(f">{row.proteinID}\n{seq}\n")
        
        # Save FASTA text
        with open(fasta_file, "w") as f:
            f.write("".join(fasta_parts))

