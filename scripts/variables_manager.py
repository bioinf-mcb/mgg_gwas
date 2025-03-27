
class VariablesManager:
    def __init__(self, structure, params):
        
        self.structure = structure
        self.params = params


    def get_mmseqs_vars(self):
        
        # params
        mmseqs_binary = self.params['mmseqs']['mmseqs_binary']
        min_identity = self.params['mmseqs']['min_identity']
        min_coverage = self.params['mmseqs']['min_coverage']
        max_evalue = self.params['mmseqs']['max_evalue']
        sensitivity = self.params['mmseqs']['sensitivity']
        
        
        # paths
        input_fasta =  self.structure['input_processed']['proteins_file']
        mmseqs_dir =  self.structure['mmseqs']['mmseqs_dir']
        raw_clusters =  self.structure['mmseqs']['raw_clusters']
        clusters =  self.structure['mmseqs']['clusters']
        clusters_matrix =  self.structure['mmseqs']['clusters_matrix']
        alignments_dir =  self.structure['mmseqs']['alignments_dir']

        # list
        mmseqs_vars = [input_fasta, mmseqs_binary, mmseqs_dir, \
                       raw_clusters, clusters, clusters_matrix, alignments_dir, \
                       min_identity, min_coverage, max_evalue, sensitivity]

        return mmseqs_vars

    def get_process_clustering_vars(self):

        # paths
        prophages_tsv = self.structure['input_processed']['prophages_tsv']

        # list
        proc_clus_vars = [prophages_tsv]

        return proc_clus_vars


    def get_matrix_vars(self):
        

        # params
        min_pc_freq = self.params['mmseqs']['min_pc_freq']
        max_pc_freq = self.params['mmseqs']['max_pc_freq']
        min_pc_freq_in_SCs = self.params['mmseqs']['min_pc_freq_in_SCs']

        # paths
        bacteria_tsv = self.structure['input_processed']['bacteria_tsv']
        matrix = self.structure['mmseqs']['clusters_matrix']
        
        # list
        matrix_vars = [bacteria_tsv, matrix, min_pc_freq, max_pc_freq, min_pc_freq_in_SCs]

        return matrix_vars


    def get_map_functions_vars(self):

        # params
        mmseqs_version = self.params['mmseqs']['version']
        report_topology_map = self.params['report_topology_map']

        # paths
        pc80_map_tsv = self.structure['input_processed']['pc80_map_tsv']
        pc80_functions_best_tsv = self.structure['input_processed']['pc80_functions_best_tsv']

        clusters_map = self.structure['functions']['clusters_map']
        clusters_functions = self.structure['functions']['clusters_functions']

        # list
        map_functions_vars = [pc80_map_tsv, pc80_functions_best_tsv, clusters_map, clusters_functions, mmseqs_version, report_topology_map]

        return map_functions_vars
    


    def get_gwas_vars(self, mode):

        # select dict
        if mode == 'lasso': mode_dict = self.params['gwas']['lasso']
        elif mode == 'elastic_net': mode_dict = self.params['gwas']['elastic_net']
        elif mode == 'incorrect': mode_dict = self.params['gwas']['incorrect']
        else: 
            print('Error in get_gwas_vars!')
            exit()

        # params
        correction = mode_dict['correction']
        alpha = mode_dict['alpha']
        cor_filter = self.params['gwas']['cor_filter']
        min_af =  self.params['mmseqs']['min_pc_freq']
        max_af = self.params['mmseqs']['max_pc_freq']
    
        gwas_vars = [mode, correction, alpha, cor_filter, min_af, max_af]

        return gwas_vars

