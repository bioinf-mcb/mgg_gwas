import os
import shutil
from pathlib import Path

class FolderStructure:
    def __init__(self, input_dir, output_dir):
        
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)


    def get_params_dict(self):

        # input params
        input_processing_dict = {'ecod_prob_threshold': 0.7,
                                'phrogs_prob_threshold': 0.9,
                                'pfam_prob_threshold': 0.7,
                                'alan_prob_threshold': 0.7,
                                'annotation_tool': 'hhsearch', # DO NOT CHANGE. Another, unused option: hhblits
                                'min_protein_length': 300,
                                'min_ncontigs': 500,
                                'retain_species': ['KPN', 'KVV', 'KQQ', 'KQS']}

        # mmseqs params
        mmseqs_dict = {'mmseqs_binary': "/opt/homebrew/bin/mmseqs",
                    'min_identity': 0.8,
                    'min_coverage': 0.8,
                    'max_evalue': 10**-3,
                    'sensitivity': 7.5,
                    'min_pc_freq': 0.001, # equivalent of min_af
                    'max_pc_freq': 0.70,  # equivalent of max_af
                    'min_pc_freq_in_SCs': 2}  # min number of SCs for PC

        min_identity = (str(int(float(mmseqs_dict['min_identity']) * 100))).zfill(2)
        min_coverage = (str(int(float(mmseqs_dict['min_coverage']) * 100))).zfill(2)
        mmseqs_dict['version'] = f'PCI{min_identity}C{min_coverage}'

        report_topology_map = [
            (0, 'Pectin lyase-like', 'Pectin lyase-like'),
            (1, 'SGNH hydrolase', 'SGNH hydrolase'),
            (2, 'Alanine racemase-C', 'Alanine racemase-C'),
            (3, 'Intramolecular chaperone', 'Intramolecular chaperone'),
            (4, 'Concanavalin A-like', 'Concanavalin A-like'),
            (5, 'bladed', 'x-bladed'),
            (6, 'spike', 'tail spike')
            ]


        # gwas
        elastic_net_dict = {'correction': True, 'alpha': 0.0069}
        lasso_dict = {'correction': True, 'alpha': 0.8}
        incorrect_dict = {'correction': False, 'alpha': 0.2}
    

        gwas_dict = {'phenotypes': ['KL1', 'KL2','KL3', 'KL6','KL7','KL8', \
                                    'KL10','KL11','KL13','KL14','KL15','KL16','KL21','KL22', \
                                    'KL24','KL25','KL28','KL30','KL35','KL38','KL39','KL47', \
                                    'KL48','KL55','KL57','KL60','KL61','KL62','KL64','KL103', \
                                    'KL111','KL122','KL125','KL127','KL142'],
                     'bootstrap_sample_size': 2478,
                     'nbootstrap': 20, 
                     'cor_filter': 0.25,
                     'elastic_net': elastic_net_dict,
                     'lasso': lasso_dict,
                     'incorrect': incorrect_dict}

        colors_and_shapes_dict = {'ecod_colors': {'Pectin lyase-like': '#0E470E',
                                        'SGNH hydrolase': '#FFD700',
                                        'Concanavalin A-like': '#0000FF',
                                        'x-bladed': '#E63900',
                                        'Alanine racemase-C': '#8A2BE2',
                                        'Intramolecular chaperone': '#4682B4', 
                                        'tail spike': '#00ffff',
                                        'other': '#A9A4A4', 
                                        'no hit': '#BDBDF8'
                                        },
                                          
                        'recombinant_depos_colors': {'lysogenic_zdk_produced_active': '#EF3B2C', 
                                                    'lysogenic_zdk_not_produced': '#BDBDBD',
                                                    'lysogenic_zdk_produced_inactive': '#FC9272',
                                                    'lysogenic_genscript_produced_active': '#0570B0',
                                                    'lysogenic_genscript_not_produced': '#A6BDDB', 
                                                    'lytic':'#DF65B0'},
                        'collections_colors': {'Common SCs': '#3CB371', 
                                              'KASPAH unique SCs': '#ADD8E6', 
                                              'KLEBPAVIA unique SCs': '#DDA0DD',
                                              'KASPAH99': '#D22B2B'},
                        'min_precision': 0.5, 
                        'min_recall': 0.5,
                        'high_precision': 0.8
                        }

        # change pyseer alpha to 0.6 to avoid error (lasso)
        manual_alpha_adjustment = {'PCI00C50': ['KL11', 'KL15', 'KL48', 'KL57'],
                                   'PCI00C80': ['KL2', 'KL28', 'KL48', 'KL57'],
                                   'PCI50C50': ['KL48'],
                                   'PCI50C80': ['KL28', 'KL48']}



        # params
        params_dict = {}
        params_dict['input']  = input_processing_dict
        params_dict['mmseqs'] = mmseqs_dict
        params_dict['gwas'] = gwas_dict
        params_dict['colors_and_shapes'] = colors_and_shapes_dict
        params_dict['manual_alpha_adjustment'] = manual_alpha_adjustment
        params_dict['report_topology_map'] = report_topology_map

        return params_dict



    def create_structure(self):
        """
        Create a folder structure based on hardcoded directory definitions.
        Returns a dictionary containing paths to various directories and files.
        """

        ### versioning

        # params
        params = self.get_params_dict()

        min_protein_length = str(params['input']['min_protein_length'])
        species = '-'.join(params['input']['retain_species'])
        min_ncontigs = str(params['input']['min_ncontigs'])
        min_identity, min_coverage = (str(int(float(params['mmseqs']['min_identity']) * 100))).zfill(2), (str(int(float(params['mmseqs']['min_coverage']) * 100))).zfill(2)
        annotation_tool = str(params['input']['annotation_tool']).upper()
        nboostrap = str(params['gwas']['nbootstrap'])

        input_processed_version = f'PROTMINLEN{min_protein_length}_SPECIES-{species}_MINNCONTIGS{min_ncontigs}'
        mmseqs_version = f'PCI{min_identity}C{min_coverage}'
        pyseer_version = f'{input_processed_version}_{mmseqs_version}_BTSP{nboostrap}X'
        output_processed_version = f'{pyseer_version}_{annotation_tool}'


        ### define folder structure

        # main directories
        input_dir = self.input_dir
        output_dir = self.output_dir

        # INPUT DATA
        bacteria_dir = Path(input_dir, '1_BACTERIA')
        prophages_dir = Path(input_dir, '2_PROPHAGES')
        gwas_dir = Path(input_dir, '3_GWAS')

        bacteria_tsv = Path(bacteria_dir, 'bacteria_metadata.tsv')
        bacteria_tree_nwk = Path(bacteria_dir, 'bacteria_iqtree.nwk')
        prophages_tsv = Path(prophages_dir, 'prophages_metadata.tsv')
        pcs2proteins_tsv = Path(prophages_dir, 'pcs2proteins.tsv')
        search_tsv = Path(prophages_dir, 'raw_hhsuite.tsv')
        proteins_dir = Path(prophages_dir, '4_FASTA_CDS_AA')
        proteins_files = list(Path(proteins_dir).glob('*.faa'))

        lytic_table = Path(gwas_dir, 'DEPOLYMERASES_RECOMBINANT/LITERATURE_SEARCH/2025-01-22_LITERATURE_SEARCH.tsv')
        lysogenic_table = Path(gwas_dir, 'DEPOLYMERASES_RECOMBINANT/PROPHAGES/2025-03-21_PROPHAGE_DEPOLYMERASES.xlsx')
        predictions_table = Path(gwas_dir, 'DEPOLYMERASES_PREDICTIONS/BEST_DEPO_HITS.xlsx')


        # dict
        input_dict = {
            "input_files_dir": input_dir,
            "bacteria_tsv": bacteria_tsv,
            "bacteria_tree": bacteria_tree_nwk,
            "prophages_tsv": prophages_tsv,
            "search_tsv": search_tsv,
            'pcs2proteins_tsv': pcs2proteins_tsv,
            "proteins_dir": proteins_dir,
            "proteins_files": proteins_files,
            "lytic_table": lytic_table,
            "lysogenic_table": lysogenic_table,
            "predictions_table": predictions_table
        }


        ## intermediate
        intermediate_dir = Path(output_dir, '1_INTERMEDIATE')

        # processed
        processed_dir = Path(intermediate_dir, '1_PROCESSED_INPUT', input_processed_version)
        processed_proteins_dir = Path(processed_dir, '1_PROTEINS')
        processed_functions_dir = Path(processed_dir, '2_FUNCTIONS')

        bacteria_tsv = Path(processed_dir, 'bacteria.tsv')
        prophages_tsv = Path(processed_dir, 'prophages.tsv')
        processed_bacteria_tree_nwk = Path(processed_dir, 'bacteria_iqtree.nwk')

        proteins_file = Path(processed_proteins_dir, 'prophage_proteins.fasta')
        protein_length_tsv = Path(processed_proteins_dir, 'protein_length.tsv')
        protein_sequences_tsv = Path(processed_proteins_dir, 'protein_seq.tsv')
        recombinant_depos_table = Path(processed_proteins_dir, 'enzymes.tsv')
        pc80_map_tsv = Path(processed_functions_dir, 'pc80_map.tsv')
        pc80_functions_tsv  = Path(processed_functions_dir, 'pc80_functions.tsv')
        pc80_functions_best_tsv  = Path(processed_functions_dir, 'pc80_functions_best.tsv')

        # dict
        input_processed_dict = {
                    "processed_dir": processed_dir,
                    "processed_proteins_dir": processed_proteins_dir,
                    "processed_functions_dir": processed_functions_dir,
                    "bacteria_tsv": bacteria_tsv,
                    "prophages_tsv": prophages_tsv,
                    "processed_bacteria_tree_nwk": processed_bacteria_tree_nwk,
                    "proteins_file": proteins_file,
                    "protein_length_tsv": protein_length_tsv,
                    "protein_sequences_tsv": protein_sequences_tsv,
                    "recombinant_depos_table": recombinant_depos_table,
                    "pc80_map_tsv": pc80_map_tsv,
                    "pc80_functions_tsv": pc80_functions_tsv,
                    "pc80_functions_best_tsv": pc80_functions_best_tsv
                    }
        

        # mmseqs (clustering)
        mmseqs_dir = Path(intermediate_dir, '2_MMSEQS', mmseqs_version)
        raw_clusters = Path(mmseqs_dir, f'1_raw_clusters.tsv')
        clusters = Path(mmseqs_dir, f'2_clusters.tsv')
        clusters_matrix = Path(mmseqs_dir, f'3_binary_matrix.tsv')
        alignments_dir = Path(mmseqs_dir, 'alignments')

        # dict
        mmseqs_dict = {
                    "mmseqs_dir": mmseqs_dir,
                    "raw_clusters": raw_clusters,
                    "clusters": clusters,
                    "clusters_matrix": clusters_matrix,
                    "alignments_dir": alignments_dir
        }


        # functions
        functions_dir = Path(intermediate_dir, '3_FUNCTIONS')        

        versions_functions_dir = Path(functions_dir, annotation_tool, mmseqs_version)
        clusters_map = Path(versions_functions_dir, 'clusters_map.tsv')
        clusters_functions = Path(versions_functions_dir, 'clusters_functions.tsv')
        clusters_functions_best = Path(versions_functions_dir, 'clusters_functions_best.tsv')
        clusters_functions_best_all = Path(functions_dir, 'clusters_functions_best_all.tsv')
        
        # dict
        functions_dict = {
                    "functions_dir": functions_dir,
                    "versions_functions_dir": versions_functions_dir,
                    "clusters_map": clusters_map,
                    "clusters_functions": clusters_functions,
                    "clusters_functions_best": clusters_functions_best,
                    "clusters_functions_best_all": clusters_functions_best_all
                    }



        ## gwas
        pyseer_dir = Path(output_dir, '2_PYSEER', pyseer_version)

        phenotypes_dir = Path(pyseer_dir, '1_PHENOTYPES')
        variants_dir   = Path(pyseer_dir, '2_VARIANTS')
        covariates_dir = Path(pyseer_dir, '3_COVARIATES')
        scripts_dir    = Path(pyseer_dir, '4_SCRIPTS')
        results_dir    = Path(pyseer_dir, '5_RESULTS')
        
        phenotypes_vectors_dir   = Path(phenotypes_dir, '1_PHENOTYPES')
        bootstrap_phenotypes_dir = Path(phenotypes_dir, '2_BOOTSTRAP')
        
        phenotypes_matrix = Path(phenotypes_dir, 'matrix.tsv')
        bootstrap_phenotypes_matrix = Path(phenotypes_dir, 'bootstrap_matrix.tsv')

        n_variants_tsv   = Path(variants_dir, 'n_variants.tsv')
        variants_tsv   = Path(variants_dir, 'variants.tsv')
        covariates_tsv = Path(covariates_dir, 'covariates.tsv')
        
        lasso_sh   = Path(scripts_dir, 'lasso.sh')
        elastic_net_sh   = Path(scripts_dir, 'elastic_net.sh')
        incorrect_sh = Path(scripts_dir, 'incorrect.sh')

        lasso_dir   = Path(results_dir, '1_LASSO')
        elastic_net_dir = Path(results_dir, '2_ELASTIC_NET')
        incorrect_dir = Path(results_dir, '3_INCORRECT')


        # dict
        gwas_dict = {
            "pyseer_dir": pyseer_dir,
            "phenotypes_dir": phenotypes_dir,
            "variants_dir"   : variants_dir,
            "covariates_dir" : covariates_dir,
            "scripts_dir": scripts_dir,
            "results_dir": results_dir,
            "phenotypes_vectors_dir": phenotypes_vectors_dir,
            "bootstrap_phenotypes_dir": bootstrap_phenotypes_dir,
            "phenotypes_matrix": phenotypes_matrix,
            "bootstrap_phenotypes_matrix": bootstrap_phenotypes_matrix,
            "n_variants_tsv": n_variants_tsv,
            "variants_tsv": variants_tsv,
            "covariates_tsv": covariates_tsv,
            "lasso_dir": lasso_dir,
            "elastic_net_dir": elastic_net_dir,
            "incorrect_dir": incorrect_dir,
            "lasso_sh": lasso_sh,
            "elastic_net_sh": elastic_net_sh,
            "incorrect_sh": incorrect_sh
        }


        ## postprocessing
        processed_dir = Path(output_dir, '3_PROCESSING')
        versions_dir = Path(processed_dir, 'VERSIONS', output_processed_version)
        intermediate_dir = Path(versions_dir, 'INTERMEDIATE')

        raw_pyseer_table = Path(intermediate_dir, '1_raw_pyseer.tsv')
        clean_pyseer_table = Path(intermediate_dir, '2_clean_pyseer.tsv')
        pcs_metrics_table = Path(intermediate_dir, '3_pcs_metrics.tsv')
        pyseer_bootstrap_table = Path(intermediate_dir, '4_pyseer_bootstrap.tsv')
        
        pyseer_hits_version_table_all = Path(versions_dir, 'pyseer_hits_all.tsv')
        pyseer_hits_version_table_filtered = Path(versions_dir, 'pyseer_hits_filtered.tsv')
        
        pyseer_hits_final_table_all = Path(processed_dir, 'pyseer_hits_all.tsv')
        pyseer_hits_final_table_filtered = Path(processed_dir, 'pyseer_hits_filtered.tsv')


        # dict
        processed_dict = {
            "processing_dir": processed_dir,
            "versions_dir": versions_dir,
            "intermediate_dir": intermediate_dir,
            "raw_pyseer_table": raw_pyseer_table,
            "clean_pyseer_table": clean_pyseer_table,
            "pcs_metrics_table": pcs_metrics_table,
            "pyseer_bootstrap_table": pyseer_bootstrap_table,
            "pyseer_hits_version_table_all": pyseer_hits_version_table_all,
            "pyseer_hits_version_table_filtered": pyseer_hits_version_table_filtered,
            "pyseer_hits_final_table_all": pyseer_hits_final_table_all,
            "pyseer_hits_final_table_filtered": pyseer_hits_final_table_filtered,
        }


        ## analze
        analyze_dir = Path(output_dir, '4_ANALYZE')        
        per_locus_dir = Path(analyze_dir, '1_PER_LOCUS')
        aggregated_data_dir = Path(analyze_dir, '2_AGGREGATED_DATA')
        predictions_and_enzymes_table = Path(aggregated_data_dir, 'predictions_and_enzymes.tsv')
        tmp_blast_dir = '/Users/januszkoszucki/tmp_blast_dir'
        

        # analyze dict
        analyze_dict = {"analyze_dir": analyze_dir,
                        "per_locus_dir": per_locus_dir,
                        "aggregated_data_dir": aggregated_data_dir,
                        "predictions_and_enzymes_table": predictions_and_enzymes_table, 
                        "tmp_blast_dir": tmp_blast_dir,
                        }

        ## figures
        figures_dir = Path(output_dir, '5_FIGURES')
        figure4_dir = Path(figures_dir, 'FIG4')
        figure4B_dir = Path(figures_dir, 'FIG4', 'FIG4B')
        tmp_blast_dir = '/Users/januszkoszucki/tmp_blast_dir'

        FIG4A = Path(figure4_dir, 'FIG4A.pdf')
        FIG4B_NODES = Path(figure4B_dir, 'FIG4B_NODES.tsv')
        FIG4B_EDGES = Path(figure4B_dir, 'FIG4B_EDGES.tsv')
        FIG4C = Path(figure4_dir, 'FIG4C.pdf')
        
        # figures dict 
        figures_dict = {"figures_dir": figures_dir,
                        "figure4_dir": figure4_dir,
                        "figure4B_dir": figure4B_dir,
                        "FIG4A": FIG4A,
                        "FIG4B_NODES": FIG4B_NODES,
                        "FIG4B_EDGES": FIG4B_EDGES,
                        "FIG4C": FIG4C,
                        "tmp_blast_dir": tmp_blast_dir
                        } 

        #### STRUCTURE ####

        structure = {
            "input_processed": input_processed_dict,
            "mmseqs": mmseqs_dict,
            "functions": functions_dict,
            "gwas": gwas_dict,
            "processed": processed_dict,
            "analyze": analyze_dict,
            "figures": figures_dict
            }


        # clean tmp folder
        self._remove_tmp_dir(tmp_blast_dir)

        # create folders
        self._create_folders(structure)

        # append input folder
        structure["input"] = input_dict

        return structure


    def _create_folders(self, structure):
        """
        Private method to recursively create directories from a dictionary.
        Folders are created based on keys ending in '_dir'.
        """
        for key, value in structure.items():
            if isinstance(value, dict):
                self._create_folders(value)  # Recursively handle nested dictionaries
            elif key.endswith('_dir'):
                # Create the directory if it does not exist
                if not Path(value).exists():
                    Path(value).mkdir(exist_ok=True, parents=True)
                    print(f"Created directory: {value}")

    def _remove_tmp_dir(self, tmp_dir):
        if Path(tmp_dir).exists():
            shutil.rmtree(tmp_dir)

    def get_paths_dict(self):
        """
        Method to return the folder structure without creating directories.
        Useful for loading paths in other classes.
        """
        return self.create_structure()

