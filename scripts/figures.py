
def get_k_loci(structure, params):

    import pandas as pd

    # Retrieve file paths and parameters.
    custom_gwas_kloci = params['gwas']['phenotypes']
    bacteria_table = structure['input_processed']['bacteria_tsv']
    predictions_and_enzymes = structure['analyze']['predictions_and_enzymes_table']
    output_pdf_path = structure['figures']['FIG4A']
    
    # Read input tables.
    bacteria_df = pd.read_csv(bacteria_table, sep='\t')
    predictions_and_enzymes_df = pd.read_csv(predictions_and_enzymes, sep='\t')


    # Filter rows based on the source.
    zdk_df = predictions_and_enzymes_df.query('source == "PROPHAGE_ZDKLAB"')
    lytic_df = predictions_and_enzymes_df.query('source == "LITERATURE_SEARCH"')
    predicted_df = predictions_and_enzymes_df.query('source == "PREDICTION"')

    zdk_specificity = zdk_df['specificity'].unique()
    zdk_host = zdk_df['K_locus_host'].unique()
    literature = lytic_df['specificity'].unique()

    k_loci_dict = {'gwas': custom_gwas_kloci,
                    'zdk_specificity': zdk_specificity,
                    'zdk_hosts': zdk_host,
                    'literature': literature}

    return k_loci_dict


def FIGURE4_PANELA(structure, params, merge_expr=True, cell_number_color="#000000", 
                   target_fig_height=3.5, target_fig_width=8, labels_fontsize=6, cells_fontsize=4, 
                   margin=0.3, header_space=0.1):
    """
    Generates FIG4 PANEL A visualization with an extra final row for LITERATURE_SEARCH DEPOLYMERASES.
    
    In the non-merged version, the rows are:
      0. KASPAH_REF isolates
      1. GWAS_KSC counts
      2. Tested expressed proteinsnota
      3. Tested not expressed proteins
      4. Predicted proteins
      5. LITERATURE_SEARCH DEPOLYMERASES
      
    In the merged version, the rows are:
      0. KASPAH_REF isolates
      1. GWAS_KSC counts
      2. Tested proteins (merged expressed)
      3. Predicted proteins
      4. LITERATURE_SEARCH DEPOLYMERASES
    """
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import re


    def generate_pdf_visualization(data, output_pdf_path, custom_gwas_kloci, merge_expr=True, 
                                   cell_number_color="#000000", target_fig_height=3.5, target_fig_width=8, 
                                   labels_fontsize=6, cells_fontsize=4, margin=0.3, header_space=0.1):

        def sort_key(key):
            # Extract the two-letter prefix (if available)
            prefix = key[:2]
            # Define the order: "KL" should come before "KN"; all others get a higher value.
            if prefix == "KL":
                prefix_order = 0
            elif prefix == "KN":
                prefix_order = 1
            else:
                prefix_order = 2
            # Extract the numeric portion from the key.
            match = re.search(r'\d+', key)
            num = int(match.group()) if match else float('inf')
            return (prefix_order, num)


        # Set number of rows and row labels based on merge_expr flag.
        if merge_expr:
            num_rows = 5
            row_labels = [
                "# kaspah ref isolates", 
                "# Kp complex SCs", 
                "Tested proteins", 
                "Predicted proteins",
                "Literature search proteins"
            ]
        else:
            num_rows = 6
            row_labels = [
                "# kaspah ref isolates", 
                "# Kp complex SCs", 
                "Tested expressed proteins", 
                "Tested not expressed proteins", 
                "Predicted proteins",
                "Literature search proteins"
            ]
        
        # Ensure every K locus from the custom list is in data.
        for kl in custom_gwas_kloci:
            if kl not in data:
                data[kl] = {
                    "KASPAH_REF": 0,
                    "GWAS_KSC": 0,
                    "EXPRESSED": {
                        "ACTIVE": 0,
                        "INACTIVE": 0,
                        "NOT_PRODUCED": 0
                    },
                    "PREDICTION": {},
                    "LITERATURE_SEARCH_DEPOLYMERASES": 0
                }
        
        # Sort keys (e.g. KL1, KL2, …)
        sorted_keys = sorted(data.keys(), key=sort_key)
        num_cols = len(sorted_keys)

        # Adjust cell_width based on target_fig_width.
        cell_width = (target_fig_width - 2 * margin) / num_cols
        cell_height = (target_fig_height - 2 * margin) / num_rows
        fig_width = target_fig_width
        fig_height = target_fig_height
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.set_xlim(0, fig_width)
        ax.set_ylim(0, fig_height)
        ax.axis('off')
        
        # Draw column headers (K locus labels).
        header_y = fig_height - margin + header_space
        for i, k_locus in enumerate(sorted_keys):
            x = margin + i * cell_width + cell_width / 2
            ax.text(x, header_y, k_locus, ha='center', va='bottom',
                    fontsize=labels_fontsize, fontweight='bold', rotation=60)
        
        # Define color mapping.
        color_mapping = {
            "KASPAH_REF": "#f06c6c",
            "ACTIVE": "#fa61bd",
            "INACTIVE": "#f5bfe0",
            "NOT_PRODUCED": "#D0D0D0",
            "GWAS_KSC": "#FFBD4D",
            "PREDICTION_PERFECT": "#85C187",
            "PREDICTION_GOOD": "#DBEFDC",
            "LITERATURE_SEARCH_DEPOLYMERASES": "#bc96d9"
        }
        
        for row in range(num_rows):
            for col, k_locus in enumerate(sorted_keys):
                x = margin + col * cell_width
                y = margin + (num_rows - 1 - row) * cell_height
                
                # Row 0: KASPAH_REF isolates.
                if row == 0:
                    count = data[k_locus].get("KASPAH_REF", 0)
                    cell_color = color_mapping["KASPAH_REF"] if count > 0 else "white"
                    rect = patches.Rectangle((x, y), cell_width, cell_height,
                                             facecolor=cell_color, edgecolor="black")
                    ax.add_patch(rect)
                    if count:
                        ax.text(x + cell_width/2, y + cell_height/2, str(count),
                                ha='center', va='center', color=cell_number_color,
                                fontsize=cells_fontsize, fontweight='bold')
                
                # Row 1: GWAS_KSC counts.
                elif row == 1:
                    count = data[k_locus].get("GWAS_KSC", 0)
                    cell_color = color_mapping["GWAS_KSC"] if (k_locus in custom_gwas_kloci and count > 0) else "white"
                    rect = patches.Rectangle((x, y), cell_width, cell_height,
                                             facecolor=cell_color, edgecolor="black")
                    ax.add_patch(rect)
                    if count:
                        ax.text(x + cell_width/2, y + cell_height/2, str(count),
                                ha='center', va='center', color=cell_number_color,
                                fontsize=cells_fontsize, fontweight='bold')
                
                # Expressed / Not Produced proteins row(s)
                if merge_expr:
                    # Merged version: row 2.
                    if row == 2:
                        expr = data[k_locus].get("EXPRESSED", {})
                        active_count = expr.get("ACTIVE", 0)
                        inactive_count = expr.get("INACTIVE", 0)
                        not_prod_count = expr.get("NOT_PRODUCED", 0)
                        num_categories = sum(1 for c in [active_count, inactive_count, not_prod_count] if c > 0)
                        if num_categories <= 1:
                            if active_count > 0:
                                cell_color = color_mapping["ACTIVE"]
                                count = active_count
                            elif inactive_count > 0:
                                cell_color = color_mapping["INACTIVE"]
                                count = inactive_count
                            elif not_prod_count > 0:
                                cell_color = color_mapping["NOT_PRODUCED"]
                                count = not_prod_count
                            else:
                                cell_color = "white"
                                count = 0
                            rect = patches.Rectangle((x, y), cell_width, cell_height,
                                                     facecolor=cell_color, edgecolor="black")
                            ax.add_patch(rect)
                            if count:
                                ax.text(x + cell_width/2, y + cell_height/2, str(count),
                                        ha='center', va='center', color=cell_number_color,
                                        fontsize=cells_fontsize, fontweight='bold')
                        else:
                            # More than one category present.
                            produced = []
                            if active_count > 0:
                                produced.append(("ACTIVE", active_count))
                            if inactive_count > 0:
                                produced.append(("INACTIVE", inactive_count))
                            np_present = (not_prod_count > 0)
                            num_subcells = len(produced) + (1 if np_present else 0)
                            subcell_width = cell_width / num_subcells
                            # Draw produced subcells on the left.
                            for i, (cat, count) in enumerate(produced):
                                sub_x = x + i * subcell_width
                                sub_color = color_mapping[cat]
                                sub_rect = patches.Rectangle((sub_x, y), subcell_width, cell_height,
                                                             facecolor=sub_color, edgecolor="black")
                                ax.add_patch(sub_rect)
                                ax.text(sub_x + subcell_width/2, y + cell_height/2, str(count),
                                        ha='center', va='center', color=cell_number_color,
                                        fontsize=cells_fontsize, fontweight='bold')
                            # Draw NOT_PRODUCED subcell on the right if present.
                            if np_present:
                                sub_x = x + len(produced) * subcell_width
                                sub_color = color_mapping["NOT_PRODUCED"]
                                sub_rect = patches.Rectangle((sub_x, y), subcell_width, cell_height,
                                                             facecolor=sub_color, edgecolor="black")
                                ax.add_patch(sub_rect)
                                ax.text(sub_x + subcell_width/2, y + cell_height/2, str(not_prod_count),
                                        ha='center', va='center', color=cell_number_color,
                                        fontsize=cells_fontsize, fontweight='bold')
                else:
                    # Not merged version: separate rows.
                    if row == 2:
                        # Row 2: Expressed proteins (only ACTIVE and INACTIVE).
                        expr = data[k_locus].get("EXPRESSED", {})
                        active_count = expr.get("ACTIVE", 0)
                        inactive_count = expr.get("INACTIVE", 0)
                        num_categories = sum(1 for c in [active_count, inactive_count] if c > 0)
                        if num_categories <= 1:
                            if active_count > 0:
                                cell_color = color_mapping["ACTIVE"]
                                count = active_count
                            elif inactive_count > 0:
                                cell_color = color_mapping["INACTIVE"]
                                count = inactive_count
                            else:
                                cell_color = "white"
                                count = 0
                            rect = patches.Rectangle((x, y), cell_width, cell_height,
                                                     facecolor=cell_color, edgecolor="black")
                            ax.add_patch(rect)
                            if count:
                                ax.text(x + cell_width/2, y + cell_height/2, str(count),
                                        ha='center', va='center', color=cell_number_color,
                                        fontsize=cells_fontsize, fontweight='bold')
                        else:
                            subcell_width = cell_width / 2
                            for i, (cat, count) in enumerate([("ACTIVE", active_count), ("INACTIVE", inactive_count)]):
                                if count:
                                    sub_x = x + i * subcell_width
                                    sub_color = color_mapping[cat]
                                    sub_rect = patches.Rectangle((sub_x, y), subcell_width, cell_height,
                                                                 facecolor=sub_color, edgecolor="black")
                                    ax.add_patch(sub_rect)
                                    ax.text(sub_x + subcell_width/2, y + cell_height/2, str(count),
                                            ha='center', va='center', color=cell_number_color,
                                            fontsize=cells_fontsize, fontweight='bold')
                    if not merge_expr and row == 3:
                        # Row 3: Not produced proteins.
                        expr = data[k_locus].get("EXPRESSED", {})
                        count = expr.get("NOT_PRODUCED", 0)
                        cell_color = color_mapping["NOT_PRODUCED"] if count > 0 else "white"
                        rect = patches.Rectangle((x, y), cell_width, cell_height,
                                                 facecolor=cell_color, edgecolor="black")
                        ax.add_patch(rect)
                        if count:
                            ax.text(x + cell_width/2, y + cell_height/2, str(count),
                                    ha='center', va='center', color=cell_number_color,
                                    fontsize=cells_fontsize, fontweight='bold')
                
                # Prediction row.
                if (merge_expr and row == 3) or (not merge_expr and row == 4):
                    pred_counts = data[k_locus].get("PREDICTION", {})
                    if not pred_counts:
                        rect = patches.Rectangle((x, y), cell_width, cell_height,
                                                 facecolor="white", edgecolor="black")
                        ax.add_patch(rect)
                    elif len(pred_counts) == 1:
                        pred_strength, count = list(pred_counts.items())[0]
                        mapped_category = "PREDICTION_PERFECT" if pred_strength.lower() == "yes" else "PREDICTION_GOOD"
                        cell_color = color_mapping[mapped_category]
                        rect = patches.Rectangle((x, y), cell_width, cell_height,
                                                 facecolor=cell_color, edgecolor="black")
                        ax.add_patch(rect)
                        ax.text(x + cell_width/2, y + cell_height/2, str(count),
                                ha='center', va='center', color=cell_number_color,
                                fontsize=cells_fontsize, fontweight='bold')
                    else:
                        n_categories = len(pred_counts)
                        subcell_width = cell_width / n_categories
                        sorted_pred = sorted(pred_counts.items(), key=lambda x: 0 if x[0].lower() == "yes" else 1)
                        for i, (pred_strength, count) in enumerate(sorted_pred):
                            mapped_category = "PREDICTION_PERFECT" if pred_strength.lower() == "yes" else "PREDICTION_GOOD"
                            sub_x = x + i * subcell_width
                            sub_rect = patches.Rectangle((sub_x, y), subcell_width, cell_height,
                                                         facecolor=color_mapping[mapped_category],
                                                         edgecolor="black")
                            ax.add_patch(sub_rect)
                            ax.text(sub_x + subcell_width/2, y + cell_height/2, str(count),
                                    ha='center', va='center', color=cell_number_color,
                                    fontsize=cells_fontsize, fontweight='bold')
                
                # Final row: LITERATURE_SEARCH DEPOLYMERASES.
                if (merge_expr and row == 4) or (not merge_expr and row == 5):
                    count = data[k_locus].get("LITERATURE_SEARCH_DEPOLYMERASES", 0)
                    cell_color = color_mapping["LITERATURE_SEARCH_DEPOLYMERASES"] if count > 0 else "white"
                    rect = patches.Rectangle((x, y), cell_width, cell_height,
                                             facecolor=cell_color, edgecolor="black")
                    ax.add_patch(rect)
                    if count:
                        ax.text(x + cell_width/2, y + cell_height/2, str(count),
                                ha='center', va='center', color=cell_number_color,
                                fontsize=cells_fontsize, fontweight='bold')
        
        # Draw row labels.
        for r in range(num_rows):
            y = margin + (num_rows - 1 - r) * cell_height + cell_height/2
            ax.text(margin - 0.2, y, row_labels[r], ha='right', va='center',
                    fontsize=labels_fontsize, fontweight='bold')
        
        plt.savefig(output_pdf_path, format='pdf', bbox_inches='tight')
        plt.close(fig)
    
    # Retrieve file paths and parameters.
    custom_gwas_kloci = set(params['gwas']['phenotypes'])
    bacteria_table = structure['input_processed']['bacteria_tsv']
    predictions_and_enzymes = structure['analyze']['predictions_and_enzymes_table']
    output_pdf_path = structure['figures']['FIG4A']
    
    # Read input tables.
    bacteria_df = pd.read_csv(bacteria_table, sep='\t')
    predictions_and_enzymes_df = pd.read_csv(predictions_and_enzymes, sep='\t')
    
    # Filter rows based on the source.
    zdk_df = predictions_and_enzymes_df.query('source == "PROPHAGE_ZDKLAB"')
    lytic_df = predictions_and_enzymes_df.query('source == "LITERATURE_SEARCH"')
    predicted_df = predictions_and_enzymes_df.query('source == "PREDICTION"')
    
    # --- Process dual specificity for predicted proteins ---
    if 'specificity' in predicted_df.columns:
        predicted_df = predicted_df.copy()
        predicted_df['assing_specificity'] = predicted_df.apply(
            lambda row: [s.strip() for s in row['specificity'].split('/')]
            if pd.notnull(row['specificity']) and "/" in row['specificity'] 
            else [row['assing_K_locus_host_when_specificity_missing']],
            axis=1)
        predicted_df = predicted_df.explode('assing_specificity')
    else:
        predicted_df['assing_specificity'] = predicted_df['assing_K_locus_host_when_specificity_missing']
    
    # --- Process dual specificity for literature search proteins ---
    if 'specificity' in lytic_df.columns:
        lytic_df = lytic_df.copy()
        lytic_df['assing_specificity'] = lytic_df.apply(
            lambda row: [s.strip() for s in row['specificity'].split('/')]
            if pd.notnull(row['specificity']) and "/" in row['specificity'] 
            else [row['assing_K_locus_host_when_specificity_missing']],
            axis=1)
        lytic_df = lytic_df.explode('assing_specificity')
    else:
        lytic_df['assing_specificity'] = lytic_df['assing_K_locus_host_when_specificity_missing']
    
    # --- Generate counts for each category ---
    kaspa_ref_series = pd.concat([zdk_df, predicted_df]).groupby('assing_K_locus_host_when_specificity_missing')["# KASPAH-REF"].first().fillna(0)
    
    gwas_series = bacteria_df.groupby("MGG_K_locus")["MGG_SC"].nunique().fillna(0)
    
    active_series = zdk_df.query('expression == "PRODUCED" and activity == "ACTIVE"') \
                            .groupby('assing_K_locus_host_when_specificity_missing').size().fillna(0)
    
    inactive_series = zdk_df.query('expression == "PRODUCED" and activity == "INACTIVE"') \
                              .groupby('assing_K_locus_host_when_specificity_missing').size().fillna(0)
    
    not_expressed_series = zdk_df.query('expression == "NOT_PRODUCED"') \
                                  .groupby('assing_K_locus_host_when_specificity_missing').size().fillna(0)
    
    # Group predicted proteins by the (possibly exploded) specificity.
    predicted_group = predicted_df.groupby('assing_specificity')
    
    # Group literature search proteins by the (possibly exploded) specificity.
    literature_depolymerases_series = lytic_df.groupby('assing_specificity').size().fillna(0)
    
    predicted_dict = {}
    for kl, group in predicted_group:
        counts = group['prediction_strength'].value_counts().to_dict()
        predicted_dict[kl] = counts
    
    all_kloci = set(kaspa_ref_series.index) | set(active_series.index) | set(custom_gwas_kloci) | \
                set(inactive_series.index) | set(not_expressed_series.index) | set(predicted_dict.keys()) | \
                set(literature_depolymerases_series.index)
    
    input_dict = {}
    for kl in all_kloci:
        input_dict[kl] = {
            "KASPAH_REF": int(kaspa_ref_series.get(kl, 0)),
            "GWAS_KSC": int(gwas_series.get(kl, 0)),
            "EXPRESSED": {
                "ACTIVE": int(active_series.get(kl, 0)),
                "INACTIVE": int(inactive_series.get(kl, 0)),
                "NOT_PRODUCED": int(not_expressed_series.get(kl, 0))
            },
            "PREDICTION": predicted_dict.get(kl, {}),
            "LITERATURE_SEARCH_DEPOLYMERASES": int(literature_depolymerases_series.get(kl, 0))
        }
    
    # Run the visualization.
    generate_pdf_visualization(input_dict, output_pdf_path, custom_gwas_kloci=custom_gwas_kloci, merge_expr=merge_expr, 
                               cell_number_color=cell_number_color, target_fig_height=target_fig_height, target_fig_width=target_fig_width, 
                               labels_fontsize=labels_fontsize, cells_fontsize=cells_fontsize, margin=margin, header_space=header_space)


def FIGURE4_PANELB(structure, params, min_cov=0.5, min_pident=0, max_evalue=1e-3, use_coverage=False):
    import pandas as pd
    import numpy as np
    import subprocess
    from pathlib import Path

    def generate_edges(
        nodes_df: pd.DataFrame,
        tmp_dir: str,
        edges_file: str,
        min_cov: float = 0.5,
        min_pident: float = 20,
        max_evalue: float = 1e-3,
        seqidcol='proteinID',
        use_coverage: bool = False,
        intervals: list = None  # Expected as [bottom, mid, top]
    ) -> pd.DataFrame:
        """
        Generate a DataFrame of edges between proteins based on BLASTp similarity.

        Parameters
        ----------
        nodes_df : pd.DataFrame
            Must include at least columns: ['proteinID', 'seq'].
        tmp_dir : str
            Directory to store intermediate BLAST files.
        edges_file : str
            Path to final edges TSV output file.
        min_cov : float, optional
            Minimum coverage threshold. Default=0.5.
        min_pident : float, optional
            Minimum percent-identity threshold (ignored if use_coverage is True). Default=20.
        max_evalue : float, optional
            Maximum e-value threshold. Default=1e-3.
        seqidcol : str, optional
            Column name for sequence IDs.
        use_coverage : bool, optional
            If True, filtering and interval assignment use coverage instead of percent identity.
        intervals : list, optional
            List of three numbers defining the interval boundaries.
            For pident: [bottom, mid, top] (e.g., [min_pident, 50, 80])
            For coverage: [bottom, mid, top] in percentage (e.g., [min_cov*100, 70, 90])
        
        Returns
        -------
        edges_df : pd.DataFrame
            The filtered edges, ready for downstream usage.
        """
        from pathlib import Path
        import subprocess
        import numpy as np

        # Ensure tmp_dir is a Path object.
        tmp_dir = Path(tmp_dir)
        tmp_dir.mkdir(exist_ok=True, parents=True)

        # Define working file paths.
        sequences_file = tmp_dir / 'sequences.fasta'
        blastp_raw     = tmp_dir / 'blastp_raw.tsv'
        blast_not_filtered = tmp_dir / 'blast_not_filtered.tsv'
        blast_processed    = tmp_dir / 'blast_processed.tsv'
        
        # Build FASTA from nodes_df.
        fasta_strs = []
        for _, row in nodes_df.iterrows():
            header = f">{row[seqidcol]}\n"
            seq    = f"{row['seq']}\n"
            fasta_strs.append(header + seq)
        with open(sequences_file, 'w') as f:
            f.write(''.join(fasta_strs))

        # Create BLAST database and run BLASTp.
        makeblastdb_cmd = [
            'makeblastdb',
            '-in', str(sequences_file),
            '-dbtype', 'prot'
        ]
        blastp_cmd = [
            'blastp',
            '-query', str(sequences_file),
            '-db', str(sequences_file),
            '-out', str(blastp_raw),
            '-outfmt', '6 qseqid sseqid evalue pident bitscore qlen qstart qend slen sstart send'
        ]
        subprocess.run(makeblastdb_cmd, check=True)
        subprocess.run(blastp_cmd, check=True)

        # Load raw BLAST results.
        edges_df = pd.read_csv(blastp_raw, sep='\t', header=None)
        cols = [
            'query', 'target', 'evalue', 'pident', 'bitscore',
            'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send'
        ]
        edges_df.columns = cols

        # Calculate coverage for query and subject.
        edges_df['qcov'] = np.round((edges_df['qend'] - edges_df['qstart'] + 1) / edges_df['qlen'], 2)
        edges_df['scov'] = np.round((edges_df['send'] - edges_df['sstart'] + 1) / edges_df['slen'], 2)

        # Deduplicate self-hits and symmetrical hits.
        edges_df['query-target'] = edges_df.apply(
            lambda row: '-'.join(sorted([row['query'], row['target']])), axis=1
        )
        edges_df = edges_df.drop_duplicates(subset='query-target')

        # Save unfiltered BLAST results.
        edges_df.to_csv(blast_not_filtered, sep='\t', index=False)

        # Apply filters.
        is_self = edges_df['query'] == edges_df['target']
        filt_cov = (edges_df['scov'] >= min_cov) & (edges_df['qcov'] >= min_cov)
        filt_evalue = edges_df['evalue'] <= max_evalue

        if use_coverage:
            # Use coverage as the metric.
            filt = filt_cov & filt_evalue & (~is_self)
        else:
            # Filter on both percent identity and coverage.
            filt_pident = edges_df['pident'] >= min_pident
            filt = filt_cov & filt_pident & filt_evalue & (~is_self)
        edges_df = edges_df.loc[filt]

        # Assign intervals based on the chosen metric.
        if use_coverage:
            # Calculate combined coverage (minimum of query and subject coverage, in percentage).
            edges_df['coverage'] = edges_df[['qcov', 'scov']].min(axis=1) * 100
            # Set default intervals if not provided.
            if intervals is None:
                intervals = [min_cov * 100, 70, 90]
            bottom, mid, top = intervals[0], intervals[1], intervals[2]
            remove_low_cov = edges_df['coverage'] >= bottom
            filt_low = (edges_df['coverage'] >= bottom) & (edges_df['coverage'] < mid)
            filt_med = (edges_df['coverage'] >= mid) & (edges_df['coverage'] < top)
            filt_high = edges_df['coverage'] >= top

            edges_df = edges_df.loc[remove_low_cov]
            edges_df.loc[filt_low, 'coverage_interval'] = f'low ({bottom}-{mid})'
            edges_df.loc[filt_med, 'coverage_interval'] = f'medium ({mid}-{top})'
            edges_df.loc[filt_high, 'coverage_interval'] = f'high ({top}-100)'
        else:
            # Use percent identity.
            if intervals is None:
                intervals = [min_pident, 50, 80]
            bottom, mid, top = intervals[0], intervals[1], intervals[2]
            remove_low_pident = edges_df['pident'] >= bottom
            filt_low = (edges_df['pident'] >= bottom) & (edges_df['pident'] < mid)
            filt_med = (edges_df['pident'] >= mid) & (edges_df['pident'] < top)
            filt_high = edges_df['pident'] >= top

            edges_df = edges_df.loc[remove_low_pident]
            edges_df.loc[filt_low, 'pident_interval'] = f'low ({bottom}-{mid})'
            edges_df.loc[filt_med, 'pident_interval'] = f'medium ({mid}-{top})'
            edges_df.loc[filt_high, 'pident_interval'] = f'high ({top}-100)'

        # Save processed BLAST results.
        edges_df.to_csv(blast_processed, sep='\t', index=False)

        # Add singletons (nodes with no edges).
        all_nodes = list(nodes_df[seqidcol].unique())
        blast_nodes = set(edges_df['query'].tolist() + edges_df['target'].tolist())
        missing_nodes = [node for node in all_nodes if node not in blast_nodes]
        if missing_nodes:
            singletons_df = pd.DataFrame({
                'query': missing_nodes,
                'target': missing_nodes,
                'pident_interval': ['high ({0}-100)'.format(top)] * len(missing_nodes) if not use_coverage 
                                    else ['high ({}-100)'.format(top)] * len(missing_nodes)
            })
            edges_df = pd.concat([edges_df, singletons_df], ignore_index=True)

        edges_df['interaction'] = 'A'
        edges_df.to_csv(edges_file, sep='\t', index=False)
        return edges_df

    

    def assign_categories_and_labels(df):
        """
        Assigns a category and final label to each row in the DataFrame based on:
          - source (e.g., PROPHAGE_ZDKLAB, LITERATURE_SEARCH, GENSCRIPT, PREDICTION)
          - expression and activity (for PROPHAGE_ZDKLAB and GENSCRIPT)
          - prediction_strength (for PREDICTION)
        
        Final label construction:
          - For categories that map to PROTEINID_K_LOCUS_HOST:
                label = proteinID + "_" + K_locus_host
          - For categories that map to PROTEINID_SPECIFICITY:
                label = proteinID + "_" + specificity
          - For LITERATURE_SEARCH, label = proteinID directly
        """
        
        # Function to determine category for each row
        def get_category(row):
            source = row.get('source', '')
            if source == "PROPHAGE_ZDKLAB":
                if row.get('expression', '') == "NOT_PRODUCED":
                    return "ZDK_NOT_PRODUCED"
                elif row.get('expression', '') == "PRODUCED":
                    if row.get('activity', '') == "ACTIVE":
                        return "ZDK_ACTIVE"
                    elif row.get('activity', '') == "INACTIVE":
                        return "ZDK_INACTIVE"
            elif source == "LITERATURE_SEARCH":
                return "LITERATURE_SEARCH"
            elif source == "GENSCRIPT":
                if row.get('activity', '') == "ACTIVE":
                    return "GENSCRIPT_ACTIVE"
                elif row.get('activity', '') == "INACTIVE":
                    return "GENSCRIPT_INACTIVE"
            elif source == "PREDICTION":
                ps = str(row.get('prediction_strength', ''))
                if ps == "yes":
                    return "PREDICTION_PERFECT"
                elif ps == "maybe (+)":
                    return "PREDICTION_GOOD"
            return None

        # Apply category assignment
        df['category'] = df.apply(get_category, axis=1)

        # Function to build the final label by concatenating columns.
        def get_final_label(row):
            cat = row.get('category', None)
            protein_id = str(row.get('proteinID', ''))
            if cat == "LITERATURE_SEARCH":
                # This category already has specificity in the proteinID.
                return f"{protein_id}_{row.get('specificity', '')}"
            elif cat in ["ZDK_NOT_PRODUCED", "ZDK_INACTIVE", "GENSCRIPT_INACTIVE"]:
                # Use K_locus_host column
                return f"{protein_id}_HOST_{row.get('K_locus_host', '')}"
            elif cat in ["ZDK_ACTIVE", "GENSCRIPT_ACTIVE", "PREDICTION_PERFECT", "PREDICTION_GOOD"]:
                # Use specificity column
                return f"{protein_id}_{row.get('specificity', '')}"
            return None

        # Apply final label construction
        df['label1'] = df.apply(get_final_label, axis=1)
        return df

    def get_label2(label):
        if 'HOST' in label:
            return '_'.join(label.split('_')[-2:])
        else: 
            return label.split('_')[-1]

    # paths
    predictions_and_enzymes = structure['analyze']['predictions_and_enzymes_table']
    nodes_outfile = structure['figures']['FIG4B_NODES']
    edges_outfile = structure['figures']['FIG4B_EDGES']
    tmp_dir = structure['figures']['tmp_blast_dir']

    # read
    predictions_and_enzymes_df = pd.read_csv(predictions_and_enzymes, sep='\t')

    # nodes: generate and save
    cols = ['query', 'proteinID', 'source', 'expression', 'category', 'label1', 'label2', 'seq']
    predictions_and_enzymes_df = assign_categories_and_labels(predictions_and_enzymes_df)

    predictions_and_enzymes_df.index = predictions_and_enzymes_df.index + 1
    predictions_and_enzymes_df['query'] = predictions_and_enzymes_df.index.astype(str) + '_' + predictions_and_enzymes_df['label1']
    predictions_and_enzymes_df['label2'] = predictions_and_enzymes_df['label1'].apply(get_label2)

    predictions_and_enzymes_df[cols].to_csv(nodes_outfile, sep='\t', index=False)

    # edges: generate and save
    edges_df = generate_edges(predictions_and_enzymes_df, tmp_dir, edges_outfile,
                            min_cov=min_cov, min_pident=min_pident, max_evalue=max_evalue,
                            seqidcol='query', use_coverage=True,
                            intervals=[00, 50, 80])


