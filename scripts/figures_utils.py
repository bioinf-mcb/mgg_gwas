import pandas as pd
from pathlib import Path
from Bio import SeqIO
from subprocess import run

import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import numpy as np
import math

from matplotlib.ticker import AutoMinorLocator, MultipleLocator


def process_pyseer_results(pyseer):
    
    # filter results
    npcs = 1690
    minus_log10_bonferroni_pvalue = -1 * math.log10(0.05 / (1 * npcs))

    # load
    results_df = pd.read_csv(pyseer, sep='\t')
    results_df['PC_locus_mode'] = results_df['PC'] + '_' + results_df['locus'] + '_' + results_df['mode']

    # lasso & elastic net
    results_df = results_df.loc[results_df['mode'].isin(['lasso', 'elastic_net'])].reset_index(drop=True)

    # filter
    filt_main_run = (results_df['nbootstrap'] == 0)
    filt_significance = (results_df['-log10(pvalue)'] >= minus_log10_bonferroni_pvalue)
    filt_beta = (results_df['beta'] > 0)

    filt = filt_main_run & filt_significance & filt_beta
    pc_locus_mode = results_df.loc[filt, 'PC_locus_mode'].to_list()

    results_df = results_df.loc[results_df['PC_locus_mode'].isin(pc_locus_mode)].reset_index(drop=True)

    return results_df

def get_quantiles(df, column, Q1=0.25, Q2=0.5, Q3 =0.75):   
    
    # precision quantiles
    quantiles = df[column].quantile([Q1, Q2, Q3])
    
    Q1_value = round(quantiles.loc[Q1], 3)
    Q2_value = round(quantiles.loc[Q2], 3)
    Q3_value = round(quantiles.loc[Q3], 3)

    return Q1_value, Q2_value, Q3_value


def get_quantiles_per_pc(results_df):
    data = []
    for (pc, locus, mode), group in results_df.groupby(['PC', 'locus', 'mode']):

        # precision | recall for main run
        not_sampled_run = group.loc[group['nbootstrap'] == 0]

        beta_value = not_sampled_run['beta'].iloc[0]
        precision_value = not_sampled_run['precision'].iloc[0]
        recall_value = not_sampled_run['recall'].iloc[0] 
        
        REPORTED_T_ECOD = not_sampled_run['REPORTED_T_ECOD'].iloc[0]
        ntopologies = not_sampled_run['ntopologies'].iloc[0]
        multitopologies = not_sampled_run['multitopologies'].iloc[0]
        PC_freq = not_sampled_run['PC-freq'].iloc[0]

        precision_Q1, precision_Q2, precision_Q3 = get_quantiles(group, column='precision')
        recall_Q1, recall_Q2, recall_Q3 = get_quantiles(group, column='recall')

        row = {'PC': pc,
            'locus': locus,
            'mode': mode,
            'beta': beta_value,
            'precision': precision_value,
            'recall': recall_value, 
            f'precision_Q1': precision_Q1, 
            f'precision_Q2': precision_Q2,
            f'precision_Q3': precision_Q3,
            f'recall_Q1': recall_Q1,
            f'recall_Q2': recall_Q2,
            f'recall_Q3': recall_Q3,
            'REPORTED_T_ECOD': REPORTED_T_ECOD,
            'ntopologies': ntopologies,
            'multitopologies': multitopologies,
            'PC-freq': PC_freq
              }

        data.append(row)
        
    predictions_df = pd.DataFrame(data)    
    predictions_df = predictions_df.sort_values(['mode', 'locus', 'precision'], ascending=[False, True, False]).reset_index(drop=True)

    return predictions_df


def process_predictions(predictions_df, min_precision_q25th_circle, min_precision_q25th_square, ecod_colors):

    # filters
    filt_circle = (predictions_df['precision_Q1'] >= min_precision_q25th_circle)
    filt_square = (predictions_df['precision_Q1'] >= min_precision_q25th_square)

    predictions_df = predictions_df.loc[filt_circle].copy()

    # shapes
    predictions_df.loc[filt_circle, 'shape'] = 'o'
    predictions_df.loc[filt_square, 'shape'] = 's'

    # colors
    predictions_df['color'] = predictions_df['REPORTED_T_ECOD'].map(ecod_colors)

    # PC - locus
    predictions_df['PC_locus'] = predictions_df['PC'] + '_' + predictions_df['locus']

    # split
    filt_lasso = (predictions_df['mode'] == 'lasso')
    filt_elastic_net = (predictions_df['mode'] == 'elastic_net')

    is_lasso = (predictions_df['mode'] == 'lasso')
    is_square = (predictions_df['precision_Q1'] >= 0.8)
    is_pectin_lyase = (predictions_df['REPORTED_T_ECOD'] == 'Pectin lyase-like')
    is_no_hit = (predictions_df['REPORTED_T_ECOD'] == 'no hit')
    
    filt_to_compare = (is_lasso & is_pectin_lyase) | (is_lasso & is_square)
    
    lasso_df = predictions_df.loc[filt_lasso]
    elastic_net_df = predictions_df.loc[filt_elastic_net]
    to_compare_with_experiments_df = predictions_df.loc[filt_to_compare]

    return lasso_df, elastic_net_df, to_compare_with_experiments_df


def sequences2clusters(mmseqs_clusters, proteins_fasta, pcs):

    # load
    clusters_df = pd.read_csv(mmseqs_clusters, sep='\t', usecols=[0,2])
    records = list(SeqIO.parse(proteins_fasta, 'fasta'))

    # only repr 
    clusters_df = clusters_df.drop_duplicates('PC')

    # select pcs
    clusters_df = clusters_df.loc[clusters_df['PC'].isin(pcs)]

    # clean
    clusters_df = clusters_df.reset_index(drop=True)

    # map sequences
    seqs = []
    for idx, row in clusters_df.iterrows():
        for r in records:
            if row['repr'] == r.id:
                seqs.append(str(r.seq))
                break

    clusters_df['seq'] = seqs

    return clusters_df




def get_alignments(mmseqs_clusters, proteins_fasta, pcs, alignments_dir):

    # load
    clusters_df = pd.read_csv(mmseqs_clusters, sep='\t')

    # filter
    filt_pcs = clusters_df['PC'].isin(pcs)
    clusters_df = clusters_df.loc[filt_pcs].reset_index(drop=True)

    # each PC
    for pc, group in clusters_df.groupby('PC'):

        # path
        fasta_file = Path(alignments_dir, f'{pc}.fasta')

        # checkpoint
        if Path(fasta_file).exists(): continue

        # alignment
        fasta = []
        for row in group.itertuples():
            for r in SeqIO.parse(proteins_fasta, 'fasta'):
                if row.proteinID == r.id:
                    fasta.append(f'>{row.proteinID}\n{r.seq}\n')
                    break
        fasta = ''.join(fasta)

        # save
        with open(fasta_file, 'w') as f:
            f.write(fasta)



def FIG2_GWAS(bacteria_df, pdf, min_sc_freq = 3, threshold_line=9.5, colors = {'Common SCs': '#3CB371','KASPAH unique SCs': '#ADD8E6', 'KLEBPAVIA unique SCs': '#DDA0DD'}):
    
    # Initialize a dictionary to store counts per K locus
    counts_per_klocus = {}
    
    # Group data by 'MGG_K_locus' and compute counts
    for klocus_id, group in bacteria_df.groupby('MGG_K_locus'):
        # Get unique 'MGG_SC' values for each collection
        kaspah_scs = set(group[group['collection'].isin(['kaspah', 'kaspah_complete'])]['MGG_SC'].unique())
        klebpavia_scs = set(group[group['collection'] == 'klebpavia']['MGG_SC'].unique())
        
        # Determine common and unique sequence clusters
        common_scs = kaspah_scs.intersection(klebpavia_scs)
        kaspah_unique_scs = kaspah_scs - common_scs
        klebpavia_unique_scs = klebpavia_scs - common_scs
    
        # Count the number of sequence clusters
        klocus_counts = {
            'Common SCs': len(common_scs),
            'KASPAH unique SCs': len(kaspah_unique_scs),
            'KLEBPAVIA unique SCs': len(klebpavia_unique_scs),
            'Total SCs': len(common_scs) + len(kaspah_unique_scs) + len(klebpavia_unique_scs)
        }
    
        # Store the counts in the dictionary
        counts_per_klocus[klocus_id] = klocus_counts
    
    # Convert the dictionary to a DataFrame
    counts_df = pd.DataFrame.from_dict(counts_per_klocus, orient='index').reset_index().rename(columns={'index': 'MGG_K_locus'})
    
    # Remove rows where 'MGG_K_locus' is 'UNK'
    counts_df = counts_df[counts_df['MGG_K_locus'] != 'UNK']
    
    # Create boolean flags based on counts
    counts_df['Is_KASPAH_Only'] = counts_df['KASPAH unique SCs'] > 0
    counts_df['Is_Common'] = counts_df['Common SCs'] > 0
    counts_df['Is_KASPAH'] = counts_df['Is_KASPAH_Only'] | counts_df['Is_Common']
    
    # Sort the DataFrame by 'Total SCs' in descending order
    counts_df = counts_df.sort_values('Total SCs', ascending=False)
    
    # Filter based on minimum frequency
    counts_df = counts_df[counts_df['Total SCs'] >= min_sc_freq]
    
    # Prepare data for plotting
    plot_df = counts_df.sort_values('Total SCs', ascending=True)
    plot_df = plot_df[['MGG_K_locus', 'Common SCs', 'KASPAH unique SCs', 'KLEBPAVIA unique SCs']]
    
    # Plotting the horizontal stacked bar plot
    fig, ax = plt.subplots(figsize=(4, 16))
    
    # Plot each category as a stacked bar with outlines
    ax.barh(
        plot_df['MGG_K_locus'], 
        plot_df['Common SCs'], 
        color=colors['Common SCs'], 
        label='Common SCs',
        edgecolor='black', 
        linewidth=0.5
    )
    ax.barh(
        plot_df['MGG_K_locus'], 
        plot_df['KASPAH unique SCs'], 
        left=plot_df['Common SCs'], 
        color=colors['KASPAH unique SCs'], 
        label='KASPAH unique SCs',
        edgecolor='black', 
        linewidth=0.5
    )
    ax.barh(
        plot_df['MGG_K_locus'], 
        plot_df['KLEBPAVIA unique SCs'], 
        left=plot_df['Common SCs'] + plot_df['KASPAH unique SCs'], 
        color=colors['KLEBPAVIA unique SCs'], 
        label='KLEBPAVIA unique SCs',
        edgecolor='black', 
        linewidth=0.5
    )

    # Customize labels and title
    ax.set_xlabel('# bacterial lineages (SCs)', fontsize=12, fontweight='bold')
    # ax.set_ylabel('K locus', fontsize=12)
    ax.set_title(f'MGG K Locus Sequence Clusters (Minimum Frequency >= {min_sc_freq})', fontsize=14)
    ax.legend(fontsize=10)
    
    # Customize ticks and grid
    plt.xticks(fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')
    plt.tight_layout()
    
    # Optional: Add a vertical line at x=10
    ax.axvline(threshold_line, color='black', linestyle='dashed', linewidth=1)
    
    # Add minor ticks to both axes
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(MultipleLocator(1))

    # Adjust grid lines to show on both major and minor ticks for both axes
    ax.grid(which='both', axis='both', linestyle='--', linewidth=0.5)

    # Adjust the y-axis grid to align with categorical data
    ax.set_axisbelow(True)

    # save
    ax.margins(y=0.001)
    plt.savefig(pdf, dpi=400)
    plt.close()



def FIG2(bacteria_df, pdf, min_sc_freq=3, threshold_line=9.5,
         colors={'Common SCs': '#3CB371', 'KASPAH unique SCs': '#ADD8E6', 'KLEBPAVIA unique SCs': '#DDA0DD'},
         kaspah_color='#D22B2B'):
    
    # Filter out 'UNK' in 'MGG_K_locus' for both plots
    bacteria_df = bacteria_df[bacteria_df['MGG_K_locus'] != 'UNK']

    # Process the complete dataset for the top plot
    counts_per_klocus = {}
    for klocus_id, group in bacteria_df.groupby('MGG_K_locus'):
        kaspah_scs = set(group[group['collection'].isin(['kaspah', 'kaspah_complete'])]['MGG_SC'].unique())
        klebpavia_scs = set(group[group['collection'] == 'klebpavia']['MGG_SC'].unique())
        common_scs = kaspah_scs.intersection(klebpavia_scs)
        kaspah_unique_scs = kaspah_scs - common_scs
        klebpavia_unique_scs = klebpavia_scs - common_scs

        klocus_counts = {
            'Common SCs': len(common_scs),
            'KASPAH unique SCs': len(kaspah_unique_scs),
            'KLEBPAVIA unique SCs': len(klebpavia_unique_scs),
            'Total SCs': len(common_scs) + len(kaspah_unique_scs) + len(klebpavia_unique_scs)
        }
        counts_per_klocus[klocus_id] = klocus_counts

    counts_df = pd.DataFrame.from_dict(counts_per_klocus, orient='index').reset_index().rename(columns={'index': 'MGG_K_locus'})
    counts_df = counts_df[counts_df['Total SCs'] >= min_sc_freq]

    # Sort by 'Total SCs' in descending order
    counts_df = counts_df.sort_values('Total SCs', ascending=False)

    # KASPAH99 subset
    kaspah_df = bacteria_df[bacteria_df['collection'] == 'kaspah_complete']
    kaspah_counts_df = kaspah_df.groupby('MGG_K_locus')['MGG_SC'].nunique().reset_index().rename(columns={'MGG_SC': 'Total SCs'})

    # Merge both datasets to align on MGG_K_locus
    merged_df = pd.merge(counts_df, kaspah_counts_df, on='MGG_K_locus', how='left', suffixes=('', '_kaspah')).fillna(0)

    # Set x-axis positions for each MGG_K_locus based on the sorted complete dataset
    x_labels = merged_df['MGG_K_locus']
    x = np.arange(len(x_labels))

    # Initialize the figure with two subplots
    fig, ax = plt.subplots(nrows=2, figsize=(16, 6), sharex=True, gridspec_kw={'height_ratios': [4, 1]})

    # Top Plot - Complete Dataset
    bar1 = ax[0].bar(x, merged_df['Common SCs'], color=colors['Common SCs'], label='Common SCs', edgecolor='black')
    bar2 = ax[0].bar(x, merged_df['KASPAH unique SCs'], bottom=merged_df['Common SCs'], color=colors['KASPAH unique SCs'], label='KASPAH unique SCs', edgecolor='black')
    bar3 = ax[0].bar(x, merged_df['KLEBPAVIA unique SCs'], bottom=merged_df['Common SCs'] + merged_df['KASPAH unique SCs'], color=colors['KLEBPAVIA unique SCs'], label='KLEBPAVIA unique SCs', edgecolor='black')
    ax[0].axhline(threshold_line, color='black', linestyle='dashed', linewidth=1)
    ax[0].grid(True, which='major', axis='both', linestyle='-', linewidth=0.2)

    ax[0].xaxis.set_minor_locator(AutoMinorLocator())
    ax[0].yaxis.set_minor_locator(MultipleLocator(1))
    ax[0].set_axisbelow(True)


    # Bottom Plot - KASPAH99 Distribution
    bar4 = ax[1].bar(x, merged_df['Total SCs_kaspah'], color=kaspah_color, label='KASPAH99 Total SCs', edgecolor='black')
    ax[1].grid(True, which='major', axis='both', linestyle='-', linewidth=0.2)
    ax[1].xaxis.set_minor_locator(AutoMinorLocator())
    ax[1].yaxis.set_minor_locator(MultipleLocator(1))
    ax[1].set_axisbelow(True)

    # y label
    fig.supylabel('# Bacterial Lineages (SCs)\n')
    

    # Set x-axis labels and tick alignment
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(x_labels, rotation=90, fontsize=8, fontdict={'fontweight': 'bold'})

    # Create a combined legend in the top-left corner
    fig.legend(handles=[bar1, bar2, bar3, bar4],
               labels=['Common SCs', 'KASPAH unique SCs', 'KLEBPAVIA unique SCs', 'KASPAH99 Total SCs'],
               bbox_to_anchor=(0.95, 0.9), fontsize=10, frameon=False)

    # Final layout adjustments
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(pdf, dpi=500)  # Saving disabled
    plt.close()  # Close disabled
