import pandas as pd
import numpy as np
import cooltools
import cooler
import os
import pyranges as pr
import re

from matplotlib.ticker import EngFormatter
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings(action='ignore', message='Mean of empty slice')
pd.options.mode.chained_assignment = None
warnings.filterwarnings("ignore")


REPLACE_DICT = {'chr1': '1', 'chr2': '2', 'chr3': '3', 'chr4': '4', 'chr5': '5', 'chr6': '6', 'chr7': '7', 'chr8': '8',
                'chr9': '9', 'chr10': '10', 'chr11': '11', 'chr12': '12', 'chr13': '13', 'chr14': '14', 'chr15': '15',
                'chr16': '16', 'chr17': '17', 'chr18': '18', 'chr19': '19', 'chr20': '20', 'chr21': '21'}

def plot_tads(tads_annot, boundaries_df_clr_filename, window):
    boundaries_df_clr = pd.read_csv(f'{boundaries_df_clr_filename}', index_col=0)
    
    if boundaries_df_clr.shape[1] != 3:
        boundaries_df_clr = boundaries_df_clr.loc[boundaries_df_clr[f'is_boundary_{window}'] == True]
        boundaries_df_clr = boundaries_df_clr.reset_index()
        boundaries_df_clr.drop('index', axis= 1 , inplace= True)   
        boundaries_df_clr = boundaries_df_clr[['chrom', 'start', 'end']]

    boundaries_df_clr['chrom'] = boundaries_df_clr['chrom'].replace(REPLACE_DICT)
    boundaries_df_clr = boundaries_df_clr.rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'})
    boundaries_df_clr = boundaries_df_clr[['Chromosome', 'Start', 'End']]

    merged_clr = tads_annot.merge(boundaries_df_clr, on='Chromosome').query('(Start_x - 5 * 100_000 <= Start_y <= End_x + 5 * 100_000) or (Start_x - 5 * 100_000 <= End_y <= End_x + 5 * 100_000)')
    return merged_clr


# change: intensity or split/merge
def visualisation(file_name_1, file_name_2, boundaries_df_clr1_filename, boundaries_df_clr2_filename,
                  resolution, binsize, window, rslt_df_name, change, save_directory):
    rslt_df = pd.read_csv(f'{rslt_df_name}', index_col=0)

    if change == 'intensity':
        df = rslt_df.sort_values('pvalue', key=abs, ascending=True).head(5)  # .dropna().tail(5)  # .head(5)
        most_diff_tads = []
        binsize = binsize * 1.5
        for index, row in df.iterrows():
            most_diff_tads.append([row['chrom'],  max(row['start_1'], row['start_2']) - binsize * 5,
                                   max(row['end_1'], row['end_2']) + binsize * 5])

    else:
        df = rslt_df.sort_values('pvalue', key=abs, ascending=True)  # .dropna().tail(5)
        if change == 'split':
            df = df.drop_duplicates(subset=['chrom', 'start_1', 'end_1'], keep='first')
        elif change == 'merge':
            df = df.drop_duplicates(subset=['chrom', 'start_2', 'end_2'], keep='first')

        df = df.head(5)
        most_diff_tads = []
        binsize = binsize * 1.5
        if change == 'merge':
            for index, row in df.iterrows():
                most_diff_tads.append([row['chrom'],  row['start_2'] - binsize * 5,
                                       row['end_2'] + binsize * 5])
        elif change == 'split':
            for index, row in df.iterrows():
                most_diff_tads.append([row['chrom'],  row['start_1'] - binsize * 5,
                                       row['end_1'] + binsize * 5])

    tads_annot = pd.DataFrame(most_diff_tads, columns=['Chromosome', 'Start', 'End'])
    tads_annot['Start'] = tads_annot['Start'].apply(lambda x: x + 5 * binsize)
    tads_annot['End'] = tads_annot['End'].apply(lambda x: x - 5 * binsize)

    genes = pd.read_csv('data/ncbi_dataset.tsv',sep='\t')[['Chromosome', 'Begin', 'End', 'Gene_name', 'Symbol', 'Orientation']]
    tads_annot['Chromosome'] = tads_annot['Chromosome'].replace(REPLACE_DICT)
    genes = genes.rename(columns={'Begin': 'Start'})
    gr1, gr2 = pr.PyRanges(genes), pr.PyRanges(tads_annot)
    gr = gr1.intersect(gr2)
    tad_annots = gr.df
    if not os.path.exists(f'{save_directory}/graphics'):
        os.makedirs(f'{save_directory}/graphics')
    tad_annots.to_csv(f"{save_directory}/graphics/annotation_{change}.csv")

    bp_formatter = EngFormatter('b')
    clr_1 = cooler.Cooler(f'{file_name_1}::resolutions/{resolution}')
    clr_2 = cooler.Cooler(f'{file_name_2}::resolutions/{resolution}')

    def format_ticks(ax, x=True, y=True, rotate=True):
        if y:
            ax.yaxis.set_major_formatter(bp_formatter)
        if x:
            ax.xaxis.set_major_formatter(bp_formatter)
            ax.xaxis.tick_bottom()
        if rotate:
            ax.tick_params(axis='x',rotation=45)
    count_top = 0
    for i in most_diff_tads:
        count_top+=1
        if i[0] == 'chrX':
            chrom = 'X'
        elif i[0] == 'chrY':
            chrom = 'Y'
        else:  # i[0] != 'chrX' or i[0] != 'chrY':
            chrom = re.findall("\d+", i[0])[0]

        tad_annot = tad_annots[tad_annots['Chromosome']==str(chrom)][tad_annots['Gene_name'].str.contains('uncharacterized|RNA|miR|pseudogene')==False][tad_annots['Symbol'].str.contains('LOC')==False] 
        tad_annot = tad_annot[tad_annot['Start'] >= i[1] - binsize][tad_annot['End'] <= i[2] + binsize]

        f, axs = plt.subplots(
            figsize=(12,4),
            ncols=2)

        ax = axs[0]
        start, end = i[1], i[2]
        region = (i[0], i[1], i[2])
        im = ax.matshow(
            clr_1.matrix(balance=False).fetch(region),
            vmax=600,
            extent=(start, end, end, start), cmap='gist_heat_r'
        );


        clr1_tads = plot_tads(tads_annot[tads_annot['Chromosome']==str(chrom)][tads_annot['Start'] >= i[1] - binsize][tads_annot['End'] <= i[2] + binsize],
                              boundaries_df_clr1_filename, window)
        tad_bgn = i[1]
        for index, row in clr1_tads.iterrows():
            tad_end = (row[3] + row[4])/2
            if tad_bgn == i[1]:
                ax.plot([tad_end, tad_end], [tad_bgn, tad_end], color='blue', linewidth=1, linestyle='-')
            else:
                ax.plot([tad_bgn, tad_end], [tad_bgn, tad_bgn], color='blue', linewidth=1, linestyle='-')
                ax.plot([tad_end, tad_end], [tad_bgn, tad_end], color='blue', linewidth=1, linestyle='-')
            tad_bgn = tad_end
        tad_end = i[2]
        ax.plot([tad_bgn, tad_end], [tad_bgn, tad_bgn], color='blue', linewidth=1, linestyle='-')


        k=1
        for index, row in tad_annot.iterrows():
            if row['Orientation'] == 'minus':
                bbox_props = dict(boxstyle="larrow,pad=0.3", fc="palegreen", ec="g", lw=1)
                ax.text(row['Start'], start - (k/1.6)* binsize, row['Symbol'], ha="center", va="top", rotation=0, size=5, bbox=bbox_props)
            else:
                bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="palegreen", ec="g", lw=1)
                ax.text(row['Start'], start - (k/1.6)* binsize, row['Symbol'], ha="center", va="top", rotation=0, size=5, bbox=bbox_props)
            k=2 if k==1 else 1

        ax.set_title(f'{change} {count_top} {i[0]}: {binsize * 5 +i[1]:,}-{i[2]- binsize * 5:,}', y=1.08)
        
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='raw counts');
        ax.set_ylabel('position, Mb')
        format_ticks(ax)

        ax = axs[1]
        start, end = i[1], i[2]
        region = (i[0], i[1], i[2])
        im = ax.matshow(
            clr_2.matrix(balance=False).fetch(region),
            vmax=600,
            extent=(start, end, end, start), cmap='gist_heat_r'
        );


        clr2_tads = plot_tads(tads_annot[tads_annot['Chromosome']==str(chrom)][tads_annot['Start'] >= i[1] - binsize][tads_annot['End'] <= i[2] + binsize],
                              boundaries_df_clr2_filename, window)
        tad_bgn = i[1]
        for index, row in clr2_tads.iterrows():
            tad_end = (row[3] + row[4])/2
            if tad_bgn == i[1]:
                ax.plot([tad_end, tad_end], [tad_bgn, tad_end], color='blue', linewidth=1, linestyle='-')
            else:
                ax.plot([tad_bgn, tad_end], [tad_bgn, tad_bgn], color='blue', linewidth=1, linestyle='-') #horizontal
                ax.plot([tad_end, tad_end], [tad_bgn, tad_end], color='blue', linewidth=1, linestyle='-')
            tad_bgn = tad_end
        tad_end = i[2]
        ax.plot([tad_bgn, tad_end], [tad_bgn, tad_bgn], color='blue', linewidth=1, linestyle='-')


        k=1
        for index, row in tad_annot.iterrows():
            if row['Orientation'] == 'minus':
                bbox_props = dict(boxstyle="larrow,pad=0.3", fc="palegreen", ec="g", lw=1)
                ax.text(row['Start'], start - (k/1.6)* binsize, row['Symbol'],
                        ha="center", va="top", rotation=0, size=5, bbox=bbox_props)
            else:
                bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="palegreen", ec="g", lw=1)
                ax.text(row['Start'], start - (k/1.6)* binsize, row['Symbol'],
                        ha="center", va="top", rotation=0, size=5, bbox=bbox_props)
            k=2 if k==1 else 1

        ax.set_title(f'{change} {count_top} {i[0]}: {binsize * 5 +i[1]:,}-{i[2]- binsize * 5:,}', y=1.08)
        
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='raw counts');
        format_ticks(ax)
        plt.savefig(f"{save_directory}/graphics/{change}_{i[0]}_{i[1]}_{i[2]}.jpg")
        plt.tight_layout()
