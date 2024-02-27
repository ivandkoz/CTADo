import pandas as pd
import numpy as np
import cooltools
import cooler
import os

from matplotlib.ticker import EngFormatter
import matplotlib.pyplot as plt


def visualisation(file_name_1, file_name_2, resolution, binsize, rslt_df_name):
    rslt_df = pd.read_csv(f'{rslt_df_name}', index_col=0)
    df = rslt_df.sort_values('difference', ascending=False).head(5)
    most_diff_tads = []
    binsize = binsize * 1.5
    for index, row in df.iterrows():
        most_diff_tads.append([row['chrom'],  max(row['start_1'], row['start_2']) - binsize * 5,
                               max(row['end_1'], row['end_2']) + binsize * 5])

    if not os.path.exists('graphics'):
        os.makedirs('graphics')

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

    for i in most_diff_tads:
        f, axs = plt.subplots(
            figsize=(12,4),
            ncols=2)

        ax = axs[0]
        start, end = i[1], i[2]
        region = (i[0], i[1], i[2])
        im = ax.matshow(
            clr_1.matrix(balance=False).fetch(region),
            vmax=800,
            extent=(start, end, end, start), cmap='gist_heat_r'
        );
        tad_bgn = start + binsize * 5
        tad_end = end - binsize * 5
        ax.plot([tad_bgn, tad_end], [tad_bgn, tad_bgn], color='blue', linewidth=2, linestyle=':')
        ax.plot([tad_end, tad_end], [tad_bgn, tad_end], color='blue', linewidth=2, linestyle=':')
        ax.set_title(f'graphics/clr_1 {i[0]}: {binsize * 5 +i[1]:,}-{i[2]- binsize * 5:,}', y=1.08)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='raw counts');
        ax.set_ylabel('position, Mb')
        format_ticks(ax)

        ax = axs[1]
        start, end = i[1], i[2]
        region = (i[0], i[1], i[2])
        im = ax.matshow(
            clr_2.matrix(balance=False).fetch(region),
            vmax=800,
            extent=(start, end, end, start), cmap='gist_heat_r'
        );
        tad_bgn = start + binsize * 5
        tad_end = end - binsize * 5
        ax.plot([tad_bgn, tad_end], [tad_bgn, tad_bgn], color='blue', linewidth=2, linestyle=':')
        ax.plot([tad_end, tad_end], [tad_bgn, tad_end], color='blue', linewidth=2, linestyle=':')
        ax.set_title(f'clr_2 {i[0]}: {binsize * 5 +i[1]:,}-{i[2]- binsize * 5:,}', y=1.08)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='raw counts');
        format_ticks(ax)
        plt.savefig(f"graphics/output_{i[0]}_{i[1]}_{i[2]}.jpg")
        plt.tight_layout()


file_name_1 = '4DNFIL6BHWZL_rs.mcool'
file_name_2 = '4DNFIHXCPUAP_rs.mcool'
resolution = 100_000
window = 500_000
flank = 200_000
saving = True
binsize = 100_000
rslt_df_name = 'result_dataframe.csv'
visualisation(file_name_1, file_name_2, resolution, binsize, rslt_df_name)