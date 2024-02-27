import pandas as pd
import numpy as np
import cooltools
import cooler
from cooltools import insulation
import pyranges as pr

import warnings
warnings.filterwarnings(action='ignore', message='Mean of empty slice')
pd.options.mode.chained_assignment = None


def get_boundaries(file_name, resolution, window, save=False):
    '''
    Sets the insulation_table for the given mcool file and window size.

    :param file_name: name of mcool file
    :param resolution: resolution of mcool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :save: if True saves file, else not
    :return: dataframe with all boundaries
    '''
    clr = cooler.Cooler(f'{file_name}::resolutions/{resolution}')
    boundaries_df = insulation(clr, [window], verbose=True)

    if save:
        boundaries_df.to_csv(f'{file_name}_{window}_boundaries.csv')
    return boundaries_df


def creation_result_dataframe(file_name, resolution, window, boundaries_df_name=None, save=False):
    '''
    Creates a dataframe containing the chromosome name and tad boundaries in the format.
    If boundaries_df_name is missing, a dataframe is created using get_boundaries function

    :param file_name: name of mcool file
    :param resolution: resolution of mcool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :param boundaries_df_name: name of dataframe with all boundaries
    :param save: if True saves file, else not
    :return: dataframe with TADs boundaries
    '''
    
    if not boundaries_df_name:
        df = get_boundaries(file_name, resolution, window, save=save)
    else:
        df = pd.read_csv(f'{boundaries_df_name}', index_col=0)
    df = df.loc[df[f'is_boundary_{window}'] == True]
    df = df.reset_index()
    df.drop('index', axis= 1 , inplace= True)

    out_list = []
    for index, row in df.iterrows():
        if row['chrom'] == df.iloc[index+1][0]:
          start = (row['start'] + row['end']) / 2
          end = (df.iloc[index+1][1] + df.iloc[index+1][2]) / 2
          out_list.append([row['chrom'], int(start), int(end)])
        if index+1 == df.shape[0]-1:
            break

    result_dataframe = pd.DataFrame(out_list, columns=['chrom', 'start', 'end'])

    if save:
        result_dataframe.to_csv(f'{file_name}_{window}_result_df.csv')
    return result_dataframe


def pileup_over_expected(file_name, resolution, window, flank, result_dataframe_name=None, save=False):
    '''
    Counting the average intensity in a given flank within the TAD boundary and normalizing the values.
    Creating a dataframe with chrom, start & end of TADs boundary and average intensity

    :param file_name: name of mcool file
    :param resolution: resolution of mcool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :param result_dataframe_name: ame of dataframe with all boundaries
    :param flank: how much to flank the center of the features by, in bp
    :param save: if True saves file, else not
    :return: dataframe with chrom, start & end of TADs boundaries and average intensity columns
    '''
    if not result_dataframe_name:
        result_dataframe = creation_result_dataframe(file_name, resolution, window, save=save)
    else:
        result_dataframe = pd.read_csv(f'{result_dataframe_name}', index_col=0)

    clr = cooler.Cooler(f'{file_name}::resolutions/{resolution}')
    chroms_view = pd.DataFrame(data={
        'chrom': clr.chromsizes.index,
        'start': [0] * len(clr.chromsizes),
        'end': clr.chromsizes.values,
        'name': clr.chromsizes.index
        })
    expected = cooltools.expected_cis(clr, view_df=chroms_view, nproc=2, chunksize=1_000_000)

    result_dataframe['mean_intensity'] = np.nan
    for index, row in result_dataframe.iterrows():
        data = pd.DataFrame([row['chrom'], row['start'], row['end']]).T
        data.columns=['chrom', 'start', 'end']
        data[['start', 'end']] = data[['start', 'end']].astype(int)
        mean_intensity = np.nanmean(cooltools.pileup(clr, data, expected_df=expected, flank=flank))
        result_dataframe['mean_intensity'][index] = mean_intensity

    if save:
        result_dataframe.to_csv(f'pileup_result_{file_name}_{window}.csv')
    return result_dataframe


def intersect_tads(file_name_1, file_name_2, resolution, window, flank, binsize,
                    result_1_name=None, result_2_name=None, save=False):
    '''
    Creating a table with boundaries intersecting by no more than 1.5 bins from mcool source file or
    from two dataframe with chrom, start & end of TADs boundaries and average intensity columns (optional).

    :param file_name_1: name of first mcool file
    :param file_name_2: name of second mcool file
    :param resolution: resolution of mcool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :param flank: how much to flank the center of the features by, in bp
    :param binsize: bin size, bp
    :param result_1_name: dataframe name with chrom, start & end of TADs boundary and average intensity columns
    :param result_2_name: dataframe name with chrom, start & end of TADs boundary and average intensity columns
    :param save: if True saves file, else not
    :return: dataframe with boundaries intersecting of two mcool files
    '''

    if not result_1_name or not result_2_name:
        result_1 = pileup_over_expected(file_name=file_name_1, resolution=resolution,
                                        window=window, flank=flank, save=saving)
        result_2 = pileup_over_expected(file_name=file_name_2, resolution=resolution,
                                        window=window, flank=flank, save=saving)
    else:
        result_1 = pd.read_csv(f'{result_1_name}', index_col=0)
        result_2 = pd.read_csv(f'{result_2_name}', index_col=0)
    binsize = 1.5 * binsize
    cols = "Chromosome Start End mean_intensity".split()
    result_1.columns = cols
    result_2.columns = cols
    gr1, gr2 = pr.PyRanges(result_1), pr.PyRanges(result_2)
    gr = gr1.join(gr2)
    df = gr.df
    rslt_df = df.loc[((df['End'] - df['End_b'] <= binsize) & (df['End'] - df['End_b'] >= -binsize)) &
                     ((df['Start'] - df['Start_b'] <= binsize) & (df['Start'] - df['Start_b'] >= -binsize))]

    rslt_df.columns = ['chrom', 'start_1', 'end_1', 'mean_intensity_1', 'start_2', 'end_2', 'mean_intensity_2']
    rslt_df['difference'] = abs(np.log2(rslt_df['mean_intensity_1'] / rslt_df['mean_intensity_2']))

    if save:
        rslt_df.to_csv('result_dataframe.csv')
    return rslt_df


file_name_1 = '4DNFIL6BHWZL_rs.mcool'
file_name_2 = '4DNFIHXCPUAP_rs.mcool'
resolution = 100_000
window = 500_000
flank = 200_000
saving = True
binsize = 100_000
result = intersect_tads(file_name_1, file_name_2, resolution, window, flank, binsize, save=saving)