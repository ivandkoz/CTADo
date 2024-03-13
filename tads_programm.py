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


def intersect_tads(file_name_1, file_name_2, resolution, window, binsize,
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
        result_1 = creation_result_dataframe(file_name=clr1_filename, resolution=resolution,
                                        window=window, save=saving)
        result_2 = creation_result_dataframe(file_name=clr2_filename, resolution=resolution,
                                        window=window, save=saving)
    else:
        result_1 = pd.read_csv(f'{result_1_name}', index_col=0)
        result_2 = pd.read_csv(f'{result_2_name}', index_col=0)
   
    binsize = binsize * 1.5
    merged = result_1.merge(result_2, on='chrom').query(f'((start_x - {binsize} <= start_y <= end_x + {binsize}) and (start_x - {binsize} <= end_y <= end_x + {binsize})) or'+
                                                        f'((start_y - {binsize} <= start_x <= end_y + {binsize}) and (start_y - {binsize} <= end_x <= end_y + {binsize}))')
    
    df = merged.loc[(abs(merged['end_x'] - merged['end_y']) <= binsize) & (abs(merged['start_x'] - merged['start_y']) <= binsize)]
    df.columns = ['chrom', 'start_1', 'end_1', 'start_2', 'end_2']
    
    if save:
        df.to_csv('intersect_result_df.csv')
    return df

def pileup_over_expected(clr1_filename, clr2_filename, resolution, window, flank, binsize,
                         result_1_name=None, result_2_name=None, result_dataframe_name=None, save=False):
    '''
    Creates a dataframe with boundaries intersecting of two mcool files, their intensity and log2FC

    :param file_name_1: name of first mcool file
    :param file_name_2: name of second mcool file
    :param resolution: resolution of mcool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :param flank: how much to flank the center of the features by, in bp
    :param binsize: bin size, bp
    :param result_1_name: dataframe name with chrom, start & end of TADs boundary and average intensity columns
    :param result_2_name: dataframe name with chrom, start & end of TADs boundary and average intensity columns
    :param result_dataframe_name: dataframe with boundaries intersecting of two mcool files
    :param save: if True saves file, else not
    :return: dataframe with boundaries intersecting of two mcool files, their intensity and log2FC
    '''

    if not result_dataframe_name:
        result_dataframe = intersect_tads(clr1_filename, clr2_filename, resolution, window, binsize,
                                          result_1_name, result_2_name, save=True)
    else:
        result_dataframe = pd.read_csv(f'{result_dataframe_name}', index_col=0)

    clr1 = cooler.Cooler(f'{clr1_filename}::resolutions/{resolution}')
    chroms_view1 = pd.DataFrame(data={
        'chrom': clr1.chromsizes.index,
        'start': [0] * len(clr1.chromsizes),
        'end': clr1.chromsizes.values,
        'name': clr1.chromsizes.index
        })
    expected1 = cooltools.expected_cis(clr1, view_df=chroms_view1, nproc=2, chunksize=1_000_000)

    clr2 = cooler.Cooler(f'{clr2_filename}::resolutions/{resolution}')
    chroms_view2 = pd.DataFrame(data={
        'chrom': clr2.chromsizes.index,
        'start': [0] * len(clr2.chromsizes),
        'end': clr2.chromsizes.values,
        'name': clr2.chromsizes.index
        })
    expected2 = cooltools.expected_cis(clr2, view_df=chroms_view2, nproc=2, chunksize=1_000_000)

    result_dataframe['mean_intensity_1'] = np.nan
    result_dataframe['mean_intensity_2'] = np.nan
    result_dataframe['log2_intensity'] = np.nan

    for index, row in result_dataframe.iterrows():
        start = min(row['start_1'], row['start_2'])
        end = min(row['end_1'], row['end_2'])

        data = pd.DataFrame([row['chrom'], start, end]).T
        data.columns=['chrom', 'start', 'end']
        data[['start', 'end']] = data[['start', 'end']].astype(int)
        mean_intensity = np.nanmean(cooltools.pileup(clr1, data, expected_df=expected1, flank=flank))
        result_dataframe['mean_intensity_1'][index] = mean_intensity

        data = pd.DataFrame([row['chrom'], start, end]).T
        data.columns=['chrom', 'start', 'end']
        data[['start', 'end']] = data[['start', 'end']].astype(int)
        mean_intensity = np.nanmean(cooltools.pileup(clr2, data, expected_df=expected2, flank=flank))
        result_dataframe['mean_intensity_2'][index] = mean_intensity

        result_dataframe['log2_intensity'] = np.log2(result_dataframe['mean_intensity_1'] / result_dataframe['mean_intensity_2'])

    if save:
        result_dataframe.to_csv(f'pileup_df.csv')
    return result_dataframe



clr1_filename = '4DNFIL6BHWZL_rs.mcool'
clr2_filename = '4DNFIHXCPUAP_rs.mcool'
resolution = 100_000
window = 400_000
flank = 200_000
binsize = 100_000
result_dataframe_name = 'intersect_result_df.csv'
saving = True
pileup_over_expected(clr1_filename, clr2_filename, resolution, window, flank, binsize,
                     result_dataframe_name = result_dataframe_name, save=saving)
