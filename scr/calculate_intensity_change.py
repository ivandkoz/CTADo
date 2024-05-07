import pandas as pd
import numpy as np
import cooltools
import cooler
from cooltools import insulation
import pyranges as pr
import scipy.stats as stats

import warnings
warnings.filterwarnings(action='ignore', message='Mean of empty slice')
pd.options.mode.chained_assignment = None

def get_boundaries(file_name, resolution, window, save=False):
    '''
    Sets the insulation_table for the given mcool file and window size.

    :param file_name: name of mcool/cool file
    :param resolution: resolution of mcool/cool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :save: if True saves file, else not
    :return: dataframe with all boundaries
    '''
    clr = cooler.Cooler(f'{file_name}::resolutions/{resolution}')
    boundaries_df = insulation(clr, [window], verbose=True)

    if save:
        boundaries_df.to_csv(f'{file_name}_{window}_boundaries.csv')
    return boundaries_df


def creation_tads_dataframe(filename, resolution, window, boundaries_df_name, save=False, save_directory):
    '''
    Creates a dataframe containing the chromosome name and tad boundaries on provided resolution.
    If boundaries_df_name is missing, a dataframe is created using get_boundaries function

    :param filename: name of mcool/cool file
    :param resolution: resolution of mcool/cool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :param boundaries_df_name: name of dataframe with all boundaries
    :param save: if True saves file, else not
    :return: dataframe with TADs
    '''
    df = pd.read_csv(f'{boundaries_df_name}', index_col=0)
    
    if df.shape[1] != 3:
        df = df.loc[df[f'is_boundary_{window}'] == True]
        df = df.reset_index()
        df.drop('index', axis= 1 , inplace= True)   
        df = df[['chrom', 'start', 'end']]

    out_list = []
    for index, row in df.iterrows():
        if row['chrom'] == df.iloc[index+1][0]:
          start = (row['start'] + row['end']) / 2
          end = (df.iloc[index+1][1] + df.iloc[index+1][2]) / 2
          out_list.append([row['chrom'], int(start), int(end)])
        if index+1 == df.shape[0]-1:
            break

    result_dataframe = pd.DataFrame(out_list, columns=['chrom', 'start', 'end'])

    result_dataframe.to_csv(f'{save_directory}/{filename}_{window}_result_df.csv')
    return result_dataframe


def intersect_tads(clr1_filename, clr2_filename, resolution, window, binsize, clr1_boundaries_name, clr2_boundaries_name, 
                   result_df_1_name=None, result_df_2_name=None, save=False, save_directory):
    '''
    Creating a table with boundaries intersecting by no more than 1.5 bins from mcool/cool source file or
    from two dataframe with chrom, start & end of TADs boundaries and average intensity columns (optional).

    :param clr1_filename: name of first mcool file
    :param clr2_filename: name of second mcool file
    :param resolution: resolution of mcool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :param flank: how much to flank the center of the features by, in bp
    :param binsize: bin size, bp
    :param clr1_boundaries_name: dataframe name with chrom, start & end of TADs boundary
    :param clr2_boundaries_name: dataframe name with chrom, start & end of TADs boundary
    :param save: if True saves file, else not
    :return: dataframe with boundaries intersecting of two mcool/cool files s TADs
    '''


    if not result_df_1_name or not result_df_2_name:
        result_1 = creation_tads_dataframe(filename=clr1_filename, resolution=resolution,
                                           window=window, boundaries_df_name=clr1_boundaries_name,
                                           save=save, save_directory)
        result_2 = creation_tads_dataframe(filename=clr2_filename, resolution=resolution,
                                           window=window, boundaries_df_name=clr2_boundaries_name,
                                           save=save, save_directory)
    else:
        result_1 = pd.read_csv(f'{result_df_1_name}', index_col=0)
        result_2 = pd.read_csv(f'{result_df_2_name}', index_col=0)


    binsize = binsize * 1.5
    merged = result_1.merge(result_2, on='chrom').query(f'((start_x - {binsize} <= start_y <= end_x + {binsize}) and (start_x - {binsize} <= end_y <= end_x + {binsize})) or'+
                                                        f'((start_y - {binsize} <= start_x <= end_y + {binsize}) and (start_y - {binsize} <= end_x <= end_y + {binsize}))')
    df = merged.loc[(abs(merged['end_x'] - merged['end_y']) <= binsize) & (abs(merged['start_x'] - merged['start_y']) <= binsize)]
    df.columns = ['chrom', 'start_1', 'end_1', 'start_2', 'end_2']

    if save:
        df.to_csv(f'{save_directory}/intersect_result_df.csv')
    return df


def get_pval(x, mean_lmi1, std_lmi1):
    z_score = (x - mean_lmi1) / std_lmi1
    p_value = 1 - stats.norm.cdf(z_score)
    return p_value

def count_pvalue(result_df):
    """
    Calculates p-value for changed intensity TADs
    """
    result_df['log2_mi_1'] = np.log2(result_df['mean_intensity_1'])
    mean_lmi1 = result_df['log2_mi_1'].mean()
    std_lmi1 = result_df['log2_mi_1'].std()
    result_df['pvalue'] = result_df['log2_intensity'].apply(lambda x: get_pval(x, mean_lmi1, std_lmi1))
    del result_df['log2_mi_1']
    return result_df

def count_tads_change_intensity(clr1_filename, clr2_filename, resolution, window, flank, binsize,
                                clr1_boundaries_name, clr2_boundaries_name,
                                result_df_1_name=None, result_df_2_name=None, result_dataframe_name=None,
                                save=False, save_directory):
    '''
    Creating a table with boundaries intersecting by no more than 1.5 bins from mcool/cool source file or
    from two dataframe with chrom, start & end of TADs boundaries and average intensity columns (optional).

    :param clr1_filename: name of first mcool file
    :param clr2_filename: name of second mcool file
    :param resolution: resolution of mcool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :param flank: how much to flank the center of the features by, in bp
    :param binsize: bin size, bp
    :param clr1_boundaries_name: The first contact matrix boundaries argument, a dataframe name with TADs 
                                 boundaries in chrom, start, end format or cooler insulation table
    :param clr2_boundaries_name: The second contact matrix boundaries argument, a dataframe name with TADs
                                 boundaries in chrom, start, end format or cooler insulation table
    :param result_df_1_name: dataframe name with chrom, start & end of TADs
    :param result_df_2_name: dataframe name with chrom, start & end of TADs
    :param save: if True saves file, else not
    :return: an output dataframe with information of two mcool/cool files s TADs that changed their intensity in format: 
            chrom, start_1, end_1, start_2, end_2, mean_intensity_1, mean_intensity_2, log2_intensity, pvalue
    '''
    if not result_dataframe_name:
        result_dataframe = intersect_tads(clr1_filename, clr2_filename, resolution, window, binsize,
                                          clr1_boundaries_name, clr2_boundaries_name, result_df_1_name,
                                          result_df_2_name, save=save, save_directory)
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

    result_dataframe = count_pvalue(result_dataframe)
    result_dataframe.to_csv(f'intensity_change_result.csv')

    return result_dataframe
