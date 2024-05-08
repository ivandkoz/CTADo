import numpy as np
import pandas as pd
import os
import warnings
import cooltools
import cooler
import pyranges as pr
from cooltools import insulation
from scipy.stats import mannwhitneyu
from src.func_condition_wrapper import wrapper_print

warnings.simplefilter(action='ignore', category=FutureWarning)
BINSIZE_COEF = 1.5


def create_tads_tables(path_tad_1: os.path, path_tad_2: os.path):
    tad1 = pd.read_csv(path_tad_1, index_col=0)
    tad2 = pd.read_csv(path_tad_2, index_col=0)
    return tad1, tad2


def get_chrom_list(tad1: pd.DataFrame, tad2: pd.DataFrame) -> np.ndarray:
    tad1_chrom_list = tad1["chrom"].unique()
    tad2_chrom_list = tad2["chrom"].unique()
    # print(tad1_chrom_list, tad2_chrom_list)
    if len(tad1_chrom_list) != len(tad2_chrom_list):
        warnings.warn("Different numbers of chromosomes were detected!", UserWarning)
    tad_chrom_list = [chrom for chrom in tad1_chrom_list if chrom in tad2_chrom_list]
    return tad_chrom_list


def get_chroms_coords(tad1: pd.DataFrame, tad2: pd.DataFrame, chrom: str):
    tad1_chr_coords = tad1.query("chrom == @chrom")
    tad2_chr_coords = tad2.query("chrom == @chrom")
    return tad1_chr_coords, tad2_chr_coords


def modify_tads_map_by_condition(tad1_chr_coords: pd.DataFrame, binsize: int, length_flexibility: float):
    tad1_search_regs = pd.DataFrame()
    tad1_search_regs['chrom'] = tad1_chr_coords['chrom']
    tad1_search_regs['start'] = tad1_chr_coords['start'] - BINSIZE_COEF * binsize  # нет проверки на нулевые координаты
    tad1_search_regs['end'] = tad1_chr_coords['end'] + BINSIZE_COEF * binsize
    tad1_search_regs['size'] = (tad1_chr_coords['end'] - tad1_chr_coords['start']) * length_flexibility
    return tad1_search_regs


def find_min_and_max_tad_coords(tad1_2_regions: pd.DataFrame):
    tad2_regions = pd.DataFrame
    tad2_regions['start_tad2'] = tad1_2_regions.start_tad2.min()
    tad2_regions['end_tad2'] = tad1_2_regions.end_tad2.max()
    tad2_regions['size_tad2'] = tad1_2_regions.size_tad2.sum()
    return tad2_regions


def add_size_column(tad2_chr_coords: pd.DataFrame):
    tad2_comp_regs = tad2_chr_coords
    tad2_comp_regs = tad2_comp_regs.assign(size=tad2_comp_regs['end'] - tad2_comp_regs['start'])

    return tad2_comp_regs


def demodify_tads_map(tads_region_intersect: pd.DataFrame, binsize: int, length_flexibility: float):
    tads_region_intersect['start_tad1'] = tads_region_intersect['start_tad1'] + BINSIZE_COEF * binsize
    tads_region_intersect['end_tad1'] = tads_region_intersect['end_tad1'] - BINSIZE_COEF * binsize
    tads_region_intersect['size_tad1'] = tads_region_intersect['end_tad1'] - tads_region_intersect['start_tad1']
    return tads_region_intersect


def find_split(tad1_chr_coords: pd.DataFrame, tad2_chr_coords: pd.DataFrame,
               binsize: int = 100_000, length_flexibility: float = 1.1):
    tad1_search_regs = modify_tads_map_by_condition(tad1_chr_coords, binsize, length_flexibility)
    tad2_comp_regs = add_size_column(tad2_chr_coords)
    tads_region_intersect = pd.merge(tad1_search_regs, tad2_comp_regs, on='chrom', how='outer',
                                     suffixes=('_tad1', '_tad2'))
    tads_region_intersect = tads_region_intersect.loc[
        (tads_region_intersect.start_tad1 <= tads_region_intersect.start_tad2) &
        (tads_region_intersect.end_tad1 >= tads_region_intersect.end_tad2) &
        (tads_region_intersect.size_tad1 >= tads_region_intersect.size_tad2)]

    tads_region_intersect = tads_region_intersect[tads_region_intersect.duplicated(subset='start_tad1', keep=False)]
    # display(tads_region_intersect)
    tads_region_intersect_size = tads_region_intersect.groupby(['chrom', 'start_tad1', 'end_tad1', 'size_tad1']).agg(
        {'start_tad2': 'min', 'end_tad2': 'max', 'size_tad2': 'sum'})
    tads_region_intersect_size = tads_region_intersect_size.reset_index()
    tads_region_intersect_size = tads_region_intersect_size.query('size_tad1 >= size_tad2')
    tads_region_intersect = tads_region_intersect[
        tads_region_intersect['start_tad1'].isin(tads_region_intersect_size['start_tad1'])]
    # tads_region_intersect = demodify_tads_map(tads_region_intersect, binsize, length_flexibility)
    return tads_region_intersect


def choose_region(df, clr_1, clr_2, binsize):
    intensity_dict = {}
    small_tads_choords = []
    df_with_value = df
    df_with_value['pvalue'] = None
    for index, row in df.iterrows():
        if index == 0:
            main_tad_choords = row[['chrom', 'start_tad1', 'end_tad1']].to_list()
        # print(main_tad_choords)
        # print(small_tads_choords)
        if main_tad_choords != row[['chrom', 'start_tad1', 'end_tad1']].to_list():
            df_with_value.loc[index, 'pvalue'] = create_diff_matrix(main_tad_choords, small_tads_choords, clr_1, clr_2,
                                                                    binsize)
            # intensity_dict[main_tad_choords] = (square_mean, square_var, hill_mean, hill_var)
            main_tad_choords = row[['chrom', 'start_tad1', 'end_tad1']].to_list()
            small_tads_choords = []
        small_tads_choords.append(row[['start_tad2', 'end_tad2']].to_list())
    return df_with_value


def find_region(main_tad_choords, small_tads_choords):
    main_start = main_tad_choords[1]
    main_end = main_tad_choords[2]
    small_start = max([max(pair) for pair in small_tads_choords])
    small_end = min([min(pair) for pair in small_tads_choords])
    chrom = main_tad_choords[0]
    start = min([main_start, small_start])
    end = max([main_end, small_end])
    region = (chrom, start, end)
    return region


def create_diff_matrix(main_tad_choords, small_tads_choords, clr_1, clr_2, binsize):
    # region = main_tad_choords
    region = find_region(main_tad_choords, small_tads_choords)
    try:
        matrix1 = clr_1.matrix(balance=False).fetch(region)
        matrix2 = clr_2.matrix(balance=False).fetch(region)
    except ValueError:
        print('ALARM')
        return None
    bins = clr_1.bins().fetch(region)
    coords = list(bins[['start', 'end']].itertuples(index=False, name=None))

    # return clr_df
    diff_matrix = np.log(matrix1 + 1) - np.log(matrix2 + 1)
    diff_matrix_df = pd.DataFrame(diff_matrix, columns=coords, index=coords)
    return calculate_intensity(diff_matrix_df, main_tad_choords, small_tads_choords, binsize, coords, region)


def find_coords(position, coords):
    for i, (first, second) in enumerate(coords):
        if first <= position and second >= position:
            return i


def calculate_pvalue(square_mean, hill_mean, square_var, hill_var):
    sample_size = 1000
    square = np.random.normal(square_mean, np.sqrt(square_var), size=sample_size)
    hill = np.random.normal(hill_mean, np.sqrt(hill_var), size=sample_size)

    stat, pvalue = mannwhitneyu(square, hill)
    return pvalue


def calculate_intensity(diff_matrix, main_tad_choords, small_tads_choords, binsize, coords, region):
    start, end = region[1], region[2]
    square_intensity = []
    hill_intensity = []
    square_intensity_var = []
    hill_intensity_var = []
    intensity_big_tad = np.mean(diff_matrix)
    for tad_id, small_tad in enumerate(small_tads_choords):
        if tad_id == (len(small_tads_choords) - 1):
            hill_intensity.append(
                diff_matrix.iloc[start_2_corrected:end_2_corrected + 1,
                start_2_corrected:end_2_corrected + 1].mean().mean())
            hill_intensity_var.append(
                diff_matrix.iloc[start_2_corrected:end_2_corrected + 1,
                start_2_corrected:end_2_corrected + 1].var().var())
            continue
        start1, end1 = small_tad[0], small_tad[1]
        start2, end2 = small_tads_choords[tad_id + 1][0], small_tads_choords[tad_id + 1][1]
        start_1_corrected, end_1_corrected = find_coords(start1, coords), find_coords(end1, coords)
        start_2_corrected, end_2_corrected = find_coords(start2, coords), find_coords(end2, coords)
        # print(diff_matrix.iloc[start_1_corrected:end_1_corrected, start_2_corrected + 1:end_2_corrected + 1].mean().mean())
        square_intensity.append(
            diff_matrix.iloc[start_1_corrected:end_1_corrected,
            start_2_corrected + 1:end_2_corrected + 1].mean().mean())
        square_intensity_var.append(
            diff_matrix.iloc[start_1_corrected:end_1_corrected, start_2_corrected + 1:end_2_corrected + 1].var().var())
        hill_intensity.append(
            diff_matrix.iloc[start_1_corrected:end_1_corrected + 1,
            start_1_corrected:end_1_corrected + 1].mean().mean())
        hill_intensity_var.append(
            diff_matrix.iloc[start_1_corrected:end_1_corrected + 1, start_1_corrected:end_1_corrected + 1].var().var())

    square_mean = np.mean(square_intensity)
    square_var = np.mean(square_intensity_var)
    hill_mean = np.mean(hill_intensity)
    hill_var = np.mean(hill_intensity_var)
    # print(square_mean, hill_mean)
    # diff = square_mean - hill_mean
    pvalue = calculate_pvalue(square_mean, hill_mean, square_var, hill_var)
    return pvalue


def save_frame(path_save: os.path, option: str, saving_dataframe: pd.DataFrame) -> None:
    # saving_dataframe = saving_dataframe.reset_index()
    if option == "merge":
        saving_dataframe = saving_dataframe.rename(columns={"start_tad1": "start_2", "start_tad2": "start_1",
                                                            "end_tad1": "end_2", "end_tad2": "end_1"})

    else:
        saving_dataframe = saving_dataframe.rename(columns={"start_tad1": "start_1", "start_tad2": "start_2",
                                                            "end_tad1": "end_1", "end_tad2": "end_2"})
    # display(saving_dataframe)
    saving_dataframe = saving_dataframe.drop(columns=['size_tad1', 'size_tad2'])
    # saving_dataframe = saving_dataframe.dropna()
    # display(saving_dataframe)
    save_name = f'{option}_coords.csv'
    save_path_df = os.path.join(path_save, save_name)
    saving_dataframe.to_csv(save_path_df)
    return

@wrapper_print
def main_split_merge_detection(clr1_filename, clr2_filename, resolution, binsize,
                               path_tad_1: os.path, path_tad_2: os.path, path_save: os.path = './'):
    clr_1 = cooler.Cooler(f'{clr1_filename}::resolutions/{resolution}')
    clr_2 = cooler.Cooler(f'{clr2_filename}::resolutions/{resolution}')
    split_merge_episodes = []
    for option in ['split', 'merge']:
        tad_split_table = pd.DataFrame()
        tad1, tad2 = create_tads_tables(path_tad_1, path_tad_2)
        if option == 'merge':
            tad1, tad2 = tad2, tad1

        tad_chrom_list = get_chrom_list(tad1, tad2)
        for chrom in tad_chrom_list:
            tad1_chr_coords, tad2_chr_coords = get_chroms_coords(tad1, tad2, chrom)
            if tad_split_table.empty is True:
                tad_split_table = find_split(tad1_chr_coords, tad2_chr_coords)
            else:
                tad_split_table = pd.concat([tad_split_table, find_split(tad1_chr_coords, tad2_chr_coords)],
                                            ignore_index=True)

        final_table = choose_region(tad_split_table, clr_1, clr_2, binsize)
        save_frame(path_save, option, final_table)
        split_merge_episodes.append(len(final_table.drop_duplicates()))
    # print(tad_split_table)
    return tuple(split_merge_episodes)

