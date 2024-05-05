import numpy as np
import pandas as pd
import os
import warnings
import cooltools
import cooler
from cooltools import insulation
import pyranges as pr

BIN_SIZE_COEF = 1.5


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


def modify_tads_map_by_condition(tad1_chr_coords: pd.DataFrame, bin_size: int, length_flexibility: float):
    tad1_search_regs = pd.DataFrame()
    tad1_search_regs['chrom'] = tad1_chr_coords['chrom']
    tad1_search_regs['start'] = tad1_chr_coords[
                                    'start'] - BIN_SIZE_COEF * bin_size  # нет проверки на нулевые координаты
    tad1_search_regs['end'] = tad1_chr_coords['end'] + BIN_SIZE_COEF * bin_size
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


def find_split(tad1_chr_coords: pd.DataFrame, tad2_chr_coords: pd.DataFrame,
               bin_size: int = 100_000, length_flexibility: float = 1.1):
    tad1_search_regs = modify_tads_map_by_condition(tad1_chr_coords, bin_size, length_flexibility)
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
    # display(tads_region_intersect)
    return tads_region_intersect


def save_frame(path_save: os.path, option: str, saving_dataframe: pd.DataFrame) -> None:
    saving_dataframe = saving_dataframe.reset_index()
    if option == "merge":
        saving_dataframe = saving_dataframe.rename(columns={"start_tad1": "start_2", "start_tad2": "start_1",
                                                            "end_tad1": "end_2", "end_tad2": "end_1"})

    else:
        saving_dataframe = saving_dataframe.rename(columns={"start_tad1": "start_1", "start_tad2": "start_2",
                                                            "end_tad1": "end_1", "end_tad2": "end_2"})
    print(saving_dataframe)
    saving_dataframe = saving_dataframe.drop(columns=['size_tad1', 'size_tad2', 'index'])
    saving_dataframe = saving_dataframe.dropna()
    save_name = f'{option}_coords.csv'
    save_path_df = os.path.join(path_save, save_name)
    saving_dataframe.to_csv(save_path_df)
    return


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


def create_diff_matrix(main_tad_choords, small_tads_choords, clr_1, clr_2, binsize):
    region = main_tad_choords
    try:
        matrix1 = clr_1.matrix(balance=False).fetch(region)
        matrix2 = clr_2.matrix(balance=False).fetch(region)
    except ValueError:
        return None

    diff_matrix = np.log(matrix1 + 1) - np.log(matrix2 + 1)
    return calculate_intensity(diff_matrix, main_tad_choords, small_tads_choords, binsize)


def calculate_intensity(diff_matrix, main_tad_choords, small_tads_choords, binsize):
    start, end = main_tad_choords[1], main_tad_choords[2]
    square_intensity = []
    hill_intensity = []
    square_intensity_var = []
    hill_intensity_var = []
    intensity_big_tad = np.mean(diff_matrix)
    for tad_id, small_tad in enumerate(small_tads_choords):
        if tad_id == len(small_tads_choords) - 1:
            hill_intensity_tad.append(
                np.mean(diff_matrix[start_1_corrected:end_1_corrected1 + 1, start_1_corrected:end_1_corrected1 + 1]))
            continue
        start1, end1 = small_tad[0], small_tad[1]
        start2, end2 = small_tads_choords[tad_id + 1][0], small_tads_choords[tad_id + 1][1]
        start_1_corrected, end_1_corrected1 = int((start1 - start) / binsize), int((end1 - start) / binsize)
        start_2_corrected, end_2_corrected1 = int((start2 - start) / binsize), int((end2 - start) / binsize)

        square_intensity.append(
            np.mean(diff_matrix[start_1_corrected:end_1_corrected1, start_2_corrected + 1:end_2_corrected1 + 1]))
        square_intensity_var.append(
            np.var(diff_matrix[start_1_corrected:end_1_corrected1, start_2_corrected + 1:end_2_corrected1 + 1]))
        hill_intensity.append(
            np.mean(diff_matrix[start_1_corrected:end_1_corrected1 + 1, start_1_corrected:end_1_corrected1 + 1]))
        hill_intensity_var.append(
            np.var(diff_matrix[start_1_corrected:end_1_corrected1 + 1, start_1_corrected:end_1_corrected1 + 1]))
    square_mean = np.mean(square_intensity)
    square_var = np.mean(square_intensity_var)
    hill_mean = np.mean(hill_intensity)
    hill_var = np.mean(hill_intensity_var)
    diff = square_mean - hill_mean
    return diff


def main_split_merge_detection(clr1_filename, clr2_filename, resolution, binsize,
                               path_tad_1: os.path, path_tad_2: os.path,
                               option: str, path_save: os.path = './'):
    clr_1 = cooler.Cooler(f'{clr1_filename}::resolutions/{resolution}')
    clr_2 = cooler.Cooler(f'{clr2_filename}::resolutions/{resolution}')
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
    # print(tad_split_table)
    return

