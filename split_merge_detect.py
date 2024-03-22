import numpy as np
import pandas as pd
import os
import warnings


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
    tads_region_intersect = tads_region_intersect.groupby(['chrom', 'start_tad1', 'end_tad1', 'size_tad1']).agg(
        {'start_tad2': 'min', 'end_tad2': 'max', 'size_tad2': 'sum'})
    tads_region_intersect = tads_region_intersect.reset_index()
    tads_region_intersect = tads_region_intersect.query('size_tad1 >= size_tad2')
    # display(tads_region_intersect)
    return tads_region_intersect


def save_frame(path_save: os.path, option: str, saving_dataframe: pd.DataFrame) -> None:
    save_name = f'{option}_coords.csv'
    save_path_df = os.path.join(path_save, save_name)
    saving_dataframe.to_csv(save_path_df)
    return


def main_split_merge_detection(path_tad_1: os.path, path_tad_2: os.path,
                               option: str, path_save: os.path = './'):
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
            tad_split_table = pd.concat([tad_split_table, find_split(tad1_chr_coords, tad2_chr_coords)])
    save_frame(path_save, option, tad_split_table)
    return tad_split_table

##########################
# Пример использования:
# p1 = "./4DNFIHXCPUAP_rs.mcool_400000_result_df.csv"
# p2 = "./4DNFIL6BHWZL_rs.mcool_400000_result_df.csv"
# split_map  = main_split_merge_detection(p1, p2, 'split')
