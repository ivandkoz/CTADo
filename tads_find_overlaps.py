import numpy as np
import pandas as pd

def find_overlaps(table1: pd.DataFrame, table2: pd.DataFrame, bins: int = 100_000) -> pd.DataFrame:
    """
    Finds overlaps between two tables of TADs with a tolerance of 1.5 bins.

    Args:
        table1 (pd.DataFrame): DataFrame representing the first table of TADs. Must include specified columns:
            - 'chrom' with chromosome number (int ot str)
            - 'start' - coordinate of TAD start
            - 'end' - coordinate of TAD end
        table2 (pd.DataFrame): DataFrame representing the second table of TADs. Must include same specified columns
        bins (int, optional): Size of bins. Defaults to 100'000.

    Returns:
        pd.DataFrame: DataFrame containing overlaps between the two tables. Columns include:
            - Region1 - chromosome number from first table
            - Region2 - chromosome number from second table
            - Table1_length - length for each TAD from first table
            - Table2_length - length for each TAD from second table
            - Table1_coordinates - coordinates for each TAD from first table
            - Table2_coordinates - coordinates for each TAD from second table
    Example:
        results = find_overlaps(table1, table2)
    """
    bins_tolerance = bins*1.5
    overlaps = []
    for idx1, row1 in table1.iterrows():
        for idx2, row2 in table2.iterrows():
            if row1['chrom'] != row2['chrom']: #проверка по хромосоме
              continue
            if (
                (abs(row1['start'] - row2['start']) <= bins_tolerance)  and (row2['end'] < row1['end']) #первый микротад
                ) or (
                (row1['start'] < row2['start']) and (row2['end'] < row1['end']) #средний микротад
                )  or (
               (row1['start'] < row2['start']) and (abs(row1['end'] - row2['end']) <= bins_tolerance)  #последний микротад
               ):
                overlaps.append((row1['chrom'], row2['chrom'],
                                 abs(row1['end'] - row1['start']),
                                 abs(row2['end'] - row2['start']),
                                 (row1['start'], row1['end']), #tuple  или pandas.interval (?)
                                 (row2['start'], row2['end'])))
    return pd.DataFrame(overlaps, columns=['Region1', 'Region2', 'Table1_length', 'Table2_length', 'Table1_coordinates', 'Table2_coordinates'])
