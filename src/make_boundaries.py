import os
import cooler
from cooltools import insulation


def get_boundaries(file_name: os.Path, resolution: int, window: int) -> pd.DataFrame:
    '''
    Sets the insulation_table for the given mcool file and window size.

    :param file_name: name of mcool/cool file
    :param resolution: resolution of mcool/cool file
    :param window: the size of the sliding diamond window used to calculate the insulation score
    :return: dataframe with all boundaries
    '''
    clr = cooler.Cooler(f'{file_name}::resolutions/{resolution}')
    boundaries_df = insulation(clr, [window], verbose=True)

    boundaries_df.to_csv(f'{file_name}_{window}_boundaries.csv')
    return boundaries_df
