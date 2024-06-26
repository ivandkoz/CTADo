import argparse
import logging
import os

from src.calculate_intensity_change import count_tads_change_intensity
from src.func_condition_wrapper import parser_wrapper
from src.split_merge_detect import main_split_merge_detection
from src.tads_plot import visualisation

INTENSITY = 'intensity_change_result.csv'
SPLIT = 'split_coords.csv'
MERGE = 'merge_coords.csv'


def counting_tads(file: os.path) -> int:
    """
    Count the number of TADs in a file.

    :param file: The path to the file containing TAD information.
    :return int: The number of TADs in the file.
    """
    with open(file, 'r') as fp:
        tads = sum(1 for row in fp)
        return tads-1


@parser_wrapper
def parse() -> os.path:
    """
    Perform parsing of command line arguments and execute the analysis.

    :return os.path: The path to the output directory and the count of TADs for each map.
    """
    parser = argparse.ArgumentParser(
        prog='Differential analysis of interacting domains between two contact matrices',
        description='This tool is needed to find four types of changes in TADs\
                           between two contact matrices',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='Good luck!   (∿°○°)∿ .・。.・゜✭・.・。.・゜✭・.・。.・゜✭')

    parser.add_argument('clr1_filename', type=str,
                        help='Name of first contact matrix in mcool/cool format')
    parser.add_argument('clr2_filename', type=str,
                        help='Name of second contact matrix in mcool/cool format')
    parser.add_argument('resolution', type=int, help='Resolution of mcool file')
    parser.add_argument('window', type=int, help='Size of the sliding diamond window')
    parser.add_argument('flank', type=int, help='Flank size in bp')
    parser.add_argument('binsize', type=int, help='Bin size in bp')
    parser.add_argument('clr1_boundaries_name', type=str,
                        help='The first contact matrix boundaries argument,\
                             a dataframe name with TADs boundaries in chrom, start, end format or cooler insulation table')
    parser.add_argument('clr2_boundaries_name', type=str,
                        help='The second contact matrix boundaries argument,\
                             a dataframe name with TADs boundaries in chrom, start, end format or cooler insulation table')
    parser.add_argument('-r1', '--result_df_1_name',
                        default=None, type=str, help='The first contact matrix dataframe name with chrom,\
                             start & end of TADs')
    parser.add_argument('-r2', '--result_df_2_name',
                        default=None, type=str, help='The second contact matrix dataframe name with chrome,\
                             start & end of TADs')
    parser.add_argument('-df', '--result_dataframe_name',
                        default=None, type=str, help='Dataframe name with intersecting TADs of two contact matrices')
    parser.add_argument('-od', '--output_directory', type=str,
                        default=f'{os.getcwd()}', help='The path to the save directory')
    parser.add_argument('-nc', '--number_of_charts', type=int,
                        default=5, help='The number of output charts for each type of change.\
                             If the specified number is greater than the number of events, then all of them will be output.\
                             If number is -1 than all of them will be output.')
    parser.add_argument('-lg', '--logging', type=bool, choices=(True, False),
                        default=False, help='Enables logging')
    parser.add_argument('-t', '--threads', type=int,
                        default=1, help='Parameter for specifying the number of threads')
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.disabled = not args.logging

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    count_tads_change_intensity(args.clr1_filename, args.clr2_filename, args.resolution, args.window,
                                args.flank, args.binsize, args.clr1_boundaries_name, args.clr2_boundaries_name,
                                args.result_df_1_name, args.result_df_2_name, args.result_dataframe_name,
                                args.output_directory, args.threads)
    main_split_merge_detection(args.clr1_filename, args.clr2_filename, args.resolution, args.binsize,
                               f'{args.output_directory}/{args.clr1_filename}_{args.window}_result_df.csv',
                               f'{args.output_directory}/{args.clr2_filename}_{args.window}_result_df.csv',
                               args.output_directory)
    for file in [INTENSITY, SPLIT, MERGE]:
        type_of_change = file[:file.find('_')]
        visualisation(args.clr1_filename, args.clr2_filename, args.clr1_boundaries_name, args.clr2_boundaries_name,
                      args.resolution, args.binsize, args.window, f'{args.output_directory}/{file}',
                      type_of_change, args.output_directory, args.number_of_charts)
    output_dir = os.path.abspath(args.output_directory)
    map1_tad_count = counting_tads(f'{args.output_directory}/{args.clr1_filename}_{args.window}_result_df.csv',)
    map2_tad_count = counting_tads(f'{args.output_directory}/{args.clr2_filename}_{args.window}_result_df.csv',)
    return output_dir, map1_tad_count, map2_tad_count


if __name__ == "__main__":
    parse()
