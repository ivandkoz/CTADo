import argparse
from calculate_intensity_change import count_tads_change_intensity
from tads_plot import visualisation


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                      prog='TADs intensity', 
                      description='This tool is needed to find changes in TADs intensity', 
                      epilog='Good luck!   (∿°○°)∿ .・。.・゜✭・.・。.・゜✭・.・。.・゜✭')
    parser.add_argument('clr1_filename', type=str, help='Name of first mcool file')
    parser.add_argument('clr2_filename', type=str, help='Name of second mcool file')
    parser.add_argument('resolution', type=int, help='Resolution of mcool file')
    parser.add_argument('window', type=int, help='Size of the sliding diamond window')
    parser.add_argument('flank', type=int, help='Flank size in bp')
    parser.add_argument('binsize', type=int, help='Bin size in bp')
    parser.add_argument('clr1_boundaries_name', type=str, help='The first contact matrix boundaries argument, a dataframe name with TADs boundaries in chrom, start, end format or cooler insulation table')
    parser.add_argument('clr2_boundaries_name', type=str, help='The second contact matrix boundaries argument, a dataframe name with TADs boundaries in chrom, start, end format or cooler insulation table')
    parser.add_argument('-r1', '--result_df_1_name', default=None, type=str, help='The first contact matrix dataframe name with chrom, start & end of TADs boundary')
    parser.add_argument('-r2', '--result_df_2_name', default=None, type=str, help='The second contact matrix dataframe name with chrome, start & end of TADs boundary')
    parser.add_argument('-df', '--result_dataframe_name', default=None, type=str, help='Dataframe name with intersecting TADs of two contact matrixes')
    parser.add_argument('-s', '--save', type=bool, choices=(True, False), default=False, help='True if all result files should be saved, else False')

    args = parser.parse_args()

    count_tads_change_intensity(args.clr1_filename, args.clr2_filename, args.resolution, args.window, args.flank, args.binsize, args.clr1_boundaries_name, args.clr2_boundaries_name, args.result_df_1_name, args.result_df_2_name, args.result_dataframe_name, args.save)

    visualisation(args.clr1_filename, args.clr2_filename, args.clr1_boundaries_name, args.clr2_boundaries_name, args.resolution, args.binsize, args.window, 'intensity_change_result.csv', 'intensity')
