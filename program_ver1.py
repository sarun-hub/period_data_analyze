#!/usr/bin/env python3
from utils import *
import fire

# Version1: without subdirectory 
def main1():
    # Set directory
    directory = 'data_0_01'
    conc = getConcentrationValue(directory)
    conc_string = getConcentrationString(directory)

    all_compounds = getAllFilesInDirectory(directory)

    data_dict = {}
    period_dict = {}
    for compound in all_compounds:
        # get name
        compound_name = getCompoundName(compound)

        # get path
        compound_path = os.path.join(directory,compound)
        # use path to get DataFrame
        print(f'Get data of {compound_name}')
        compound_df = getValidDataFrame(compound_path)

        # Resample data (interpolation)
        print(f'Interpolate {compound_name}')
        avg_resampled, t_resampled  = getInterpolatedData(compound_df)
        
        # Remove noise (use Running Average with window = 200)
        print(f'Running Average {compound_name}')
        avg_good = run_ave(avg_resampled, nave =200)

        # find the first period (two peak for illustration)
        print(f'Finding peak {compound_name}')
        first_period, two_peak = findFirstPeriod2(avg_good, t_resampled, window_size=100, start_time=12)

        # Collect Data in Dataframe dict (for plotting)
        data_dict[compound_name] = {}
        data_dict[compound_name]['Time (h)'] = t_resampled
        data_dict[compound_name]['Average'] = avg_resampled
        data_dict[compound_name]['Running Average'] = avg_good

        # Collect results for plotting
        period_dict[compound_name] = {}
        period_dict[compound_name]['period'] = round(float(first_period), 2)
        period_dict[compound_name]['peak'] = two_peak
        print(f'Done {compound_name}')

    result_period = pd.DataFrame.from_dict(period_dict, orient='index').reset_index()
    result_period.columns = ['Compound Name', 'Period', 'Peak']

    output_dir = createOutputDirectoryVer1(directory)
    period_csv_name = f'period_{conc_string}_YP.csv' if checkYP(directory) else f'period_{conc_string}.csv'
    result_period.to_csv(os.path.join(output_dir,period_csv_name),index=False)

    print('Start Plotting!')
    plotFittedGraphOnebyOne(data_dict, conc, output_dir)
    plotRawGraphAllInOne(data_dict, conc, output_dir)
    plotCleanGraphAllInOne(data_dict, period_dict, conc, output_dir, spot_peak= True)
    plotCleanGraphAllInOne(data_dict, period_dict, conc, output_dir, spot_peak= False)
    print('All process done!')

def main2():
    print('VERSION 2 IS NOT READY, PLEASE SELECT VERSION 1.')

def versionSelection(version = 1):
    if version == 1:
        main1()
    elif version == 2:
        main2()
    else :
        print('Selected version is not found')

if __name__ == '__main__':
    fire.Fire(versionSelection)