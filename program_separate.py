#!/usr/bin/env python3
from utils import *
import fire

# Version1: without subdirectory 
def main():
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
        target_rep = ['rep1','rep2','rep3']
        print(f'Get data of {compound_name}')
        compound_df = getValidDataFrame(compound_path,target_rep)

        # Resample data (interpolation)
        print(f'Interpolate {compound_name}')
        rep_samples, t_resampled  = getInterpolatedDataSeparate(compound_df,target_rep)
        
        # Remove noise (use Running Average with window = 200)
        print(f'Running Average {compound_name}')
        period_dict[compound_name] = {}

        for rep_sample,rep_name in zip(rep_samples,target_rep):
            rep_good = run_ave(rep_sample, nave =200)

            # find the first period (two peak for illustration)
            # print(f'Finding peak {compound_name}')
            first_period, two_peak = findFirstPeriod3(rep_good, t_resampled, window_size=100, start_time=12)

            # Collect results for plotting
            period_dict[compound_name][f'{rep_name} period'] = round(float(first_period), 2)
        
        print(f'Done {compound_name}')

    result_period = pd.DataFrame.from_dict(period_dict, orient='index').reset_index()
    result_period.columns = ['Compound Name', 'Rep1 Period', 'Rep2 Period', 'Rep3 Period']
    output_dir = createOutputDirectorySeparate(directory)
    period_csv_name = f'period_{conc_string}_YP.csv' if checkYP(directory) else f'period_{conc_string}.csv'
    result_period.to_csv(os.path.join(output_dir,period_csv_name),index=False)



if __name__ == '__main__':
    main()