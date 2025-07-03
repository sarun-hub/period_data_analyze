import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# For Interpolation
from scipy.interpolate import interp1d

# =================================== Data Related ==================================================

# from filename, get the dataframe with Average value (only valid data)
def getValidDataFrame(filename):
    df = pd.read_csv(filename)
    df['Average'] = df.loc[:,['rep1','rep2','rep3']].mean(axis=1)
    # Check the row that valid both Timea and Average
    valid_mask = (~np.isnan(df['Time (h)'])) & (~np.isnan(df['Average'])) & (~np.isinf(df['Time (h)'])) & (~np.isinf(df['Average']))
    df = df[valid_mask]
    return df

def getInterpolatedData(data_dict):
    t = data_dict['Time (h)'].to_list()
    avg = data_dict['Average'].to_list()

    dt_new = 0.01
    t_resampled = np.arange(t[0], t[-1], dt_new)

    interp_func = interp1d(t,avg, kind='cubic')
    avg_resampled = interp_func(t_resampled)

    return  avg_resampled, t_resampled

def run_ave(data,nave = 200):
    nn= len(data)
    y = np.array(data)
    for i in range(nn):
        if i < nave: # (for start) nave is possible 
            y[i] = np.mean(data[:i+nave])
        elif i > nn - nave: # (for end) nn - nave is possible
            y[i] = np.mean(data[i-nave:])
        else: # for middle 
            y[i] = np.mean(data[i-nave: i+nave])

    return y

def find_two_peaks(x, window_size = 100):
    x = np.asarray(x)
    # First peak: global maximum
    first_peak = int(np.argmax(x))

    # Second peak: first local maximum after the first peak
    second_peak = None
    # Check if the center is the maximum in the window
    for i in range(first_peak + (window_size), len(x) - (window_size)):  # shift range to avoid index error
        window = x[i - (window_size):i + (window_size+1)]  # includes x[i-2], x[i-1], x[i], x[i+1], x[i+2]
        center = window[window_size]
        neighbors = np.delete(window, window_size)

        if np.all(center > neighbors):
            second_peak = i
            break

    return (first_peak, second_peak)

def findFirstPeriod2(avg_value, time_value,  window_size = 100):
    first_peak, second_peak = find_two_peaks(avg_value, window_size)
    first_period = time_value[second_peak] - time_value[first_peak]
    return first_period, (first_peak, second_peak)

# ===================================== Directory ==================================================

# For running directory number
def makeNumberedDir(base_path):
    if not os.path.exists(base_path):
        os.makedirs(base_path)
        return base_path

    i = 1
    while True:
        new_path = f"{base_path}({i})"
        if not os.path.exists(new_path):
            os.makedirs(new_path)
            return new_path
        i += 1

def getAllSubdirectories(dir):
    return [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]

# Format of Subdirectory (NOT YET DONE TBD)
# data_0_01_yp -> concentration 0.01, yellow pea compound
def createOutputDirectory(subdirs):
    output_dir = makeNumberedDir('output')
    for subdir in subdirs:
        tag = subdir.split('_')
        if len(tag) == 4 :
            # don't use 
            data_tag = tag[0]
            zero_tag = tag[1]

            conc_tag = tag[2]
            yellow_pea_tag = tag[3]

def createOutputDirectoryVer1(dir):
    return makeNumberedDir('output')
    
def getAllFilesInDirectory(dir):
    return [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir,f))]

def getCompoundName(filename):
    return filename.split('.')[0]

def getConcentrationFromDir(dir):
    if '_' in dir:
        tag = dir.split('_')
        conc_tag = tag[2]
        return float('0.'+conc_tag)
    else :
        return 0
        
# ===================================== Plotting ===================================================
 
# Plot all compounds graph one by one
def plotFittedGraphOnebyOne(data_dicts, conc = 0.01 ,output_path = None):
    if output_path == None or not os.path.isdir(output_path):
        os.makedirs(f'Temp_output_{conc}',exist_ok=True)
        print('No output directory exist, adding output directory or change it.')

    for compound ,data_dict in data_dicts.items():
        t = data_dict['Time (h)']
        avg = data_dict['Average']
        avg_good = data_dict.get('Running Average')
        if avg_good is None:
            print(f"'Running Average' is missing for {compound}!")

        fig = plt.figure()

        plt.plot(t, avg, color='black',label='Raw Data')
        plt.plot(t, avg_good, color='red',label='Run Graph')

        plt.grid(True)
        custom_yticks = [500,1000,1500]
        plt.ylim(0,1500)
        plt.yticks(custom_yticks)
        custom_xticks = [0,24,48,72,96,120,144,168]
        plt.xticks(custom_xticks)
        plt.xlim(0,168)
        plt.xlabel('Time (h)')
        plt.ylabel('Luminescence (count/sec)')
        plt.grid(axis='y', visible=False)
        plt.legend()
        plt.title(f'Compounds {compound}({conc} mg/mL)')
        # plt.show()

        os.makedirs(os.path.join(output_path,'smooth_check'), exist_ok= True)

        fig.savefig(os.path.join(output_path, 'smooth_check', f'graph_{compound}.jpg'))

def plotRawGraphAllInOne(data_dicts, conc = 0.01 ,output_path = None):
    if output_path == None or not os.path.isdir(output_path):
        os.makedirs(f'Temp_output_{conc}',exist_ok=True)
        print('No output directory exist, adding output directory or change it.')

    fig = plt.figure()
    for compound ,data_dict in data_dicts.items():
        t = data_dict['Time (h)']
        avg = data_dict['Average']
        
        plt.plot(t, avg,label=f'{compound}')

    plt.grid(True)
    custom_yticks = [500,1000,1500]
    plt.ylim(0,1500)
    plt.yticks(custom_yticks)
    custom_xticks = [0,24,48,72,96,120,144,168]
    plt.xticks(custom_xticks)
    plt.xlim(0,168)
    plt.xlabel('Time (h)')
    plt.ylabel('Luminescence (count/sec)')
    plt.grid(axis='y', visible=False)
    plt.legend()
    plt.title(f'Compounds ({conc} mg/mL)')
    # plt.show()

    fig.savefig(os.path.join(output_path,f'graph_all_raw.png'))


def plotCleanGraphAllInOne(data_dicts, period_dict,  conc = 0.01 ,output_path = None, spot_peak = False):
    if output_path == None or not os.path.isdir(output_path):
        os.makedirs(f'Temp_output_{conc}',exist_ok=True)
        print('No output directory exist, adding output directory or change it.')

    
    fig = plt.figure()
    for compound ,data_dict in data_dicts.items():
        t = data_dict['Time (h)']
        avg_good = data_dict.get('Running Average')
        if avg_good is None:
            print(f"'Running Average' is missing for {compound}!")

        peaks = period_dict[compound]['peak']

        if spot_peak:
            for peak in peaks:
                plt.scatter(t[peak],avg_good[peak],color = 'red')

        plt.plot(t, avg_good,label=f'{compound}')

    plt.grid(True)
    custom_yticks = [500,1000,1500]
    plt.ylim(0,1500)
    plt.yticks(custom_yticks)
    custom_xticks = [0,24,48,72,96,120,144,168]
    plt.xticks(custom_xticks)
    plt.xlim(0,168)
    plt.xlabel('Time (h)')
    plt.ylabel('Luminescence (count/sec)')
    plt.grid(axis='y', visible=False)
    plt.legend()
    plt.title(f'Compounds ({conc} mg/mL)')
    # plt.show()

    graph_name = 'graph_all_smooth_peak.png' if spot_peak else 'graph_all_smooth.png'
    fig.savefig(os.path.join(output_path,graph_name))