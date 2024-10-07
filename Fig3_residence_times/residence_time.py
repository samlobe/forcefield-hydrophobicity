#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def process_file(file, color, fig_name, show_plot=False):
    y_data = pd.read_csv(file, index_col=0)['COM_y'].values

    # Lists to store the start and end times of docking events
    docking_events = []

    # Define the threshold values
    start_threshold = 63
    end_threshold = 70

    # Initialize the index
    i = 0
    while i < len(y_data):
        # Look for the start of a docking event
        if y_data[i] < start_threshold:
            start_time = i
            # March through indices until the end of the docking event
            while i < len(y_data) and y_data[i] <= end_threshold:
                i += 1
            end_time = i
            # Store the start and end times as a tuple
            docking_events.append((start_time, end_time))
        i += 1

    # Convert frames to time in nanoseconds
    times = np.arange(len(y_data)) * 2 / 1000  # Convert frames to ns

    # Plot the time series data with docking events highlighted
    plt.figure()
    plt.plot(times[::200], y_data[::200] - 60.75, linewidth=1, label='Center-of-mass distance', color=color)
    # for event in docking_events:
    #     plt.axvspan(times[event[0]], times[event[1]-1], color='red', alpha=0.3)
    plt.xlabel('time (ns)', fontsize=14)
    plt.ylabel('peptide center-of-mass height (Å)', fontsize=14)
    plt.title(file[10:-4], fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim(500, 800)
    plt.ylim(top=20)
    # plt.legend()
    plt.savefig(fig_name, dpi=300)
    if show_plot:
        plt.show()
    # plt.close()

    # Calculate residence times
    residence_times = []
    for event in docking_events:
        residence_times.append(times[event[1] - 1] - times[event[0]])

    # Plot histogram of residence times
    plt.figure()
    plt.hist(residence_times, bins=20)
    plt.xlabel('residence time (ns)', fontsize=14)
    plt.ylabel('frequency', fontsize=14)
    plt.title(file[10:-4], fontsize=14)
    plt.savefig(f"{fig_name}_hist.png", dpi=300)
    if show_plot:
        plt.show()
    plt.close()

    return residence_times

# List of files and their corresponding colors and figure names
files = [
    ('com_y_VLG_a03ws.csv', 'red', 'a03ws_timeseries.png'),
    ('com_y_VLG_a99SBdisp.csv', 'blue', 'a99SBdisp_timeseries.png'),
    ('com_y_VLG_C36m.csv', 'green', 'C36m_timeseries.png')
]

# Loop through the files and process each one
all_residence_times = {}
for file, color, fig_name in files:
    residence_times = process_file(file, color, fig_name)
    all_residence_times[file] = residence_times

# output each file's residence times to a csv
for file, residence_times in all_residence_times.items():
    pd.DataFrame(residence_times).to_csv(f"{file[10:-4]}_residence_times.csv", index=False)


# Print summary of residence times
for file, residence_times in all_residence_times.items():
    print(f"File: {file}")
    print(f"Mean residence time: {np.mean(residence_times):.2f} ns")
    print(f"First quartile: {np.percentile(residence_times, 25):.2f} ns")
    print(f"Third quartile: {np.percentile(residence_times, 75):.2f} ns")
    print()

# %%
