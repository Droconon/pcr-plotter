import argparse as ap

# Global variable for the default path of the label file
LABEL_FILE_PATH = "" #'Labels\\'
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def parse_input():
    "Parse the file name from the command line input using argparse"
    
    parser = ap.ArgumentParser(description="Specify the experiments to run")
    parser.add_argument('-p', '--path', type=str, help="Path and folder with the run data", default="Runs\\")
    parser.add_argument('-o', '--output-path', type=str, help="Path and folder to output the plots", default="Plots\\")
    parser.add_argument('-s', '--start', type=int, required=True, help="Starting experiment number")
    parser.add_argument('-e', '--end', type=int, required=False, help="Ending experiment number")
    parser.add_argument('-n', '--normalize-together', action='store_true', help='Flag to normalize all replicates together instead of individually')
    parser.add_argument('-b', '--error-bars', action='store_true', help='Flag to add error bars to the plot')
    parser.add_argument('-r', '--row', type=str, required=True, help='Rows to be used (e.g., "ABCG")')
    parser.add_argument('-w', '--well', type=str, required=True, help='Pair of numbers for well range (e.g., "3-10")')
    parser.add_argument('-l', '--labels-file', type=str, help='Name of CSV file containing labels for the graph legend', default=None)

    return parser.parse_args()


def normalize_data(data, scaler=None, well_letters=None, well_range=None, together=False):
    """
    Normalize the PCR data using max min formula for each group of replicates.

    Parameters:
    - data (DataFrame): The PCR data to be normalized.
    - scaler (MinMaxScaler, optional): Scaler instance to use for normalization.
    - well_letters (list of str): List of well letters to normalize.

    Returns:
    - normalized_data (DataFrame): Normalized data.
    """
    normalized_data = data.copy()
    # Normalize all specified columns together
    cols_to_normalize = [col for letter in well_letters for col in data.columns if col != 'Cycle' and col.startswith(letter) and col[1:].isdigit() and well_range[0] <= int(col[1:]) <= well_range[1]]
    if len(cols_to_normalize) > 0:

        if together is True:
            # Calculate global min and max across all specified columns
            global_min = data[cols_to_normalize].min().min()
            global_max = data[cols_to_normalize].max().max()
            print(f"Global Max value before normalization: {global_max}")
            print(f"Global Min value before normalization: {global_min}")
            # Normalize each column using the global min and max
            for col in cols_to_normalize:
                normalized_data[col] = (data[col] - global_min) / (global_max - global_min)
        else:
        # Calculate global min and max across all specified columns
        # Normalize each column using the global min and max
            for letter in well_letters:
                sub_normalize = [col for col in cols_to_normalize if col.startswith(letter)]
                print(sub_normalize)
                local_min = data[sub_normalize].min().min()
                local_max = data[sub_normalize].max().max()
                print(f"Global Max value before normalization: {local_max}")
                print(f"Global Min value before normalization: {local_min}")
                for col in cols_to_normalize:
                    if str(col).startswith(letter):
                        normalized_data[col] = (data[col] - local_min) / (local_max - local_min)
            
    return normalized_data


def import_file_data(file_name, well_letters, well_range):
    "Open the csv folder, pull the data into a dataframe, and return it"
    data = pd.read_csv(file_name)
    # Drop any unnamed columns that may have been created due to extra commas
    data = data.loc[:, ~data.columns.str.contains('^Unnamed')]
    # Ensure cycle column is numeric
    data[data.columns[0]] = pd.to_numeric(data[data.columns[0]], errors='coerce')
    data.dropna(subset=[data.columns[0]], inplace=True)  # Drop rows where cycle column couldn't be converted
    # Filter columns based on specified rows and well range
    valid_columns = [data.columns[0]]  # Keep the Cycle column
    for letter in well_letters:
        valid_columns.extend([col for col in data.columns if col.startswith(letter) and col[1:].isdigit() and well_range[0] <= int(col[1:]) <= well_range[1]])
    data = data[valid_columns]  # Drop rows where cycle column couldn't be converted
    return data


def average_replicates(data, well_letters, well_range):
    """
    Average the replicate wells and calculate their standard deviation if they exist in the dataset.

    Parameters:
    - data (DataFrame): The PCR data with potential replicates.
    - well_letters (list of str): List of well letters to average.
    - well_range (tuple): Tuple specifying the range of well numbers.

    Returns:
    - averaged_data (DataFrame): DataFrame with averaged replicates.
    - std_data (DataFrame): DataFrame with standard deviations for each replicate group.
    """
    averaged_data = data.copy()
    std_data = pd.DataFrame(index=data.index)
    replicate_groups = {}

    # Identify columns belonging to each replicate group
    for col in data.columns[1:]:
        for letter in well_letters:
            if col != 'Cycle' and col.startswith(letter) and col[1:].isdigit() and well_range[0] <= int(col[1:]) <= well_range[1]:
                replicate_number = col[1:]
                if replicate_number not in replicate_groups:
                    replicate_groups[replicate_number] = []
                replicate_groups[replicate_number].append(col)
    
    # Debug: Print replicate groups identified
    print("Replicate Groups Identified:")
    for replicate_number, cols in replicate_groups.items():
        print(f"Replicate {replicate_number}: Columns {cols}")
    
    # Average the replicates and place standard deviation in a separate DataFrame
    for replicate_number, cols in replicate_groups.items():
        # Only average if we have more than one column (i.e., replicates exist)
        if len(cols) > 1:
            averaged_data[cols[0]] = data[cols].mean(axis=1)
            std_data[cols[0]] = data[cols].std(axis=1)
            averaged_data.drop(columns=cols[1:], inplace=True)
    
    return averaged_data, std_data


def check_files(path, start, end=None):
    """
    Process files from the starting experiment number to the ending experiment number (if provided).

    Parameters:
    - path (str): Path to the directory containing the run data files.
    - start (int): Starting experiment number.
    - end (int, optional): Ending experiment number.
    - multiple (bool): Flag to indicate whether to handle experiments with multiple variants.
    """
    if end is None:
        end = start

    file_list = []
    for experiment_number in range(start, end + 1):
        prefix = f"exp{experiment_number}"
        for file_name in os.listdir(path):
            if file_name.startswith(prefix) and file_name.endswith("_quantification_amplification_results.csv"):
                if '-' not in file_name or file_name.count('-') == 1:
                    file_path = os.path.join(path, file_name)
                    file_list.append(file_path)
    
    return file_list


def plot_pcr_samples(args, data, output_path, file_name, std_data, labels=None, graph_title='PCR Sample Amplification Curves', error_bars=False):
    """
    Plot the PCR samples using Seaborn, including error bars from the passed std_data.

    Parameters:
    - args: Arguments containing row and well information.
    - data (DataFrame): The PCR data to be plotted.
    - std_data (DataFrame): DataFrame containing standard deviations for each column.
    - output_path (str): Path to save the plot.
    - file_name (str): Name of the file to include in the plot output name.
    - labels (list of str, optional): Custom labels for the plot legend.
    - graph_title (str, optional): Title for the graph.
    - error_bars (bool, optional): Whether to plot error bars using std_data.
    """
    sns.set_theme(style="whitegrid")  # Set a nicer theme
    sns.set_palette("bright") # Set vibrant palette
    
    # Ensure data is a DataFrame with columns to be used for plotting
    plt.figure(figsize=(12, 8))

    # Plot each column's mean and standard deviation as separate lines
    lines = []
    for col in data.columns[1:]:
        # Plot mean values for each sample/column
        line, = plt.plot(data[data.columns[0]], data[col], linewidth=2, label=col)
        lines.append(line)

    # Plot shaded error bands for each column's mean ± standard deviation using std_data
    if error_bars:
        for col in data.columns[1:]:
            if col in std_data.columns:
                plt.fill_between(
                    data[data.columns[0]], 
                    data[col] - std_data[col], 
                    data[col] + std_data[col], 
                    alpha=0.1  # Reduced alpha for better visibility
                )
    
    plt.title(graph_title, fontsize=18, weight='bold')
    plt.xlabel("Cycle Number", fontsize=14)
    plt.ylabel("Normalized Fluorescence", fontsize=14)
    plt.ylim(-0.1, 1.1)  # Set y-axis limits to better visualize normalized data
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    sns.despine()
    plot_filename = os.path.join(output_path, f"{os.path.splitext(file_name)[0]}_{args.row}_{args.well}_curves.png")
    plt.tight_layout()

    # Add legend if labels are provided
    if labels is not None:
        if len(labels) > len(lines):
            labels = labels[:len(lines)]
        stripped_labels = [s.strip() for s in labels]
        plt.legend(handles=lines, labels=stripped_labels, title='Sample', title_fontsize='13', fontsize='12', loc='upper left')
    else:
        plt.legend(title='Sample', title_fontsize='13', fontsize='12', loc='upper left')

    # Save figure and close
    plt.savefig(plot_filename, dpi=300)
    plt.close()


def main():
    "Main function"
    args = parse_input()
    well_letters = list(args.row.upper())
    well_range = list(map(int, args.well.split('-')))
    file_list = check_files(args.path, args.start, args.end)
    
    # Debug: Print the list of files found
    print("Files found:", file_list)
    
    if not file_list:
        print("No files found matching the criteria. Exiting.")
        return

    labels = None
    if args.labels_file is not None:
        labels_df = pd.read_csv(os.path.join(LABEL_FILE_PATH, args.labels_file), header=None)
        graph_title = labels_df.iloc[0, 0]  # The first cell contains the graph title
        labels = labels_df.iloc[0, 1:].tolist()
    else:
        graph_title = 'PCR Sample Amplification Curves'
    
    # Normalize each file individually
    for file in file_list:
        data = import_file_data(file, well_letters, well_range)
        if data is not None:
            print("Raw Data Head:")
            print(data.head())
            normalized_data = normalize_data(data, well_letters=well_letters, well_range=well_range, together=args.normalize_together)
            print("Normalized Data Head:")
            print(normalized_data.head())
            averaged_data, std_data = average_replicates(normalized_data, well_letters, well_range)
            print("Averaged Data Head:")
            print(averaged_data.head())
            print("Averaged Standard Deviation Head:")
            print(std_data.head())
            plot_pcr_samples(args, averaged_data, args.output_path, os.path.basename(file), std_data, labels, graph_title, error_bars=args.error_bars)


if __name__ == "__main__":
    main()