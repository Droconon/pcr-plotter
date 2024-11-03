# PCR Plotter Script

This script is designed to process and visualize PCR amplification data from multiple experimental runs. The data is read from CSV files, normalized, and plotted, with options for averaging replicate wells and adding custom labels to the legend. 

## Features

- Processes PCR data from multiple experiments.
- Normalizes data either individually or across all experiments together.
- Averages replicates from the same experiment.
- Customizable legends from a CSV file.
- Generates amplification plots for easy analysis.

## Requirements

- Python 3.7+
- Required Python libraries:
  - pandas
  - matplotlib
  - seaborn

You can install the required libraries by running:

```sh
pip install pandas matplotlib seaborn
```

## Usage

## Directory

The default directory for data files is the directory `.\Runs` however this can be changed within the arguments

The directory `Labels` is for the csv files containing the title and legend data for the runs.

### Command-line Arguments

The script uses the following command-line arguments:

- `-p`, `--path`: Path to the folder containing the run data files. Default is `Runs\`.
- `-o`, `--output-path`: Path to the folder where output plots will be saved. Default is `Plots\`.
- `-s`, `--start`: Starting experiment number (required).
- `-e`, `--end`: Ending experiment number (optional).
- `-n`, `--normalize-together`: Flag to normalize all runs together instead of individually.
- `-b`, `--error-bars`: Flag to add error bars/regions to the plot.
- `-r`, `--row`: Rows to be used (e.g., "ABCG") (required).
- `-w`, `--well`: Pair of numbers for well range (e.g., "3-10") (required).
- `-l`, `--labels-file`: Name of the CSV file containing labels for the graph legend.

### Example Command

```sh
python pcr_plotter.py -s 1 -e 5 -r ABC -w 3-10 -l pcr_labels.csv
```

### Arguments Explanation

- `-s 1 -e 5`: Process experiments 1 through 5.
- `-r ABC`: Use rows A, B, and C.
- `-w 3-10`: Use wells 3 through 10.
- `-l pcr_labels.csv`: Use `pcr_labels.csv` for legend labels.

## File Requirements

- **Run Data Files**: Files must be in CSV format and named as `exp<experiment_number>_quantification_amplification_results.csv`.
- **Label File**: If provided, the CSV file (`pcr_labels.csv`) must contain a single row of labels to be used in the graph legend, in the same order as the wells.

## Output

The script generates amplification plots and saves them to the specified output folder. Each plot is named according to the experiment, rows, and well range used, e.g., `exp9_A-D_3-10_amplification_curves.png`.

## Notes

- Ensure the run data CSV files contain well data in columns with names starting with row letters and followed by well numbers (e.g., C3, D4).
- The "Cycle" column must be present in the data, and should contain numeric values representing the PCR cycle numbers.

## Contributors

- Connor Reintjes [@Droconon] - Author
