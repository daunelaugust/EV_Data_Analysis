import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Function to clean xtick labels by removing text after 'S' followed by numbers
def clean_xtick_label(label):
    # Remove leading underscores
    label = re.sub(r'^_', '', label)
    # Remove everything after 'S' followed by digits
    label = re.sub(r'_S\d+.*', '', label)
    return label

# Function to process the spreadsheet and plot the sorted histogram
def plot_sorted_histogram(file_path):
    # Load the spreadsheet
    spreadsheet = pd.read_excel(file_path, engine='openpyxl')

    # Sum each column (ignoring the 'miRNA' column)
    column_sums = spreadsheet.iloc[:, 1:].sum()

    # Define the ranges for the different runs based on columns
    run1_columns = spreadsheet.columns[1:31]  # B-AE (Run1)
    run2_columns = spreadsheet.columns[31:122]  # AE-DU (Run2)
    run3_columns = spreadsheet.columns[124:]  # Remaining columns (Run3)

    # Classify columns into runs
    run1_sums = column_sums[run1_columns]
    run2_sums = column_sums[run2_columns]
    run3_sums = column_sums[run3_columns]

    # Combine all columns and sort them
    all_sums = pd.concat([run1_sums, run2_sums, run3_sums]).sort_values()

    # Create colors based on the run the column belongs to
    colors = ['blue' if col in run1_columns else 'green' if col in run2_columns else 'red' for col in all_sums.index]

    # # Plot the histogram
    # plt.figure(figsize=(14, 8))
    # plt.bar(all_sums.index, all_sums.values, color=colors)

    # # Log scale for the y-axis
    # plt.yscale('log')

    # # Add labels and title
    # plt.xlabel('Sample')
    # plt.ylabel('Sum of Expression Values (log scale)')
    # plt.title('Sorted Histogram of Expression Values by Run (Log Scale)')

    # # Clean the xtick labels
    # cleaned_xticks = [clean_xtick_label(label) for label in all_sums.index]
    # plt.xticks(np.arange(len(cleaned_xticks)), cleaned_xticks, rotation=90, fontsize=8)

    # # Create a legend
    # plt.legend(handles=[
    #     plt.Line2D([0], [0], color='blue', lw=4, label='Run1'),
    #     plt.Line2D([0], [0], color='green', lw=4, label='Run2'),
    #     plt.Line2D([0], [0], color='red', lw=4, label='Run3')
    # ])

    # # Show the plot
    # plt.tight_layout()
    # plt.show()

# Example usage:
file_path = 'All WOBEX RUNS 090424.xlsx'  # Replace with the actual path to your file
plot_sorted_histogram(file_path)

