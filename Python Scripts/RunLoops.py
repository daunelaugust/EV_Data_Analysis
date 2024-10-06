import os
import pandas as pd

def process_folder(folder_path):
    # List to store dataframes and file names for each file in the folder
    dfs = []
    filenames = []

    # Get all Excel files in the folder
    for file in os.listdir(folder_path):
        if file.endswith(".xlsx"):
            file_path = os.path.join(folder_path, file)
            
            # Read the file
            df = pd.read_excel(file_path)
            
            # Check if "Name" and "Mature" columns exist
            if 'Name' not in df.columns or 'Mature' not in df.columns:
                print(f"{file} does not contain the necessary columns.")
                continue
            
            # Append dataframe and filename
            dfs.append(df)
            filenames.append(file)
    
    # Check if we have files to process
    if len(dfs) == 0:
        print(f"No valid Excel files found in {folder_path}.")
        return
    
    # Check if "Name" column matches across all files
    first_name_col = dfs[0]['Name']
    for i, df in enumerate(dfs[1:], start=1):
        if not first_name_col.equals(df['Name']):
            print(f"Mismatch found in 'Name' column for file: {filenames[i]}")
            return

    # Create a new DataFrame for the result
    result_df = pd.DataFrame()
    result_df['miRNA'] = first_name_col

    # Loop through each file and add the "Mature" column to the result
    for df, file in zip(dfs, filenames):
        # Rename the "Mature" column
        file_renamed = file.replace("Meany", "").replace("expression values", "").replace(".xlsx", "").strip()
        result_df[file_renamed] = df['Mature']

    # Create a new filename based on the folder name
    folder_name = os.path.basename(folder_path)
    output_file = os.path.join(folder_path, f"{folder_name}.xlsx")

    # Save the result DataFrame to a new Excel file
    result_df.to_excel(output_file, index=False)

    print(f"Processed folder '{folder_name}' successfully. Saved to {output_file}.")

# Example usage:
# Define the path to the folder you want to process
folder_to_process = 'WOBEX Run3'  # Change this to the actual folder path
process_folder(folder_to_process)
