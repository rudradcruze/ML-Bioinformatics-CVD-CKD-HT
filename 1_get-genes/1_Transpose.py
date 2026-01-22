import pandas as pd

# Load the CSV file
file_path = "filtered_data_GSE22356.csv"  
data = pd.read_csv(file_path)

# Transpose the DataFrame
transposed_data = data.T

# Change the output file path to a directory with write access
output_file_path = 'GSE22356_matrix_Transpose.csv'

# Save the transposed DataFrame to a new CSV file
transposed_data.to_csv(output_file_path, header=False, index=True)

print(f"Transposed data saved to {output_file_path}")
