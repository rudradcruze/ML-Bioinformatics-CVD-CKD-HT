import pandas as pd
import numpy as np

# Load the CSV file
file_path = 'GSE22356_series_matrix.csv'  # Replace with your CSV file path
data = pd.read_csv(file_path)

# Step 1: Get the number of columns and rows
num_columns = data.shape[1]
num_rows = data.shape[0]

# Step 2: Print column names
column_names = data.columns.tolist()
print("Column Names:")
for col in column_names:
    print(f"- {col}")

# Step 3: Take column name and threshold as input from the user
column_name = input("Enter the column name: ")

# Step 4: Calculate mean, median, mode, standard deviation, min, and max for the column
mean_value = data[column_name].mean()
median_value = data[column_name].median()
mode_value = data[column_name].mode()[0]  # Mode might return multiple values; take the first
std_dev = data[column_name].std()
min_value = data[column_name].min()
max_value = data[column_name].max()

print(f"Mean: {mean_value}\nMax: {max_value}")

threshold = float(input("Enter the threshold value: "))

# Step 5: Filter the data based on the input threshold
filtered_data = data[data[column_name] < threshold]
filtered_num_rows = filtered_data.shape[0]

# Step 6: Identify rows with any blank (NaN) values and remove them
blank_rows = filtered_data[filtered_data.isna().any(axis=1)].shape[0]

# Remove rows with any blank values
filtered_data_cleaned = filtered_data.dropna()

# Save the cleaned filtered data to a new CSV file
cleaned_file_path = 'filtered_data_GSE22356.csv'  # New CSV file for cleaned filtered data
filtered_data_cleaned.to_csv(cleaned_file_path, index=False)

# Step 7: Print the number of rows after removing blank rows
cleaned_num_rows = filtered_data_cleaned.shape[0]

# Step 8: Generate a report and save it to a text file
report = f"""
CSV File Analysis Report
------------------------

Number of Columns: {num_columns}
Number of Rows: {num_rows}

Column Names:
{', '.join(column_names)}

Selected Column: {column_name}

Statistics for {column_name}:
- Mean: {mean_value}
- Median: {median_value}
- Mode: {mode_value}
- Standard Deviation: {std_dev}
- Minimum: {min_value}
- Maximum: {max_value}

Filtering Data:
- Threshold: {threshold}
- Rows after filtering: {filtered_num_rows}

Blank Rows:
- Blank Row Numbers: {blank_rows}
- Total Rows after removing blank rows: {cleaned_num_rows}

Cleaned Filtered Data saved to: {cleaned_file_path}
"""

report_file_path = 'analysis_report_GSE22356.txt'
with open(report_file_path, 'w') as report_file:
    report_file.write(report)

print(f"Report generated and saved to {report_file_path}")


print(f"Cleaned data saved to {cleaned_file_path}")