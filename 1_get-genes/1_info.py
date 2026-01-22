import pandas as pd

csv_file = 'GSE44076_matrix_Transpose_labeled.csv'

df = pd.read_csv(csv_file)

print("Dataset Information:")
print(f"Number of Rows: {df.shape[0]}")
print(f"Number of Columns: {df.shape[1]}")

print("\nColumn Names and Data Types:")
print(df.dtypes)

print("\nColumn Names:")
print(df.columns.tolist())

print("\nFirst Few Rows of the Dataset:")
print(df.head())
print("\nLast Few Rows of the Dataset:")
print(df.tail())

print("\nSummary Statistics for Numerical Columns:")
print(df.describe(include=[int, float]))

print("\nMissing Values in Each Column:")
print(df.isnull().sum())

print("\nData Types of Each Column:")
print(df.dtypes)
