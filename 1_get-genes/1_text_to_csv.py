import csv

# Specify the input and output file paths
input_file = 'GSE113439_series_matrix.txt'
output_file = 'GSE113439_series_matrix.csv'

# Open the input file and read it
with open(input_file, 'r') as infile:
    # Use csv.reader with a tab delimiter
    reader = csv.reader(infile, delimiter='\t')
    
    # Open the output CSV file and write the data
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        
        # Write each row from the reader to the writer
        for row in reader:
            writer.writerow(row)

print(f"File has been successfully converted to {output_file}")
