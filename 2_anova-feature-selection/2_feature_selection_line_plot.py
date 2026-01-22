import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def plot_combined_line_plot(folder_path, output_folder):
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Get all CSV files in the folder
    csv_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) 
                 if f.endswith('.csv') and not f.startswith('._')]
    
    # Initialize the plot
    plt.figure(figsize=(6, 6))
    plt.rcParams.update({'font.size': 8})
    # plot background color
    plt.gca().set_facecolor('whitesmoke')
    
    # Define a color map
    # ['Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
    #                   'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b',
    #                   'tab20c'])
    colors = plt.cm.tab20b(np.linspace(0, 1, len(csv_files)))
    
    for i, csv_file in enumerate(csv_files):
        try:
            # Extract the GSE name from the filename
            file_name = os.path.basename(csv_file)
            title_part = file_name.split('_')[0]  # Take the part before the first underscore

            # Load the CSV file
            df = pd.read_csv(csv_file)
            # Visualizing Column Number (adjust if necessary)
            column_number = 1
            df_sorted = df.sort_values(by=df.columns[column_number], ascending=False).reset_index(drop=True)
            
            # Plot the sorted data as a line plot
            plt.plot(df_sorted[df_sorted.columns[column_number]], 
                     label=title_part, color=colors[i], linewidth=2)
        
        except Exception as e:
            print(f"Error processing file '{csv_file}': {e}. Skipping this file.")
            continue

    # Configure the plot
    # plt.title('Combined Line Plot of ANOVA F-Values')
    plt.xlabel('Number of Genes')
    plt.ylabel('ANOVA F-Value')
    plt.grid(True, linestyle='--', color='gray', linewidth=0.5)
    plt.legend(title='Datasets', fontsize=10, title_fontsize=8)
    
    # Save the figure as an image file in the specified output folder
    output_path_base = os.path.join(output_folder, 'line_plot_anova_f_scores')
    plt.savefig(f'{output_path_base}.png', format='png', dpi=1000)
    plt.savefig(f'{output_path_base}.jpg', format='jpg', dpi=1000)
    plt.savefig(f'{output_path_base}.pdf', format='pdf', dpi=1000)
    plt.savefig(f'{output_path_base}.svg', format='svg', dpi=1000)
    
    # Uncomment the next line to display the plot when running the script interactively
    # plt.show()
    
    print("Combined line plot has been saved successfully.")

# Example usage
folder_path = './Processed_Datasets/Score'  # Folder containing the CSV files
output_folder = './Processed_Datasets/Plots'  # Folder to save the output plot
plot_combined_line_plot(folder_path, output_folder)
