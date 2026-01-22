import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

Set_1_Label = 'Chronic Kidney - GSE62792'
Set_2_Label = 'HyperTension - GSE22356'
Set_3_Label = 'Cardiovascular - GSE66360'

only_set1 = 1735
only_set1_percentage = '(31.4%)'

only_set2 = 1677
only_set2_percentage = '(30.4%)'

only_set3 = 1657
only_set3_percentage = '(30.0%)'

set1_set2 = 112
set1_set2_percentage = '(2.0%)'

set1_set3 = 132
set1_set3_percentage = '(2.4%)'

set2_set3 = 190
set2_set3_percentage = '(2.4%)'

all_three = 21
all_three_percentage = '(0.4%)'

# Create the diagram
fig, ax = plt.subplots(figsize=(12, 8))

# Shift the diagram downwards by adjusting the ellipse centers
ellipse1 = Ellipse((0.35, 0.3), 0.6, 0.35, angle=-30, edgecolor='green', facecolor='green', alpha=0.3, label=Set_1_Label)
ellipse2 = Ellipse((0.65, 0.3), 0.6, 0.35, angle=30, edgecolor='red', facecolor='red', alpha=0.3, label=Set_2_Label)
ellipse3 = Ellipse((0.5, 0.5), 0.6, 0.35, angle=90, edgecolor='blue', facecolor='blue', alpha=0.3, label=Set_3_Label)

# Add ellipses to the plot
ax.add_patch(ellipse1)
ax.add_patch(ellipse2)
ax.add_patch(ellipse3)

# Add labels for each region
ax.text(0.22, 0.35, str(only_set1), fontsize=18, ha='center', va='center')
ax.text(0.22, 0.32, str(only_set1_percentage), fontsize=8, ha='center', va='center')

ax.text(0.78, 0.35, str(only_set2), fontsize=18, ha='center', va='center')
ax.text(0.78, 0.32, str(only_set2_percentage), fontsize=8, ha='center', va='center')

ax.text(0.5, 0.65, str(only_set3), fontsize=18, ha='center', va='center')
ax.text(0.5, 0.62, str(only_set3_percentage), fontsize=8, ha='center', va='center')

ax.text(0.5, 0.16, str(set1_set2), fontsize=18, ha='center', va='center')
ax.text(0.5, 0.13, str(set1_set2_percentage), fontsize=8, ha='center', va='center')

ax.text(0.6, 0.4, str(set2_set3), fontsize=18, ha='center', va='center')
ax.text(0.6, 0.37, str(set1_set3_percentage), fontsize=8, ha='center', va='center')

ax.text(0.4, 0.4, str(set1_set3), fontsize=18, ha='center', va='center')
ax.text(0.4, 0.37, str(set2_set3_percentage), fontsize=8, ha='center', va='center')

ax.text(0.5, 0.28, str(all_three), fontsize=18, ha='center', va='center', weight='bold')
ax.text(0.5, 0.25, str(all_three_percentage), fontsize=8, ha='center', va='center', weight='bold')


# Add title and legend
# ax.set_title("Venn Diagram Adjusted Downwards", fontsize=20)
ax.legend(loc='upper right')

# Adjust axis and save the plot
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.axis('off')

# Save and show the diagram
output_path = "Common_Genes_Code_Percentage.png"
plt.savefig(output_path, dpi=1000, bbox_inches='tight')
# plt.show()
