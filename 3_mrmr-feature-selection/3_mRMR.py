import pandas as pd
from mrmr import mrmr_classif  # Ensure this library is installed
import os  # For directory handling

# Load the dataset
dataset_name = './data/GSE66360_matrix_Transpose_labeled_Selected_Top_5000_Features.csv'
df = pd.read_csv(dataset_name)

# Extract the first word from the dataset name by splitting on '_'
dataset_prefix = os.path.basename(dataset_name).split('_')[0]  # Get "GSE66360"

# Top features to select
top_features = 2000

# Separate features and target
X = df.drop(columns=['target'])  # Features
y = df['target']  # Target

# Apply mRMR feature selection
selected_features = mrmr_classif(X=X, y=y, K=top_features)  # Select top features

# Extract feature scores (optional: add your method for scoring if needed)
# Here, we assume scores are ranked by selection order (if actual scores are not available)
feature_scores = {feature: idx + 1 for idx, feature in enumerate(selected_features)}

# Create a DataFrame for the selected features and their scores
features_with_scores = pd.DataFrame({
    'Feature': selected_features,
    'Score': [feature_scores[feature] for feature in selected_features]
})

# Sort the features by Score in descending order
features_with_scores = features_with_scores.sort_values(by='Score', ascending=False)

# Construct the output directory and file name dynamically
output_dir = f'./{top_features}_selected_features'
output_file = f'{output_dir}/{dataset_prefix}_selected_{top_features}_features_with_scores_mRMR.csv'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Save the sorted features with scores to a CSV file
features_with_scores.to_csv(output_file, index=False)

# Print confirmation message
print(f"Selected features and their scores saved to: {output_file}")