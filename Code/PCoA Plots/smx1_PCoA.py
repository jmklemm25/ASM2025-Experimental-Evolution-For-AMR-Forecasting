# Import all necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, permdisp, DistanceMatrix

# Load the distance matrix
df = pd.read_csv('/Users/jasonklemm/Desktop/SMX1_fna/SMX1DV.csv', index_col=0)

# Load the metadata
metadata = pd.read_csv('/Users/jasonklemm/Desktop/SMX1meta.csv')

# Clean the distance matrix (remove rows/columns with all NaN)
df_clean = df.dropna(how='all', axis=0).dropna(how='all', axis=1)

# Remove P12 from the distance matrix -- failed QC
omit_sample = 'P12.fna'
df_clean_omit = df_clean.drop(index=omit_sample, columns=omit_sample)

# Remove P12 from metadata and get sample names without .fna extension
metadata_omit = metadata[metadata['Sample_ID'] != 'P12']
sample_names_omit = [name.replace('.fna', '') for name in df_clean_omit.index]

# Get MIC values for remaining samples
mic_values_omit = []
for sample in sample_names_omit:
    mic_val = metadata_omit[metadata_omit['Sample_ID'] == sample]['MIC'].iloc[0]
    mic_values_omit.append(mic_val)

# Assign MIC groups
mic_group_labels_omit = []
for mic in mic_values_omit:
    if mic < 32:
        mic_group_labels_omit.append('<32')
    elif 32 <= mic <= 128:
        mic_group_labels_omit.append('32-128')
    else:
        mic_group_labels_omit.append('>128')

# Define colors for groups
colors = {'<32': '#2E8B57', '32-128': '#FF8C00', '>128': '#DC143C'}

# Run PCoA
pcoa_results_omit = pcoa(df_clean_omit.values)

# Create DataFrame for plotting
pcoa_df_omit = pd.DataFrame(pcoa_results_omit.samples.iloc[:, :3], columns=['PC1', 'PC2', 'PC3'])
pcoa_df_omit['Sample'] = df_clean_omit.index
pcoa_df_omit['MIC_Group'] = mic_group_labels_omit

# Helper function to plot 95% CI ellipses
def plot_ellipse(ax, x, y, color, label=None):
    if len(x) < 3:
        return  # Not enough points for ellipse
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * 1.96 * np.sqrt(vals)  # 95% CI
    mean_x, mean_y = np.mean(x), np.mean(y)
    ellipse = Ellipse((mean_x, mean_y), width, height, angle=theta, edgecolor=color, facecolor='none', lw=2, label=label)
    ax.add_patch(ellipse)

# Create the plot
fig, ax = plt.subplots(figsize=(10, 8))

# Plot points and ellipses for each group
for group in ['<32', '32-128', '>128']:
    group_data = pcoa_df_omit[pcoa_df_omit['MIC_Group'] == group]
    ax.scatter(group_data['PC1'], group_data['PC2'],
               c=colors[group], label='MIC ' + group, s=100, alpha=0.7, edgecolors='black', linewidth=0.5)
    plot_ellipse(ax, group_data['PC1'], group_data['PC2'], colors[group])

# Add sample labels
#for i, row in pcoa_df_omit.iterrows():
    #ax.annotate(row['Sample'], (row['PC1'], row['PC2']),
               # xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.7)

# Format the plot
ax.set_xlabel('PC1 (42.2% variance)', fontsize = 14)
ax.set_ylabel('PC2 (16.4% variance)', fontsize = 14)
ax.set_title('SMX1 2D PCoA with 95% CI Ellipses', fontsize=16, fontweight='bold')
ax.legend(title='MIC Groups (ug/mL)', title_fontsize=11, fontsize=10, frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('/Users/jasonklemm/Desktop/SMX1_fna/pcoa_2d_dv_matrix_ellipses.png', dpi=300, bbox_inches='tight')
plt.show()



# Print summary statistics
print('Analysis Summary:')
print('- Samples analyzed:', len(sample_names_omit))
print('- Variance explained by PC1:', round(pcoa_results_omit.proportion_explained[0] * 100, 1), '%')
print('- Variance explained by PC2:', round(pcoa_results_omit.proportion_explained[1] * 100, 1), '%')
print('- Group counts:', pd.Series(mic_group_labels_omit).value_counts().to_dict())
