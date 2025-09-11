# Load the distance matrix file
# Install required packages
%pip install pandas numpy
import pandas as pd
import numpy as np

# Load the DV.csv file
dv_matrix = pd.read_csv('/Users/jasonklemm/Desktop/DV.csv')

print('DV.csv file loaded successfully')
print('Shape:', dv_matrix.shape)
print('\
First few rows and columns:')
print(dv_matrix.head())
print('\
Column names:')
print(dv_matrix.columns.tolist())

# Prepare the distance matrix for scikit-bio PCoA
# Set the index to the first column and drop the 'Unnamed: 0' column
dv_matrix.set_index('Unnamed: 0', inplace=True)
# Ensure all values are float
for col in dv_matrix.columns:
    dv_matrix[col] = pd.to_numeric(dv_matrix[col], errors='coerce')

# Check the shape and head to confirm
print('Distance matrix shape:', dv_matrix.shape)
print(dv_matrix.head())

# Import required libraries for PCoA analysis

%pip install matplotlib scikit-bio

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix


# Load metadata and assign MIC groups (using the same metadata as before)
metadata = pd.read_csv('/Users/jasonklemm/Desktop/MRP1_metadata_with_mutations.csv')

def assign_mic_group(mic_value):
    if mic_value < 0.125:
        return 'MIC <0.125'
    elif mic_value <= 0.38:
        return 'MIC 0.125-0.38'
    else:
        return 'MIC >0.38'

metadata['MIC_Group'] = metadata['MIC'].apply(assign_mic_group)
metadata['Sample_ID_clean'] = metadata['Sample_ID'] + '.fna'

# Create mapping from sample ID to MIC group
mic_group_map = dict(zip(metadata['Sample_ID_clean'], metadata['MIC_Group']))

print('Metadata loaded and MIC groups assigned')
print('MIC group distribution:')
print(metadata['MIC_Group'].value_counts())

# Prepare the distance matrix for scikit-bio PCoA
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix

# Ensure the order of samples matches between distance matrix and metadata
ordered_samples = [s for s in dv_matrix.index if s in mic_group_map]
dv_matrix_ordered = dv_matrix.loc[ordered_samples, ordered_samples]

# Create a scikit-bio DistanceMatrix object
skbio_dm = DistanceMatrix(dv_matrix_ordered.values, ids=ordered_samples)

# Run PCoA
pcoa_results = pcoa(skbio_dm)
coords_pcoa_2d = pcoa_results.samples.values[:, :2]
explained_var_2d = pcoa_results.proportion_explained[:2]

# Prepare color mapping and group assignments
color_map = {
    'MIC <0.125': '#2E8B57',
    'MIC 0.125-0.38': '#FF8C00',
    'MIC >0.38': '#DC143C'
}
mic_groups = [mic_group_map.get(g, 'Unknown') for g in ordered_samples]
labels = [g.replace('.fna', '') for g in ordered_samples]

print('PCoA completed using scikit-bio')
print('Variance explained by first two axes:')
print('PCoA1:', round(explained_var_2d[0]*100, 2), '%')
print('PCoA2:', round(explained_var_2d[1]*100, 2), '%')

# Define function for 95% confidence ellipses
def plot_confidence_ellipse(ax, x, y, n_std=1.96, facecolor='none', **kwargs):
    """Create a plot of the covariance confidence ellipse of *x* and *y*."""
    if x.size <= 1:
        return
    
    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x and y
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_x = np.mean(x)
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# Create the 2D PCoA plot with 95% CI ellipses using the DV.csv distance matrix
fig, ax = plt.subplots(figsize=(10, 8))

# Plot each MIC group with different colors and add 95% CI ellipses
for group, color in color_map.items():
    # Get indices for this group
    group_indices = [i for i, g in enumerate(mic_groups) if g == group]
    
    if group_indices:
        # Extract coordinates for this group
        group_x = coords_pcoa_2d[group_indices, 0]
        group_y = coords_pcoa_2d[group_indices, 1]
        
        # Plot scatter points
        ax.scatter(group_x, group_y, c=color, label=group, 
                  s=100, edgecolor='black', alpha=0.7, linewidth=0.5)
        
        # Add 95% confidence ellipse if we have enough points
        if len(group_indices) > 1:
            plot_confidence_ellipse(ax, group_x, group_y, n_std=1.96, 
                                  edgecolor=color, linewidth=2) #alpha=0.7)

# Add passage labels to each point
#for i, (x, y) in enumerate(coords_pcoa_2d):
    #ax.text(x, y, labels[i], fontsize=9, weight='bold',
            #ha='center', va='bottom', color='black')

# Customize the plot
ax.set_xlabel('PC1 (' + str(round(explained_var_2d[0]*100, 1)) + '% variance)', fontsize=14)
ax.set_ylabel('PC2 (' + str(round(explained_var_2d[1]*100, 1)) + '% variance)', fontsize=14)
ax.set_title('MRP1 2D PCoA with 95% CI Ellipses', fontsize=16, fontweight='bold')

# Add legend
ax.legend(title='MIC Groups (ug/mL)', title_fontsize=11, fontsize=10, 
          loc='upper right', frameon=True, fancybox=True, shadow=True)

# Add grid and improve appearance
ax.grid(True, alpha=0.3)
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('/Users/jasonklemm/Desktop/MRP1_fna/pcoa_2d_dv_matrix_ellipses.png', dpi=300, bbox_inches='tight')
plt.show()

print('2D PCoA plot with 95% CI ellipses using DV.csv matrix saved as pcoa_2d_dv_matrix_ellipses.png')
print('Total variance explained by first 2 axes: ' + str(round(sum(explained_var_2d)*100, 2)) + '%')