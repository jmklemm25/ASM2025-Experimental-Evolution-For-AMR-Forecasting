# Load the SMX metadata file  
import pandas as pd  
import numpy as np  
import matplotlib.pyplot as plt  
import seaborn as sns  
  
# Load the data 
smx_df = pd.read_csv('/Users/jasonklemm/Downloads/SMX_Metadata_consolidated_genes_position_specific_p.28L.csv', encoding='ascii')  
  
# Check the basic structure  
print("Dataset shape:", smx_df.shape)  
print("\nColumn names:")  
print(smx_df.columns.tolist())  
print("\nFirst few rows:")  
print(smx_df.head())  

# Add a project-specific folP status column for plotting   
def assign_folP_status(row):  
    if row['Project_code'] == 'SMX1':  
        if row['folP_p.F28L'] == 1:  
            return 'folP p.F28L'  
        else:  
            return 'folP WT'  
    elif row['Project_code'] == 'SMX2':  
        if row['folP_p.F28L'] == 1:  
            return 'folP p.F28L'  
        elif row['folP_p.P64S'] == 1:  
            return 'folP p.P64S'  
        else:  
            return 'folP WT'  
    elif row['Project_code'] == 'SMX3':  
        if row['folP_p.P64S'] == 1:  
            return 'folP p.P64S'  
        else:  
            return 'folP WT'  
    else:  
        return 'Other'  
  
smx_df['project_folP_status'] = smx_df.apply(assign_folP_status, axis=1)  
  
# Show the new column for the first few rows  
display_cols = ['Sample_ID', 'Project_code', 'MIC', 'folP_p.F28L', 'folP_p.P64S', 'project_folP_status']  
print("Updated folP status assignments:")  
print(smx_df[display_cols].head(15))  
  
# Check the distribution of folP status across projects  
print("\nDistribution of folP status by project:")  
status_counts = smx_df.groupby(['Project_code', 'project_folP_status']).size().unstack(fill_value=0)  
print(status_counts)  

# Prepare data and plot the grouped bar chart with updated logic for SMX2  
projects = ['SMX1', 'SMX2', 'SMX3']  
wt_means = []  
wt_stds = []  
wt_counts = []  
f28l_means = []  
f28l_stds = []  
f28l_counts = []  
p64s_means = []  
p64s_stds = []  
p64s_counts = []  
  
for project in projects:  
    project_data = smx_df[smx_df['Project_code'] == project]  
      
    # WT data  
    wt_data = project_data[project_data['project_folP_status'] == 'folP WT']  
    wt_means.append(wt_data['MIC'].mean())  
    wt_stds.append(wt_data['MIC'].std() if len(wt_data) > 1 else 0)  
    wt_counts.append(len(wt_data))  
      
    # F28L data  
    f28l_data = project_data[project_data['project_folP_status'] == 'folP p.F28L']  
    f28l_means.append(f28l_data['MIC'].mean() if len(f28l_data) > 0 else 0)  
    f28l_stds.append(f28l_data['MIC'].std() if len(f28l_data) > 1 else 0)  
    f28l_counts.append(len(f28l_data))  
      
    # P64S data  
    p64s_data = project_data[project_data['project_folP_status'] == 'folP p.P64S']  
    p64s_means.append(p64s_data['MIC'].mean() if len(p64s_data) > 0 else 0)  
    p64s_stds.append(p64s_data['MIC'].std() if len(p64s_data) > 1 else 0)  
    p64s_counts.append(len(p64s_data))  
  
print("\nSummary statistics:")  
for i, project in enumerate(projects):  
    print(f"{project}:")  
    print(f"  WT: n={wt_counts[i]}, mean={wt_means[i]:.1f}")  
    print(f"  F28L: n={f28l_counts[i]}, mean={f28l_means[i]:.1f}")  
    print(f"  P64S: n={p64s_counts[i]}, mean={p64s_means[i]:.1f}")  

    # Create the updated grouped bar chart with adjusted spacing for SMX3  
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)  
  
# Define positions for each project  
x_smx1 = 0  
x_smx2 = 1.5  # Normal spacing  
x_smx3_wt = 3.0  # Start of SMX3 group  
x_smx3_p64s = 3.25  # Close together for SMX3  
  
x_positions = {  
    'SMX1': {'WT': x_smx1 - 0.125, 'F28L': x_smx1 + 0.125, 'P64S': x_smx1 + 0.375},  
    'SMX2': {'WT': x_smx2 - 0.125, 'F28L': x_smx2 + 0.125, 'P64S': x_smx2 + 0.375},  
    'SMX3': {'WT': x_smx3_wt, 'P64S': x_smx3_p64s}  
}  
  
width = 0.25  

# Plot bars for each project individually  
# SMX1 bars  
i = 0  
if wt_counts[i] > 0:  
    ax.bar(x_positions['SMX1']['WT'], wt_means[i], width, yerr=wt_stds[i], capsize=5,   
           color='#87CEEB', alpha=0.8, label='folP WT')  
if f28l_counts[i] > 0:  
    ax.bar(x_positions['SMX1']['F28L'], f28l_means[i], width, yerr=f28l_stds[i], capsize=5,   
           color='#FF6B6B', alpha=0.8, label='folP p.F28L')  
  
# SMX2 bars  
i = 1  
if wt_counts[i] > 0:  
    ax.bar(x_positions['SMX2']['WT'], wt_means[i], width, yerr=wt_stds[i], capsize=5,   
           color='#87CEEB', alpha=0.8)  
if f28l_counts[i] > 0:  
    ax.bar(x_positions['SMX2']['F28L'], f28l_means[i], width, yerr=f28l_stds[i], capsize=5,   
           color='#FF6B6B', alpha=0.8)  
if p64s_counts[i] > 0:  
    ax.bar(x_positions['SMX2']['P64S'], p64s_means[i], width, yerr=p64s_stds[i], capsize=5,   
           color='#4ECDC4', alpha=0.8, label='folP p.P64S')  
  
# SMX3 bars 
i = 2  
if wt_counts[i] > 0:  
    ax.bar(x_positions['SMX3']['WT'], wt_means[i], width, yerr=wt_stds[i], capsize=5,   
           color='#87CEEB', alpha=0.8)  
if p64s_counts[i] > 0:  
    ax.bar(x_positions['SMX3']['P64S'], p64s_means[i], width, yerr=p64s_stds[i], capsize=5,   
           color='#4ECDC4', alpha=0.8)  

    # Overlay individual data points  
for i, project in enumerate(projects):  
    project_data = smx_df[smx_df['Project_code'] == project]  
      
    if project == 'SMX1':  
        # WT points  
        wt_data = project_data[project_data['project_folP_status'] == 'folP WT']  
        if len(wt_data) > 0:  
            x_points = np.random.normal(x_positions['SMX1']['WT'], 0.02, size=len(wt_data))  
            ax.scatter(x_points, wt_data['MIC'], color='black', alpha=0.6, s=30, zorder=3)  
          
        # F28L points  
        f28l_data = project_data[project_data['project_folP_status'] == 'folP p.F28L']  
        if len(f28l_data) > 0:  
            x_points = np.random.normal(x_positions['SMX1']['F28L'], 0.02, size=len(f28l_data))  
            ax.scatter(x_points, f28l_data['MIC'], color='black', alpha=0.6, s=30, zorder=3)  
      
    elif project == 'SMX2':  
        # WT points  
        wt_data = project_data[project_data['project_folP_status'] == 'folP WT']  
        if len(wt_data) > 0:  
            x_points = np.random.normal(x_positions['SMX2']['WT'], 0.02, size=len(wt_data))  
            ax.scatter(x_points, wt_data['MIC'], color='black', alpha=0.6, s=30, zorder=3)  
          
        # F28L points  
        f28l_data = project_data[project_data['project_folP_status'] == 'folP p.F28L']  
        if len(f28l_data) > 0:  
            x_points = np.random.normal(x_positions['SMX2']['F28L'], 0.02, size=len(f28l_data))  
            ax.scatter(x_points, f28l_data['MIC'], color='black', alpha=0.6, s=30, zorder=3)  
          
        # P64S points  
        p64s_data = project_data[project_data['project_folP_status'] == 'folP p.P64S']  
        if len(p64s_data) > 0:  
            x_points = np.random.normal(x_positions['SMX2']['P64S'], 0.02, size=len(p64s_data))  
            ax.scatter(x_points, p64s_data['MIC'], color='black', alpha=0.6, s=30, zorder=3)  
      
    elif project == 'SMX3':  
        # WT points  
        wt_data = project_data[project_data['project_folP_status'] == 'folP WT']  
        if len(wt_data) > 0:  
            x_points = np.random.normal(x_positions['SMX3']['WT'], 0.02, size=len(wt_data))  
            ax.scatter(x_points, wt_data['MIC'], color='black', alpha=0.6, s=30, zorder=3)  
          
        # P64S points  
        p64s_data = project_data[project_data['project_folP_status'] == 'folP p.P64S']  
        if len(p64s_data) > 0:  
            x_points = np.random.normal(x_positions['SMX3']['P64S'], 0.02, size=len(p64s_data))  
            ax.scatter(x_points, p64s_data['MIC'], color='black', alpha=0.6, s=30, zorder=3)  

 # Add sample size annotations
# SMX1
if wt_counts[0] > 0:
    ax.text(x_positions['SMX1']['WT'], wt_means[0] + wt_stds[0] + 10, 'n=' + str(wt_counts[0]), 
            ha='center', va='bottom', fontsize=10, fontweight='bold')
if f28l_counts[0] > 0:
    ax.text(x_positions['SMX1']['F28L'], f28l_means[0] + f28l_stds[0] + 10, 'n=' + str(f28l_counts[0]), 
            ha='center', va='bottom', fontsize=10, fontweight='bold')

# SMX2
if wt_counts[1] > 0:
    ax.text(x_positions['SMX2']['WT'], wt_means[1] + wt_stds[1] + 10, 'n=' + str(wt_counts[1]), 
            ha='center', va='bottom', fontsize=10, fontweight='bold')
if f28l_counts[1] > 0:
    ax.text(x_positions['SMX2']['F28L'], f28l_means[1] + f28l_stds[1] + 10, 'n=' + str(f28l_counts[1]), 
            ha='center', va='bottom', fontsize=10, fontweight='bold')
if p64s_counts[1] > 0:
    ax.text(x_positions['SMX2']['P64S'], p64s_means[1] + p64s_stds[1] + 10, 'n=' + str(p64s_counts[1]), 
            ha='center', va='bottom', fontsize=10, fontweight='bold')

# SMX3
if wt_counts[2] > 0:
    ax.text(x_positions['SMX3']['WT'], wt_means[2] + wt_stds[2] + 10, 'n=' + str(wt_counts[2]), 
            ha='center', va='bottom', fontsize=10, fontweight='bold')
if p64s_counts[2] > 0:
    ax.text(x_positions['SMX3']['P64S'], p64s_means[2] + p64s_stds[2] + 10, 'n=' + str(p64s_counts[2]), 
            ha='center', va='bottom', fontsize=10, fontweight='bold')

# Formatting
ax.set_xlabel('Project', fontsize=14)
ax.set_ylabel('MIC (ug/mL)', fontsize=14)
ax.set_title(r'Mean MIC by $\boldsymbol{folP}$ Mutation Status', fontsize=16, fontweight='bold')

# Set custom x-tick positions and labels
ax.set_xticks([x_smx1, x_smx2, (x_smx3_wt + x_smx3_p64s)/2])
ax.set_xticklabels(['SMX1', 'SMX2', 'SMX3'], fontsize=12)

import matplotlib.patches as mpatches
legend_elements = [
    mpatches.Patch(color='#87CEEB', label='$\it{folP}$ WT'),
    mpatches.Patch(color='#FF6B6B', label='$\it{folP}$ p.F28L'),
    mpatches.Patch(color='#4ECDC4', label='$\it{folP}$ p.P64S')
]
ax.legend(handles=legend_elements,loc='upper left', frameon=True, fancybox=True, shadow=True, fontsize=11)
ax.grid(True, alpha=0.3, axis='y')

# Set y-axis limit
max_val = max([max(wt_means), max(f28l_means), max(p64s_means)]) + max([max(wt_stds), max(f28l_stds), max(p64s_stds)])
ax.set_ylim(0, max_val + 50)

plt.tight_layout()
plt.savefig('/Users/jasonklemm/Downloads/publication_ready_SMX_folP_final_adjusted.png', dpi=300, bbox_inches='tight')
plt.show()

print('Saved updated chart with SMX3 bars side-by-side as publication_ready_SMX_folP_final_adjusted.png')
