# Load the MRP1 metadata file
import pandas as pd

mrp1_df_old = pd.read_csv('/Users/jasonklemm/Downloads/MRP1_metadata_with_mutations (1).csv')

# Show the shape, columns, and first few rows to understand the data
print('Shape:', mrp1_df_old.shape)
print('Columns:', mrp1_df_old.columns.tolist())
print(mrp1_df_old.head())

# Create the updated chart with explicit WT labels
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

bar_colors = [
    '#87CEEB',  # WT (light blue)
    '#FF6B6B',  # mrdA p.T547I (red)
    '#4ECDC4',  # mrdA p.T547I + ftsI p.A513S (teal)
    '#FFD166',  # mrdA p.T547I, p.G505R + ftsI p.A513S (yellow)
    '#6A4C93'   # mrdA p.T547I, p.G505R + ftsI p.A513S, p.R167L (purple)
]

combos = [
    (0, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 1, 0, 1),
    (1, 1, 0, 1),
    (1, 1, 1, 1)
]

# Updated legend labels to be more explicit about WT status
legend_labels = [
    '$\it{mrdA}$ WT + $\it{ftsI}$ WT',
    '$\it{mrdA}$ p.T547I + $\it{ftsI}$ WT',
    '$\it{mrdA}$ p.T547I + $\it{ftsI}$ p.A513S',
    '$\it{mrdA}$ p.T547I, p.G505R + $\it{ftsI}$ p.A513S',
    '$\it{mrdA}$ p.T547I, p.G505R + $\it{ftsI}$ p.A513S, p.R167L'
]

means = []
stds = []
counts = []
for combo in combos:
    subset = mrp1_df_old[
        (mrp1_df_old['mrdA_499237'] == combo[0]) &
        (mrp1_df_old['mrdA_499364'] == combo[1]) &
        (mrp1_df_old['ftsI_24331'] == combo[2]) &
        (mrp1_df_old['ftsI_25368'] == combo[3])
    ]
    means.append(subset['MIC'].mean())
    stds.append(subset['MIC'].std())
    counts.append(len(subset))

print('Group counts:', counts)
print('Group means:', [round(m, 3) for m in means])

# Keep the chart size matched to the SMX chart for consistency
fig, ax = plt.subplots(figsize=(10, 7), dpi=300)
x = np.arange(len(legend_labels))
width = 0.6

bars = ax.bar(x, means, width, yerr=stds, capsize=5, 
              color=bar_colors, alpha=0.85)

# Overlay individual data points for each group
for i, combo in enumerate(combos):
    subset = mrp1_df_old[
        (mrp1_df_old['mrdA_499237'] == combo[0]) &
        (mrp1_df_old['mrdA_499364'] == combo[1]) &
        (mrp1_df_old['ftsI_24331'] == combo[2]) &
        (mrp1_df_old['ftsI_25368'] == combo[3])
    ]
    y = subset['MIC']
    x_points = np.random.normal(i, 0.04, size=len(y))
    ax.scatter(x_points, y, color='black', alpha=0.7, s=32, zorder=3)

# Annotate sample sizes above each bar
for i, n in enumerate(counts):
    ax.text(i, means[i] + (stds[i] if not np.isnan(stds[i]) else 0) + 0.02, 
            'n=' + str(n), ha='center', va='bottom', fontsize=10, fontweight='bold')

# Descriptive legend with explicit WT status
legend_handles = [mpatches.Patch(color=bar_colors[i], label=legend_labels[i]) for i in range(len(legend_labels))]
ax.legend(handles=legend_handles, loc='upper left', frameon=True, fancybox=True, shadow=True, 
          fontsize=10, borderpad=0.4, labelspacing=0.4, handlelength=1.3, handletextpad=0.5)

ax.set_xlabel('Project', fontsize=14)
ax.set_ylabel('MIC (ug/mL)', fontsize=14)
ax.set_title(r'Mean MIC by $\boldsymbol{ftsI}$ and $\boldsymbol{mrdA}$ Mutation Status', fontsize=16, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(['MRP1'] * len(legend_labels), fontsize=12)
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim(0, 0.8)

plt.tight_layout()
plt.savefig('/Users/jasonklemm/Desktop/MRP1_MIC_barplot_explicit_WT_labels.png', dpi=300, bbox_inches='tight')
plt.show()

print('Saved chart with explicit WT labels as MRP1_MIC_barplot_explicit_WT_labels.png')