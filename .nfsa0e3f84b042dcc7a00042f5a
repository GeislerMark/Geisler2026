#!/usr/bin/env python3
"""Comprehensive genomic bin visualization"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

filepath = "/path/to/input/multiBigWigMergeBinCounts.tsv"

print("="*60)
print("GENOMIC BIN VISUALIZATION")
print("="*60)
print("\nLoading data...")
df = pd.read_csv(filepath, sep='\t')
df.columns = df.columns.str.strip("'#")
df['midpoint'] = (df['start'] + df['end']) / 2

# Order chromosomes (chrHis last)
chrom_order = sorted([c for c in df['chr'].unique() if c != 'chrHis'])
if 'chrHis' in df['chr'].unique():
    chrom_order.append('chrHis')

df['chr'] = pd.Categorical(df['chr'], categories=chrom_order, ordered=True)
df = df.sort_values(['chr', 'start'])

chr_lengths = df.groupby('chr', observed=True)['end'].max()
chr_cumsum = chr_lengths.cumsum().shift(fill_value=0)

# Convert to numeric explicitly before adding
chr_offset = df['chr'].map(chr_cumsum).astype(float)
df['cumulative_pos'] = chr_offset + df['midpoint']

sample_cols = [col for col in df.columns if col.endswith('.bw')]
df['mean_signal'] = df[sample_cols].mean(axis=1)

# Calculate condition means for MA plot
mute_cols = [col for col in sample_cols if 'mute' in col]
mxc_cols = [col for col in sample_cols if 'mxc' in col]

df['mute_mean'] = df[mute_cols].mean(axis=1)
df['mxc_mean'] = df[mxc_cols].mean(axis=1)

print(f"Loaded: {len(df)} bins, {len(chrom_order)} chromosomes, {len(sample_cols)} samples")
print(f"Chromosomes: {', '.join(chrom_order)}")

chr_centers = [df[df['chr'] == chrom]['cumulative_pos'].mean() for chrom in df['chr'].cat.categories]

# PLOT 1: MA Plot (Scatter plot: Mute vs MXC)
print("\n1. Creating MA plot (Mute vs MXC comparison)...")

fig, ax = plt.subplots(figsize=(10, 10))

# Separate chrHis from other chromosomes
is_chrhis = df['chr'] == 'chrHis'
other_data = df[~is_chrhis]
chrhis_data = df[is_chrhis]

# Plot non-chrHis points first (background)
ax.scatter(other_data['mute_mean'], other_data['mxc_mean'], 
          c='#CCCCCC', s=20, alpha=0.4, label='Other chromosomes', edgecolors='none', 
          rasterized=True)

# Plot chrHis points on top (highlighted)
ax.scatter(chrhis_data['mute_mean'], chrhis_data['mxc_mean'], 
          c='#FF6B6B', s=30, alpha=0.8, label='chrHis', edgecolors='black', linewidth=0.5, 
          rasterized=True)

# Add 45-degree reference line (y=x, equal signal)
max_val = max(df['mute_mean'].max(), df['mxc_mean'].max())
min_val = min(df['mute_mean'].min(), df['mxc_mean'].min())
ax.plot([min_val, max_val], [min_val, max_val], 
       'k--', linewidth=1.5, alpha=0.5, label='Equal signal (y=x)')

# Scale axes to emphasize chrHis in upper right
# Get percentiles to set reasonable limits while keeping chrHis visible
mute_95 = df['mute_mean'].quantile(0.95)
mxc_95 = df['mxc_mean'].quantile(0.95)

# If chrHis mean is much higher, extend limits to include it
if len(chrhis_data) > 0:
    chrhis_mute_max = chrhis_data['mute_mean'].max()
    chrhis_mxc_max = chrhis_data['mxc_mean'].max()
    x_max = max(mute_95, chrhis_mute_max) * 1.1
    y_max = max(mxc_95, chrhis_mxc_max) * 1.1
else:
    x_max = mute_95 * 1.1
    y_max = mxc_95 * 1.1

ax.set_xlim(0, x_max)
ax.set_ylim(0, y_max)

# Label top 5 bins by total signal
df['total_signal'] = df['mute_mean'] + df['mxc_mean']
top5 = df.nlargest(5, 'total_signal')

for _, row in top5.iterrows():
    label = f"{row['chr']}:{int(row['start'])}-{int(row['end'])}"
    ax.annotate(label, 
               xy=(row['mute_mean'], row['mxc_mean']),
               xytext=(8, 8), textcoords='offset points',
               fontsize=8, 
               bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8, 
                        edgecolor='black', linewidth=0.5),
               arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2', lw=0.5))

ax.set_xlabel('Mute Signal (CPM)', fontsize=12, fontweight='bold')
ax.set_ylabel('MXC Signal (CPM)', fontsize=12, fontweight='bold')
ax.set_title('Mute vs MXC Signal Comparison', fontsize=14, fontweight='bold', pad=20)
ax.legend(loc='upper left', fontsize=10, frameon=True, fancybox=True)
ax.grid(True, alpha=0.3)
ax.set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.savefig('MA_plot_mute_vs_mxc.pdf', dpi=300, bbox_inches='tight')
print("   Saved: MA_plot_mute_vs_mxc.pdf")
plt.savefig('MA_plot_mute_vs_mxc.svg', bbox_inches='tight')
print("   Saved: MA_plot_mute_vs_mxc.svg")
plt.close()

# PLOT 2: Manhattan by Condition (with top 5 labeled)
print("\n2. Creating condition comparison with top 5 peaks labeled...")
conditions = {
    'IgG Control': [col for col in sample_cols if 'IgG' in col],
    'Mute': [col for col in sample_cols if 'mute' in col],
    'MXC': [col for col in sample_cols if 'mxc' in col]
}

fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
for idx, (condition, samples) in enumerate(conditions.items()):
    ax = axes[idx]
    if not samples:
        continue
    df['cond_mean'] = df[samples].mean(axis=1)
    
    # Plot points
    for i, chrom in enumerate(df['chr'].cat.categories):
        chrom_data = df[df['chr'] == chrom]
        color = '#FF6B6B' if chrom == 'chrHis' else ('#4ECDC4' if i % 2 == 0 else '#45B7D1')
        ax.scatter(chrom_data['cumulative_pos'], chrom_data['cond_mean'],
                  c=color, s=8, alpha=0.6, rasterized=True)
    
    # Find and label top 5 peaks
    top5 = df.nlargest(5, 'cond_mean')
    for _, row in top5.iterrows():
        label = f"{row['chr']}:{int(row['start'])}-{int(row['end'])}"
        ax.annotate(label, 
                   xy=(row['cumulative_pos'], row['cond_mean']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, 
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7, edgecolor='black', linewidth=0.5),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=0.5))
    
    # Labels only on bottom
    if idx == 2:
        ax.set_xticks(chr_centers)
        ax.set_xticklabels(df['chr'].cat.categories, rotation=45, ha='right')
        ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    
    # Chromosome boundaries
    for pos in chr_cumsum.values[1:]:
        ax.axvline(x=pos, color='gray', linestyle='--', alpha=0.3, linewidth=0.5)
    
    ax.set_ylabel('Signal (CPM)', fontsize=11, fontweight='bold')
    ax.set_title(f'{condition}', fontsize=12, fontweight='bold', loc='left')
    ax.grid(True, alpha=0.3)

fig.suptitle('Genome-wide Signal Distribution by Condition', 
            fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout()

plt.savefig('manhattan_by_condition.pdf', dpi=300, bbox_inches='tight')
print("   Saved: manhattan_by_condition.pdf")
plt.savefig('manhattan_by_condition.svg', bbox_inches='tight')
print("   Saved: manhattan_by_condition.svg")
plt.close()

# PLOT 3: chrHis Focus
if 'chrHis' in df['chr'].cat.categories:
    print("3. Creating chrHis-focused plots...")
    chrhis = df[df['chr'] == 'chrHis'].copy()
    other = df[df['chr'] != 'chrHis'].copy()
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Detailed tracks
    ax = axes[0]
    for sample in sample_cols:
        name = sample.replace('_Merged_chrHis_div100.bw', '').replace('_', ' ')
        ax.plot(chrhis['midpoint'], chrhis[sample], label=name, linewidth=1.5, alpha=0.7)
    ax.set_xlabel('Position on chrHis (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Signal (CPM)', fontsize=12, fontweight='bold')
    ax.set_title('Detailed Signal on chrHis', fontsize=13, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Distribution
    ax = axes[1]
    bins = np.linspace(0, max(df['mean_signal'].max(), 5), 50)
    ax.hist(other['mean_signal'], bins=bins, alpha=0.6, label='Rest of Genome', 
           color='#4ECDC4', edgecolor='black')
    ax.hist(chrhis['mean_signal'], bins=bins, alpha=0.6, label='chrHis', 
           color='#FF6B6B', edgecolor='black')
    
    genome_mean = other['mean_signal'].mean()
    chrhis_mean = chrhis['mean_signal'].mean()
    ax.axvline(genome_mean, color='#4ECDC4', linestyle='--', linewidth=2,
              label=f'Genome mean: {genome_mean:.2f}')
    ax.axvline(chrhis_mean, color='#FF6B6B', linestyle='--', linewidth=2,
              label=f'chrHis mean: {chrhis_mean:.2f}')
    
    ax.set_xlabel('Mean Signal (CPM)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title('Signal Distribution Comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('chrHis_focus.pdf', dpi=300, bbox_inches='tight')
    print("   Saved: chrHis_focus.pdf")
    plt.savefig('chrHis_focus.svg', bbox_inches='tight')
    print("   Saved: chrHis_focus.svg")
    plt.close()

# PLOT 4: Correlation
print("4. Creating correlation heatmap...")
corr = df[sample_cols].corr()
names = [col.replace('_Merged_chrHis_div100.bw', '').replace('_', ' ') for col in sample_cols]
corr.index = names
corr.columns = names

fig, ax = plt.subplots(figsize=(10, 9))
im = ax.imshow(corr, cmap='RdYlBu_r', vmin=0, vmax=1, aspect='auto')
ax.set_xticks(range(len(names)))
ax.set_yticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha='right')
ax.set_yticklabels(names)

for i in range(len(names)):
    for j in range(len(names)):
        ax.text(j, i, f'{corr.iloc[i, j]:.2f}', ha="center", va="center", 
               color="black", fontsize=9)

cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Correlation', rotation=270, labelpad=20, fontweight='bold')
ax.set_title('Sample Correlation', fontsize=14, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig('sample_correlation.pdf', dpi=300, bbox_inches='tight')
print("   Saved: sample_correlation.pdf")
plt.savefig('sample_correlation.svg', bbox_inches='tight')
print("   Saved: sample_correlation.svg")
plt.close()

# SUMMARY
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
for condition, samples in conditions.items():
    if samples:
        mean_val = df[samples].mean(axis=1).mean()
        print(f"{condition:15s}: {mean_val:.3f} CPM")

if 'chrHis' in df['chr'].cat.categories:
    enrichment = chrhis_mean / genome_mean
    print(f"\nchrHis enrichment: {enrichment:.2f}x")

# MA plot statistics
print("\nMute vs MXC comparison:")
print(f"  Pearson correlation: {df['mute_mean'].corr(df['mxc_mean']):.3f}")
if len(chrhis_data) > 0:
    print(f"  chrHis - Mute mean: {chrhis_data['mute_mean'].mean():.2f} CPM")
    print(f"  chrHis - MXC mean: {chrhis_data['mxc_mean'].mean():.2f} CPM")

print("\n" + "="*60)
print("OUTPUT FILES (vector format for Illustrator):")
print("="*60)
print("  MA_plot_mute_vs_mxc.pdf/.svg")
print("  manhattan_by_condition.pdf/.svg")
print("  chrHis_focus.pdf/.svg")
print("  sample_correlation.pdf/.svg")
print("\n✓ Complete!")
