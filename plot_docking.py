import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw

# Define paths
HTVS_PATH = "/home/ankit/Thesis/SRCVS/out/03_dock/htvs_scores.csv"
XP_PATH = "/home/ankit/Thesis/SRCVS/out/03_dock/xp_scores.csv"
LEADS_PATH = "/home/ankit/Thesis/SRCVS/out/05_select/final_leads.csv"
OUT_DIR = "/home/ankit/Thesis/SRCVS/out/thesis_plots"

os.makedirs(OUT_DIR, exist_ok=True)

# 1. Docking Score Distribution (HTVS vs XP)
print("Plotting Docking Score Distributions...")
htvs_df = pd.read_csv(HTVS_PATH)
xp_df = pd.read_csv(XP_PATH)
leads_df = pd.read_csv(LEADS_PATH)

plt.figure(figsize=(10, 6))
sns.kdeplot(data=htvs_df['score'], fill=True, label='HTVS (Initial Screen)', color='gray', alpha=0.5)
sns.kdeplot(data=xp_df['score'], fill=True, label='XP (Extra Precision)', color='teal', alpha=0.7)

# Add vertical lines for final leads
colors = ['#d73027', '#4575b4'] # Red for 44259, Blue for 72307
for i, row in leads_df.iterrows():
    c = '#35978f' if '44259' in row['name'] else '#f46d43'
    plt.axvline(x=row['score'], color=c, linestyle='--', linewidth=2, 
                label=f"Lead: {row['name']} ({row['score']:.2f})")

plt.xlabel('Docking Score (kcal/mol)', fontsize=14)
plt.ylabel('Density', fontsize=14)
plt.title('Distribution of Docking Scores Across Screening Phases', fontsize=16)
plt.legend(fontsize=12)
plt.grid(axis='x', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "11_docking_score_distribution.png"), dpi=300)
plt.close()

# 2. 2D Chemical Structures
print("Generating 2D Chemical Structures...")
mols = []
legends = []
for i, row in leads_df.iterrows():
    mol = Chem.MolFromSmiles(row['smiles'])
    if mol:
        mols.append(mol)
        legends.append(f"{row['name']}\nScore: {row['score']:.2f} kcal/mol")

if mols:
    img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(500, 500), legends=legends, useSVG=False)
    img.save(os.path.join(OUT_DIR, "12_lead_structures_2D.png"))
    
# 3. Admet radar chart or table
print("Generating ADMET profile comparison...")
# Normalizing values for a radar chart
# mw (ideal < 500), logp (ideal < 5), hbd (<5), hba (<10), tpsa (<140)
categories = ['MW/100', 'LogP', 'HBD', 'HBA', 'TPSA/20', 'LogBB + 1']
import numpy as np

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, polar=True)

angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
angles += angles[:1]

for i, row in leads_df.iterrows():
    values = [
        row['mw']/100,
        row['logp'],
        row['hbd'],
        row['hba']/2, # scale hba to be comparable
        row['tpsa']/20,
        row['logbb'] + 1 # shift logbb to be positive for visualization
    ]
    values += values[:1]
    
    c = '#35978f' if '44259' in row['name'] else '#f46d43'
    ax.plot(angles, values, linewidth=2, label=row['name'], color=c)
    ax.fill(angles, values, alpha=0.25, color=c)

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.set_thetagrids(np.degrees(angles[:-1]), categories, fontsize=12)
plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
plt.title('ADMET Property Radar Chart', y=1.1, fontsize=16)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "13_admet_radar_chart.png"), dpi=300)
plt.close()

# 4. Virtual Screening Funnel
print("Generating Virtual Screening Funnel...")
try:
    htvs_count = len(pd.read_csv(HTVS_PATH))
    sp_count = len(pd.read_csv("/home/ankit/Thesis/SRCVS/out/03_dock/sp_scores.csv"))
    xp_count = len(pd.read_csv(XP_PATH))
    admet_count = len(pd.read_csv("/home/ankit/Thesis/SRCVS/out/04_admet/admet_passed.csv"))
    leads_count = len(leads_df)
    
    stages = ['HTVS', 'SP', 'XP', 'ADMET\nPassed', 'Final\nLeads']
    counts = [htvs_count, sp_count, xp_count, admet_count, leads_count]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    x = range(len(stages))
    bar_colors = ['#cccccc', '#969696', '#636363', '#252525', '#e05c5c']
    bars = ax.bar(x, counts, color=bar_colors, width=0.6)
    
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval + (max(counts)*0.02), 
                int(yval), ha='center', va='bottom', fontsize=12, fontweight='bold')
        
    ax.plot(x, counts, color='black', linestyle='--', marker='o', alpha=0.5)
    
    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=12)
    ax.set_ylabel('Number of Compounds', fontsize=14)
    ax.set_title('Virtual Screening Pipeline Attrition', fontsize=16)
    ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "14_virtual_screening_funnel.png"), dpi=300)
    plt.close()
except Exception as e:
    print(f"Could not generate funnel plot: {e}")

# 5. XP Candidates Docking Scores
print("Generating XP Candidates Docking Scores Bar Chart...")
try:
    admet_df = pd.read_csv("/home/ankit/Thesis/SRCVS/out/04_admet/admet_descriptors.csv")
    admet_df = admet_df.sort_values(by='score', ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    names = admet_df['name'].tolist()
    scores = admet_df['score'].tolist()
    
    lead_names = leads_df['name'].tolist()
    colors_list = []
    for name in names:
        if name in lead_names:
            c = '#35978f' if '44259' in name else '#f46d43'
            colors_list.append(c)
        else:
            colors_list.append('#cccccc')
            
    bars = ax.bar(range(len(names)), scores, color=colors_list)
    
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha='right')
    ax.set_ylabel('Docking Score (kcal/mol)', fontsize=14)
    ax.set_title('Top 15 XP Candidates Docking Scores', fontsize=16)
    
    import matplotlib.patches as mpatches
    legend_elements = [
        mpatches.Patch(facecolor='#cccccc', label='Filtered Candidates'),
        mpatches.Patch(facecolor='#35978f', label='Lead 44259_01'),
        mpatches.Patch(facecolor='#f46d43', label='Lead 72307_01')
    ]
    ax.legend(handles=legend_elements, loc='best')
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "15_xp_candidates_scores.png"), dpi=300)
    plt.close()
except Exception as e:
    print(f"Could not generate XP candidates bar chart: {e}")

print("Done generating intuitive plots.")
