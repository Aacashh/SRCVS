# Results and Discussion

## Molecular Dynamics Simulation and Binding Affinity Analysis

To investigate the dynamic behavior and binding stability of the top two lead compounds (44259_01 and 72307_01) within the target protein binding pocket, we performed 2.0 ns molecular dynamics (MD) simulations in an implicit solvent (OBC2) environment. The structural stability was assessed through Root Mean Square Deviation (RMSD) and Root Mean Square Fluctuation (RMSF) analysis, while the thermodynamic binding energy was estimated using the MM-GBSA method. 

### Structural Stability and Flexibility (RMSD & RMSF)

The structural integrity of the protein-ligand complexes was evaluated by tracking the RMSD of the protein backbone and the heavy atoms of the ligands over the trajectory (Figure X, `09_combined_rmsd_panel.png`). 

During the production phase (1.0 - 2.0 ns), the **44259_01** complex demonstrated superior structural stability. The protein backbone RMSD for 44259_01 equilibrated at a mean of **3.72 ± 0.43 Å**, which was noticeably lower and less variable than the **4.73 ± 0.39 Å** observed for 72307_01. This suggests that the binding of 44259_01 induces a more stable overall conformation of the target protein.

More importantly, the positional stability of the ligand itself within the binding pocket strongly favors 44259_01. The ligand RMSD for 44259_01 during the production phase was **8.67 ± 0.82 Å**, whereas 72307_01 exhibited a higher deviation and larger fluctuations with a mean RMSD of **10.86 ± 0.64 Å**. The increased translational and rotational movement of 72307_01 indicates a looser, less optimized fit within the binding site compared to 44259_01.

Residue-level flexibility, as measured by RMSF, further supports this conclusion (Figure X, `05_rmsf_difference.png`). The differential RMSF profile reveals that the 72307_01 complex exhibits higher flexibility across several key loop regions and binding site residues, indicating that 44259_01 more effectively anchors the surrounding protein structure, likely through stronger or more persistent intermolecular interactions.

### Thermodynamic Binding Energy (MM-GBSA)

To quantify the energetic favorability of the interactions, MM-GBSA binding free energies ($\Delta G_{bind}$) were calculated from trajectory frames.

Compound **44259_01** exhibited a highly favorable mean binding energy of **-48.19 ± 7.21 kcal/mol**. In stark contrast, **72307_01** displayed a much weaker binding affinity, with a mean $\Delta G_{bind}$ of **-22.45 ± 3.51 kcal/mol** (Figure X, `08_mmgbsa_boxplot.png`). 

This represents a significant thermodynamic advantage for 44259_01, with a mean binding free energy difference ($\Delta\Delta G$) of **25.74 kcal/mol**. A Welch's t-test confirms that this difference is highly statistically significant ($t = -9.635, p < 0.0001$). The substantial gap in binding energy suggests that 44259_01 forms a more robust and complementary interaction network with the target protein.

### Correlation between Stability and Affinity

The relationship between ligand positional stability (RMSD) and binding affinity (MM-GBSA) reveals a compelling negative correlation (Figure X, `10_energy_vs_rmsd.png`). As the ligand RMSD increases (indicating a loss of stable binding pose), the binding energy becomes less favorable (less negative). 

Compound 44259_01 occupies the optimal region of this phase space: lower RMSD and highly negative binding energy. Conversely, 72307_01 is shifted towards higher RMSD and less favorable binding energy. This analysis directly correlates the superior structural stability of 44259_01 to its vastly improved thermodynamic binding profile, validating it as the definitively superior lead compound.
