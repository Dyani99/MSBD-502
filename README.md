# MSBD-502
Computational Biology Final Project 
Computational Assessment of the Structural Stability and Flexibility of Human TNNT3 Using Molecular Dynamics Simulations
Course: MSBD 502 — Computational Structural Biology
Project scope: Wild-type human TNNT3 only (no variant comparison)
Method: Explicit-solvent molecular dynamics (MD) with OpenMM; exploratory RMSD and RMSF analyses
Overview
This project explores the short-timescale structural behavior of wild-type human fast skeletal troponin T (TNNT3). An AlphaFold-predicted protein structure was prepared, solvated, energy-minimized, and simulated in OpenMM. Backbone root-mean-square deviation (RMSD) and Cα root-mean-square fluctuation (RMSF) were calculated to characterize conformational changes and residue-level motion.

Interpretation boundary: The ~50 ps trajectory is a preliminary feasibility/relaxation experiment. It does not establish thermodynamic equilibration, long-term structural stability, or the effect of any TNNT3 mutation.
Starting structure and simulation settings
Item	Recorded setting
Structure	AlphaFold model AF-P45378-2-F1-model_v6.pdb (human TNNT3 isoform 2)
Prepared structure	TNNT3_WT_small_minimized.pdb
MD engine	OpenMM, executed in Google Colab
Platform	OpenCL on NVIDIA Tesla T4; mixed precision
Protein force field	amber19-all.xml
Water model	OPC (amber19/opc.xml)
Solvent preparation	Dodecahedral periodic box; 0.8 nm padding; 0.15 M ionic strength
Nonbonded interactions	PME; 0.8 nm cutoff
Constraints	Hydrogen bonds constrained; rigid water
Integrator	Langevin Middle; 298 K; 5 ps⁻¹ friction; 0.0005 ps (0.5 fs) timestep in the successful run
Ensemble	NVT (fixed volume and target temperature)
Total duration	Approximately 50.05 ps
Solvated atom count	470,389


The initial prepared/minimized structure and the exact OpenMM system/topology/checkpoint files should be archived with the notebook. Do not regenerate the starting structure and assume the coordinates or atom ordering are identical to the original run.
Preliminary results
- Final temperature: 298.18 K.
- Final potential energy: −6,282,298.54 kJ/mol (whole solvated system).
- Final backbone RMSD relative to the minimized structure: 0.1911 nm (1.911 Å).
- During the final 34 ps extension, 68 frames were saved at 0.5 ps intervals; mean backbone RMSD relative to the first extension frame: 0.920 Å; final: 1.231 Å.
- Preliminary Cα RMSF: mean 0.543 Å; maximum 1.218 Å at residue 1; residue 63 0.380 Å. These values were computed against the first frame and should be distinguished from mean-centered RMSF.
- Minimum protein-heavy-atom separation from a neighboring periodic image at one snapshot: 3.479 nm (not a guarantee for all frames).
The RMSD reference differs between the full-run endpoint and extension-trajectory analyses; the two values are not directly interchangeable. The RMSD trend still showed drift over this short interval.
Suggested repository layout
TNNT3_WT_MD/
├── README.md
├── requirements.txt                 # Export and pin versions from the working Colab environment
├── notebooks/
│   └── TNNT3_WT_OpenMM.ipynb         # Add the actual, cleaned simulation/analysis notebook
├── data/
│   ├── README.md                     # Structure source, download link, and file provenance
│   └── TNNT3_WT_small_minimized.pdb  # Include if permitted and size allows
├── results/
│   ├── TNNT3_WT_equilibration.csv
│   ├── TNNT3_WT_equilibration.dcd
│   ├── rmsd.png
│   └── rmsf.png
└── report.pdf
The directory tree above is a proposed submission structure, not a claim that all listed files are present in this download. Copy your original notebook, trajectory, analysis plots, and report into the corresponding folders. Large trajectories/checkpoints may be provided through an accessible archive with documented download instructions if GitHub size limits apply.
Reproducing the analysis
1. Open the original simulation notebook in Google Colab with a GPU runtime. Keep the original minimized PDB, system XML, and checkpoint together; the checkpoint is tied to its corresponding system and environment.
2. Install the versions used in the successful environment. In Colab, record them before pinning requirements.txt:
   import importlib.metadata as md
   for package in ["openmm", "mdtraj", "numpy", "scipy", "pandas", "matplotlib"]:
       try:
           print(package, md.version(package))
       except md.PackageNotFoundError:
           print(package, "not installed")
3. To reproduce the reported results exactly, use the original trajectory and topology, rather than rerunning stochastic MD. The DCD trajectory TNNT3_WT_equilibration.dcd is paired with TNNT3_WT_small_minimized.pdb; verify that their atom counts and ordering match.
4. For RMSD, load the DCD with MDTraj, select protein backbone atoms, align frames, and calculate RMSD relative to the first saved frame. Convert nanometers to ångströms by multiplying by 10.
5. For the reported preliminary RMSF, align frames using the backbone and calculate Cα fluctuations relative to the first saved frame. If calculating conventional mean-centered RMSF, label it separately and do not substitute it for the reported preliminary values.
6. Regenerate RMSD/RMSF plots with Matplotlib, including labeled axes, units, and figure captions. Use the original OpenMM CSV for temperature and energy plots.
Minimal analysis example
import mdtraj as md
import numpy as np

traj = md.load("results/TNNT3_WT_equilibration.dcd",
               top="data/TNNT3_WT_small_minimized.pdb")
backbone = traj.topology.select("protein and backbone")
ca = traj.topology.select("protein and name CA")
traj.superpose(traj, 0, atom_indices=backbone)
rmsd_A = 10 * md.rmsd(traj, traj, 0, atom_indices=backbone)
rmsf_first_frame_A = 10 * md.rmsf(traj, traj, 0, atom_indices=ca)
print("Frames:", traj.n_frames)
print("Mean RMSD (Å):", rmsd_A.mean())
print("Final RMSD (Å):", rmsd_A[-1])
print("Mean Cα RMSF (Å):", rmsf_first_frame_A.mean())
Reproducibility notes and limitations
- This is a single short trajectory, not a replicated or converged MD study. RMSD/RMSF estimates are sensitive to alignment, reference choice, periodic boundary handling, and sampling time.
- The 0.8 nm solvent padding and periodic-image separation were screened at one time point, not validated across an equilibrated trajectory.
- The isolated predicted TNNT3 structure does not represent the full troponin complex. AlphaFold confidence and isoform/sequence provenance should be checked before biological conclusions.
- Stochastic dynamics may differ on rerun. Record the random seed if explicitly set; otherwise state that it was not controlled.
- Do not upload an incomplete requirements.txt as if it captured the original environment. Pin actual package versions from the notebook runtime.
- Confirm that any included PDB, DCD, or other data files may be redistributed and document their sources.
References and source attribution
- Jumper, J., et al. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596, 583–589. https://doi.org/10.1038/s41586-021-03819-2
- Wei, B., & Jin, J.-P. (2016). TNNT1, TNNT2, and TNNT3: Isoform genes, regulation, and structure-function relationships. Gene, 582(1), 1–13. https://doi.org/10.1016/j.gene.2016.01.006
- AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk/
- OpenMM documentation: https://docs.openmm.org/
Submission checklist
- [ ] Add cleaned and commented Colab notebook(s)
- [ ] Add requirements.txt with actual tested versions
- [ ] Add structure provenance and raw-data/download instructions
- [ ] Add trajectory, analysis CSV, and Python-generated figures or links to large files
- [ ] Add the final PDF report
- [ ] Confirm that the repository reproduces reported plots and statistics
- [ ] Share the jointly owned GitHub repository with the instructor before the deadline
Course submission note: The assignment calls for a shared GitHub repository containing notebooks, a README, requirements, and relevant data, plus an individual written report and a group presentation.
