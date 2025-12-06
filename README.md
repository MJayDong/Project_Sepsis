# Project_Sepsis
A reproducible, modular single-cell analysis framework for dissecting cellular dynamics, regulatory networks, and intercellular communication in sepsis-induced cardiomyopathy (SICM).

Temporal Single-Cell Landscape and Regulatory Network Remodeling in Sepsis-Induced Cardiac Dysfunction

This repository provides a fully modular and reproducible analysis framework for our single-cell study of sepsis-induced cardiomyopathy (SICM).
The project integrates multi-layer computational workflows, including preprocessing, cell-type–resolved analysis, trajectory inference, gene regulatory network modeling, cell–cell communication analysis, and interpretable machine learning.
All scripts are organized into independent, plug-and-play modules for transparent, review-friendly, and easily replicable computational research.

This repository contains:
	•	Preprocessing pipelines for scRNA-seq data
	•	Cell-type–specific analyses covering immune, cardiomyocyte, endothelial, mesenchymal, and niche populations
	•	Cross-cell analysis, including trajectory inference, Milo/OR differential neighborhood testing, GeneSwitches pseudo-R² gene switching
	•	Network inference modules using SCENIC and hdWGCNA
	•	Intercellular communication modeling via CellChat
	•	Machine-learning interpretation using SHAP values
	•	Figure-generation scripts used in the manuscript

The repository is designed for clarity, reproducibility, and ease of navigation, enabling researchers and reviewers to follow each step from raw data processing to final biological interpretation.
