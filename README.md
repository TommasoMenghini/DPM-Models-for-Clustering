# DPM-Models-for-Clustering

This repository provides detailed materials for implementing the work my friend (and colleague) Tommaso Pozzi and I conducted for our Non-Parametric Statistics exam project.

We studied Dirichlet Process Mixture (DPM) models for clustering, taking our first steps from two key articles:

[Wade and Ghahramani (2018). Bayesian cluster analysis: Point estimation and credible balls]
[Dahl, Johnson, and Müller (2022). Search algorithms and loss functions for Bayesian clustering]
Our exploration focused on two different loss functions: Binder's loss and the Variation of Information loss.

After a brief digression into lattice theory, we compared the performance of two algorithms for Bayesian clustering optimization via decision theory:

The greedy search algorithm proposed by Wade and Ghahramani
SALSO
To evaluate these algorithms, we conducted a simulation study. Lastly, we applied our findings to a real-world context using the Country-Data dataset (Kaggle, HELP International). This dataset includes socio-economic and health-related features for 167 countries in 2010, with the goal of helping HELP International allocate donations strategically by identifying the nations in greatest need.

Furthermore, we extended our analysis by comparing these results to more recent data (2022) obtained from the World Bank Group's website.

All the analyses are performed with a Dell G5 15 (Windows 10, 64-bit), using a R version 4.3.3.
