# Thermodynamic Analysis of a Lattice Polymer (HP Model)

This project implements a numerical and statistical-mechanical analysis of the 2D HP lattice polymer model.

The workflow includes:

- Scaling analysis of self-avoiding random walks
- Density-of-states based thermodynamic calculations
- Heat capacity evaluation from energy fluctuations
- Native-state probability analysis
- Monte Carlo simulation data validation

The implementation reproduces key statistical mechanics relations directly from the partition function and compares exact density-of-states results with Monte Carlo sampling.
## Model description
The HP model is a simplified representation of proteins and is used to capture the basics of protein folding.
In this model a polymer is represented as a self-avoiding walk on a 2D lattice, where each monomer is classified as either
hydrophobic (H) or polar (P). The energy of a configuration is determined by the number of H-H contacts, with an energy contribution of -1 for each contact. The partition function is computed by summing over all possible configurations, weighted by their Boltzmann factors.


## Scaling Behavior of Self-Avoiding Walks

## Thermodynamics from Density of States

## Monte Carlo Validation

## Results

## How to Run

