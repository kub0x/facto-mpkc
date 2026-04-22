# Code and Experimental Material for Multivariate Factorisation-Based Schemes

This repository contains the scripts, sample outputs, and auxiliary material used in the study of factorisation-based multivariate schemes.

## Repository structure

### `analysed-samples/`
This folder contains the results used in Section 5, **Experimental results**.  
These files are included so that anyone can inspect the generated systems and validate the reported experiments.

### `unpublished/`
This folder contains older and unpublished work dating from 2020, related to the factorisation idea later seen in the X-1 algorithm. These draft contain several erroneous claims and math errors. 

## Main scripts

### `facto-dsa.py`
Generates a system of equations with:

- `2n - 1` equations
- `2n` variables

This script implements the Facto-DSA construction.

### `edf.py`
Generates a system of equations with:

- `2n` equations
- `2n` variables

This script implements the EDF construction.

### `x-1.py`
Generates a system of equations with:

- `2n - 1` equations
- `2n` variables

This script implements the X-1 construction.

### `limit.sh`
Utility script to run commands with resource control, including swap limits.  
It is mainly used to launch heavy MAGMA computations.

### `run.sh`
Runs all MAGMA `.txt` files found in the `systems/` folder in order to compute a Gröbner basis of each system with the `F4` algorithm.

The generated outputs are stored in:

- `systems/results/`

## Generated MAGMA files

Each Python script generates two MAGMA `.txt` files per run:

1. one with the original polynomial system,
2. one with the same system augmented with field equations.

This is done to compare Gröbner basis computations with and without field equations.

## Validation and correctness checks

All scripts include tests for invertibility of the corresponding scheme.  
That is, besides generating the public system, they also verify that the construction can recover or validate the expected internal relation.

## Purpose of the repository

The goal of this repository is to provide:

- the code used to generate the polynomial systems,
- the MAGMA inputs used in Gröbner basis experiments,
- the output samples discussed in the paper,
- and previous unpublished material related to the same research line.

## Notes

The MAGMA runs are intended to study the algebraic behaviour of the generated systems through Gröbner basis computation via `F4`.  
The included sample results allow independent verification of the experimental section.
