#!/bin/bash

# Create a new conda environment from the environment.yml file
mamba env create -f environment.yml

# Init conda
conda init bash
