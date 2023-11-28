#!/bin/bash

# Change to the /workspace directory to persist our BAER data
cd /workspace

# Create a new conda environment from the environment.yml file
mamba env create -f environment.yml

# Create the data directory if it doesn't exist
mkdir -p data

# Download the data
curl -L -o data/baer_classifications.zip https://storage.googleapis.com/national_park_service/joshua_tree/BAER_Classifications/Fire_data_bundles_tiAWDyiXROUl0fjFRr8r.zip

# Unzip the data
unzip data/baer_classifications.zip -d data

# Activate the conda environment
source activate burn-severity