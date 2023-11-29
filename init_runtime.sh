# Check if the data/baer directory exists
if [ ! -d "/data" ]
then
    # Create the data directory
    mkdir -p data

    # Download the data
    curl -L -o data/baer_classifications.zip https://storage.googleapis.com/national_park_service/joshua_tree/BAER_Classifications/Fire_data_bundles_tiAWDyiXROUl0fjFRr8r.zip

    # Unzip the data
    unzip data/baer_classifications.zip -d data
fi

# Init the conda environment
conda init bash

# Need to source bashrc to get the init changes in the current shell
source ~/.bashrc

# Activate our environment
conda activate burn-severity