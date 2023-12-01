echo "Hello from runtime!"

# Check if the data/baer directory exists, from a previous initialization of this container
if [ ! -d "/workspace/burn-severity-mapping-poc/data" ]
then
    # Create the data directory
    mkdir -p /workspace/burn-severity-mapping-poc/data

    # Download the data
    curl -L -o /workspace/burn-severity-mapping-poc/data/baer_classifications.zip https://storage.googleapis.com/national_park_service/joshua_tree/BAER_Classifications/Fire_data_bundles_tiAWDyiXROUl0fjFRr8r.zip

    # Unzip the data
    unzip /workspace/burn-severity-mapping-poc/data/baer_classifications.zip -d data
fi