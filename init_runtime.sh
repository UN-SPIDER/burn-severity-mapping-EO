# Create the data directory if it doesn't exist
mkdir -p workspace/data

# Download the data
curl -L -o workspace/data/baer_classifications.zip https://storage.googleapis.com/national_park_service/joshua_tree/BAER_Classifications/Fire_data_bundles_tiAWDyiXROUl0fjFRr8r.zip

# Unzip the data
unzip workspace/data/baer_classifications.zip -d workspace/data
