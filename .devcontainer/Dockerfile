FROM condaforge/mambaforge

# Get necessary utils
RUN apt-get update && apt-get install -y bash unzip curl

# Copy repo into 
COPY . .

# Create a new conda environment from the environment.yml file 
RUN mamba env create -f environment.yml

# Install nb_conda_kernels in base env to allow for env discovery in jupyter
RUN mamba install -n base nb_conda_kernels