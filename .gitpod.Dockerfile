FROM condaforge/mambaforge

RUN apt-get update && apt-get install -y bash unzip curl