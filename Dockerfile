# RotaHost sc-eQTL Docker image
# Hosts the limix_qtl pipeline for Terra/Cromwell execution
#
# Build: docker build -t ghcr.io/shaleklab/rotahost-eqtl:latest .
# Push:  docker push ghcr.io/shaleklab/rotahost-eqtl:latest
# Dockstore: https://dockstore.org/organizations/ShalekLab

FROM continuumio/miniconda3:23.10.0-1

LABEL maintainer="ShalekLab <strianas@mit.edu>"
LABEL description="RotaHost sc-eQTL pipeline (limix_qtl + dependencies)"
LABEL version="1.0.0"

# System dependencies
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
COPY environment_eqtl.yml /tmp/
RUN conda env create -f /tmp/environment_eqtl.yml && \
    conda clean -afy

# Install limix_qtl (apply numpy 2.x compatibility patches)
RUN git clone https://github.com/single-cell-genetics/limix_qtl.git /limix_qtl && \
    # Fix numpy 2.x compatibility in qtl_utilities.py
    sed -i 's/np\.in1d(\(.*\))/np.isin(\1)/g' /limix_qtl/Limix_QTL/qtl_utilities.py && \
    sed -i 's/snpC\.index\[np\.where(snpC==1)\]\.values/snpC[snpC == 1].index.tolist()/g' \
        /limix_qtl/Limix_QTL/qtl_utilities.py && \
    conda run -n eqtl pip install tables==3.11.1

# Add limix_qtl to path
ENV LIMIX_QTL=/limix_qtl/Limix_QTL
ENV CONDA_ENV=eqtl

# Entry point wrapper
COPY docker_entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker_entrypoint.sh

ENTRYPOINT ["/usr/local/bin/docker_entrypoint.sh"]
