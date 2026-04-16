# RotaHost sc-eQTL Docker image
# Hosts the limix_qtl pipeline for Terra/Cromwell execution
#
# Build: docker build -t ghcr.io/shts2123/rotahost-eqtl:latest .
# Push:  docker push ghcr.io/shts2123/rotahost-eqtl:latest
# Dockstore: https://dockstore.org/organizations/ShalekLab

FROM continuumio/miniconda3:23.10.0-1

LABEL maintainer="ShalekLab <strianas@mit.edu>"
LABEL description="RotaHost sc-eQTL pipeline (limix_qtl + dependencies)"
LABEL version="1.0.0"

# System dependencies (no gsutil needed — Terra PAPI localizes File inputs automatically)
RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python 3.11 via conda (base), then pip install all dependencies
# This avoids conda/pip resolver conflicts by using pip exclusively for packages
RUN conda install -y python=3.11 pip && conda clean -afy

RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    statsmodels \
    h5py \
    tables \
    bgen-reader \
    cbgen \
    glimix-core \
    ndarray-listener \
    pandas-plink \
    limix \
    google-cloud-storage

# Install limix_qtl (apply numpy 2.x compatibility patches)
RUN git clone https://github.com/single-cell-genetics/limix_qtl.git /limix_qtl && \
    sed -i 's/np\.in1d(\(.*\))/np.isin(\1)/g' /limix_qtl/Limix_QTL/qtl_utilities.py && \
    sed -i 's/snpC\.index\[np\.where(snpC==1)\]\.values/snpC[snpC == 1].index.tolist()/g' \
        /limix_qtl/Limix_QTL/qtl_utilities.py

# NOTE: checkpoint resume patches (patches/checkpoint_resume.py) are NOT applied here.
# checkpointFile runtime attribute is not supported on GCPBatch (Terra backend since June 2025).
# Patches are archived in patches/ for future use if Terra adds GCPBatch checkpoint support.

# Add limix_qtl to path
ENV LIMIX_QTL=/limix_qtl/Limix_QTL

# Entry point wrapper
COPY docker_entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker_entrypoint.sh

ENTRYPOINT ["/usr/local/bin/docker_entrypoint.sh"]
