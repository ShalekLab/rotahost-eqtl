#!/bin/bash
# Activate conda environment and run command
source /opt/conda/etc/profile.d/conda.sh
conda activate eqtl
exec "$@"
