#!/bin/bash
# Apply the patch and run Snakemake

# Apply the monkey patch to the Snakemake logger
cd /usr/src/planter
python patches/monkey_patch.py

# Run the command that was passed to this script
exec "$@"