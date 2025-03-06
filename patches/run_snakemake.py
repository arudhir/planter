#!/usr/bin/env python3
"""
Wrapper script to run Snakemake with the necessary patches.
"""
import sys
import os
import importlib.util
import subprocess

# Add the parent directory to the path so we can import the patch module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import and apply the patch
from patches.snakemake_logger_patch import patch_snakemake_logger

# Apply the patch
patch_snakemake_logger()

# Now run the original snakemake executable with all arguments
snakemake_args = sys.argv[1:] if len(sys.argv) > 1 else []
cmd = ["snakemake"] + snakemake_args

# Execute Snakemake
try:
    result = subprocess.run(cmd, check=True)
    sys.exit(result.returncode)
except subprocess.CalledProcessError as e:
    sys.exit(e.returncode)
except Exception as e:
    print(f"Error running Snakemake: {e}", file=sys.stderr)
    sys.exit(1)