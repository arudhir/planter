# Troubleshooting Guide

This document contains solutions for common issues encountered when using the Planter pipeline.

## Table of Contents
- [Snakemake Logger Issues](#snakemake-logger-issues)
- [S3 Storage Plugin Issues](#s3-storage-plugin-issues)

## Snakemake Logger Issues

### Problem: AttributeError for Logger Methods

**Error messages like:**
```
AttributeError: 'Logger' object has no attribute 'run_info'
AttributeError: 'Logger' object has no attribute 'resources_info'
AttributeError: 'Logger' object has no attribute 'get_logfile'
AttributeError: 'Logger' object has no attribute 'logfile_hint'
```

**Cause:**
This issue occurs when using Snakemake 8.12.0 with certain package managers (like uv) that resolve dependencies differently than pip. The logger methods exist in the Snakemake codebase, but aren't accessible at runtime due to how the modules are initialized or imported.

Strangely, these methods are visible when inspecting the Logger class directly:
```python
import snakemake
from snakemake.logging import logger
dir(logger)  # Shows all methods including run_info, etc.
```

But they fail when Snakemake tries to access them during workflow execution.

**Solution:**
The solution is to patch the Snakemake workflow.py file to replace calls to non-standard logger methods with standard ones. We've created a patch script that does this automatically:

1. Create a file named `direct_patch.py`:
```python
#!/usr/bin/env python3

"""
This script directly patches Snakemake's workflow.py file to replace all occurrences
of logger.run_info() with logger.info() and handles other methods similarly.
"""

# Patch workflow.py directly
workflow_path = '/path/to/your/venv/lib/python3.x/site-packages/snakemake/workflow.py'

# Read the workflow.py file
with open(workflow_path, 'r') as f:
    content = f.read()

# Replace all problematic methods
replacements = {
    'logger.run_info(': 'logger.info(',
    'logger.resources_info(': 'logger.info(',
    'logger.logfile_hint()': 'None',
    'logger.get_logfile()': 'None'
}

for old, new in replacements.items():
    content = content.replace(old, new)

# Write the patched content back
with open(workflow_path, 'w') as f:
    f.write(content)

print("Successfully patched Snakemake workflow.py")
```

2. Run the script to patch your Snakemake installation:
```bash
python direct_patch.py
```

3. For Docker-based execution, include this in your Dockerfile:
```dockerfile
# Copy the patch script
COPY direct_patch.py /tmp/direct_patch.py

# Apply the patch
RUN python /tmp/direct_patch.py
```

**Alternative Solution:**
If the above patch doesn't resolve all logger issues, you might need a more comprehensive approach that creates a SafeLogger wrapper to intercept all attribute access. See the project's GitHub repository for the more advanced patch script.

## S3 Storage Plugin Issues

### Problem: Missing S3 Storage Plugin

**Error message:**
```
InvalidPluginException in file .../finalize.smk, line XX:
Error loading Snakemake plugin s3: The package snakemake-storage-plugin-s3 is not installed.
```

**Cause:**
The workflow is configured to use S3 storage (`storage: provider = "s3"`) but the required plugin isn't installed.

**Solution:**
1. Add the S3 storage plugin to your dependencies in `pyproject.toml`:
```toml
dependencies = [
    # Other dependencies...
    "snakemake-storage-plugin-s3",
]
```

2. Install the plugin:
```bash
pip install snakemake-storage-plugin-s3
# or
uv pip install snakemake-storage-plugin-s3
```

3. For Docker-based execution, ensure the plugin is installed in your Dockerfile.