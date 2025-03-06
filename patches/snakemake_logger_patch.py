#!/usr/bin/env python3
"""
Patch for Snakemake Logger to add missing methods in Snakemake 7.x.
This patch should be applied before running the workflow.
"""
import logging
import sys

def patch_snakemake_logger():
    """Add missing methods to the Snakemake Logger class."""
    try:
        import snakemake.logging
        
        # Only apply if the method is missing
        if not hasattr(snakemake.logging.Logger, 'resources_info'):
            # Add the missing method
            def resources_info(self, msg):
                self.info(msg)
            
            # Attach the method to the Logger class
            snakemake.logging.Logger.resources_info = resources_info
            
            # Add host_info if it's also missing
            if not hasattr(snakemake.logging.Logger, 'host_info'):
                def host_info(self):
                    self.info("Host information: Running on local host")
                
                snakemake.logging.Logger.host_info = host_info
            
            print("Successfully patched Snakemake Logger", file=sys.stderr)
        else:
            print("Snakemake Logger already has resources_info method", file=sys.stderr)
    
    except ImportError:
        print("Could not import snakemake.logging module", file=sys.stderr)
    except Exception as e:
        print(f"Error patching Snakemake Logger: {e}", file=sys.stderr)

if __name__ == "__main__":
    patch_snakemake_logger()