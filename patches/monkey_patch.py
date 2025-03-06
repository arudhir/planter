import sys
import os
import importlib.util

# Find the snakemake package
snakemake_spec = importlib.util.find_spec("snakemake")
if snakemake_spec:
    snakemake_path = os.path.dirname(snakemake_spec.origin)
    logging_path = os.path.join(snakemake_path, "logging.py")
    
    # Check if the file exists
    if os.path.exists(logging_path):
        print(f"Found Snakemake logging module at: {logging_path}")
        
        # Read the file content
        with open(logging_path, 'r') as f:
            content = f.read()
        
        # Check if resources_info method is missing
        if "def resources_info" not in content:
            print("Adding missing resources_info method...")
            
            # Add methods at the end of the Logger class
            logger_class_end = content.find("def setup_logfile")
            if logger_class_end > -1:
                # Insert our new methods before setup_logfile
                new_methods = """
    def resources_info(self, msg):
        self.info(msg)
    
    def host_info(self):
        self.info("Host information: Running on local host")
    
"""
                patched_content = content[:logger_class_end] + new_methods + content[logger_class_end:]
                
                # Write the patched file
                with open(logging_path, 'w') as f:
                    f.write(patched_content)
                    
                print("Successfully patched Snakemake logging.py file")
            else:
                print("Could not find the right location to insert the patch")
        else:
            print("resources_info method already exists, no patching needed")
    else:
        print(f"Could not find Snakemake logging module at: {logging_path}")
else:
    print("Could not find Snakemake package")