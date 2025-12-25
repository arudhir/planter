"""
Pipeline Execution Page
"""
import streamlit as st
import subprocess
import json
import time
from datetime import datetime

st.set_page_config(page_title="Run Pipeline", page_icon="⚙️", layout="wide")

st.title("⚙️ Run Planter Pipeline")

st.markdown("""
Run the Planter pipeline to process RNA-seq samples. Enter one or more sample IDs (SRA accessions) to process.
""")

# Sample ID input
st.subheader("Sample IDs")

col1, col2 = st.columns([4, 1])

with col1:
    sample_ids_text = st.text_area(
        "Sample IDs (one per line)",
        height=150,
        placeholder="ERR13770160\nSRR12068547\nSRR12068548\n...",
        help="Enter SRA sample IDs, one per line"
    )

with col2:
    st.write("")  # Spacing
    st.write("")  # Spacing
    if st.button("Load Example IDs", use_container_width=True):
        example_ids = "ERR13770160\nSRR12068547\nSRR12068548\nSRR12068549\nSRR12068550\nSRR12068551"
        st.session_state.sample_ids_text = example_ids
        st.rerun()

# Use session state if exists
if 'sample_ids_text' in st.session_state:
    sample_ids_text = st.session_state.sample_ids_text

# Parse sample IDs
sample_ids = [s.strip() for s in sample_ids_text.split('\n') if s.strip()]

if sample_ids:
    st.info(f"📋 {len(sample_ids)} sample(s) entered")

# Pipeline configuration
st.subheader("Configuration")

col1, col2, col3 = st.columns(3)

with col1:
    cores = st.number_input("CPU Cores", min_value=1, max_value=128, value=16)

with col2:
    outdir = st.text_input("Output Directory", value="outputs")

with col3:
    s3_bucket = st.text_input("S3 Bucket", value="recombia.planter")

# Options
col1, col2 = st.columns(2)

with col1:
    dry_run = st.checkbox("Dry Run", value=False, help="Simulate pipeline execution without running tasks")

with col2:
    show_output = st.checkbox("Show Live Output", value=True, help="Display pipeline output in real-time")

# Run pipeline button
st.divider()

if not sample_ids:
    st.warning("⚠️ Please enter at least one sample ID")
    st.button("Run Pipeline", disabled=True, use_container_width=True)
else:
    if st.button("▶️ Run Pipeline", type="primary", use_container_width=True):
        # Create command
        samples_json = json.dumps(sample_ids)

        cmd = [
            'docker-compose', 'run', '--rm', 'planter',
            'snakemake', '--snakefile', 'planter/workflow/Snakefile',
            '--cores', str(cores),
            '--rerun-incomplete',
            '--config',
            f"outdir={outdir}",
            f"s3_bucket={s3_bucket}",
            f"samples={samples_json}"
        ]

        if dry_run:
            cmd.insert(-4, '--dry-run')  # Insert before --config

        # Store job in session state
        st.session_state.pipeline_job = {
            'command': ' '.join(cmd),
            'start_time': datetime.now(),
            'status': 'running',
            'samples': sample_ids,
            'dry_run': dry_run
        }

        # Display command
        with st.expander("📝 Command Being Executed"):
            st.code(' '.join(cmd), language='bash')

        # Run pipeline
        if show_output:
            st.subheader("Pipeline Output")

            output_container = st.empty()
            status_container = st.empty()

            try:
                # Run process
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    bufsize=1
                )

                # Capture output
                output_lines = []

                for line in iter(process.stdout.readline, ''):
                    if line:
                        output_lines.append(line.rstrip())
                        # Show last 100 lines
                        display_lines = output_lines[-100:]
                        output_container.text_area(
                            "Console Output",
                            value='\n'.join(display_lines),
                            height=400,
                            disabled=True
                        )

                # Wait for completion
                return_code = process.wait()

                # Update job status
                st.session_state.pipeline_job['status'] = 'completed' if return_code == 0 else 'failed'
                st.session_state.pipeline_job['return_code'] = return_code
                st.session_state.pipeline_job['output'] = output_lines
                st.session_state.pipeline_job['end_time'] = datetime.now()

                if return_code == 0:
                    status_container.success("✅ Pipeline completed successfully!")
                else:
                    status_container.error(f"❌ Pipeline failed with return code {return_code}")

            except Exception as e:
                st.error(f"Error running pipeline: {str(e)}")
                st.session_state.pipeline_job['status'] = 'failed'
                st.session_state.pipeline_job['error'] = str(e)

        else:
            # Run in background (no live output)
            with st.spinner("Running pipeline in background..."):
                try:
                    process = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=3600  # 1 hour timeout
                    )

                    st.session_state.pipeline_job['status'] = 'completed' if process.returncode == 0 else 'failed'
                    st.session_state.pipeline_job['return_code'] = process.returncode
                    st.session_state.pipeline_job['output'] = process.stdout.split('\n')
                    st.session_state.pipeline_job['stderr'] = process.stderr
                    st.session_state.pipeline_job['end_time'] = datetime.now()

                    if process.returncode == 0:
                        st.success("✅ Pipeline completed successfully!")
                    else:
                        st.error(f"❌ Pipeline failed with return code {process.returncode}")

                        with st.expander("Error Output"):
                            st.text(process.stderr)

                except subprocess.TimeoutExpired:
                    st.error("❌ Pipeline timed out after 1 hour")
                    st.session_state.pipeline_job['status'] = 'timeout'
                except Exception as e:
                    st.error(f"Error running pipeline: {str(e)}")
                    st.session_state.pipeline_job['status'] = 'failed'
                    st.session_state.pipeline_job['error'] = str(e)

# Show last job results
if 'pipeline_job' in st.session_state:
    st.divider()
    st.subheader("Last Pipeline Job")

    job = st.session_state.pipeline_job

    col1, col2, col3 = st.columns(3)

    with col1:
        status_emoji = {
            'running': '🔄',
            'completed': '✅',
            'failed': '❌',
            'timeout': '⏱️'
        }
        st.metric("Status", f"{status_emoji.get(job['status'], '❓')} {job['status'].title()}")

    with col2:
        st.metric("Samples", len(job['samples']))

    with col3:
        if 'end_time' in job:
            duration = job['end_time'] - job['start_time']
            st.metric("Duration", f"{duration.total_seconds():.1f}s")
        else:
            st.metric("Duration", "Running...")

    # Show output
    if 'output' in job and job['output']:
        with st.expander("📄 Pipeline Output"):
            st.text_area(
                "Output",
                value='\n'.join(job['output']),
                height=400,
                disabled=True
            )

    # Show command
    with st.expander("📝 Command Executed"):
        st.code(job['command'], language='bash')

# Help section
with st.expander("❓ Help"):
    st.markdown("""
    ### Running the Pipeline

    **Sample IDs**
    - Enter SRA accession IDs (e.g., ERR13770160, SRR12068547)
    - One sample per line
    - The pipeline will download, assemble, and annotate each sample

    **Configuration**
    - **CPU Cores**: Number of cores to use for parallel processing
    - **Output Directory**: Where to save pipeline results
    - **S3 Bucket**: S3 bucket for storing/retrieving data
    - **Dry Run**: Test the pipeline without actually running it

    **Notes**
    - Large samples may take hours to process
    - Enable "Show Live Output" to see progress in real-time
    - Pipeline runs inside Docker container
    - Requires docker-compose to be configured
    """)
