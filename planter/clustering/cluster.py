#!/usr/bin/env python3
"""
Pure-function clustering module using MMseqs2.

This module provides a clean interface for sequence clustering that:
- Takes sequences in, returns cluster assignments out
- Has no side effects on databases
- Is easily testable and reproducible
"""

import hashlib
import logging
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ClusteringParams:
    """Parameters for MMseqs2 clustering."""

    min_seq_id: float = 0.3  # --min-seq-id
    coverage: float = 0.8  # -c
    coverage_mode: int = 0  # --cov-mode
    cluster_mode: int = 0  # --cluster-mode
    threads: int = 4

    def to_mmseqs_args(self) -> List[str]:
        """Convert to MMseqs2 command line arguments."""
        return [
            "--min-seq-id", str(self.min_seq_id),
            "-c", str(self.coverage),
            "--cov-mode", str(self.coverage_mode),
            "--cluster-mode", str(self.cluster_mode),
            "--threads", str(self.threads),
        ]

    def to_dict(self) -> Dict:
        """Convert to dictionary for storage."""
        return {
            "min_seq_id": self.min_seq_id,
            "coverage": self.coverage,
            "coverage_mode": self.coverage_mode,
            "cluster_mode": self.cluster_mode,
            "threads": self.threads,
        }


@dataclass
class ClusteringResult:
    """Result of a clustering operation."""

    assignments: pd.DataFrame  # columns: seqhash_id, cluster_id, is_representative
    sequence_count: int
    cluster_count: int
    params: ClusteringParams
    success: bool = True
    error_message: Optional[str] = None

    @property
    def representatives(self) -> pd.DataFrame:
        """Get only representative sequences."""
        return self.assignments[self.assignments["is_representative"]]


def _generate_cluster_id(representative_seqhash_id: str) -> str:
    """
    Generate a stable cluster ID from the representative sequence ID.

    Uses a hash to create a shorter, stable identifier that won't change
    unless the representative changes.
    """
    # Use first 12 chars of MD5 hash for reasonable uniqueness
    hash_val = hashlib.md5(representative_seqhash_id.encode()).hexdigest()[:12]
    return f"CLU_{hash_val}"


def _run_mmseqs_command(cmd: List[str], log_file: Optional[Path] = None) -> bool:
    """Run an MMseqs2 command with proper error handling."""
    logger.info(f"Running: {' '.join(cmd)}")

    try:
        if log_file:
            with open(log_file, "a") as lf:
                result = subprocess.run(
                    cmd,
                    stdout=lf,
                    stderr=subprocess.STDOUT,
                    check=True
                )
        else:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        logger.error(f"Return code: {e.returncode}")
        if e.stderr:
            logger.error(f"Stderr: {e.stderr}")
        return False


def run_clustering(
    fasta_path: Path,
    params: Optional[ClusteringParams] = None,
    work_dir: Optional[Path] = None,
    keep_temp: bool = False,
) -> ClusteringResult:
    """
    Run MMseqs2 clustering on a FASTA file.

    This is a pure function: sequences in → cluster assignments out.
    No database interaction, no side effects.

    Args:
        fasta_path: Path to input FASTA file
        params: Clustering parameters (uses defaults if None)
        work_dir: Working directory for temp files (uses system temp if None)
        keep_temp: Whether to keep temporary files after completion

    Returns:
        ClusteringResult with assignments DataFrame
    """
    if params is None:
        params = ClusteringParams()

    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        return ClusteringResult(
            assignments=pd.DataFrame(),
            sequence_count=0,
            cluster_count=0,
            params=params,
            success=False,
            error_message=f"Input FASTA not found: {fasta_path}"
        )

    # Create working directory
    if work_dir:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = None
    else:
        temp_dir = tempfile.TemporaryDirectory()
        work_dir = Path(temp_dir.name)

    try:
        # Define paths
        seq_db = work_dir / "seqDB"
        cluster_db = work_dir / "clusterDB"
        tsv_path = work_dir / "clusters.tsv"
        tmp_dir = work_dir / "tmp"
        log_file = work_dir / "mmseqs.log"

        tmp_dir.mkdir(exist_ok=True)

        # Step 1: Create sequence database
        if not _run_mmseqs_command(
            ["mmseqs", "createdb", str(fasta_path), str(seq_db)],
            log_file
        ):
            return ClusteringResult(
                assignments=pd.DataFrame(),
                sequence_count=0,
                cluster_count=0,
                params=params,
                success=False,
                error_message="Failed to create sequence database"
            )

        # Step 2: Run clustering
        cluster_cmd = [
            "mmseqs", "cluster",
            str(seq_db),
            str(cluster_db),
            str(tmp_dir),
        ] + params.to_mmseqs_args()

        if not _run_mmseqs_command(cluster_cmd, log_file):
            return ClusteringResult(
                assignments=pd.DataFrame(),
                sequence_count=0,
                cluster_count=0,
                params=params,
                success=False,
                error_message="Clustering failed"
            )

        # Step 3: Create TSV output
        if not _run_mmseqs_command(
            ["mmseqs", "createtsv", str(seq_db), str(seq_db), str(cluster_db), str(tsv_path)],
            log_file
        ):
            return ClusteringResult(
                assignments=pd.DataFrame(),
                sequence_count=0,
                cluster_count=0,
                params=params,
                success=False,
                error_message="Failed to create TSV output"
            )

        # Step 4: Parse results into DataFrame
        assignments = _parse_cluster_tsv(tsv_path)

        sequence_count = len(assignments)
        cluster_count = assignments["cluster_id"].nunique()

        logger.info(f"Clustering complete: {sequence_count} sequences in {cluster_count} clusters")

        return ClusteringResult(
            assignments=assignments,
            sequence_count=sequence_count,
            cluster_count=cluster_count,
            params=params,
            success=True
        )

    except Exception as e:
        logger.error(f"Clustering failed with exception: {e}")
        return ClusteringResult(
            assignments=pd.DataFrame(),
            sequence_count=0,
            cluster_count=0,
            params=params,
            success=False,
            error_message=str(e)
        )

    finally:
        if temp_dir and not keep_temp:
            temp_dir.cleanup()


def _parse_cluster_tsv(tsv_path: Path) -> pd.DataFrame:
    """
    Parse MMseqs2 cluster TSV output into a clean DataFrame.

    MMseqs2 TSV format: representative_id<TAB>member_id

    Returns DataFrame with columns:
        - seqhash_id: sequence identifier
        - cluster_id: stable cluster identifier (hash-based)
        - representative_seqhash_id: the representative of this cluster
        - is_representative: boolean flag
    """
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        header=None,
        names=["representative_seqhash_id", "seqhash_id"]
    )

    # Generate stable cluster IDs from representative sequences
    cluster_ids = {
        rep: _generate_cluster_id(rep)
        for rep in df["representative_seqhash_id"].unique()
    }

    df["cluster_id"] = df["representative_seqhash_id"].map(cluster_ids)
    df["is_representative"] = df["seqhash_id"] == df["representative_seqhash_id"]

    return df[["seqhash_id", "cluster_id", "representative_seqhash_id", "is_representative"]]


def export_sequences_to_fasta(
    sequences: pd.DataFrame,
    output_path: Path,
    id_column: str = "seqhash_id",
    seq_column: str = "sequence"
) -> int:
    """
    Export sequences DataFrame to FASTA format.

    Args:
        sequences: DataFrame with sequence data
        output_path: Path to write FASTA file
        id_column: Column name for sequence IDs
        seq_column: Column name for sequences

    Returns:
        Number of sequences written
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with open(output_path, "w") as f:
        for _, row in sequences.iterrows():
            f.write(f">{row[id_column]}\n")
            f.write(f"{row[seq_column]}\n")
            count += 1

    logger.info(f"Exported {count} sequences to {output_path}")
    return count
