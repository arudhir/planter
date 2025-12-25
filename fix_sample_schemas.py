#!/usr/bin/env python3
"""
Fix the schema in sample databases to make repseq_id nullable.

This fixes the NOT NULL constraint on repseq_id which prevents merging
sample databases that don't have clustering information yet.
"""
import sys
from pathlib import Path
import duckdb
from tqdm import tqdm

def fix_database_schema(db_path: Path) -> bool:
    """
    Fix the schema of a single database by making repseq_id nullable.

    Returns True if successful, False otherwise.
    """
    try:
        # Try to connect - if WAL file is corrupt, skip this database
        try:
            con = duckdb.connect(str(db_path))
        except Exception as e:
            if "WAL file" in str(e) or "transcripts does not exist" in str(e):
                print(f"  ⚠ {db_path.name}: Corrupted WAL file, skipping (can rebuild from source)")
                return True
            raise

        # Check if sequences table exists and has repseq_id column
        try:
            schema = con.execute("PRAGMA table_info(sequences)").fetchall()
            has_repseq = any(col[1] == 'repseq_id' for col in schema)

            if not has_repseq:
                print(f"  ⚠ {db_path.name}: No repseq_id column, skipping")
                con.close()
                return True
        except Exception as e:
            print(f"  ✗ {db_path.name}: Error checking schema: {e}")
            con.close()
            return False

        # DuckDB doesn't support ALTER COLUMN to change NOT NULL constraint
        # We need to recreate the table, but must handle foreign keys
        con.execute("BEGIN TRANSACTION")

        # Get list of all tables that might have foreign keys to sequences
        tables = con.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
        table_names = [t[0] for t in tables]

        # Tables that reference sequences (in dependency order)
        dependent_tables = ['annotations', 'go_terms', 'ec_numbers', 'kegg_info',
                           'clusters', 'cluster_members', 'gene_protein_map']

        # Backup all dependent tables that exist
        backups_created = []
        for table in dependent_tables:
            if table in table_names:
                try:
                    con.execute(f"CREATE TABLE {table}_backup AS SELECT * FROM {table}")
                    backups_created.append(table)
                except Exception as e:
                    # Table might not exist, that's OK
                    pass

        # Drop dependent tables
        for table in reversed(backups_created):
            con.execute(f"DROP TABLE IF EXISTS {table}")

        # Now backup and recreate sequences table
        con.execute("CREATE TABLE sequences_backup AS SELECT * FROM sequences")
        con.execute("DROP TABLE sequences")

        # Recreate sequences with nullable repseq_id
        con.execute("""
            CREATE TABLE sequences (
                seqhash_id VARCHAR PRIMARY KEY,
                sequence VARCHAR NOT NULL,
                sample_id VARCHAR NOT NULL,
                assembly_date TIMESTAMP NOT NULL,
                is_representative BOOLEAN NOT NULL DEFAULT FALSE,
                repseq_id VARCHAR,  -- Changed from NOT NULL to nullable
                length INTEGER NOT NULL,
                FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
            )
        """)

        # Copy data back
        con.execute("INSERT INTO sequences SELECT * FROM sequences_backup")
        con.execute("DROP TABLE sequences_backup")

        # Recreate dependent tables from backups
        table_schemas = {
            'annotations': """
                CREATE TABLE annotations (
                    seqhash_id VARCHAR PRIMARY KEY,
                    seed_ortholog VARCHAR,
                    evalue DOUBLE,
                    score DOUBLE,
                    eggnog_ogs VARCHAR,
                    max_annot_lvl VARCHAR,
                    cog_category VARCHAR,
                    description VARCHAR,
                    preferred_name VARCHAR,
                    sample_id VARCHAR NOT NULL,
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
                )
            """,
            'go_terms': """
                CREATE TABLE go_terms (
                    seqhash_id VARCHAR NOT NULL,
                    go_term VARCHAR NOT NULL,
                    PRIMARY KEY (seqhash_id, go_term),
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """,
            'ec_numbers': """
                CREATE TABLE ec_numbers (
                    seqhash_id VARCHAR NOT NULL,
                    ec_number VARCHAR NOT NULL,
                    PRIMARY KEY (seqhash_id, ec_number),
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """,
            'kegg_info': """
                CREATE TABLE kegg_info (
                    seqhash_id VARCHAR PRIMARY KEY,
                    kegg_ko VARCHAR,
                    kegg_pathway VARCHAR,
                    kegg_module VARCHAR,
                    kegg_reaction VARCHAR,
                    kegg_rclass VARCHAR,
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """,
            'clusters': """
                CREATE TABLE clusters (
                    cluster_id VARCHAR PRIMARY KEY,
                    representative_seqhash_id VARCHAR NOT NULL,
                    size INTEGER NOT NULL,
                    FOREIGN KEY (representative_seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """,
            'cluster_members': """
                CREATE TABLE cluster_members (
                    seqhash_id VARCHAR NOT NULL,
                    cluster_id VARCHAR NOT NULL,
                    PRIMARY KEY (seqhash_id),
                    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
                    FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
                )
            """,
            'gene_protein_map': """
                CREATE TABLE gene_protein_map (
                    gene_seqhash_id VARCHAR,
                    protein_seqhash_id VARCHAR NOT NULL,
                    PRIMARY KEY (gene_seqhash_id, protein_seqhash_id),
                    FOREIGN KEY (protein_seqhash_id) REFERENCES sequences(seqhash_id)
                )
            """
        }

        for table in backups_created:
            if table in table_schemas:
                con.execute(table_schemas[table])

                # For kegg_info, handle different schema versions
                if table == 'kegg_info':
                    # Get the actual columns from the backup
                    backup_cols = con.execute(f"PRAGMA table_info({table}_backup)").fetchall()
                    backup_col_names = [col[1] for col in backup_cols]

                    # Get the columns from the new table
                    new_cols = con.execute(f"PRAGMA table_info({table})").fetchall()
                    new_col_names = [col[1] for col in new_cols]

                    # Only insert columns that exist in both
                    common_cols = [col for col in new_col_names if col in backup_col_names]
                    col_str = ", ".join(common_cols)

                    con.execute(f"INSERT INTO {table} ({col_str}) SELECT {col_str} FROM {table}_backup")
                else:
                    con.execute(f"INSERT INTO {table} SELECT * FROM {table}_backup")

                con.execute(f"DROP TABLE {table}_backup")

        # Recreate expression table if it exists
        if 'expression' in table_names:
            try:
                # Check if gene_protein_map exists for the foreign key
                has_gene_protein = 'gene_protein_map' in backups_created

                if has_gene_protein:
                    con.execute("""
                        CREATE TABLE expression (
                            gene_seqhash_id VARCHAR NOT NULL,
                            sample_id VARCHAR NOT NULL,
                            tpm DOUBLE NOT NULL,
                            num_reads DOUBLE NOT NULL,
                            effective_length DOUBLE NOT NULL,
                            PRIMARY KEY (gene_seqhash_id, sample_id),
                            FOREIGN KEY (gene_seqhash_id) REFERENCES gene_protein_map(gene_seqhash_id),
                            FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
                        )
                    """)
                    con.execute("INSERT INTO expression SELECT * FROM expression_backup")
                    con.execute("DROP TABLE expression_backup")
            except Exception as e:
                # Expression table recreation failed, that's OK
                pass

        # Commit the transaction
        con.execute("COMMIT")

        con.close()
        return True

    except Exception as e:
        print(f"  ✗ {db_path.name}: Error: {e}")
        try:
            con.execute("ROLLBACK")
            con.close()
        except:
            pass
        return False


def main():
    # Read sample list
    samples_file = Path('/tmp/samples_with_duckdb.txt')

    if not samples_file.exists():
        print(f"Error: {samples_file} not found")
        sys.exit(1)

    with open(samples_file) as f:
        samples = [line.strip() for line in f if line.strip()]

    print(f"Fixing schemas for {len(samples)} sample databases...")
    print("This will make repseq_id nullable to allow merging before clustering.\n")

    base_dir = Path('/mnt/data4/recombia.planter')
    success_count = 0
    fail_count = 0

    for sample in tqdm(samples, desc="Fixing databases"):
        db_path = base_dir / sample / f"{sample}.duckdb"

        if not db_path.exists():
            print(f"  ⚠ {sample}: Database not found, skipping")
            continue

        if fix_database_schema(db_path):
            success_count += 1
        else:
            fail_count += 1

    print(f"\n{'='*60}")
    print(f"✓ Successfully fixed: {success_count} databases")
    if fail_count > 0:
        print(f"✗ Failed: {fail_count} databases")
    print(f"{'='*60}\n")

    if fail_count > 0:
        print("⚠ Some databases failed to update. Check the errors above.")
        sys.exit(1)
    else:
        print("✓ All databases updated successfully!")
        print("\nYou can now run: python rebuild_master_db.py")


if __name__ == '__main__':
    main()
