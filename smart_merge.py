#!/usr/bin/env python3
"""
Smart merge that handles schema differences
"""
import duckdb
from pathlib import Path

# Read samples
with open('/tmp/samples_with_duckdb_clean.txt') as f:
    samples = [line.strip() for line in f if line.strip()]

print(f"Merging {len(samples)} sample databases...")

# Create new database
output_db = '/mnt/data4/master_simple.duckdb'
if Path(output_db).exists():
    Path(output_db).unlink()

con = duckdb.connect(output_db)

# Create tables based on first database
first_sample = samples[0]
first_db = f'/mnt/data4/recombia.planter/{first_sample}/{first_sample}.duckdb'

print(f"Using {first_db} as template...")
con.execute(f"ATTACH '{first_db}' AS template")

# Get table definitions - copy the full schema
tables = con.execute("SELECT name FROM template.sqlite_master WHERE type='table' ORDER BY name").fetchall()
table_names = [t[0] for t in tables]

print(f"Found tables: {', '.join(table_names)}")

# Copy schema from template (this includes all columns and types)
for table in table_names:
    print(f"Creating table {table}...")
    # Create table by selecting from template with WHERE 1=0 (no data)
    con.execute(f"CREATE TABLE {table} AS SELECT * FROM template.{table} WHERE 1=0")

con.execute("DETACH template")

# Now merge all samples
merged_count = 0
failed_samples = []

for i, sample in enumerate(samples, 1):
    db_path = f'/mnt/data4/recombia.planter/{sample}/{sample}.duckdb'

    if not Path(db_path).exists():
        print(f"Skipping {sample} - database not found")
        failed_samples.append((sample, "not found"))
        continue

    print(f"[{i}/{len(samples)}] Merging {sample}...")

    try:
        alias = f"db{i}"
        con.execute(f"ATTACH '{db_path}' AS {alias}")

        # Copy data from each table with column matching
        for table in table_names:
            try:
                # Get columns from master table
                master_cols = con.execute(f"PRAGMA table_info({table})").fetchall()
                master_col_names = [col[1] for col in master_cols]

                # Get columns from source table
                try:
                    source_cols = con.execute(f"PRAGMA table_info({alias}.{table})").fetchall()
                    source_col_names = [col[1] for col in source_cols]
                except:
                    # Table doesn't exist in source, skip
                    continue

                # Find common columns
                common_cols = [col for col in master_col_names if col in source_col_names]

                if common_cols:
                    col_str = ", ".join(common_cols)
                    con.execute(f"INSERT INTO {table} ({col_str}) SELECT {col_str} FROM {alias}.{table}")
                else:
                    print(f"  Warning: No common columns in {table}")

            except Exception as e:
                print(f"  Warning: Error merging {table}: {e}")

        con.execute(f"DETACH {alias}")
        merged_count += 1

    except Exception as e:
        print(f"  ERROR: Failed to attach {sample}: {e}")
        failed_samples.append((sample, str(e)))
        continue

# Check results
result = con.execute("SELECT COUNT(*) as total, COUNT(DISTINCT sample_id) as samples FROM sequences").fetchone()
print(f"\n✓ Merge complete!")
print(f"  Successfully merged: {merged_count}/{len(samples)} samples")
print(f"  Total sequences: {result[0]:,}")
print(f"  Total samples: {result[1]}")

if failed_samples:
    print(f"\n⚠ Failed samples ({len(failed_samples)}):")
    for sample, reason in failed_samples[:10]:
        print(f"  - {sample}: {reason}")
    if len(failed_samples) > 10:
        print(f"  ... and {len(failed_samples) - 10} more")

con.close()

print(f"\nDatabase created: {output_db}")
print(f"\nNext: Load clusters with:")
print(f"  python scripts/load_clusters.py -d {output_db} -t /mnt/data4/planter_outputs/repseq_v3/output130/newClusterDB.tsv")
