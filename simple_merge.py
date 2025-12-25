#!/usr/bin/env python3
"""
Simple direct merge of DuckDB files
"""
import duckdb
from pathlib import Path

# Read samples
with open('/tmp/samples_with_duckdb.txt') as f:
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

# Get table definitions
tables = con.execute("SELECT name FROM template.sqlite_master WHERE type='table' ORDER BY name").fetchall()
table_names = [t[0] for t in tables]

print(f"Found tables: {', '.join(table_names)}")

# Create tables in master
for table in table_names:
    print(f"Creating table {table}...")
    con.execute(f"CREATE TABLE {table} AS SELECT * FROM template.{table} WHERE 1=0")

con.execute("DETACH template")

# Now merge all samples
for i, sample in enumerate(samples, 1):
    db_path = f'/mnt/data4/recombia.planter/{sample}/{sample}.duckdb'

    if not Path(db_path).exists():
        print(f"Skipping {sample} - database not found")
        continue

    print(f"[{i}/{len(samples)}] Merging {sample}...")

    alias = f"db{i}"
    con.execute(f"ATTACH '{db_path}' AS {alias}")

    # Copy data from each table
    for table in table_names:
        try:
            con.execute(f"INSERT INTO {table} SELECT * FROM {alias}.{table}")
        except Exception as e:
            print(f"  Warning: Error merging {table}: {e}")

    con.execute(f"DETACH {alias}")

# Check results
result = con.execute("SELECT COUNT(*) as total, COUNT(DISTINCT sample_id) as samples FROM sequences").fetchone()
print(f"\n✓ Merge complete!")
print(f"  Total sequences: {result[0]:,}")
print(f"  Total samples: {result[1]}")

con.close()

print(f"\nDatabase created: {output_db}")
print(f"\nNext: Load clusters with:")
print(f"  python scripts/load_clusters.py -d {output_db} -t /mnt/data4/planter_outputs/repseq_v3/output130/newClusterDB.tsv")
