# Database Update Process

This document explains the database update process in the Planter project, particularly focusing on how databases are merged and updated with new clustering information.

## Overview

The database update process consists of several key steps:

1. **Merging databases**: Combining multiple sample databases into a master database
2. **Clustering sequences**: Grouping similar sequences to identify representatives
3. **Updating cluster information**: Applying clustering results to the master database

All of these operations are now schema-aware and support forward/backward compatibility between different schema versions.

## Update Database Workflow

The main workflow for updating the database is implemented in the `update_database` rule in `planter/workflow/rules/finalize.smk`. This rule:

1. Downloads the master database from S3
2. Ensures schema compatibility
3. Extracts representative sequences
4. Runs MMSeqs2 clustering
5. Merges sample databases into the master database
6. Updates cluster information
7. Verifies database integrity
8. Uploads the updated database back to S3

### Schema Compatibility

The update process automatically handles schema compatibility:

- If the master database has an older schema, it will be upgraded
- Sample databases with different schema versions will be correctly merged
- Schema-specific operations (like updating `is_representative` flags) are only performed when appropriate

## Database Merge Process

The database merge process is handled by the `merge_duckdbs` function in `planter/database/utils/duckdb_utils.py`:

```python
merged_db_path = merge_duckdbs(
    duckdb_paths=[sample1.duckdb, sample2.duckdb],
    master_db_path=master.duckdb,
    schema_sql_path=schema.sql,
    upgrade_schema=True,  # Automatically upgrade schema if needed
    target_schema_version=None  # Use the latest available version
)
```

This function:

1. Creates or opens the master database
2. Ensures schema compatibility
3. Merges each source database according to their schema version
4. Adapts column mappings for schema differences
5. Reports on the merge process

## Cluster Update Process

The cluster update process is handled by the `update_duckdb_with_cluster_info` function:

```python
update_duckdb_with_cluster_info(
    db_path=master.duckdb,
    tsv_path=clusters.tsv,
    upgrade_schema=True  # Automatically upgrade schema if needed
)
```

This function:

1. Detects the database schema version
2. Loads clustering information from the TSV file
3. Updates sequence representative pointers
4. Rebuilds cluster information
5. Performs schema-specific operations (like setting `is_representative` flags in v2+ schemas)

## Manual Database Updates

To manually merge and update databases:

```python
from planter.database.utils.duckdb_utils import merge_duckdbs, update_duckdb_with_cluster_info
from pathlib import Path

# Merge databases
master_path = merge_duckdbs(
    duckdb_paths=["sample1.duckdb", "sample2.duckdb"],
    master_db_path="master.duckdb",
    schema_sql_path="schema.sql"
)

# Update clustering
update_duckdb_with_cluster_info(
    db_path=master_path,
    tsv_path="clusters.tsv"
)
```

## Handling Different Schema Versions

The database update process handles different schema versions by:

1. **Dynamic column detection**: Identifying available columns in each database
2. **Common column mapping**: Working with the intersection of available columns
3. **Schema-specific operations**: Running operations only when appropriate for the schema version
4. **Automatic upgrades**: Upgrading schemas when required

## Best Practices for Database Updates

1. **Test schema changes**: Test schema changes thoroughly before deploying
2. **Monitor logs**: Watch the update logs for schema adaptation messages
3. **Verify integrity**: Always verify database integrity after updates
4. **Backup before upgrading**: Create backups before performing schema upgrades
5. **Consider versioning**: Use schema version control for significant changes

## Troubleshooting

If issues occur during database updates, check:

1. **Schema compatibility**: Ensure the schema versions are compatible
2. **Column mappings**: Check that column names match between versions
3. **Clustering format**: Verify MMSeqs2 clustering output format
4. **Transaction errors**: Look for transaction errors in the logs

For more detailed information on schema versioning, see [schema_versioning.md](schema_versioning.md).