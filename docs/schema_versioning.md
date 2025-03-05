# Database Schema Versioning Guide

This document explains the database schema versioning system used in the Planter project, including how to update the schema and ensure compatibility between different versions.

## Schema Version Overview

Planter uses a version-tracked schema system to manage changes to the database structure over time. Each schema version is documented and migrations between versions are automated.

### Current Schema Versions

- **Version 1**: Basic schema with sequences, clusters, and related tables
- **Version 2**: Added `is_representative` column to sequences table

## How Schema Versioning Works

The schema versioning system is implemented in `planter/database/schema/schema_version.py` and provides:

1. **Version Detection**: Automatically detects the schema version of a database
2. **Automated Upgrades**: Upgrades older schemas to newer versions when needed
3. **Forward Compatibility**: Newer code can work with older database schemas
4. **Backward Compatibility**: Older code can (within limits) work with newer schemas

## Updating the Schema

To update the database schema, follow these steps:

### 1. Define the Schema Changes

Add your schema changes to the appropriate migration file in `planter/database/schema/migrations/`:

```sql
-- Example: Adding a new column to an existing table
ALTER TABLE sequences ADD COLUMN new_column VARCHAR;
```

### 2. Update the Schema Version Registry

In `planter/database/schema/schema_version.py`, update the `SCHEMA_VERSIONS` dictionary to include your new version:

```python
SCHEMA_VERSIONS = {
    1: {
        # Version 1 schema definition
    },
    2: {
        # Version 2 schema definition
    },
    3: {  # New version
        "sequences": {
            # Include the updated schema with your new column
            "seqhash_id": "VARCHAR",
            "sequence": "VARCHAR",
            "sample_id": "VARCHAR",
            "assembly_date": "TIMESTAMP",
            "is_representative": "BOOLEAN",
            "repseq_id": "VARCHAR",
            "length": "INTEGER",
            "new_column": "VARCHAR",  # New column added
        },
        # Include other tables as needed
    },
}
```

### 3. Add Upgrade Logic

Add the SQL statements needed to upgrade from the previous version to your new version:

```python
def get_schema_upgrade_sql(from_version: int, to_version: int) -> List[str]:
    # Existing upgrade paths...
    
    # Version 2 to 3: Add new_column to sequences table
    if from_version == 2 and to_version >= 3:
        sql_statements.append("""
            ALTER TABLE sequences 
            ADD COLUMN new_column VARCHAR DEFAULT NULL
        """)
    
    return sql_statements
```

### 4. Test the Schema Update

Create tests to verify that the schema upgrade works correctly. You can model these after the existing tests in `tests/database/test_schema_version.py`.

### 5. Apply the Updates

The schema updates will be automatically applied when:

- Running the workflow with `snakemake`: The `update_database` rule will detect and apply schema updates
- Using database utility functions: Functions like `merge_duckdbs` and `update_duckdb_with_cluster_info` will automatically apply schema updates

## Manual Schema Updates

To manually upgrade a database to a specific schema version:

```python
from planter.database.schema.schema_version import upgrade_schema

# Upgrade to the latest schema version
upgrade_schema("path/to/database.duckdb")

# Upgrade to a specific version
upgrade_schema("path/to/database.duckdb", target_version=3)
```

## Forward/Backward Compatibility Considerations

When making schema changes, consider these best practices:

1. **Add, don't remove or rename**: Add new columns or tables rather than removing or renaming existing ones
2. **Use nullable columns**: Make new columns nullable with defaults when possible
3. **Update utility functions**: Ensure database utility functions handle both old and new schemas
4. **Test with both schemas**: Test code changes with both old and new database schemas

## Schema Compatibility API

The schema compatibility system provides these main functions:

- `get_db_schema_version(db_path)`: Detects the schema version of a database
- `upgrade_schema(db_path, target_version)`: Upgrades a database to a target version
- `ensure_compatibility(db_path, required_version)`: Ensures a database is compatible with code

Example usage:

```python
from planter.database.schema.schema_version import ensure_compatibility

# Check and upgrade if needed
version, was_upgraded = ensure_compatibility(
    "path/to/database.duckdb",
    required_version=2
)

if was_upgraded:
    print(f"Database was upgraded to version {version}")
else:
    print(f"Database is already at version {version}")
```