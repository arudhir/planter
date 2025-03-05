# Planter Documentation

Welcome to the Planter project documentation. This documentation provides detailed information about the Planter database system and workflow.

## Database Documentation

- [Database Update Process](database_update_process.md): Overview of how databases are merged and updated
- [Schema Versioning Guide](schema_versioning.md): Guide to the schema versioning system and update process

## Architecture

The Planter database system consists of several key components:

1. **Core Database**: Built on DuckDB for efficient sequence and annotation storage
2. **Clustering System**: Using MMSeqs2 for sequence clustering and representative selection
3. **Schema Versioning**: Supporting forward/backward compatibility between database versions
4. **Workflow Integration**: Snakemake rules for database updates and processing

## Key Features

- **Schema Version Tracking**: Automatic detection and upgrades of database schemas
- **Forward/Backward Compatibility**: Support for different schema versions
- **Error Handling**: Robust error management during database operations
- **Data Integrity**: Database integrity verification after updates
- **Performance**: Efficient database merging and updating

## Getting Started

To get started with Planter, see the main [README.md](../README.md) file at the root of the repository. For more specialized topics, refer to the documentation pages linked above.

## Developer Documentation

The database system has several key components for developers to understand:

- **Database Builder**: `planter/database/builder.py` - Main database creation and management
- **Database Utils**: `planter/database/utils/duckdb_utils.py` - Utilities for database operations
- **Schema Version**: `planter/database/schema/schema_version.py` - Schema versioning system
- **Workflow Rules**: `planter/workflow/rules/finalize.smk` - Snakemake rules for database updates