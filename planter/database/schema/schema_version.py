"""
Schema version management system for database schema compatibility.

This module provides utilities for tracking schema versions and ensuring compatibility
between different versions of the database schema.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import duckdb

logger = logging.getLogger(__name__)

# Define known schema versions and their features
# Format: version_number: {table_name: {column_name: type}}
SCHEMA_VERSIONS = {
    1: {
        "sequences": {
            "seqhash_id": "VARCHAR",
            "sequence": "VARCHAR",
            "sample_id": "VARCHAR",
            "assembly_date": "TIMESTAMP",
            "repseq_id": "VARCHAR",
            "length": "INTEGER",
        },
        "clusters": {
            "cluster_id": "VARCHAR",
            "representative_seqhash_id": "VARCHAR",
            "size": "INTEGER",
        },
        "cluster_members": {
            "seqhash_id": "VARCHAR",
            "cluster_id": "VARCHAR",
        },
    },
    2: {
        "sequences": {
            "seqhash_id": "VARCHAR",
            "sequence": "VARCHAR",
            "sample_id": "VARCHAR",
            "assembly_date": "TIMESTAMP",
            "is_representative": "BOOLEAN",
            "repseq_id": "VARCHAR",
            "length": "INTEGER",
        },
        "clusters": {
            "cluster_id": "VARCHAR",
            "representative_seqhash_id": "VARCHAR",
            "size": "INTEGER",
        },
        "cluster_members": {
            "seqhash_id": "VARCHAR",
            "cluster_id": "VARCHAR",
        },
    },
}


def get_db_schema_version(db_path: Union[str, Path]) -> int:
    """
    Determine the schema version of a database.

    Args:
        db_path: Path to the database file

    Returns:
        int: The detected schema version
    """
    db_path = str(db_path)
    con = duckdb.connect(db_path)

    try:
        # First check if we have a schema_version table
        has_version_table = (
            con.execute(
                """
            SELECT COUNT(*) FROM sqlite_master 
            WHERE type='table' AND name='schema_version'
        """
            ).fetchone()[0]
            > 0
        )

        if has_version_table:
            # Get the latest version
            latest_version = con.execute(
                """
                SELECT MAX(version) FROM schema_version
            """
            ).fetchone()[0]

            if latest_version is not None:
                return latest_version

        # If no version table or no version, check table structure
        tables_columns = {}
        for table_name in ["sequences", "clusters", "cluster_members"]:
            has_table = (
                con.execute(
                    f"""
                SELECT COUNT(*) FROM sqlite_master 
                WHERE type='table' AND name='{table_name}'
            """
                ).fetchone()[0]
                > 0
            )

            if has_table:
                columns = con.execute(
                    f"""
                    PRAGMA table_info({table_name})
                """
                ).fetchall()
                tables_columns[table_name] = {col[1] for col in columns}

        # Check if sequences table has is_representative column (v2+)
        if (
            "sequences" in tables_columns
            and "is_representative" in tables_columns["sequences"]
        ):
            return 2
        else:
            return 1

    finally:
        con.close()

    # Default to version 1 if we can't determine the version
    return 1


def get_schema_upgrade_sql(from_version: int, to_version: int) -> List[str]:
    """
    Get the SQL statements needed to upgrade from one schema version to another.

    Args:
        from_version: The current schema version
        to_version: The target schema version

    Returns:
        List[str]: SQL statements to execute for the upgrade
    """
    if from_version >= to_version:
        return []  # No upgrade needed

    sql_statements = []

    # Version 1 to 2: Add is_representative column
    if from_version == 1 and to_version >= 2:
        sql_statements.append(
            """
            ALTER TABLE sequences 
            ADD COLUMN is_representative BOOLEAN DEFAULT FALSE
        """
        )

    return sql_statements


def upgrade_schema(
    db_path: Union[str, Path], target_version: Optional[int] = None
) -> int:
    """
    Upgrade a database schema to the specified version.

    Args:
        db_path: Path to the database file
        target_version: Target version to upgrade to (default: latest)

    Returns:
        int: The new schema version
    """
    db_path = str(db_path)
    current_version = get_db_schema_version(db_path)

    if target_version is None:
        # Get the latest version from SCHEMA_VERSIONS
        target_version = max(SCHEMA_VERSIONS.keys())

    if current_version >= target_version:
        logger.info(
            f"Database already at schema version {current_version}, no upgrade needed"
        )
        return current_version

    logger.info(
        f"Upgrading database schema from version {current_version} to {target_version}"
    )

    con = duckdb.connect(db_path)
    try:
        # Get SQL statements for upgrade
        sql_statements = get_schema_upgrade_sql(current_version, target_version)

        # Execute each statement
        con.execute("BEGIN TRANSACTION")
        for sql in sql_statements:
            con.execute(sql)

        # Update schema version table
        has_version_table = (
            con.execute(
                """
            SELECT COUNT(*) FROM sqlite_master 
            WHERE type='table' AND name='schema_version'
        """
            ).fetchone()[0]
            > 0
        )

        if not has_version_table:
            con.execute(
                """
                CREATE TABLE schema_version (
                    version INTEGER PRIMARY KEY,
                    migration_name VARCHAR NOT NULL,
                    applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """
            )

        # Insert the new version
        con.execute(
            """
            INSERT INTO schema_version (version, migration_name)
            VALUES (?, ?)
        """,
            [target_version, f"upgrade_to_v{target_version}"],
        )

        con.execute("COMMIT")
        logger.info(
            f"Database schema successfully upgraded to version {target_version}"
        )
        return target_version

    except Exception as e:
        con.execute("ROLLBACK")
        logger.error(f"Error upgrading database schema: {str(e)}")
        raise
    finally:
        con.close()


def ensure_compatibility(
    db_path: Union[str, Path], required_version: int = None
) -> Tuple[int, bool]:
    """
    Ensure a database is compatible with the required schema version.
    If the database has an older schema, it will be upgraded.

    Args:
        db_path: Path to the database file
        required_version: The schema version required by the code

    Returns:
        Tuple[int, bool]: (current_version, was_upgraded)
    """
    if required_version is None:
        required_version = max(SCHEMA_VERSIONS.keys())

    current_version = get_db_schema_version(db_path)

    if current_version < required_version:
        # Upgrade needed
        new_version = upgrade_schema(db_path, required_version)
        return new_version, True

    return current_version, False
