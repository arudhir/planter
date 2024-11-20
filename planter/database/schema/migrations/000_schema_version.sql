-- Track database schema versions
CREATE TABLE IF NOT EXISTS schema_version (
    version INTEGER PRIMARY KEY,
    migration_name VARCHAR NOT NULL,
    applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
