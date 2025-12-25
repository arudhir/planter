-- Migration: Immutable Clustering Architecture
-- This migration introduces versioned, immutable clustering to replace the mutable approach.
-- Benefits:
--   - No race conditions during updates
--   - Easy rollback to previous clustering
--   - Full audit trail
--   - Clean separation between sequence data and clustering results

-- Track clustering runs with their parameters
CREATE TABLE IF NOT EXISTS clustering_runs (
    run_id INTEGER PRIMARY KEY,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP,
    status VARCHAR DEFAULT 'pending',  -- pending, running, completed, failed
    parameters VARCHAR,  -- JSON string of MMseqs2 parameters
    sequence_count INTEGER,
    cluster_count INTEGER,
    notes VARCHAR
);

-- Clusters for each run (immutable once created)
CREATE TABLE IF NOT EXISTS versioned_clusters (
    run_id INTEGER NOT NULL,
    cluster_id VARCHAR NOT NULL,  -- stable identifier (hash of founding representative)
    representative_seqhash_id VARCHAR NOT NULL,
    size INTEGER NOT NULL,
    PRIMARY KEY (run_id, cluster_id),
    FOREIGN KEY (run_id) REFERENCES clustering_runs(run_id),
    FOREIGN KEY (representative_seqhash_id) REFERENCES sequences(seqhash_id)
);

-- Cluster memberships for each run (immutable once created)
CREATE TABLE IF NOT EXISTS versioned_cluster_members (
    run_id INTEGER NOT NULL,
    seqhash_id VARCHAR NOT NULL,
    cluster_id VARCHAR NOT NULL,
    PRIMARY KEY (run_id, seqhash_id),  -- each sequence in exactly one cluster per run
    FOREIGN KEY (run_id) REFERENCES clustering_runs(run_id),
    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
    FOREIGN KEY (run_id, cluster_id) REFERENCES versioned_clusters(run_id, cluster_id)
);

-- Index for efficient lookups
CREATE INDEX IF NOT EXISTS idx_vcm_cluster ON versioned_cluster_members(run_id, cluster_id);
CREATE INDEX IF NOT EXISTS idx_vc_representative ON versioned_clusters(representative_seqhash_id);
CREATE INDEX IF NOT EXISTS idx_clustering_runs_status ON clustering_runs(status);

-- View for convenient access to current (latest completed) clustering
CREATE OR REPLACE VIEW current_clustering AS
SELECT
    vcm.seqhash_id,
    vcm.cluster_id,
    vc.representative_seqhash_id,
    vc.size as cluster_size,
    vcm.run_id
FROM versioned_cluster_members vcm
JOIN versioned_clusters vc ON vcm.run_id = vc.run_id AND vcm.cluster_id = vc.cluster_id
WHERE vcm.run_id = (
    SELECT MAX(run_id)
    FROM clustering_runs
    WHERE status = 'completed'
);

-- View to get current representatives only
CREATE OR REPLACE VIEW current_representatives AS
SELECT DISTINCT
    vc.representative_seqhash_id as seqhash_id,
    vc.cluster_id,
    vc.size as cluster_size,
    vc.run_id
FROM versioned_clusters vc
WHERE vc.run_id = (
    SELECT MAX(run_id)
    FROM clustering_runs
    WHERE status = 'completed'
);
