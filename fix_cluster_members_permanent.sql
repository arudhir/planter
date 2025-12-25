-- Permanent fix to run after update_clusters
-- Ensures ALL sequences with repseq_id are in cluster_members table

BEGIN TRANSACTION;

-- Create missing cluster entries for any repseq_id that doesn't have a cluster yet
INSERT OR IGNORE INTO clusters (cluster_id, representative_seqhash_id, size)
SELECT DISTINCT
    repseq_id AS cluster_id,
    repseq_id AS representative_seqhash_id,
    0 AS size
FROM sequences
WHERE repseq_id IS NOT NULL
  AND repseq_id NOT IN (SELECT cluster_id FROM clusters);

-- Insert missing cluster_members entries
INSERT OR IGNORE INTO cluster_members (seqhash_id, cluster_id)
SELECT 
    seqhash_id,
    repseq_id AS cluster_id
FROM sequences
WHERE repseq_id IS NOT NULL
  AND seqhash_id NOT IN (SELECT seqhash_id FROM cluster_members);

-- Update all cluster sizes
UPDATE clusters
SET size = (
    SELECT COUNT(*) 
    FROM cluster_members cm 
    WHERE cm.cluster_id = clusters.cluster_id
);

COMMIT;
