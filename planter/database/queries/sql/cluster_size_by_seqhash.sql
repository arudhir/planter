-- Get cluster member lengths by seqhash IDs
-- Returns the sequence length of EVERY UNIQUE MEMBER in the cluster(s) for the specified seqhash ID(s)
-- Useful for verifying that clusters contain full-length sequences
-- and checking for partial or incomplete sequences
-- Example usage: Replace the seqhash IDs with your own

WITH input_clusters AS (
    SELECT DISTINCT cluster_id
    FROM cluster_members
    WHERE seqhash_id IN ('v1_DLS_38e612b34e68d0fc6d2620a80d0b20bbee0e86c1aa057e221db11bb3caafb8cb.p2')  -- Replace with your seqhash IDs
)
SELECT DISTINCT
    cm.cluster_id,
    cm.seqhash_id,
    s.length AS sequence_length,
    s.is_representative,
    s.sample_id,
    sm.organism,
    a.description,
    a.preferred_name
FROM cluster_members cm
JOIN sequences s ON cm.seqhash_id = s.seqhash_id
LEFT JOIN sra_metadata sm ON s.sample_id = sm.sample_id
LEFT JOIN annotations a ON cm.seqhash_id = a.seqhash_id
WHERE cm.cluster_id IN (SELECT cluster_id FROM input_clusters)
ORDER BY cm.cluster_id, s.length DESC, cm.seqhash_id;
