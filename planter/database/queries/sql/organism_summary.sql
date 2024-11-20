SELECT 
    m.organism,
    COUNT(DISTINCT m.sample_id) AS sample_count,
    COUNT(DISTINCT s.seqhash_id) AS total_sequences,
    ROUND(AVG(s.length), 2) AS avg_sequence_length,
    COUNT(DISTINCT CASE WHEN a.seqhash_id IS NOT NULL THEN s.seqhash_id END) AS annotated_sequences,
    ROUND(
        COUNT(DISTINCT CASE WHEN a.seqhash_id IS NOT NULL THEN s.seqhash_id END) * 100.0 / 
        COUNT(DISTINCT s.seqhash_id), 2
    ) AS percent_annotated,
    COUNT(DISTINCT CASE WHEN cm.seqhash_id IS NOT NULL THEN s.seqhash_id END) AS sequences_in_clusters,
    STRING_AGG(DISTINCT m.bioproject, '; ') AS bioprojects
FROM sra_metadata m
JOIN sequences s ON s.sample_id = m.sample_id
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
WHERE m.organism IS NOT NULL
GROUP BY m.organism
ORDER BY total_sequences DESC
