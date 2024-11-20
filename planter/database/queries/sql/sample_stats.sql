SELECT 
    s.sample_id,
    COUNT(DISTINCT s.seqhash_id) AS total_sequences,
    AVG(s.length) AS avg_length,
    MIN(s.length) AS min_length,
    MAX(s.length) AS max_length,
    COUNT(DISTINCT a.seqhash_id) AS annotated_sequences,
    COUNT(DISTINCT g.seqhash_id) AS sequences_with_go,
    COUNT(DISTINCT e.seqhash_id) AS sequences_with_ec,
    COUNT(DISTINCT CASE WHEN s.is_representative THEN s.seqhash_id END) AS representative_sequences,
    COUNT(DISTINCT cm.cluster_id) AS clusters
FROM sequences s
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
WHERE (? IS NULL OR s.sample_id = ?)
GROUP BY s.sample_id
