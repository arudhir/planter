SELECT 
    COUNT(DISTINCT s.seqhash_id) as total_sequences,
    COUNT(DISTINCT s.sample_id) as total_samples,
    COUNT(DISTINCT CASE WHEN s.is_representative THEN s.seqhash_id END) as representative_sequences,
    COUNT(DISTINCT CASE WHEN a.seqhash_id IS NOT NULL THEN s.seqhash_id END) as annotated_sequences,
    COUNT(DISTINCT CASE WHEN g.seqhash_id IS NOT NULL THEN s.seqhash_id END) as sequences_with_go,
    COUNT(DISTINCT CASE WHEN e.seqhash_id IS NOT NULL THEN s.seqhash_id END) as sequences_with_ec,
    COUNT(DISTINCT c.cluster_id) as total_clusters,
    ROUND(AVG(s.length), 2) as avg_sequence_length,
    MIN(s.length) as min_sequence_length,
    MAX(s.length) as max_sequence_length
FROM sequences s
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
LEFT JOIN clusters c ON cm.cluster_id = c.cluster_id
