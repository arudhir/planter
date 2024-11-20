SELECT 
    m.*,
    COUNT(DISTINCT s.seqhash_id) AS total_sequences,
    COUNT(DISTINCT CASE WHEN a.seqhash_id IS NOT NULL THEN s.seqhash_id END) AS annotated_sequences,
    COUNT(DISTINCT CASE WHEN g.seqhash_id IS NOT NULL THEN s.seqhash_id END) AS sequences_with_go,
    COUNT(DISTINCT CASE WHEN e.seqhash_id IS NOT NULL THEN s.seqhash_id END) AS sequences_with_ec,
    COUNT(DISTINCT CASE WHEN cm.seqhash_id IS NOT NULL THEN s.seqhash_id END) AS sequences_in_clusters
FROM sra_metadata m
JOIN sequences s ON s.sample_id = m.sample_id
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
WHERE m.sample_id = ?
GROUP BY 
    m.sample_id, m.organism, m.study_title, m.study_abstract, 
    m.bioproject, m.biosample, m.library_strategy, m.library_source,
    m.library_selection, m.library_layout, m.instrument, 
    m.run_spots, m.run_bases, m.run_published
