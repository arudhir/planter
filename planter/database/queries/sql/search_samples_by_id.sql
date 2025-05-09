-- Search samples by ID
-- Returns metadata for samples matching the provided ID(s)
-- Example usage: Replace 'SRR12068547', 'SRR12068548' with your own sample IDs

SELECT 
    m.sample_id,
    m.organism,
    m.experiment_id,
    m.run_id,
    m.library_strategy,
    m.library_source,
    m.library_selection,
    m.library_layout,
    m.platform,
    m.instrument_model,
    m.collection_date,
    m.geo_loc_name,
    m.host,
    m.isolation_source,
    COUNT(DISTINCT s.seqhash_id) AS sequence_count,
    COUNT(DISTINCT CASE WHEN s.is_representative THEN s.seqhash_id END) AS representative_count,
    COUNT(DISTINCT CASE WHEN a.seqhash_id IS NOT NULL THEN s.seqhash_id END) AS annotated_count
FROM sra_metadata m
LEFT JOIN sequences s ON m.sample_id = s.sample_id
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
WHERE m.sample_id IN ('SRR12068547', 'SRR12068548')  -- Replace these with your own sample IDs
GROUP BY 
    m.sample_id,
    m.organism,
    m.experiment_id,
    m.run_id,
    m.library_strategy,
    m.library_source,
    m.library_selection,
    m.library_layout,
    m.platform,
    m.instrument_model,
    m.collection_date,
    m.geo_loc_name,
    m.host,
    m.isolation_source
ORDER BY m.sample_id