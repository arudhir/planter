SELECT 
    s.seqhash_id,
    s.sample_id,
    s.length,
    s.is_representative,
    a.description,
    a.preferred_name,
    a.cog_category,
    STRING_AGG(DISTINCT g.go_term, '; ') AS go_terms,
    STRING_AGG(DISTINCT e.ec_number, '; ') AS ec_numbers
FROM cluster_members cm
JOIN sequences s ON cm.seqhash_id = s.seqhash_id
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
WHERE cm.cluster_id = ?
GROUP BY 
    s.seqhash_id, s.sample_id, s.length, s.is_representative,
    a.description, a.preferred_name, a.cog_category
