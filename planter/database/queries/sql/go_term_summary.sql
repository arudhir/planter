SELECT 
    g.go_term,
    COUNT(DISTINCT s.seqhash_id) AS sequence_count,
    STRING_AGG(DISTINCT s.sample_id, ', ') AS sample_ids
FROM go_terms g
JOIN sequences s ON g.seqhash_id = s.seqhash_id
WHERE (? IS NULL OR s.sample_id = ?)
GROUP BY g.go_term
HAVING COUNT(DISTINCT s.seqhash_id) >= 1
ORDER BY sequence_count DESC
