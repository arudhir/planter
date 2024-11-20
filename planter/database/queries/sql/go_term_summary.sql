SELECT 
    g.go_term,
    COUNT(DISTINCT s.seqhash_id) AS sequence_count
FROM sequences s
JOIN go_terms g ON s.seqhash_id = g.seqhash_id
GROUP BY g.go_term
HAVING COUNT(DISTINCT s.seqhash_id) >= ?
ORDER BY sequence_count DESC
