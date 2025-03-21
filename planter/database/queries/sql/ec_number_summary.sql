SELECT 
    e.ec_number,
    COUNT(DISTINCT s.seqhash_id) AS sequence_count
FROM sequences s
JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
GROUP BY e.ec_number
HAVING COUNT(DISTINCT s.seqhash_id) >= 1
ORDER BY sequence_count DESC
