-- Top expressed sequences by TPM
-- Returns the most highly expressed sequences in the database
-- Can be filtered by sample if sample_id is provided
-- Limit controls how many sequences to return (default: 10)

SELECT 
    e.gene_seqhash_id,
    e.sample_id,
    s.sample_id as original_sample,
    e.tpm,
    e.num_reads,
    e.effective_length,
    s.length as sequence_length,
    s.is_representative,
    a.description,
    a.preferred_name
FROM expression e
JOIN gene_protein_map gpm ON e.gene_seqhash_id = gpm.gene_seqhash_id
JOIN sequences s ON gpm.protein_seqhash_id = s.seqhash_id
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
WHERE
    e.sample_id = COALESCE(?, e.sample_id)
ORDER BY e.tpm DESC
LIMIT {{ limit }}