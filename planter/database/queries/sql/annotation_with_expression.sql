-- Get annotations with expression levels for proteins
-- This query joins across the gene-protein map to connect annotations with expression data
-- Can be filtered by sample_id (optional)
-- Ordered by TPM (highest to lowest)

SELECT 
    a.seqhash_id as protein_id,
    gpm.gene_seqhash_id as gene_id,
    a.preferred_name,
    a.description,
    e.tpm,
    e.num_reads,
    s.length as protein_length,
    a.sample_id,
    sm.organism
FROM annotations a
JOIN gene_protein_map gpm ON a.seqhash_id = gpm.protein_seqhash_id
JOIN expression e ON gpm.gene_seqhash_id = e.gene_seqhash_id
JOIN sequences s ON a.seqhash_id = s.seqhash_id
JOIN sra_metadata sm ON a.sample_id = sm.sample_id
WHERE a.sample_id = COALESCE(?, a.sample_id)
ORDER BY e.tpm DESC
LIMIT {{ limit }}