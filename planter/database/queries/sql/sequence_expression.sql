-- Expression profile for a specific sequence
-- Returns expression data for a sequence across all samples
-- Sorted by TPM to show highest expression samples first

WITH sequence_params AS (
    SELECT
        CASE
            -- When input is a gene_seqhash_id (i.e., doesn't have a .p suffix)
            WHEN ? NOT LIKE '%.p%' THEN ?
            -- When input is a protein_seqhash_id (i.e., has a .p suffix)
            ELSE (SELECT gene_seqhash_id FROM gene_protein_map WHERE protein_seqhash_id = ?)
        END AS target_seqhash_id
)

SELECT 
    e.sample_id,
    e.tpm,
    e.num_reads,
    e.effective_length,
    sm.organism,
    sm.study_title
FROM expression e
JOIN sra_metadata sm ON e.sample_id = sm.sample_id
WHERE e.gene_seqhash_id = (SELECT target_seqhash_id FROM sequence_params)
ORDER BY e.tpm DESC