-- Expression profile for a specific sequence
-- Returns expression data for a sequence across all samples
-- Sorted by TPM to show highest expression samples first

SELECT 
    e.sample_id,
    e.tpm,
    e.num_reads,
    e.effective_length,
    sm.organism,
    sm.study_title
FROM expression e
JOIN sra_metadata sm ON e.sample_id = sm.sample_id
WHERE e.seqhash_id = ?
ORDER BY e.tpm DESC