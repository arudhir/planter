-- Expression data summary statistics
-- Provides overview metrics for expression data (TPM values)
-- Can be filtered by sample_id if provided

SELECT 
    CASE 
        WHEN ? IS NOT NULL THEN ?
        ELSE 'all'
    END as sample_id,
    COUNT(*) as expression_records,
    ROUND(MIN(tpm), 2) as min_tpm,
    ROUND(MAX(tpm), 2) as max_tpm,
    ROUND(AVG(tpm), 2) as mean_tpm,
    ROUND(MEDIAN(tpm), 2) as median_tpm,
    ROUND(STDDEV(tpm), 2) as std_tpm,
    ROUND(AVG(num_reads), 2) as mean_reads,
    ROUND(SUM(num_reads), 0) as total_reads
FROM expression
WHERE
    sample_id = COALESCE(?, sample_id)