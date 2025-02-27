-- Expression distribution breakdown
-- Categorizes sequences by expression level
-- Provides counts for each expression bucket (low, medium, high, very high)

SELECT
    CASE 
        WHEN tpm < 1 THEN 'Low (<1 TPM)'
        WHEN tpm < 10 THEN 'Medium (1-10 TPM)'
        WHEN tpm < 100 THEN 'High (10-100 TPM)'
        ELSE 'Very High (>100 TPM)'
    END as expression_level,
    COUNT(*) as count,
    ROUND(COUNT(*) * 100.0 / SUM(COUNT(*)) OVER(), 2) as percentage
FROM expression
WHERE
    sample_id = COALESCE(?, sample_id)
GROUP BY expression_level
ORDER BY 
    CASE 
        WHEN expression_level = 'Low (<1 TPM)' THEN 1
        WHEN expression_level = 'Medium (1-10 TPM)' THEN 2
        WHEN expression_level = 'High (10-100 TPM)' THEN 3
        ELSE 4
    END