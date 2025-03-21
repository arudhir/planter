-- Search sequence by seqhash ID
-- Returns comprehensive information about the specified sequence(s)
-- The {SEQHASH_IDS} placeholder will be replaced with the actual list of IDs

WITH seq_info AS (
    SELECT 
        s.seqhash_id,
        s.sample_id,
        s.is_representative,
        s.length,
        s.sequence,
        a.description,
        a.preferred_name,
        a.cog_category,
        a.eggnog_seed,
        sm.organism
    FROM sequences s
    LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
    LEFT JOIN sra_metadata sm ON s.sample_id = sm.sample_id
    WHERE s.seqhash_id IN ({SEQHASH_IDS})
),
cluster_info AS (
    SELECT 
        cm.seqhash_id,
        cm.cluster_id,
        COUNT(*) OVER (PARTITION BY cm.cluster_id) AS cluster_size,
        ARRAY_AGG(cm2.seqhash_id) OVER (PARTITION BY cm.cluster_id) AS cluster_members
    FROM cluster_members cm
    JOIN cluster_members cm2 ON cm.cluster_id = cm2.cluster_id
    WHERE cm.seqhash_id IN ({SEQHASH_IDS})
),
go_terms_info AS (
    SELECT 
        g.seqhash_id,
        ARRAY_AGG(g.go_term) AS go_terms
    FROM go_terms g
    WHERE g.seqhash_id IN ({SEQHASH_IDS})
    GROUP BY g.seqhash_id
),
ec_numbers_info AS (
    SELECT 
        e.seqhash_id,
        ARRAY_AGG(e.ec_number) AS ec_numbers
    FROM ec_numbers e
    WHERE e.seqhash_id IN ({SEQHASH_IDS})
    GROUP BY e.seqhash_id
),
expression_info AS (
    SELECT 
        gpm.protein_seqhash_id AS seqhash_id,
        AVG(expr.tpm) AS avg_tpm,
        MAX(expr.tpm) AS max_tpm,
        SUM(expr.num_reads) AS total_reads
    FROM gene_protein_map gpm
    JOIN expression expr ON gpm.gene_seqhash_id = expr.gene_seqhash_id
    WHERE gpm.protein_seqhash_id IN ({SEQHASH_IDS})
    GROUP BY gpm.protein_seqhash_id
)

SELECT 
    s.seqhash_id,
    s.sample_id,
    s.organism,
    s.is_representative,
    s.length,
    s.sequence,
    s.description,
    s.preferred_name,
    s.cog_category,
    s.eggnog_seed,
    c.cluster_id,
    c.cluster_size,
    c.cluster_members,
    g.go_terms,
    e.ec_numbers,
    expr.avg_tpm,
    expr.max_tpm,
    expr.total_reads
FROM seq_info s
LEFT JOIN cluster_info c ON s.seqhash_id = c.seqhash_id
LEFT JOIN go_terms_info g ON s.seqhash_id = g.seqhash_id
LEFT JOIN ec_numbers_info e ON s.seqhash_id = e.seqhash_id
LEFT JOIN expression_info expr ON s.seqhash_id = expr.seqhash_id
ORDER BY s.seqhash_id