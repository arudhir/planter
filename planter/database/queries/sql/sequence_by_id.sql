SELECT 
    s.*,
    a.seed_ortholog,
    a.evalue,
    a.score,
    a.eggnog_ogs,
    a.description,
    a.preferred_name,
    a.cog_category,
    STRING_AGG(DISTINCT g.go_term, '; ') AS go_terms,
    STRING_AGG(DISTINCT e.ec_number, '; ') AS ec_numbers,
    k.kegg_ko,
    k.kegg_pathway,
    k.kegg_module,
    c.cluster_id,
    c.size AS cluster_size
FROM sequences s
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
LEFT JOIN kegg_info k ON s.seqhash_id = k.seqhash_id
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
LEFT JOIN clusters c ON cm.cluster_id = c.cluster_id
WHERE s.seqhash_id = ?
GROUP BY 
    s.seqhash_id, s.sequence, s.sample_id, s.assembly_date, 
    s.is_representative, s.length, a.seed_ortholog, a.evalue, 
    a.score, a.eggnog_ogs, a.description, a.preferred_name,
    a.cog_category, k.kegg_ko, k.kegg_pathway, k.kegg_module,
    c.cluster_id, c.size
