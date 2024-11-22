SELECT 
    s.seqhash_id,
    s.sample_id,
    m.organism,
    a.description,
    a.preferred_name,
    a.cog_category,
    STRING_AGG(DISTINCT g.go_term, '; ') as go_terms,
    STRING_AGG(DISTINCT e.ec_number, '; ') as ec_numbers,
    k.kegg_ko,
    k.kegg_pathway
FROM sequences s
LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
LEFT JOIN ec_numbers e ON s.seqhash_id = e.seqhash_id
LEFT JOIN kegg_info k ON s.seqhash_id = k.seqhash_id
WHERE s.seqhash_id = ANY(?)
GROUP BY 
    s.seqhash_id, s.sample_id, m.organism,
    a.description, a.preferred_name, a.cog_category,
    k.kegg_ko, k.kegg_pathway
