SELECT DISTINCT
    s.seqhash_id,
    s.sample_id,
    s.length,
    s.is_representative,
    a.description,
    a.preferred_name,
    a.cog_category,
    c.cluster_id,
    c.size AS cluster_size
FROM sequences s
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
LEFT JOIN clusters c ON cm.cluster_id = c.cluster_id
WHERE 
    1=1
    {{sample_id_condition}}
    {{min_length_condition}}
    {{max_length_condition}}
    {{annotation_condition}}
    {{description_condition}}
    {{go_term_condition}}
    {{ec_number_condition}}
    {{representative_condition}}
    {{min_cluster_size_condition}}
LIMIT ?
