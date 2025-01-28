SELECT DISTINCT
    s.seqhash_id,
    s.sample_id,
    s.length,
    s.is_representative,
    a.description,
    a.preferred_name,
    a.cog_category,
    a.evalue,
    a.seed_ortholog,
    a.max_annot_lvl,
    a.eggnog_ogs,
    c.cluster_id,
    c.size AS cluster_size,
    m.organism,
    m.study_title,
    m.bioproject,
    m.biosample
FROM sequences s
LEFT JOIN annotations a ON s.seqhash_id = a.seqhash_id
LEFT JOIN cluster_members cm ON s.seqhash_id = cm.seqhash_id
LEFT JOIN clusters c ON cm.cluster_id = c.cluster_id
LEFT JOIN sra_metadata m ON s.sample_id = m.sample_id
LEFT JOIN go_terms g ON s.seqhash_id = g.seqhash_id
WHERE 
    1=1
    {% if sample_id_condition %}
        AND s.sample_id IN (
            {% for sample in sample_id_condition %}
                '{{ sample }}'
                {% if not loop.last %}, {% endif %}
            {% endfor %}
        )
    {% endif %}
    {% if min_length_condition %} AND {{ min_length_condition }} {% endif %}
    {% if max_length_condition %} AND {{ max_length_condition }} {% endif %}
    {% if description_condition %} AND {{ description_condition }} {% endif %}
    {% if organism_condition %} AND {{ organism_condition }} {% endif %}
    {% if cog_category_condition %}
        AND (
            {% for cog in cog_category_condition %}
                a.cog_category LIKE '%{{ cog }}%'
                {% if not loop.last %} OR {% endif %}
            {% endfor %}
        )
    {% endif %}
    {% if go_term_condition %}
        AND g.go_term IN (
            {% for go in go_term_condition %}
                '{{ go }}'
                {% if not loop.last %}, {% endif %}
            {% endfor %}
        )
    {% endif %}
-- LIMIT {{ limit }}
