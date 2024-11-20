-- Add indexes for common query patterns

-- Sequence indexes
CREATE INDEX IF NOT EXISTS idx_sequences_sample_id ON sequences(sample_id);
CREATE INDEX IF NOT EXISTS idx_sequences_is_representative ON sequences(is_representative);
CREATE INDEX IF NOT EXISTS idx_sequences_length ON sequences(length);

-- Annotation indexes
CREATE INDEX IF NOT EXISTS idx_annotations_sample_id ON annotations(sample_id);
CREATE INDEX IF NOT EXISTS idx_annotations_description ON annotations(description);
CREATE INDEX IF NOT EXISTS idx_annotations_cog_category ON annotations(cog_category);

-- Term and number indexes for frequent lookups
CREATE INDEX IF NOT EXISTS idx_go_terms_term ON go_terms(go_term);
CREATE INDEX IF NOT EXISTS idx_ec_numbers_number ON ec_numbers(ec_number);

-- KEGG pathway lookup optimization
CREATE INDEX IF NOT EXISTS idx_kegg_ko ON kegg_info(kegg_ko);
CREATE INDEX IF NOT EXISTS idx_kegg_pathway ON kegg_info(kegg_pathway);

-- Cluster search optimization
CREATE INDEX IF NOT EXISTS idx_clusters_size ON clusters(size);

-- SRA metadata search optimization (not unique indexes)
CREATE INDEX IF NOT EXISTS idx_sra_metadata_bioproject ON sra_metadata(bioproject);
