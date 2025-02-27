-- Track database schema versions
CREATE TABLE IF NOT EXISTS schema_version (
    version INTEGER PRIMARY KEY,
    migration_name VARCHAR NOT NULL,
    applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- SRA metadata (this table must be created first as it's referenced by others)
CREATE TABLE IF NOT EXISTS sra_metadata (
    sample_id VARCHAR PRIMARY KEY,
    organism VARCHAR NULL,  -- Remove any unique constraints
    study_title VARCHAR NULL,
    study_abstract VARCHAR NULL,
    bioproject VARCHAR NULL,
    biosample VARCHAR NULL,
    library_strategy VARCHAR NULL,
    library_source VARCHAR NULL,
    library_selection VARCHAR NULL,
    library_layout VARCHAR NULL,
    instrument VARCHAR NULL,
    run_spots VARCHAR NULL,
    run_bases VARCHAR NULL,
    run_published VARCHAR NULL
);

-- Primary sequence storage
CREATE TABLE IF NOT EXISTS sequences (
    seqhash_id VARCHAR PRIMARY KEY,
    sequence VARCHAR NOT NULL,
    sample_id VARCHAR NOT NULL,
    assembly_date TIMESTAMP NOT NULL,
    is_representative BOOLEAN NOT NULL DEFAULT FALSE,
    repseq_id VARCHAR NOT NULL,  -- new column to track the representative sequence
    length INTEGER NOT NULL,
    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
);

-- Sequence annotations from eggNOG
CREATE TABLE IF NOT EXISTS annotations (
    seqhash_id VARCHAR PRIMARY KEY,
    seed_ortholog VARCHAR,
    evalue DOUBLE,
    score DOUBLE,
    eggnog_ogs VARCHAR,
    max_annot_lvl VARCHAR,
    cog_category VARCHAR,
    description VARCHAR,
    preferred_name VARCHAR,
    sample_id VARCHAR NOT NULL,
    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
);

-- GO term associations (many-to-many)
CREATE TABLE IF NOT EXISTS go_terms (
    seqhash_id VARCHAR NOT NULL,
    go_term VARCHAR NOT NULL,
    PRIMARY KEY (seqhash_id, go_term),
    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
);

-- EC number associations (many-to-many)
CREATE TABLE IF NOT EXISTS ec_numbers (
    seqhash_id VARCHAR NOT NULL,
    ec_number VARCHAR NOT NULL,
    PRIMARY KEY (seqhash_id, ec_number),
    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
);

-- KEGG pathway annotations
CREATE TABLE IF NOT EXISTS kegg_info (
    seqhash_id VARCHAR PRIMARY KEY,
    kegg_ko VARCHAR,
    kegg_pathway VARCHAR,
    kegg_module VARCHAR,
    kegg_reaction VARCHAR,
    kegg_rclass VARCHAR,
    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id)
);

-- Sequence clustering information
CREATE TABLE IF NOT EXISTS clusters (
    cluster_id VARCHAR PRIMARY KEY,
    representative_seqhash_id VARCHAR NOT NULL,
    size INTEGER NOT NULL,
    FOREIGN KEY (representative_seqhash_id) REFERENCES sequences(seqhash_id)
);

CREATE TABLE IF NOT EXISTS cluster_members (
    seqhash_id VARCHAR NOT NULL,
    cluster_id VARCHAR NOT NULL,
    PRIMARY KEY (seqhash_id),
    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
    FOREIGN KEY (cluster_id) REFERENCES clusters(cluster_id)
);