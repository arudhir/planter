-- Add gene_protein_map table to link gene and protein sequence hashes
CREATE TABLE IF NOT EXISTS gene_protein_map (
    gene_seqhash_id VARCHAR PRIMARY KEY,
    protein_seqhash_id VARCHAR NOT NULL,
    FOREIGN KEY (protein_seqhash_id) REFERENCES sequences(seqhash_id)
);

-- Add index to improve query performance for protein lookups
CREATE INDEX IF NOT EXISTS idx_gene_protein_protein ON gene_protein_map(protein_seqhash_id);

-- Create a new expression table with the updated schema
-- First, back up the old data
CREATE TABLE expression_backup AS SELECT * FROM expression;

-- Drop the original table
DROP TABLE IF EXISTS expression;

-- Create the new table with gene_seqhash_id
CREATE TABLE expression (
    gene_seqhash_id VARCHAR NOT NULL,
    sample_id VARCHAR NOT NULL,
    tpm DOUBLE NOT NULL,
    num_reads DOUBLE NOT NULL,
    effective_length DOUBLE NOT NULL,
    PRIMARY KEY (gene_seqhash_id, sample_id),
    FOREIGN KEY (gene_seqhash_id) REFERENCES gene_protein_map(gene_seqhash_id),
    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
);

-- Copy data from backup, mapping seqhash_id to gene_seqhash_id
-- (This is only possible if the gene-protein relationships exist)
-- We'll handle the data migration in the loader code