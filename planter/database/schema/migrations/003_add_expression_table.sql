-- Add expression table to store RNA quantification data from Salmon
CREATE TABLE IF NOT EXISTS expression (
    seqhash_id VARCHAR NOT NULL,
    sample_id VARCHAR NOT NULL,
    tpm DOUBLE NOT NULL,
    num_reads DOUBLE NOT NULL,
    effective_length DOUBLE NOT NULL,
    PRIMARY KEY (seqhash_id, sample_id),
    FOREIGN KEY (seqhash_id) REFERENCES sequences(seqhash_id),
    FOREIGN KEY (sample_id) REFERENCES sra_metadata(sample_id)
);