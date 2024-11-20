SELECT 
    COUNT(DISTINCT cluster_id) AS total_clusters,
    AVG(size) AS avg_cluster_size,
    MIN(size) AS min_cluster_size,
    MAX(size) AS max_cluster_size,
    APPROX_QUANTILE(size, 0.5) AS median_cluster_size
FROM clusters
