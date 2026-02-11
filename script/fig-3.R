################################################################################
# Figure 3: DTW Clustering Analysis - Regional Heterogeneity
#
# This script generates Figure 3 for the manuscript, showing:
#   - 3A: Map of province clusters (all 34 provinces, k=12)
#   - 3B: Map of western province clusters with Kalimantan (k=10)
#   - 3C: Map of western province clusters without Kalimantan (k=4)
#   - 3D: Map of province dissimilarity relative to region
#
# Dependencies: Requires data-processing.R to be run first
################################################################################

#### Load dependencies (all packages loaded in data-processing.R)
# source("script/data-processing.R")

################################################################################
# SECTION 1: Prepare Cluster Data for All Three Analyses
################################################################################

# All provinces (k=12) - from cases_admin1_cluster
province_clusters_all <- cases_admin1_cluster %>%

  mutate(
    idadmin1 = as.numeric(idadmin1),
    cluster_all = factor(group, levels = as.character(1:12))
  ) %>%
  dplyr::select(idadmin1, admin1, cluster_all)

# Western with Kalimantan (k=10) - from cases_admin1_cluster_west
province_clusters_west <- cases_admin1_cluster_west %>%
  mutate(
    idadmin1 = as.numeric(idadmin1),
    cluster_west = factor(group, levels = as.character(1:10))
  ) %>%
  dplyr::select(idadmin1, admin1, cluster_west)

# Western without Kalimantan (k=4) - from cases_admin1_cluster_west2
province_clusters_west2 <- cases_admin1_cluster_west2 %>%
  mutate(
    idadmin1 = as.numeric(idadmin1),
    cluster_west2 = factor(group, levels = as.character(1:4))
  ) %>%
  dplyr::select(idadmin1, admin1, cluster_west2)

################################################################################
# SECTION 2: Prepare Spatial Data with Clusters
################################################################################

# Join all cluster assignments to shapefile
admin1_shp_clusters <- admin1_shp %>%
  left_join(province_clusters_all, by = "idadmin1") %>%
  left_join(province_clusters_west, by = c("idadmin1", "admin1")) %>%
  left_join(province_clusters_west2, by = c("idadmin1", "admin1"))

# Prepare dissimilarity data for mapping
province_dissimilarity_map <- province_dissimilarity %>%
  dplyr::select(idadmin1, admin1, region, mean_dist_to_region, z_score, relative_dissimilarity)

admin1_shp_dissimilarity <- admin1_shp %>%
  left_join(province_dissimilarity_map, by = "idadmin1")

# Define western provinces (Sumatra + Java-Bali + Kalimantan)
western_ids_with_kal <- idadmin1_indonesia_unique[c(1:17, 20:24)]
western_ids_no_kal <- idadmin1_indonesia_unique[1:17]

# Create masks for graying out non-western provinces
admin1_shp_clusters <- admin1_shp_clusters %>%
  mutate(
    is_western_with_kal = idadmin1 %in% western_ids_with_kal,
    is_western_no_kal = idadmin1 %in% western_ids_no_kal
  )

################################################################################
# SECTION 3: Figure 3A - All Provinces Cluster Map (k=12)
################################################################################

fig_3A <- ggplot() +
  geom_sf(data = other_shp, fill = "gray90", size = 0.2) +
  geom_sf(
    data = admin1_shp_clusters,
    aes(fill = cluster_all),
    color = "black", size = 0.2
  ) +
  theme_Publication(base_size = 10) +
  labs(tag = "A", x = NULL, y = NULL, fill = "Cluster",
       title = "All provinces (k=12)") +
  geom_shadowtext(
    data = admin1_ir_centroids,
    aes(x = X, y = Y, label = admin1_no),
    size = 1.8,
    color = "white",
    bg.color = "grey10",
    fontface = "bold",
    check_overlap = FALSE
  ) +
  coord_sf(
    xlim = c(95.01098, 141.0194),
    ylim = c(-11.00759, 5.906897)
  ) +
  scale_fill_manual(
    values = PairedColor12Steps,
    na.value = "gray80"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5)
  )

################################################################################
# SECTION 4: Figure 3B - Western Provinces with Kalimantan (k=10)
################################################################################

fig_3B <- ggplot() +
  geom_sf(data = other_shp, fill = "gray90", size = 0.2) +
  # Gray out non-western provinces

  geom_sf(
    data = admin1_shp_clusters %>% filter(!is_western_with_kal),
    fill = "gray80", color = "black", size = 0.2
  ) +
  # Show western provinces with clusters

  geom_sf(
    data = admin1_shp_clusters %>% filter(is_western_with_kal),
    aes(fill = cluster_west),
    color = "black", size = 0.2
  ) +
  theme_Publication(base_size = 10) +
  labs(tag = "B", x = NULL, y = NULL, fill = "Cluster",
       title = "Western (with Kalimantan) (k=10)") +
  geom_shadowtext(
    data = admin1_ir_centroids %>% filter(idadmin1 %in% western_ids_with_kal),
    aes(x = X, y = Y, label = admin1_no),
    size = 1.8,
    color = "white",
    bg.color = "grey10",
    fontface = "bold",
    check_overlap = FALSE
  ) +
  coord_sf(
    xlim = c(95.01098, 118.9873),
    ylim = c(-8.849153, 5.906897)
  ) +
  scale_fill_manual(
    values = PairedColor12Steps[1:10],
    na.value = "gray80"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5)
  )

################################################################################
# SECTION 5: Figure 3C - Western Provinces without Kalimantan (k=4)
################################################################################

# Use a cleaner 4-color palette
cluster_colors_4 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")

fig_3C <- ggplot() +
  geom_sf(data = other_shp, fill = "gray90", size = 0.2) +
  # Gray out non-western provinces
  geom_sf(
    data = admin1_shp_clusters %>% filter(!is_western_no_kal),
    fill = "gray80", color = "black", size = 0.2
  ) +
  # Show western provinces (no Kalimantan) with clusters
  geom_sf(
    data = admin1_shp_clusters %>% filter(is_western_no_kal),
    aes(fill = cluster_west2),
    color = "black", size = 0.2
  ) +
  theme_Publication(base_size = 10) +
  labs(tag = "C", x = NULL, y = NULL, fill = "Cluster",
       title = "Western (without Kalimantan) (k=4)") +
  geom_shadowtext(
    data = admin1_ir_centroids %>% filter(idadmin1 %in% western_ids_no_kal),
    aes(x = X, y = Y, label = admin1_no),
    size = 1.8,
    color = "white",
    bg.color = "grey10",
    fontface = "bold",
    check_overlap = FALSE
  ) +
  coord_sf(
    xlim = c(95.01098, 118.9873),
    ylim = c(-8.849153, 5.906897)
  ) +
  scale_fill_manual(
    values = cluster_colors_4,
    na.value = "gray80"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5)
  )

################################################################################
# SECTION 6: Figure 3D - Province Dissimilarity Map
################################################################################

fig_3D <- ggplot() +
  geom_sf(data = other_shp, fill = "gray90", size = 0.2) +
  geom_sf(
    data = admin1_shp_dissimilarity,
    aes(fill = z_score),
    color = "black", size = 0.2
  ) +
  theme_Publication(base_size = 10) +
  labs(tag = "D", x = NULL, y = NULL,
       fill = "Dissimilarity\n(z-score)",
       title = "Within-region dissimilarity") +
  geom_shadowtext(
    data = admin1_ir_centroids,
    aes(x = X, y = Y, label = admin1_no),
    size = 1.8,
    color = "white",
    bg.color = "grey10",
    fontface = "bold",
    check_overlap = FALSE
  ) +
  coord_sf(
    xlim = c(95.01098, 141.0194),
    ylim = c(-11.00759, 5.906897)
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    na.value = "gray80"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5)
  ) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5))

################################################################################
# SECTION 7: Combine and Save Figure 3
################################################################################

fig_3 <- (fig_3A) / (fig_3B + fig_3C) / (fig_3D) +
  plot_layout(heights = c(1, 1, 1))

ggsave("output/fig_3.jpg", fig_3, height = 28, width = 17, unit = "cm", dpi = 600)
ggsave("output/fig_3.pdf", fig_3, height = 28, width = 17, unit = "cm")

print("Figure 3 saved to output/fig_3.jpg and output/fig_3.pdf")

################################################################################
# SECTION 8: Summary Statistics for Figure Caption
################################################################################

# Print silhouette scores for each clustering
cat("\n=== Clustering Summary ===\n")
cat("All provinces (k=12): Silhouette =", round(res_cvi$Sil[res_cvi$k == "k_12"], 3), "\n")
cat("Western + Kalimantan (k=10): Silhouette =", round(res_cvi_west$Sil[res_cvi_west$k == "k_10"], 3), "\n")
cat("Western only (k=4): Silhouette =", round(res_cvi_west2$Sil[res_cvi_west2$k == "k_4"], 3), "\n")

# Cluster sizes
cat("\nCluster sizes (all provinces):\n")
print(table(province_clusters_all$cluster_all))

cat("\nCluster sizes (western + Kalimantan):\n")
print(table(province_clusters_west$cluster_west))

cat("\nCluster sizes (western only):\n")
print(table(province_clusters_west2$cluster_west2))

# Outlier provinces
cat("\nProvinces with high dissimilarity (z > 1.5):\n")
print(province_dissimilarity %>%
        filter(z_score > 1.65) %>%
        dplyr::select(region, admin1, z_score, relative_dissimilarity) %>%
        arrange(desc(z_score)))
