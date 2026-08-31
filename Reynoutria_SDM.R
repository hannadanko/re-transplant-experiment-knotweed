
# Reynoutria japonica - 4-panel SDM figure
#   A: PCA biplot (variable contributions, cos2 colour scale)
#   B: PCA climate space - background vs Europe occurrences vs Japan occurrences
#   C: MaxEnt predicted probability of presence - 6-class map (Europe)
#   D: violin/beeswarm plots of presence probability vs a bioclim variable,
#      with 2.5%/97.5% quantile lines and max/min labels
library(biomod2)
library(readr)
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(lubridate)
library(stringr)
library(rnaturalearth)
library(grid)
library(sf)
library(terra)
library(raster)
library(dismo)
library(ggcorrplot)
library(ENMeval)
library(CoordinateCleaner)
library(megaSDM)
library(viridis)
library(usdm)
library(countrycode)
library(FactoMineR)
library(factoextra)
library(readxl)
library(mgcv)

set.seed(12345)

# FOLDERS -
reynoutria_folder <- "/Users/handa915/Hanna_files/Proj_Reynoutria"
sdm_folder        <- "/Users/handa915/Hanna_files/Proj_Reynoutria"
stand_bioRasters  <- "/Users/handa915/Hanna_files/Proj_Reynoutria"

# 1. GBIF OCCURRENCE CLEANING (Europe) ####

setwd(reynoutria_folder)

world <- ne_countries(scale = "medium", returnclass = "sf")

species_1 <- read.csv("Reynoutria_japonica_occurrences_00_25.csv")

dat <- species_1 %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode,
                individualCount, gbifID, family, taxonRank,
                coordinateUncertaintyInMeters, year, basisOfRecord,
                institutionCode, datasetKey, occurrenceStatus, issue) %>%
  dplyr::filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

dat$countryCode[dat$countryCode %in% c("", "ZZ")] <- NA
dat$countryCode <- countrycode(dat$countryCode, origin = "iso2c", destination = "iso3c")

dat <- data.frame(dat)
flags <- clean_coordinates(
  x = dat, lon = "decimalLongitude", lat = "decimalLatitude",
  countries = "countryCode", species = "species",
  tests = c("capitals", "centroids", "equal", "zeros")
)
dat_cl <- dat[flags$.summary, ]

dat_cl2 <- dat_cl %>%
  filter(coordinateUncertaintyInMeters <= 250 | is.na(coordinateUncertaintyInMeters)) %>%
  filter(occurrenceStatus != "ABSENT") %>%
  filter(year >= 2000, year <= 2016) %>%
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN"))

dat_cl4 <- dat_cl2[!duplicated(dat_cl2[c("decimalLongitude", "decimalLatitude")]), ]

out.round <- cd_round(
  dat_cl4, lon = "decimalLongitude", lat = "decimalLatitude",
  ds = "species", value = "flagged", T1 = 7, graphs = FALSE
)
cleaned_data <- dat_cl4[out.round, ]

occ_unique <- cleaned_data[c("decimalLongitude", "decimalLatitude")]
colnames(occ_unique) <- c("lon", "lat")

# iNaturalist supplement
inat <- read.csv("observations-658464.csv")
inat <- inat %>%
  filter(num_identification_agreements >= 2,
         num_identification_disagreements == 0,
         quality_grade == "research",
         positional_accuracy <= 250, !is.na(positional_accuracy),
         coordinates_obscured == "false") %>%
  mutate(observed_on = as.Date(observed_on)) %>%
  filter(observed_on >= as.Date("2000-01-01"), observed_on <= as.Date("2016-12-31"))
inat <- inat[!duplicated(inat[c("latitude", "longitude")]), ]
inat_short <- inat[c("longitude", "latitude")]
names(inat_short) <- c("lon", "lat")

occ_unique <- rbind(occ_unique, inat_short)
occ_unique <- occ_unique[!duplicated(occ_unique[c("lat", "lon")]), ]

# Europe subset
occ_unique_1 <- occ_unique[
  occ_unique$lon >= -15.00014 & occ_unique$lon <= 69.99986 &
  occ_unique$lat >=  29.99986 & occ_unique$lat <= 83.99986, ]

# 2. ENVIRONMENTAL LAYERS ####
tif_files <- list.files(pattern = "\\.tif$", full.names = TRUE)
studyArea_EURO_2.5 <- terra::rast(tif_files)

# one occurrence per raster cell (Europe)
cells <- terra::cellFromXY(studyArea_EURO_2.5[[1]], as.matrix(occ_unique_1[, c("lon", "lat")]))
occ_final_1 <- occ_unique_1[!duplicated(cells), ]

occ_final_vect_1 <- terra::vect(occ_final_1, geom = c("lon", "lat"), crs = "EPSG:4326")
occs_matrix_1 <- terra::extract(studyArea_EURO_2.5, occ_final_vect_1)[-1]
has_data_1 <- complete.cases(occs_matrix_1)
occ_final_1 <- occ_final_1[has_data_1, ]
occ_final_vect_1 <- terra::vect(occ_final_1, geom = c("lon", "lat"), crs = "EPSG:4326")
occs_matrix_1 <- terra::extract(studyArea_EURO_2.5, occ_final_vect_1)[-1]

# 3. VARIABLE SELECTION (VIF) ####
vif_result_cor_1 <- vifcor(occs_matrix_1, th = 0.7, method = "spearman")
vif_vars_1 <- vif_result_cor_1@results$Variables
studyArea_Euro <- studyArea_EURO_2.5[[vif_vars_1]]

# 4. PANEL A - PCA biplot (variable contributions / cos2) ####
# PCA on the climatic background (M = buffered Europe occurrence area)

buf_eu <- terra::buffer(occ_final_vect_1, width = 300000)
M_eu   <- terra::aggregate(buf_eu)

env_M_eu  <- terra::mask(terra::crop(studyArea_EURO_2.5, M_eu), M_eu)
clim_M_eu <- as.data.frame(terra::values(env_M_eu))
clim_M_eu <- na.omit(clim_M_eu)

pca_eu <- FactoMineR::PCA(clim_M_eu, scale.unit = TRUE, graph = FALSE)

panel_A <- factoextra::fviz_pca_var(
  pca_eu,
  col.var      = "cos2",
  gradient.cols = c("#2DC7C4", "#FDB863", "#D7191C"),
  repel        = TRUE,
  arrowsize    = 0.9
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_equal() +
  ggtitle("Variables - PCA")

panel_A

# 5. PANEL B - PCA climate space: background vs Europe occ vs Japan occ (8x only) ####
# varsenv <- c("bio1", "bio3", "bio7", "bio8", "bio9", "bio14", "bio15", "bio18")
env_M  <- terra::mask(terra::crop(studyArea_EURO_2.5, M_eu), M_eu)
clim_M <- as.data.frame(terra::values(env_M))
clim_M <- na.omit(clim_M)
clim_M <- clim_M[vif_vars_1]

pca_b <- prcomp(clim_M, center = TRUE, scale. = TRUE)

clim_occ_eu <- occs_matrix_1[vif_vars_1]

# --- Japan: keep only octoploid (8x) occurrences ----------------------------#
clim_occ_jap <- read.csv("PloidyLevel_R.japo_Japan.csv")[-1]

scores_bg  <- predict(pca_b, newdata = clim_M)[, 1:2]
scores_occ <- predict(pca_b, newdata = clim_occ_eu)[, 1:2]
scores_jap <- predict(pca_b, newdata = clim_occ_jap)[, 1:2]

colnames(scores_bg)  <- c("PC1", "PC2")
colnames(scores_occ) <- c("PC1", "PC2")
colnames(scores_jap) <- c("PC1", "PC2")

pc1_pct <- round(summary(pca_b)$importance[2, 1] * 100, 1)
pc2_pct <- round(summary(pca_b)$importance[2, 2] * 100, 1)

df_bg  <- data.frame(scores_bg,  group = "background - Europe")
df_occ <- data.frame(scores_occ, group = "occurrence - Europe")
df_jap <- data.frame(scores_jap, group = "occurrence - Japan")
df_pca_b <- bind_rows(df_bg, df_occ, df_jap)
df_pca_b$group <- factor(df_pca_b$group,
                         levels = c("background - Europe",
                                   "occurrence - Europe",
                                   "occurrence - Japan"))

# variable loadings as arrows
load_sc <- as.data.frame(pca_b$rotation[, 1:2] * 5)
load_sc$var <- rownames(load_sc)

panel_B <- ggplot() +
  geom_point(
    data = df_pca_b[df_pca_b$group == "background - Europe", ],
    aes(PC1, PC2), colour = "grey70", size = 0.4, alpha = 0.3
  ) +
  geom_point(
    data = df_pca_b[df_pca_b$group == "occurrence - Europe", ],
    aes(PC1, PC2, colour = group), size = 0.9, shape = 1
  ) +
  geom_point(
    data = df_pca_b[df_pca_b$group == "occurrence - Japan", ],
    aes(PC1, PC2, colour = group), size = 1.2
  ) +
  geom_segment(
    data = load_sc,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.08, "inches")),
    colour = "black", linewidth = 0.6
  ) +
  geom_text(
    data = load_sc,
    aes(x = PC1, y = PC2, label = var),
    size = 3, vjust = -0.4
  ) +
  scale_colour_manual(
    name   = NULL,
    values = c("background - Europe" = "grey70",
              "occurrence - Europe" = "red3",
              "occurrence - Japan"  = "blue3")
  ) +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  ) +
  theme_pubr() +
  theme(legend.position = "top")

panel_B

# 6. PANEL C - MaxEnt probability of presence (map, Europe) ####
# setwd(sdm_folder)
occ_final_forbuff <- occ_final_1 %>% dplyr::rename(x = lon, y = lat)
occ_final_forbuff <- occ_final_forbuff[c("x", "y")]
dir.create("output/buffer", recursive = TRUE, showWarnings = FALSE)
write.csv(occ_final_forbuff, "output/buffer/Reynoutria_japonica.csv", row.names = FALSE)

BackgroundBuffers(
  occlist = "output/buffer/Reynoutria_japonica.csv",
  studyArea_Euro,
  output  = "output/buffer",
  buff_distance = NA,
  ncores  = 7
)

BackgroundBuffers(
  occlist = "output/buffer/Reynoutria_japonica.csv",
  envdata = studyArea_Euro,
  output  = "output/buffer_2",
  buff_distance = 300000,
  ncores  = 1
)

occ_buff_2      <- st_read("output/buffer_2/Reynoutria_japonica.shp")
occ_buff_2_vect <- terra::vect(occ_buff_2)
terra::crs(occ_buff_2_vect) <- terra::crs(studyArea_Euro)
studyArea_Euro_masked <- terra::mask(studyArea_Euro, occ_buff_2_vect)

nbg <- 5000
BackgroundPoints(
  spplist = "Reynoutria_japonica",
  envdata = studyArea_Euro_masked,
  output  = "output/bg",
  nbg     = nbg,
  spatial_weights = 0.75,
  buffers = list.files("output/buffer", pattern = "Reynoutria_japonica.shp", full.names = TRUE),
  method  = "Varela",
  ncores  = 7
)

bg_full <- read.csv("output/bg/Reynoutria_japonica_background.csv")
bg <- bg_full[c("x", "y")]
colnames(bg) <- c("lon", "lat")

# note, that set.seed does not guarantee full reproducibility across the parallelized steps 
# (in BackgroundPoints, ENMevaluate with parallel = TRUE) and has no effect on the external Java process 
# invoked by dismo::maxent(). we do not expect this to affect our overall conclusions, 
# however, as it concerns only minor variation in background point selection and parameter-search order 
# rather than the logic underlying final model selection 

enmeval_results <- ENMevaluate(
  occ_final_1[c("lon", "lat")], studyArea_Euro, bg = bg,
  tune.args = list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), rm = 1:5),
  partitions = "checkerboard",
  algorithm  = "maxnet",
  parallel   = TRUE,
  numCores   = 7,
  partition.settings = list(aggregation.factor = c(10, 10))
)

tmpres <- na.omit(enmeval_results@results) %>%
  dplyr::filter(auc.val.avg == max(auc.val.avg)) %>%
  dplyr::filter(or.10p.avg  == min(or.10p.avg))

fc <- tmpres$fc
rm <- tmpres$rm

# fit final model with dismo::maxent (jar) using the selected fc/rm
studyArea_Euro_r <- raster::stack(studyArea_Euro)
fc_args <- c(
  if (grepl("L", fc)) "linear=true"    else "linear=false",
  if (grepl("Q", fc)) "quadratic=true" else "quadratic=false",
  if (grepl("H", fc)) "hinge=true"     else "hinge=false",
  if (grepl("P", fc)) "product=true"   else "product=false",
  if (grepl("T", fc)) "threshold=true" else "threshold=false",
  "jackknife=true"
)

model <- dismo::maxent(
  x = studyArea_Euro_r,
  p = occ_final_1[c("lon", "lat")],
  a = bg[c("lon", "lat")],
  args = c(paste0("betamultiplier=", rm), fc_args)
)

pred <- predict(model, studyArea_Euro_r)

r_eu  <- terra::rast(pred)
cl_eu <- terra::classify(r_eu, rbind(
  c(0, 0.1, 1), c(0.1, 0.2, 2), c(0.2, 0.3, 3),
  c(0.3, 0.4, 4), c(0.4, 0.5, 5), c(0.5, 1.0, 6)
))

cl_eu_df <- as.data.frame(cl_eu, xy = TRUE, na.rm = TRUE)
colnames(cl_eu_df) <- c("x", "y", "class")
cl_eu_df$class <- factor(
  cl_eu_df$class, levels = 1:6,
  labels = c("0 - 0.1", "0.1 - 0.2", "0.2 - 0.3",
            "0.3 - 0.4", "0.4 - 0.5", "0.5 - 1.0")
)

panel_C <- ggplot(cl_eu_df) +
  geom_tile(aes(x, y, fill = class)) +
  coord_equal(ratio = 1.4, xlim = c(-15, 70), ylim = c(30, 84), expand = FALSE) +
  scale_fill_manual(
    name   = "Predicted probability\nof presence",
  values = c("0 - 0.1"   = "gray90",
             "0.1 - 0.2" = "forestgreen",
             "0.2 - 0.3" = "yellow",
             "0.3 - 0.4" = "orange", 
             "0.4 - 0.5" = "salmon",
             "0.5 - 1.0" = "red")
) +
  labs(x = "Longitude", y = "Latitude") +
  theme_pubr() +
  theme(legend.position = "right", legend.title.align = 0.5)

panel_C

# 7. PANEL D - violins ####

bioclims_explained <- read_xlsx(file.path(reynoutria_folder, "bioclims.xlsx"))

make_violin_panel <- function(var, pred, studyArea_Euro, occs_matrix,
                              bioclims_explained, prob_threshold = 0.5,
                              sample_frac = 0.05, seed = 3215) {

  pred_terra <- terra::rast(pred)        
  names(pred_terra) <- "prob"

  mask_r <- pred_terra >= prob_threshold

  env_plot_masked <- terra::mask(studyArea_Euro[[var]], mask_r, maskvalues = 0)
  env_vals <- terra::values(env_plot_masked)

  pred_plot <- terra::mask(pred_terra, mask_r, maskvalues = 0)
  pred_vals <- terra::values(pred_plot)

  ok <- !is.na(env_vals) & !is.na(pred_vals)
  env_vals  <- env_vals[ok]
  pred_vals <- pred_vals[ok]

  euro_vals <- data.frame(env_vals, pred_vals)

  set.seed(seed)
  n   <- nrow(euro_vals)
  idx <- sample(seq_len(n), size = floor(sample_frac * n))
  euro_vals_sample <- euro_vals[idx, ]

  q_vals   <- quantile(euro_vals$env_vals, c(0.025, 0.975), na.rm = TRUE)
  min_val  <- min(euro_vals$env_vals, na.rm = TRUE)
  max_val  <- max(euro_vals$env_vals, na.rm = TRUE)
  var_label <- bioclims_explained$Chelsa_Name[bioclims_explained$abbr == var]

  ggplot(euro_vals_sample, aes(x = "", y = env_vals)) +
    geom_violin(fill = "#F6D0BA", alpha = 0.25, colour = NA) +
    geom_beeswarm(cex = 0.3, alpha = 0.15, size = 1, colour = "#F6D0BA") +
    geom_boxplot(width = 0.15, notch = TRUE, outlier.shape = NA,
                fill = "white", colour = "#C30D1E", linewidth = 1) +
    geom_hline(yintercept = q_vals, linetype = "dashed",
              linewidth = 0.8, colour = "#C30D1E") +
    annotate("text", x = 1.45, y = max_val,
            label = sprintf("max %.1f", max_val), hjust = 0, size = 3.2,
            colour = "#C30D1E") +
    annotate("text", x = 1.45, y = min_val,
            label = sprintf("min %.1f", min_val), hjust = 0, size = 3.2,
            colour = "#C30D1E") +
    labs(x = NULL, y = var_label) +
    theme_pubr() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

panel_D1 <- make_violin_panel("bio1", pred, studyArea_Euro, occs_matrix_1, bioclims_explained)
panel_D2 <- make_violin_panel("bio7", pred, studyArea_Euro, occs_matrix_1, bioclims_explained)

panel_D1
panel_D2

# 8. FIGURE ####
library(patchwork)

final_fig <- (panel_A | panel_B) /
             (panel_C | (panel_D1 | panel_D2)) +
  plot_annotation(tag_levels = "A")

# final_fig

ggsave(file.path(reynoutria_folder, "Reynoutria_4panel_figure.pdf"),
      final_fig, width = 14, height = 12, device = "pdf", bg = "white")
