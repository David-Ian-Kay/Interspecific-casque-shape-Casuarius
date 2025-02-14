### << CASSOWARY INTERSPECIFIC ANALYSIS >> ###

### < INITIALIZATION > ###
library(Momocs)
library(tibble)
library(ggplot2)
library(MASS)
library(vegan)
setwd("FILE-PATH/Interspecific-casque-disparity-Casuarius")

### FUNCTION ###
euclidean <- function(a, b) sqrt(sum((a - b)^2))

### < LATERAL DATA > ###
load("./cass_all_lat_ldk_data.RData")
sp.list <- cass_all_lat_ldk_bk_slid$species
ben <- which(sp.list %in% "bennetti")
cas <- which(sp.list %in% "casuarius")
una <- which(sp.list %in% "unappendiculatus") # smallest sample #
coord.data <- cass_all_lat_ldk_bk_slid$coo

cal_cass_lat_all_ldk_bk_slid <- calibrate_harmonicpower_efourier(cass_all_lat_ldk_bk_slid, nb.h = 20, plot = T)
# cal_cass_lat_all_ldk_bk_slid
cal_cass_lat_all_ldk_bk_slid$gg + theme_minimal() + coord_cartesian(xlim = c(0.5, 17), ylim = c(0, 100)) + ggtitle("Cassowary lateral landmark Bookstein slid Harmonic calibration")
cal_cass_lat_all_ldk_bk_slid_recon <- calibrate_reconstructions_efourier(cass_all_lat_ldk_bk_slid, range = 1:16)
cal_cass_lat_all_ldk_bk_slid_recon
cass_lat_all_ldk_bk_slid_ef <- efourier(cass_all_lat_ldk_bk_slid, nb.h = 16,norm = F)

## LDA ##
cass_lat_all_ldk_bk_slid_MDS_MOMOCS <- MDS(cass_lat_all_ldk_bk_slid_ef, k = 39)
str(cass_lat_all_ldk_bk_slid_MDS_MOMOCS)
summary(cass_lat_all_ldk_bk_slid_MDS_MOMOCS)
plot_MDS(cass_lat_all_ldk_bk_slid_MDS_MOMOCS)
cass_lat_all_ldk_bk_slid_PCO_MOMOCS <- capscale(cass_lat_all_ldk_bk_slid_ef ~ 1, distance = "euclidean")


cass_test_ldk_bk_slid_PCO_MOMOCS_lda <- lda(fac.side ~ MDS1 + MDS2 + MDS3 + MDS4 + MDS5 + MDS6 + MDS7 + MDS8, data = cass_test_ldk_bk_slid_PCO_MOMOCS_dat_fac)

## PCA ##
cass_lat_all_ldk_bk_slid_PCA <- PCA(cass_lat_all_ldk_bk_slid_ef)
cass_lat_all_ldk_bk_slid_PCA
summary(cass_lat_all_ldk_bk_slid_PCA) # shows to keep 8 PCs to capture >>99% of variation
pc.shape <- cass_lat_all_ldk_bk_slid_PCA$x

## LDA ##
cass_lat_all_ldk_bk_slid_LDA <- LDA(cass_lat_all_ldk_bk_slid_PCA, "species", retain = 8) # 8 PC axes capure 99% variance
cass_lat_all_ldk_bk_slid_LDA


### PERMUTATION ANALYSIS ###
pc.ben <- pc.shape[ben, 1:8]
mean.ben <- apply(pc.ben, 2, mean)
pc.cas <- pc.shape[cas, 1:8]
mean.cas <- apply(pc.cas, 2, mean)
pc.una <- pc.shape[una, 1:8]
mean.una <- apply(pc.una, 2, mean)

 ## bennetti vs. casuarius ###
a <- pc.cas
b <- pc.una
rep <- 1000
diff <- NULL
p.list <- NULL
n = 1
for(i in 1:rep) {
  subsamp.size <- min(nrow(a), nrow(b))
  sub.a <- a[sample(nrow(a), size=subsamp.size),]
  sub.b <- b[sample(nrow(b), size=subsamp.size),]
  mean.a <- apply(sub.a, 2, mean)
  mean.b <- apply(sub.b, 2, mean)
  mean.diff <- euclidean(mean.a, mean.b)
  #print(mean.diff)
  comb <- rbind(sub.a, sub.b)
  for(j in 1:1000) {
    rand.seq <- sample(nrow(comb))
    rand.mean.1 <- apply(comb[rand.seq[1:subsamp.size],], 2, mean)
    rand.mean.2 <- apply(comb[rand.seq[(subsamp.size+1):(subsamp.size*2)],], 2, mean)
    rand.diff <- euclidean(rand.mean.1, rand.mean.2)
    #print(rand.diff)
    if(rand.diff > mean.diff) {
      n <- n + 1
      }
  }
  p.rep <- n/1000
  p.list[i] <- p.rep
  n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)

# PCA plot without labels #
plot_PCA(cass_lat_all_ldk_bk_slid_PCA, "species")
svg(filename = "./Images/cass_lat_all_ldk_bk_slid_PCA_001.svg")
plot_PCA(cass_lat_all_ldk_bk_slid_PCA)
dev.off()

### < ROSTRAL > ###
load("./cass_all_ros_ldk_data.RData")
sp.list <- cass_all_ros_ldk_bk_slid$species
ben <- which(sp.list %in% "bennetti")
cas <- which(sp.list %in% "casuarius")
una <- which(sp.list %in% "unappendiculatus") # smallest sample #
coord.data <- cass_all_ros_ldk_bk_slid$coo

cal_cass_ros_all_ldk_bk_slid <- calibrate_harmonicpower_efourier(cass_all_ros_ldk_bk_slid, nb.h = 20, plot = T)
# cal_cass_ros_all_ldk_bk_slid
cal_cass_ros_all_ldk_bk_slid$gg + theme_minimal() + coord_cartesian(xlim = c(0.5, 15), ylim = c(0, 100)) + ggtitle("Cassowary roseral landmark Bookstein slid Harmonic calibration")
cal_cass_ros_all_ldk_bk_slid_recon <- calibrate_reconstructions_efourier(cass_all_ros_ldk_bk_slid, range = 1:14)
cal_cass_ros_all_ldk_bk_slid_recon
cass_ros_all_ldk_bk_slid_ef <- efourier(cass_all_ros_ldk_bk_slid, nb.h = 14,norm = F)

## PCA ##
cass_ros_all_ldk_bk_slid_PCA <- PCA(cass_ros_all_ldk_bk_slid_ef)
cass_ros_all_ldk_bk_slid_PCA
summary(cass_ros_all_ldk_bk_slid_PCA) # shows to keep 8 PCs to capture >>99% of variation
pc.shape <- cass_ros_all_ldk_bk_slid_PCA$x

## LDA ##
cass_ros_all_ldk_bk_slid_LDA <- LDA(cass_ros_all_ldk_bk_slid_PCA, "species", retain = 9) # 9 PC axes captures 99% variance
cass_ros_all_ldk_bk_slid_LDA


### PERMUTATION ANALYSIS ###
pc.ben <- pc.shape[ben, 1:9]
mean.ben <- apply(pc.ben, 2, mean)
pc.cas <- pc.shape[cas, 1:9]
mean.cas <- apply(pc.cas, 2, mean)
pc.una <- pc.shape[una, 1:9]
mean.una <- apply(pc.una, 2, mean)

## bennetti vs. casuarius ###
a <- pc.cas
b <- pc.una
rep <- 1000
diff <- NULL
p.list <- NULL
n <- 1
for(i in 1:rep) {
  subsamp.size <- min(nrow(a), nrow(b))
  sub.a <- a[sample(nrow(a), size=subsamp.size),]
  sub.b <- b[sample(nrow(b), size=subsamp.size),]
  mean.a <- apply(sub.a, 2, mean)
  mean.b <- apply(sub.b, 2, mean)
  mean.diff <- euclidean(mean.a, mean.b)
  #print(mean.diff)
  comb <- rbind(sub.a, sub.b)
  for(j in 1:1000) {
    rand.seq <- sample(nrow(comb))
    rand.mean.1 <- apply(comb[rand.seq[1:subsamp.size],], 2, mean)
    rand.mean.2 <- apply(comb[rand.seq[(subsamp.size+1):(subsamp.size*2)],], 2, mean)
    rand.diff <- euclidean(rand.mean.1, rand.mean.2)
    #print(rand.diff)
    if(rand.diff > mean.diff) {
      n <- n + 1
    }
  }
  p.rep <- n/1000
  p.list[i] <- p.rep
  n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)


plot_PCA(cass_ros_all_ldk_bk_slid_PCA, "species")
svg(filename = "./Images/cass_lat_all_ldk_bk_slid_PCA_001.svg")
plot_PCA(cass_ros_all_ldk_bk_slid_PCA)
dev.off()

# PCA plot with point labels #
plot_PCA(cass_ros_all_ldk_bk_slid_PCA, labelpoints = T)
svg(filename = "./Images/cass_lat_all_ldk_bk_slid_PCA_002.svg")
plot_PCA(cass_ros_all_ldk_bk_slid_PCA, labelpoints = T)
dev.off()
