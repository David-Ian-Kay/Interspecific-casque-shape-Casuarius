################################################################################
##############outlier detection Interspecific cassowary#########################
################################################################################
# 2024 David Ian Kay
# Package read-in
library(Momocs)
library(tibble)
library(ggplot2)
library(MASS)
library(vegan)

setwd("FILE-PATH/Interspecific-casque-disparity-Casuarius")

## LATERAL
load("./cass_all_lat_ldk_data.RData")

# Shape quantification from outline data
cal_cass_lat_all_ldk_bk_slid <- calibrate_harmonicpower_efourier(cass_all_lat_ldk_bk_slid, nb.h = 20, plot = T)
cal_cass_lat_all_ldk_bk_slid$gg + theme_minimal() + coord_cartesian(xlim = c(0.5, 17), ylim = c(0, 100)) + ggtitle("Cassowary lateral landmark Bookstein slid Harmonic calibration")
cal_cass_lat_all_ldk_bk_slid_recon <- calibrate_reconstructions_efourier(cass_all_lat_ldk_bk_slid, range = 1:16)
cal_cass_lat_all_ldk_bk_slid_recon
cass_lat_all_ldk_bk_slid_ef <- efourier(cass_all_lat_ldk_bk_slid, nb.h = 16,norm = F)

## PCA ##
cass_lat_all_ldk_bk_slid_PCA <- PCA(cass_lat_all_ldk_bk_slid_ef)
cass_lat_all_ldk_bk_slid_PCA
summary(cass_lat_all_ldk_bk_slid_PCA) # shows to keep 8 PCs to capture >>99% of variation

## outlier Identification of combined interspecies lateral data
cass_lat_outlier<-which_out(cass_lat_all_ldk_bk_slid_PCA$x[,1:8])
cass_lat_outlier
cass_lat_all_ldk_bk_slid_PCA$x[c(8,17,18,22,28,85,93,105,123,127),1:2]
cols<-rep("black",nrow(cass_lat_all_ldk_bk_slid_PCA$x))
outliers2<-c(8,17,18,22,28,85,93,105,123,127)
cols[outliers2]<-"red"
svg('lat_outlier_all.svg')
plot(cass_lat_all_ldk_bk_slid_PCA, col=cols, cex=0.9, zoom=0.9)
dev.off()
svg('lat_all.svg')
plot_PCA(cass_lat_all_ldk_bk_slid_PCA,'species')
dev.off()

# outlier test of data separated by species
# first create speceis-specific datasets from the PCA space
lat_cas_bk_slid_PCA<-filter(cass_lat_all_ldk_bk_slid_PCA,species=='casuarius')
lat_cas_bk_slid_PCA$species
lat_cas_bk_slid_PCA
# outlier identification
lat_cas_bk_slid_PCA_out<-which_out(lat_cass_bk_slid_PCA$x[,1:8],conf=0.001)
lat_cas_bk_slid_PCA_out
# list the species ID and PCs 1 and 2
lat_cas_bk_slid_PCA$x[c(35,39,51,59,71,89,93),1:2]

lat_cas_bk_slid_PCA$x

# repeat for bennetti
lat_ben_bk_slid_PCA<-filter(cass_lat_all_ldk_bk_slid_PCA,species=='bennetti')
lat_ben_bk_slid_PCA$species
lat_ben_bk_slid_PCA
lat_ben_bk_slid_PCA_out<-which_out(lat_ben_bk_slid_PCA$x[,1:8],conf=0.001)
lat_ben_bk_slid_PCA_out
lat_ben_bk_slid_PCA$x[8:9,1:2]

lat_ben_bk_slid_PCA$x

# repeat for unnapendiculatus
lat_una_bk_slid_PCA<-filter(cass_lat_all_ldk_bk_slid_PCA,species=='unappendiculatus')
lat_una_bk_slid_PCA$species
lat_una_bk_slid_PCA
lat_una_bk_slid_PCA_out<-which_out(lat_una_bk_slid_PCA$x[,1:8],conf=0.001)
lat_una_bk_slid_PCA_out

lat_ben_bk_slid_PCA$x

### < ROSTRAL > ###
load("./cass_all_ros_ldk_data.RData")

cal_cass_ros_all_ldk_bk_slid <- calibrate_harmonicpower_efourier(cass_all_ros_ldk_bk_slid, nb.h = 20, plot = T)
cal_cass_ros_all_ldk_bk_slid$gg + theme_minimal() + coord_cartesian(xlim = c(0.5, 15), ylim = c(0, 100)) + ggtitle("Cassowary rostral landmark Bookstein slid Harmonic calibration")
cal_cass_ros_all_ldk_bk_slid_recon <- calibrate_reconstructions_efourier(cass_all_ros_ldk_bk_slid, range = 1:14)
cal_cass_ros_all_ldk_bk_slid_recon
cass_ros_all_ldk_bk_slid_ef <- efourier(cass_all_ros_ldk_bk_slid, nb.h = 14,norm = F)

## PCA ##
cass_ros_all_ldk_bk_slid_PCA <- PCA(cass_ros_all_ldk_bk_slid_ef)
cass_ros_all_ldk_bk_slid_PCA
summary(cass_ros_all_ldk_bk_slid_PCA) # shows to keep 8 PCs to capture >>99% of variation

## outlier Identification of combined interspecies lateral data
cass_ros_outlier<-which_out(cass_ros_all_ldk_bk_slid_PCA$x[,1:8])
cass_ros_outlier
cass_ros_all_ldk_bk_slid_PCA$x[c(9,12,75,82,90,93,116,120,133),1:2]
cols<-rep("black",nrow(cass_ros_all_ldk_bk_slid_PCA$x))
outliers2<-c(9,12,75,82,90,93,116,120,133,175)
cols[outliers2]<-"red"
svg('ros_outlier_all.svg')
plot(cass_ros_all_ldk_bk_slid_PCA, col=cols, cex=0.9, zoom=0.9)
dev.off()
svg('ros_all.svg')
plot_PCA(cass_ros_all_ldk_bk_slid_PCA,'species')
dev.off()


# outlier test of data separated by species
# first create species-specific datasets from the PCA space
ros_cas_bk_slid_PCA<-filter(cass_ros_all_ldk_bk_slid_PCA,species=='casuarius')
ros_cas_bk_slid_PCA$species
ros_cas_bk_slid_PCA
# outlier identification
ros_cas_bk_slid_PCA_out<-which_out(ros_cass_bk_slid_PCA$x[,1:8],conf=0.001)
ros_cas_bk_slid_PCA_out
# list the species ID and PCs 1 and 2
ros_cas_bk_slid_PCA$x[c(41,48,56,59,82),1:2]

# repeat for bennetti
ros_ben_bk_slid_PCA<-filter(cass_ros_all_ldk_bk_slid_PCA,species=='bennetti')
ros_ben_bk_slid_PCA$species
ros_ben_bk_slid_PCA
ros_ben_bk_slid_PCA_out<-which_out(ros_ben_bk_slid_PCA$x[,1:8],conf=0.001)
ros_ben_bk_slid_PCA_out
ros_ben_bk_slid_PCA$x[c(9,12),1:2]

ros_ben_bk_slid_PCA$x

# repeat for unappendiculatus
ros_una_bk_slid_PCA<-filter(cass_ros_all_ldk_bk_slid_PCA,species=='unappendiculatus')
ros_una_bk_slid_PCA$species
ros_una_bk_slid_PCA
ros_una_bk_slid_PCA_out<-which_out(ros_una_bk_slid_PCA$x[,1:8],conf=0.001)
ros_una_bk_slid_PCA_out
ros_una_bk_slid_PCA$x[c(8,12),1:2]
