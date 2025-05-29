### << CASSOWARY INTERSPECIFIC ANALYSIS >> ###

### < INITIALIZATION > ###
library(Momocs)
library(tibble)
library(ggplot2)
library(MASS)
library(vegan)
setwd("/CassowaryIntersp/")

### FUNCTION ###
euclidean <- function(a, b) sqrt(sum((a - b)^2))

### < LATERAL DATA > ###
load("./cass_all_lat_ldk_data.RData")
sp.list <- cass_all_lat_ldk_bk_slid$species
ben <- which(sp.list %in% "bennetti")
cas <- which(sp.list %in% "casuarius")
una <- which(sp.list %in% "unappendiculatus") # smallest sample #
coord.data <- cass_all_lat_ldk_bk_slid$coo

cal_cass_lat_all <- calibrate_harmonicpower_efourier(cass_all_lat_ldk_bk_slid, nb.h = 20, plot = T)
# cal_cass_lat_all
cal_cass_lat_all$gg + theme_minimal() + coord_cartesian(xlim = c(0.5, 17), ylim = c(0, 100)) + ggtitle("Cassowary Lateral View Harmonic calibration")
cal_cass_lat_all_recon <- calibrate_reconstructions_efourier(cass_all_lat_ldk_bk_slid, range = 1:16)
cal_cass_lat_all_recon
cass_lat_all_ef <- efourier(cass_all_lat_ldk_bk_slid, nb.h = 16,norm = F)

## Principal Coordinates Analysis (PCA) ##
cass_lat_all_PCA <- PCA(cass_lat_all_ef)
cass_lat_all_PCA
summary(cass_lat_all_PCA) # shows to keep 8 PCs to capture >>99% of variation
pc.shape_lat <- cass_lat_all_PCA$x

## Quadratic Discriminant Analysis (QDA) ##
# create dataframe of the 8 principal components that capture 99% of the shape
cass_lat_all_PCs <- cbind(pc.shape_lat[,1:8],cass_lat_all_PCA$fac)
cass_lat_all_QDA <- qda(species~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8, data=cass_lat_all_PCs,
prior=c(1,1,1)/3, CV=T)
cass_lat_all_QDA
# Build a table of the classification results
CV.fac_all <- cass_lat_all_QDA$class
CV.tab_all <- table(cass_lat_all_PCs[,9],CV.fac_all)
names(dimnames(CV.tab_all)) <- c('actual','classified')
CV.correct_all <- sum(diag(CV.tab_all))/sum(CV.tab_all)
tab_all <- CV.tab_all
  ce_all <- sapply(seq_along(1:nrow(tab_all)),
  function(i) 1-(sum(tab_all[i, -i])/sum(tab_all[i, ])))
  names(ce_all) <- rownames(tab_all)
# QDA overall classification rate
CV.correct_all
# QDA correct classification rate
ce_all
#Classification table
tab_all

### PERMUTATION ANALYSIS ###
pc.ben <- pc.shape_lat[ben, 1:8]
mean.ben <- apply(pc.ben, 2, mean)
pc.cas <- pc.shape_lat[cas, 1:8]
mean.cas <- apply(pc.cas, 2, mean)
pc.una <- pc.shape_lat[una, 1:8]
mean.una <- apply(pc.una, 2, mean)

 ## bennetti vs. casuarius ###
a <- pc.ben
b <- pc.cas
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
  p.rep <- n/(1000+1)
  p.list[i] <- p.rep
  n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)

## bennetti vs. unappendiculatus ###
a <- pc.ben
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
 p.rep <- n/(1000+1)
 p.list[i] <- p.rep
 n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)

## casuarius vs. unappendiculatus ###
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
 p.rep <- n/(1000+1)
 p.list[i] <- p.rep
 n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)

# PCA plot without labels #
plot_PCA(cass_lat_all_PCA, "species")
svg(filename = "./cass_lat_all_PCA_001.svg")
plot_PCA(cass_lat_all_PCA)
dev.off()

### < ROSTRAL > ###
load("./cass_all_ros_ldk_data.RData")
sp.list <- cass_all_ros_ldk_bk_slid$species
ben <- which(sp.list %in% "bennetti")
cas <- which(sp.list %in% "casuarius")
una <- which(sp.list %in% "unappendiculatus") # smallest sample #
coord.data <- cass_all_ros_ldk_bk_slid$coo

cal_cass_ros_all <- calibrate_harmonicpower_efourier(cass_all_ros_ldk_bk_slid, nb.h = 20, plot = T)
# cal_cass_ros_all
cal_cass_ros_all$gg + theme_minimal() + coord_cartesian(xlim = c(0.5, 15), ylim = c(0, 100)) + ggtitle("Cassowary Rostral View Harmonic calibration")
cal_cass_ros_all_recon <- calibrate_reconstructions_efourier(cass_all_ros_ldk_bk_slid, range = 1:14)
cal_cass_ros_all_recon
cass_ros_all_ef <- efourier(cass_all_ros_ldk_bk_slid, nb.h = 14,norm = F)

## PCA ##
cass_ros_all_PCA <- PCA(cass_ros_all_ef)
cass_ros_all_PCA
summary(cass_ros_all_PCA) # shows to keep 8 PCs to capture >>99% of variation
pc.shape_ros <- cass_ros_all_PCA$x

## QDA ##
# create dataframe of the 9 principal components that capture 99% of the shape
cass_ros_all_PCs <- cbind(pc.shape_ros[,1:9],cass_ros_all_PCA$fac)
cass_ros_all_QDA <- qda(species~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9,
data = cass_ros_all_PCs, prior=c(1,1,1)/3, CV=TRUE)
cass_ros_all_QDA

# Build a table of the classification results
CV.fac_all <- cass_ros_all_QDA$class
CV.tab_all <- table(cass_ros_all_PCs[,10],CV.fac_all)
names(dimnames(CV.tab_all)) <- c('actual','classified')
CV.correct_all <- sum(diag(CV.tab_all))/sum(CV.tab_all)
tab_all <- CV.tab_all
  ce_all <- sapply(seq_along(1:nrow(tab_all)),
  function(i) 1-(sum(tab_all[i, -i])/sum(tab_all[i, ])))
  names(ce_all) <- rownames(tab_all)
# QDA overall classification rate
CV.correct_all
# QDA Correct classification rate
ce_all
#Classification table
tab_all

### PERMUTATION ANALYSIS ###
pc.ben <- pc.shape_ros[ben, 1:9]
mean.ben <- apply(pc.ben, 2, mean)
pc.cas <- pc.shape_ros[cas, 1:9]
mean.cas <- apply(pc.cas, 2, mean)
pc.una <- pc.shape_ros[una, 1:9]
mean.una <- apply(pc.una, 2, mean)

## bennetti vs. casuarius ###
a <- pc.ben
b <- pc.cas
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
  p.rep <- n/(1000+1)
  p.list[i] <- p.rep
  n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)

## bennetti vs. unappendiculatus ###
a <- pc.ben
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
  p.rep <- n/(1000+1)
  p.list[i] <- p.rep
  n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)

## casuarius vs. unappendiculatus ###
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
  p.rep <- n/(1000+1)
  p.list[i] <- p.rep
  n <- 1
}
p.list
mean(p.list)
min(p.list)
max(p.list)

plot_PCA(cass_ros_all_PCA, "species")
svg(filename = "./cass_lat_all_PCA_001.svg")
plot_PCA(cass_ros_all_PCA)
dev.off()

# PCA plot with point labels #
plot_PCA(cass_ros_all_PCA, labelpoints = T)
svg(filename = "./cass_lat_all_PCA_002.svg")
plot_PCA(cass_ros_all_PCA, labelpoints = T)
dev.off()
