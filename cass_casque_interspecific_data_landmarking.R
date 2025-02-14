################################################################################
# Data landmarking for Interspecific analyses
# 2022 David Ian Kay
# Written in R 3.6.3
################################################################################
# PACKAGES
# Load Momocs, a morphometrics package. I am using v1.3.0
library(Momocs)
# Load tibble, a package enabling tibble data structures instead of dataframes for momocs to use v3.0.3
library(tibble)
# Command to check which versions of packages are loaded and attached in the R session
sessionInfo()
# Command to ensure that no extranneous objects exist
rm(list=ls())

getwd()
setwd("FILE-PATH/Interspecific-casque-disparity-Casuarius")

# create list of filenames of all lateral casque silhouettes
lf_lat <- list.files("./ALL_SP_LAT", full.names=FALSE)
lf_lat
# write to a .csv to be used in analyses
write.csv(lf_lat,"cass_all_sp_lat_fac.csv")


# create list of filenames of all rostral casque silhouettes
lf_ros <- list.files("./ALL_SP_ROS", full.names=FALSE)
lf_ros
# write to a csv to be used in analyses
write.csv(lf_ros,"cass_all_sp_ros_fac.csv")



# LATERAL
# Data read-in and outline creation
# Read in the factors list to categorize specimens, initialize as a tibble
cass_all_info_lat <- as_tibble(read.csv("./cass_all_sp_lat_fac.csv", header=T, row.names=1))
# Check the factors are correct
cass_all_info_lat
# Create an object containing the filenames to reference them for import
lf_lat <- list.files("./ALL_SP_LAT", full.names=TRUE)
# Double-check the filenames in the object
lf_lat
# Import the binary mask jpegs using the file list
cass_all_lat <- import_jpg(lf_lat)
# Convert them to outline, simultaneously adding factors to the outlines
cass_all_lat_out <- Out(cass_all_lat, fac=cass_all_info_lat)
# Interpolate to the average number of points stated in the object information
cass_all_lat_out
cass_all_lat_out_int <- coo_interpolate(cass_all_lat_out,n=5243)
# Double check the point numbers are correct
cass_all_lat_out_int
pile(cass_all_lat_out)
################### Landmark Assignment and Outline Processing###################
## There does not seem to be a way to go back through the specimens or edit points.
## Landmarks must be assigned in the same order for each specimen (I did anterior first, posterior second)
cass_all_lat_ldk <- def_ldk(cass_all_lat_out_int, 2)
# Inspect to make sure all specimens have 2 landmarks assigned
cass_all_lat_ldk
cass_all_lat_ldk$ldk
# Plot the outlines in their unaligned state, landmarks are red
stack(cass_all_lat_ldk)

# There are only two landmarks on the outlines at the ventral ends of the casques, so a Procrustes alignment is not an option. Bookstein coordinate alignment is though, as it works on only two landmarks. These will be used to establish a new baseline calculated from the first and last landmarks (or just first and second if there are only two), which is exactly what is needed.
cass_all_lat_ldk_bk <- coo_bookstein(cass_all_lat_ldk)
# Plot the aligned outlines to make sure the outlines are aligned along the same baseline
stack(cass_all_lat_ldk_bk) #(they are)
pile(cass_all_lat_ldk_bk)
# The start points for the outlines are inconsistent. Here I create a single origin for each outline
# Write the coordinates out to a csv to create a vector for the coo_slide function, identifying which point is closest to the origin (0,0) of the aligned outlines using the csv in Microsoft Excel
write.csv(cass_all_lat_ldk_bk$coo,file='cass_all_lat_ldk_bk_coordinates.csv')
# Write a vector for sliding the coordinate start point to the same x-y coordinate for all outlines
# Create the list of points to shift all outline starting points to the same location

# read in lateral coordinate data
data_lat <- read.csv("cass_all_lat_ldk_bk_coordinates.csv", header=T, row.names=1)
# check the data
head(data_lat)
# create blank vector for slide vector construction
slide_vector_LAT <- c()
# for loop to generate the individual dataframes
for(n in seq(1, 305, 2)){
  data_xy <- data_lat[ ,c(n:(n+1))]
  # data_xy
  # Subset the data to only include the x-y coordinates that have a y-value of less than 0.01, to avoid designating a point on the casque as the origin
  data_xy_sub <- subset(data_xy, data_xy[ ,2]<0.01)
  # data_xy_sub

  # Find the minimum x-value of the subsetted data
  min_value <- min(abs(data_xy_sub[ ,1]))
  # min_value

  which(data_xy_sub==min_value, arr.ind=T)
  name_info <- which(data_xy_sub==min_value, arr.ind=T)
  # name_info
  name_only <- row.names(name_info)
  # name_only

  # if statement to deal with negative x values that were turned positive via the absolute value function, since they produce a matrix with a length of zero.
  if(length(name_info)==0){
    min_value <- min_value*-1
    name_info <- which(data_xy_sub==min_value, arr.ind=T)
  }
  name_only <- row.names(name_info)
  name_only


  slide_vector_LAT <- c(slide_vector_LAT, name_only)
}

slide_vector_LAT
length(slide_vector_LAT)
print(slide_vector_LAT)
# coerce vector to numeric type instead of string
slide_vector_LAT <- as.numeric(slide_vector_LAT)
# write to csv in to preserve for spot-checking with exported coordinate data csv
write.csv(slide_vector_LAT,"slide_vector_LAT.csv")


# Slide the start points of all outlines to ~ (0,0) using coo_slide
cass_all_lat_ldk_bk_slid <- coo_slide(cass_all_lat_ldk_bk,id=slide_vector_LAT)
# Verify the start points are relatively the same for all specimens
stack(cass_all_lat_ldk_bk_slid)
pile(cass_all_lat_ldk_bk_slid)
# Save the data object for redundancy
save(cass_all_lat_ldk_bk_slid,file='cass_all_lat_ldk_data.RData')


# ROSTRAL
###################### Data read-in and outline creation#########################
# Read in the factors list to categorize specimens, initialize as a tibble
cass_all_info_ros <- as_tibble(read.csv("./cass_all_sp_ros_fac.csv", header=T, row.names=1))
# Check the factors are correct
cass_all_info_ros
# Create an object containing the filenames to reference them for import
lf_ros <- list.files("./ALL_SP_ros", full.names=TRUE)
# Double-check the filenames in the object
lf_ros
# Import the binary mask jpegs using the file list
cass_all_ros <- import_jpg(lf_ros)
# Convert them to outline, simultaneously adding factors to the outlines
cass_all_ros_out <- Out(cass_all_ros, fac=cass_all_info_ros)
# Interpolate to the average number of points stated in the object information
cass_all_ros_out
cass_all_ros_out_int <- coo_interpolate(cass_all_ros_out,n=3490)
# Double check the point numbers are correct
cass_all_ros_out_int
pile(cass_all_ros_out)
stack(cass_all_ros_out)


################### Landmark Assignment and Outline Processing###################
## There does not seem to be a way to go back through the specimens or edit points.
## Landmarks must be assigned in the same order for each specimen (I did right first, left second)
cass_all_ros_ldk <- def_ldk(cass_all_ros_out_int, 2)
# Inspect to make sure all specimens have 2 landmarks assigned
cass_all_ros_ldk
cass_all_ros_ldk$ldk
# Plot the outlines in their unaligned state, landmarks are red
stack(cass_all_ros_ldk)

# There are only two landmarks on the outlines at the ventral ends of the casques, so a Procrustes alignment cannot work. Bookstein coordinate alignment can, as it works on only two landmarks. These will be used to establish a new baseline calculated from the first and last landmarks (or just first and second if there are only two), which is exactly what is needed.
cass_all_ros_ldk_bk <- coo_bookstein(cass_all_ros_ldk)
# Plot the aligned outlines to make sure the outlines are aligned along the same baseline
stack(cass_all_ros_ldk_bk) #(they are)
pile(cass_all_ros_ldk_bk)
# The start points for the outlines are inconsistent. Here I create a single origin for each outline
# Write the coordinates out to a csv to create a vector for the coo_slide function, identifying which point is closest to the origin (0,0) of the aligned outlines using the csv in Microsoft Excel
write.csv(cass_all_ros_ldk_bk$coo,file='cass_all_ros_ldk_bk_coordinates.csv')
## Write a vector for sliding the coordinate start point to the same x-y coordinate for all outlines
## Create the list of points to shift all outline starting points to the same location

# slide vector construction
# read in lateral coordinate data
data_ros<-read.csv("cass_all_ros_ldk_bk_coordinates.csv",header=T,row.names=1)
# check the data
head(data_ros)
# create blank vector for slide vector construction
slide_vector_ros<-c()
# for loop to generate the individual dataframes
for(n in seq(1,277,2)){
  data_xy <- data_ros[ ,c(n:(n+1))]
  data_xy
  data_xy_sub <- subset(data_xy, data_xy[ ,2]<0.2)

  min_value <- min(abs(data_xy_sub[ ,1]))
  # min_value
  name_info <- which(data_xy_sub==min_value,arr.ind=T)
  # name_info
  # if statement in case the minimum x coordinate is negative, producing an empty matrix
  if(length(name_info)==0){
    min_value <- min_value*-1
    name_info <- which(data_xy_sub==min_value, arr.ind=T)
    # print('poop')
  }
  name_only <- row.names(name_info)
  # name_only
  name_only_num <- as.numeric(name_only)
  name_only_num
  slide_vector_ros <- c(slide_vector_ros,name_only)
}

slide_vector_ros

# Confirm that the length of the vector
length(slide_vector_ros)
print(slide_vector_ros)
# coerce the vector to be a numeric type instead of string
slide_vector_ros <- as.numeric(slide_vector_ros)
# write to csv in to preserve for spot-checking with exported coordinate data csv
write.csv(slide_vector_ros,"slide_vector_ros.csv")



# Slide the start points of all outlines to ~ (0,0) using coo_slide
cass_all_ros_ldk_bk_slid <- coo_slide(cass_all_ros_ldk_bk,id=slide_vector_ros)
# Verify the start points are relatively the same for all specimens
stack(cass_all_ros_ldk_bk_slid)
pile(cass_all_ros_ldk_bk_slid)
# Save the data object for ease of analysis
save(cass_all_ros_ldk_bk_slid,file='cass_all_ros_ldk_data.RData')




###### 09/18/2022 SEPRATE code
# NEED TO SAVE the workspace as to keep out of windows for thermal protective puropses
save.image("landmarking_space_1.RData")



load("landmarking_space_1.RData")
