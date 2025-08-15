
############################################################################################################################
# SUPPORTING INFORMATION A

# Title: Accessibility drives research efforts on Amazonian sarcosaprophagous flies
# Authors: Bruna L.B. Façanha1,2*, Raquel L. Carvalho3, Rony P. S. Almeida4,2, Filipe M. França5,6, José R. P. Sousa7, Maria C. Esposito2,6, Leandro Juen2,6

# Journal: Proceedings of the Royal Society B

#1 Advanced Research-Action Center for Conservation and Ecosystem Recovery of the Amazon (CAPACREAM), Federal University of Amapá (UNIFAP), Josmar Chaves Pinto Highway, Km2, 68903-419, Macapá, AP, Brazil
#2 Graduate Program in Zoology (PPGZOOL), Institute of Biological Sciences (ICB), Federal University of Pará (UFPA), Augusto Corrêa Street, 01, 66075-110, Belém, PA, Brazil
#3 Institute of Advanced Studies (IEA), University of São Paulo (USP), Praça do Relógio Street, 109, 05508-050, São Paulo, SP, Brazil 
#4 Laboratory of Invertebrate Biodiversity and Ecology (LABIN), Department of Biosciences, Federal University of Sergipe (UFS), Av. Vereador Olímpio Grande, s/n, 49506-036, Itabaiana, SE, Brazil
#5 School of Biological Sciences, University of Bristol, 24 Tyndall Avenue, Bristol, BS8 1TQ, UK
#6 Graduate Program in Ecology (PPGECO), Institute of Biological Sciences (ICB), Federal University of Pará (UFPA), Augusto Corrêa Street, 01, 66075-110, Belém, PA, Brazil
#7 Laboratory of Environmental Sciences and Biodiversity, Center of Agrarian Sciences, State University of Maranhão, Paulo VI University City – Lourenço Vieira da Silva Avenue, 1,000, 65.055-310, São Luís, MA, Brazil

# *Corresponding author of this script: Bruna L. B. Façanha, brubsbrr@gmail.com

# STEPS IN THIS SCRIPT:
# 1 – Random Forest model at the family level and null model, including a loop for automatic execution.
# 2 – Random Forest model at the species level, including a loop for automatic execution.
############################################################################################################################



### DIRECTORY AND PACKAGES ----
# First, clear the workspace and set your working directory::
rm(list=ls()); gc()

# Install and load the R packages required to run the analysis:
needed_packages<-c("data.table", "devtools", "dismo", "doParallel", "dplyr", 
                   "foreach", "ggspatial", "ggplot2", "glmnet", "kernlab",
                   "magrittr", "mapview",  "maxnet", "mgcv", "parallel", "plyr", "purrr",
                   "randomForest", "raster",  "rgdal", "sf", "sp", "splitTools", "stringr", "stats", "terra")

#new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(needed_packages, require, character.only = TRUE) 
removeTmpFiles(0.00001) # remove temporary files
rm(needed_packages, new.packages)

# Before you begin, set your working directory and create a basic directory structure:
dir.create("Predictors", showWarnings=F)# to store predictor layers in tif format
dir.create("Datasets", showWarnings=F) #to store csv files, including raw data and some output
# Note: other directories will be created throughout the script.




# STEP 1 - GET THE TRAINING AND TEST DATA FOLDS ----

rm(list=ls()); gc()

# Read txt file with occurrence data:
occ <- data.table::fread("Datasets/occ_AM.csv", h = T, stringsAsFactors = F)



### 1.1 - SELLING CATEGORICAL DATA: ----
# Get the number of sampling locations per family
occ$model<-paste0(occ$familia,"_")
occ_per_familia<-occ[, .(.N), by=.(model)]
#write.table(occ_per_familia, "occc_per_familia.txt")



### 1.2 - PERFORMING PCA AND TESTING FOR MULTICOLLINEARITY OF VARIABLES ----
# Load environmental variables based on orthogonal PCA axes (PCA variables):
envT <- raster::stack(list.files("Predictors/", full.names=T))
names(envT)
names(envT)<-c("Degradation", "DrySeasonLength", "LandOwnership", "DistanceResearch", "TravelTime")


# Reclassify raster values and set 65535 as NA (For traveltime):
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE)

# Convert raster to points:
envT_df<-rasterToPoints(x=envT, spatial=FALSE)
envT_df[is.nan(envT_df)]<-NA # convert NaN values to NA
envT_df<-as.data.table(envT_df)


# Rename columns (predictors) to improve interpretability:
names(envT_df)<-c("x", "y", "Degradation", "DrySeasonLength", "LandOwnership", "DistanceResearch", "TravelTime")
summary(envT_df)

# Remove pixels outside the Legal Amazon:
envT_df<-envT_df[!is.na(envT_df$TravelTime),]


# Check the multicollinearity of the background area:
Vif <- usdm::vif(envT_df[,3:ncol(envT_df)])
data.table::fwrite(Vif, file="Datasets/vif_camadas_nulo.csv")

names(envT_df)


### 1.3 - CREATE SHAPE/LAYERS FOR THE AMAZON BASED ON THE SAMPLE POINTS ----
# Create the shapefile representing the grid cells of your study area:

envT_df$CompleteObs<-rowSums(envT_df[,3:7], na.rm=F)#creating a new column 'CompleteObs', which is the sum columns (3-7) for each row.
length(which(is.na(envT_df$CompleteObs)))/nrow(envT_df) #calculating the proportion of rows with missing data
grid_cells<-envT_df[,c(1:2,8)] #creating a new 'grid_cells' object with coordinates and 'CompleteObs'
coordinates(grid_cells)<- ~ x + y   #convert to spatial points
gridded(grid_cells)<-TRUE  # gridify your point set:
plot(grid_cells)
#save(grid_cells, file="figure_defull/SampleableCells.png")

# Convert the SpatialPixelsDataFrame to a raster object and assign a coordinate reference system (CRS):
sampleable_cells <- raster(grid_cells) #converting 'grid_cells' to raster
raster::crs(sampleable_cells) <- crs(envT) #setting the CRS to that of the variables
sampleable_cells<-raster::reclassify(sampleable_cells, cbind(0, +Inf, 1), right=TRUE) #Reclassifies the raster values: any value greater than 0 will be converted to 1
save(sampleable_cells, file="RData/SampleableCells.RData")


# Mask all predictors to consider only the sampled cells:
cl <- parallel::makePSOCKcluster(parallel::detectCores()-1) #Number of cores on the computer
doParallel::registerDoParallel(cl)
getDoParWorkers()


# rebuild the projection raster to match the exact extent of the sampled cell raster:
# Extract each raster layer in parallel:
load("RData/SampleableCells.RData")
foreach(b = 1:nlayers(envT), .packages = c("raster", "sp", "terra")) %dopar% {
  
  # Load a layer for future projection:
  SingleLayer<-raster::subset(envT, subset=b)
  
  #Crop raster layer to the extent of 'sampleable_cells':
  SingleLayer<-raster::crop(SingleLayer, sampleable_cells)
  
  # Create an empty raster to receive the values from the Single Layer:
  EmptyRaster<-raster(ext=extent(sampleable_cells),
                      crs=crs(sampleable_cells),
                      nrows=dim(sampleable_cells)[1],
                      ncols=dim(sampleable_cells)[2])
  
  # Rewrite the values within the extent of the original predictors:
  RevaluedLayer <- raster::raster(vals=values(SingleLayer),
                                  ext=extent(EmptyRaster),
                                  crs=crs(EmptyRaster),
                                  nrows=dim(EmptyRaster)[1],
                                  ncols=dim(EmptyRaster)[2])
  
  # raster:
  return(writeRaster(RevaluedLayer, filename=paste0("Predictors/", names(envT)[b], ".tif"), format="GTiff", overwrite=TRUE))
  
} # End of foreach:

parallel::stopCluster(cl)


#Replace layers in the predictor folder



### 1.4 - PERFORMING K-FOLD CROSS-VALIDATION FOR THE MODEL ----

envT <- raster::stack(list.files("Predictors/", full.names=T))
names(envT)

# Create a raster of the study area, maintaining the same extent and spatial resolution as the predictor:
StudyArea <- raster::shapefile("Shapefiles/Bullock/Bullock_mask_polygon.shp") # shapefile 
StudyArea_r <- raster::rasterize(StudyArea, y = envT[[1]], background = NA)
plot(StudyArea_r)


occ <- data.table::fread("Datasets/occ_AM.csv", h = T, stringsAsFactors = F)

# Create a new column indicating the presence value:
occ$presence<-1

# Rename columns representing coordinates:
names(occ)[2:3]<-c("x", "y")

# Create an occurrence data list for each family combination:
#occ_xy <- split(occ[,-1], f= paste0(occ$familia))
occ_xy <- split(occ[,-1], f= paste0(occ$familia, "_", occ$group))
spN <- names(occ_xy) # obter nomes de organismos de habitat



### 1.5 - CREATING ABSENCE COUNT (EQUAL TO THE NUMBER OF PRESENCES) ----

Fast_Random_Points <- function(r, n, p) {
  v <- raster::getValues(r) # Get values from the raster:
  v.notNA <- which(!is.na(v)) # remove NA
  v <- v.notNA[!v.notNA%in%raster::cellFromXY(r,p)] # Select pixels not occupied by an occurrence:
  v <- sample(v, n) # Randomly sample n pixels from v:
  pts <- raster::xyFromCell(r, v) # Return as coordinates:
  return(pts)
}


# Get the random point set for each group:
absences<-list()
for(i in 1:length(spN)){
  
# Set the seed to allow reproducibility:
set.seed(123)
  

# Create the same number of pseudo-absences as presences:
new_absences <- Fast_Random_Points(r=StudyArea_r, # máscara pseudo-absences
                                     n=nrow(occ_xy[[i]]), # n pseudo-absences
                                     p=occ_xy[[i]] # list presences
                                     )
  
# Create a new column indicating that the presence is false:
new_absences<-as.data.frame(new_absences)
new_absences$familia<-unique(occ_xy[[i]]$familia)
new_absences$presence<-0
  
# Store in a list:
  absences[[i]]<-new_absences
  rm(new_absences)
  
}


names(absences)<-spN
#Did it create four lists, one for each family, with the same number of presences? If so, okay.



### 1.6 - CREATING TRAINING AND TEST PARTITIONS ----
# For each group, create training and test data partitions for use in cross-validation:
TrainingData<-list()
TestingData<-list()

for(i in 1:length(spN)){
  
# Set the seed to allow reproducibility:
set.seed(123)
  
# Create empty lists to store the training and test data for group 'i':
  
TR_Data<-list()
TS_Data<-list()
  
# Get a vector of row numbers used in each partition for group 'i':
k_index_1 <- caret::createFolds(y = 1:nrow(occ_xy[[i]]), k = 5, list=TRUE, returnTrain=TRUE) # presences
k_index_0 <- caret::createFolds(y = 1:nrow(absences[[i]]), k = 5, list=TRUE, returnTrain=TRUE) # absences

for(k in 1:length(k_index_1)){
    
# Store the training fold 'k' for group 'i':
  TR_presences<-occ_xy[[i]][k_index_1[[k]],]
  TR_absences<-absences[[i]][k_index_0[[k]],]
    
# test:
  TS_presences<-occ_xy[[i]][-k_index_1[[k]],]
  TS_absences<-absences[[i]][-k_index_0[[k]],]
    
# store each k partition:
  TR_Data[[k]]<-data.frame(rbind(TR_presences, TR_absences), Kfold=k)
  TS_Data[[k]]<-data.frame(rbind(TS_presences, TS_absences), Kfold=k)
    
# Remove unnecessary objects before the next iteration:
  rm(TR_presences, TR_absences, TS_presences, TS_absences)
  
}

# TR_Data e TS_Data:
names(TR_Data)<-c("Kfold1", "Kfold2", "Kfold3","Kfold4", "Kfold5")
names(TS_Data)<-c("Kfold1", "Kfold2", "Kfold3","Kfold4", "Kfold5")
  

TrainingData[[i]]<-TR_Data
TestingData[[i]]<-TS_Data


rm(TR_Data, TS_Data, k_index_1, k_index_0)

}
names(TrainingData)<-spN
names(TestingData)<-spN

# Export:
save(TestingData, TrainingData, file="RData/TR_and_TS_Data.RData")







# STEP 2 - BUILDING RANDOM FOREST MODELS FOR KNOWLEDGE PROBABILITY ----

rm(list=ls()); gc()

### 2.1 - Load data ----
# Load environmental variables:
envT <- raster::stack(list.files("Predictors/", full.names=T))
names(envT)
envT[[4]]<-as.factor(envT[[4]]) # converte 'LandOwnership' em fator
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) #traveltime


# Read the shapefile with study area boundaries and mask pixels outside it: 
StudyArea <- raster::shapefile("Shapefiles/Bullock/Bullock_mask_polygon.shp") # set diretory
envT <- raster::mask(x = envT, mask = StudyArea)


# Load training and test data:
load("RData/TR_and_TS_Data.RData")

# Create an empty list to store the model results:
RF_TrainingOutput<-list()
RF_TestingOutput<-list()


### 2.2 - RUNNING MODELS ----

for(i in 1:length(TrainingData)){
  
# Get the environmental data for each training fold k and test fold k of family 'i':
  TR_envData<-lapply(TrainingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  TS_envData<-lapply(TestingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  
# Create an empty list to store the model results:
  RF_training_k<-list()
  RF_testing_k<-list()
  
  for(k in 1:length(TS_envData)){
    
# Link the spatial and environmental datasets:
    MyTRData<-cbind(TrainingData[[i]][[k]], TR_envData[[k]])
    MyTSData<-cbind(TestingData[[i]][[k]], TS_envData[[k]])
    MyTRData<-MyTRData[complete.cases(MyTRData),] # remove NA's (if any)
    MyTSData<-MyTSData[complete.cases(MyTSData),] # remove NA's (if any)
    

    MyTRData$LandOwnership<-as.factor(MyTRData$LandOwnership)
    MyTSData$LandOwnership<-as.factor(MyTSData$LandOwnership)
    MyTRData$LandOwnership<-factor(MyTRData$LandOwnership, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    MyTSData$LandOwnership<-factor(MyTSData$LandOwnership, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    
# Build and store model results for partition 'k' and organism 'i':
    RF_training_k[[k]]<-randomForest::tuneRF(x = MyTRData[,c(6:10)], # predictors
                                             y = MyTRData[,4], # response
                                             trace = F,
                                             stepFactor = 1,
                                             ntreeTry = 1000,
                                             doBest = T, # repeat model using tuning parameters
                                             plot = F,
                                             importance = T, # store variable importance
    )
    
# Add a new list element containing the model outputs for each data fold:
    RF_training_k[[k]][[ (length(RF_training_k[[k]])+1) ]]<-data.frame(tree=1:length(RF_training_k[[k]]$mse),
                                                                       group=names(TrainingData)[i],
                                                                       kfold=k,
                                                                       mse=RF_training_k[[k]]$mse,
                                                                       rsq=RF_training_k[[k]]$rsq)
    names(RF_training_k[[k]])[length(RF_training_k[[k]])]<-"RF_output_k" # name the last element of the list
    
    # the variable importance calculated for each k-fold:
    RF_training_k[[k]][[ (length(RF_training_k[[k]])+1) ]]<-as.data.frame(RF_training_k[[k]]$importance)
    names(RF_training_k[[k]])[length(RF_training_k[[k]])]<-"importance_k" # name the last element of the list
    
    # Predictions for the test data:
    RF_testing_k[[k]]<-data.frame(group=names(TrainingData)[i],
                                  MyTSData, 
                                  Predicted=predict(object=RF_training_k[[k]], newdata=MyTSData[,c(6:10)], type="response"))
    
  } # finish do k for loop
  
  
  
  
# Store results for organism 'i':
  RF_TrainingOutput[[i]]<-randomForest::combine(RF_training_k[[1]], # get the original lists
                                                RF_training_k[[2]],
                                                RF_training_k[[3]], 
                                                RF_training_k[[4]],
                                                RF_training_k[[5]])
  
# Add a new list element containing the model outputs for all data folds:
  RF_TrainingOutput[[i]] [[ (length(RF_TrainingOutput[[i]])-1) ]] <-data.frame(rbind(
    RF_training_k[[1]][[(length(RF_training_k[[1]])-1)]],
    RF_training_k[[2]][[(length(RF_training_k[[2]])-1)]],
    RF_training_k[[3]][[(length(RF_training_k[[3]])-1)]],
    RF_training_k[[4]][[(length(RF_training_k[[4]])-1)]],
    RF_training_k[[5]][[(length(RF_training_k[[5]])-1)]]))
  
# Add a new list element containing variable importance values for all data folds:
  RF_TrainingOutput[[i]] [[ length(RF_TrainingOutput[[i]]) ]] <-data.frame(rbind(
    RF_training_k[[1]][[length(RF_training_k[[1]])]],
    RF_training_k[[2]][[length(RF_training_k[[2]])]],
    RF_training_k[[3]][[length(RF_training_k[[3]])]],
    RF_training_k[[4]][[length(RF_training_k[[4]])]],
    RF_training_k[[5]][[length(RF_training_k[[5]])]]))
  
  RF_TestingOutput[[i]]<-rbindlist(RF_testing_k)

  rm(RF_training_k, RF_testing_k, TR_envData, TS_envData, MyTRData, MyTSData)
  
} # finish do loop i for

# Export
names(RF_TrainingOutput)<-c("Calliphoridae_flies", "Mesembrinellidae_flies", "Model_nulo_nulo", "Sarcophagidae_flies")
names(RF_TestingOutput)<-c("Calliphoridae_flies", "Mesembrinellidae_flies", "Model_nulo_nulo", "Sarcophagidae_flies")

save(RF_TrainingOutput, RF_TestingOutput, file="RData/TR_and_TS_ModelOutputs.RData")



 
# STEP 3 - EVALUATING VARIABLE IMPORTANCE AND MODEL PERFORMANCE METRICS ----

rm(list=ls()); gc()

# Load model output for training and test data:
load("RData/TR_and_TS_ModelOutputs.RData")

# Combine into a single dataframe:
TestingOutput<-rbindlist(RF_TestingOutput)

# Results models
ModelPerformance <- lapply(RF_TrainingOutput, purrr::pluck, "RF_output_k")


# Perform model evaluation for the complete set of predicted values using the test data:
### 3.1 - PERFORMING MODEL EVALUATION AND SAVING ----

ModelEvaluation<-list()

TestingOutput$group<-as.factor(TestingOutput$group)
for(i in 1:nlevels(TestingOutput$group)){
  
  ModelEvaluation[[i]]<-ENMTML:::Eval_Jac_Sor_TMLA(
    p = TestingOutput[TestingOutput$group==levels(TestingOutput$group)[i] & TestingOutput$presence==1,]$Predicted,
    a = TestingOutput[TestingOutput$group==levels(TestingOutput$group)[i] & TestingOutput$presence==0,]$Predicted,
    thr = NULL)
  
  names(ModelEvaluation)[i]<-levels(TestingOutput$group)[i]
  
}


##Report the limitations of the models:
ModelEvaluation$Calliphoridae_flies$SorensenTHR
ModelEvaluation$Mesembrinellidae_flies$SorensenTHR
ModelEvaluation$Model_nulo_nulo$SorensenTHR
ModelEvaluation$Sarcophagidae_flies$SorensenTHR

###report  sorensenof the models
max(ModelEvaluation$upland_Calliphoridae$Sorensen, na.rm=T)
max(ModelEvaluation$upland_Mesembrinellidae$Sorensen, na.rm=T)
max(ModelEvaluation$upland_Sarcophagidae$Sorensen, na.rm=T)

max(ModelEvaluation$Calliphoridae_flies$Sorensen, na.rm=T)
max(ModelEvaluation$Mesembrinellidae_flies$Sorensen, na.rm=T)
max(ModelEvaluation$Model_nulo_nulo$Sorensen, na.rm=T)
max(ModelEvaluation$Sarcophagidae_flies$Sorensen, na.rm=T)




#TestingOutput

######exporta matrix
as.matrix(ModelEvaluation)->ModelEvaluation
write.table (ModelEvaluation, file="Datasets/ModelEvaluation.csv")
save(ModelEvaluation, file="RData/ModelEvaluation.RData")

load("RData/ModelEvaluation.RData")

ModelEvaluation<-rbindlist(ModelEvaluation)
ModelEvaluation
ModelEvaluation<-ModelEvaluation[, .(.N,
                                     avg_Sorensen = mean(Sorensen, na.rm=T),
                                     min_Sorensen = min(Sorensen, na.rm=T),
                                     max_Sorensen = max(Sorensen, na.rm=T),
                                     sd_Sorensen = sd(Sorensen, na.rm=T),
                                     avg_Jaccard= mean(Jaccard, na.rm=T),
                                     min_Jaccard = min(Jaccard, na.rm=T),
                                     max_Jaccard = max(Jaccard, na.rm=T),
                                     sd_Jaccard = sd(Jaccard, na.rm=T)),
                                 by=.(presences,absences)] 
ModelEvaluation

write.table (ModelEvaluation, file="Datasets/ModelEvaluation_medias_final.csv")



### 3.2 - EVALUATING VARIABLE IMPORTANCE ----
# Calculate the average performance metrics of the model for each group in each iteration (tree) and fold (k)::
ModelPerformance<-rbindlist(ModelPerformance)
ModelPerformance
ModelPerformance<-ModelPerformance[, .(.N,
                                        avg_mse = mean(mse, na.rm=T),
                                        min_mse = min(mse, na.rm=T),
                                        max_mse = max(mse, na.rm=T),
                                        sd_mse = sd(mse, na.rm=T),
                                        avg_rsq = mean(rsq, na.rm=T),
                                        min_rsq = min(rsq, na.rm=T),
                                        max_rsq = max(rsq, na.rm=T),
                                        sd_rsq = sd(rsq, na.rm=T)),
                                    by=.(group, tree)] 

ModelPerformance<-ModelPerformance[ModelPerformance$tree==500,]
ModelPerformance


VarImportance <- lapply(RF_TrainingOutput, purrr::pluck, "importance_k")
names(VarImportance)<-ModelPerformance$group 

# Extract variable importance metrics by group:
VarImp_Output<-data.frame()
for(i in 1:length(VarImportance)){
  
  VarImportance[[i]]$Kfold<-c(rep(1,5), rep(2,5), rep(3,5), rep(4,5), rep(5,5))
  VarImportance[[i]]$Variables<-rep(row.names(RF_TrainingOutput[[1]]$importance), 5)
  
  for(k in 1:5){
    
    MyData<-VarImportance[[i]][VarImportance[[i]]$Kfold==k,]
    MyData$group<-ModelPerformance$group[i]
    VarImp_Output<-rbind(VarImp_Output, MyData)
  }
}

write.table (VarImportance, file="Datasets/Importanciadasvariaveis_medias.csv")
rm(VarImportance)

### 3.3 - PREPARING THE DATA FRAME FOR PLOTTING: ----

VarImp_Output<-as.data.table(VarImp_Output)
VarImportance<-VarImp_Output[, .(avg_INP = mean(IncNodePurity, na.rm=T),
                                     min_INP = min(IncNodePurity, na.rm=T),
                                     max_INP = max(IncNodePurity, na.rm=T),
                                     sd_INP = sd(IncNodePurity),
                                     avg_IMSE = mean(X.IncMSE, na.rm=T),
                                     min_IMSE = min(X.IncMSE, na.rm=T),
                                     max_IMSE = max(X.IncMSE, na.rm=T),
                                     sd_IMSE = sd(X.IncMSE, na.rm=T)),
                                 by=.(group, Variables)]

# Output:
VarImportance


### 3.4 - VARIABLE IMPORTANCE FIGURE - BARS ----

names(VarImportance$group)
# Relabel levels to improve visualization:
VarImportance$Families <- factor(VarImportance$group,
                               levels = c("Calliphoridae_flies", "Mesembrinellidae_flies", "Model_nulo_nulo", "Sarcophagidae_flies"),
                               labels = c("Calliphoridae", "Mesembrinellidae ", "Modelo nulo", "Sarcophagidae"))

VarImportance$Variables1 <- factor(VarImportance$Variables, 
                                   levels=rev(c("Degradation", "DistanceResearch","DrySeasonLength", "LandOwnership", "TravelTime" )),
                                   labels=rev(c("Degradation", "Distance to\nResearch centers", "Dry season\nLength", "Land\nOwnership", "Travel Time")))


# Plot the normalized IncNodePurity for each variable and group:
MyPlot1<-
ggplot(data=VarImportance, aes(x=Families, y=avg_IMSE, shape=Families)) +
  
  # Add geopoints:
  geom_point(shape=16, aes(colour=Families, shape=Families, fill=Families), size=2) + 
  
  # Add error bars representing min and max values of variable importance:
  geom_errorbar(aes(ymin=min_IMSE, ymax=max_IMSE, col=Families), width=0.5, cex=1) +
  
  # Set the colour ramp:
  scale_color_manual(values = c("#8dd3c7", "#bebada", "black", "#fdc086"), labels = c("Calliphoridae", "Mesembrinellidae ",  "Modelo nulo", "Sarcophagidae")) +
  
  # Split geom_points according to levels of 'Variables':
  facet_wrap(~Variables1, strip.position="left", nrow=5, scales = "free_y") +
  
  # Axes titles:
  xlab("") + ylab('% Increase in MSE (Normalized)') +
  
  # Other aesthetics:
  coord_flip() +  # flip coordinates (puts labels on y axis)
  #theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold"),
        axis.title.x=element_text(size=12, colour="black", face="bold"),
        legend.position="top",
        #plot.background=element_rect(fill="transparent", colour=NA),
        strip.placement = "outside")

MyPlot1

ggsave("figure_defull/VariableImportance.png", plot=MyPlot1, width=10, height=5)
#ggsave("Figures/VariableImportance2.pdf", plot=MyPlot1, width=8, height=5, units="in", bg="transparent")

# Export:
data.table::fwrite(VarImportance, file="Datasets/VarImportance.csv")



# STEP 4 - INSPECTING THE BEHAVIOR OF RESEARCH PROBABILITY AMONG PREDICTOR VALUES (PARTIAL PLOTS) ----


rm(list=ls()); gc()

# Load training and test data for the model:
load("RData/TR_and_TS_Data.RData")

# Load the model output for training and test data:
load("RData/TR_and_TS_ModelOutputs.RData")

# Load environmental variables:
envT <- raster::stack(list.files("Predictors/", full.names=T))

envT[[4]]<-as.factor(envT[[4]]) # convert 'LandOwnership' to factor
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) # TravelTime - reclassifica os valores raster e define o valor 65535 como NA


# organism 'i':
Partial_Plots<-list()
for(i in 1:length(TestingData)){ 
  
  # Extract environmental data for all presence/pseudo-absence records:  
  All_envData<-lapply(TestingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  All_envData<-lapply(All_envData, as.data.frame)
  
  #Combine all data folds used for training:
  MyData<-cbind(rbindlist(TestingData[[i]]), rbindlist(All_envData))
  MyData<-MyData[complete.cases(MyData),] # remove NA's (if any)
  MyData$LandOwnership<-as.factor(MyData$LandOwnership)
  MyData$LandOwnership<-factor(MyData$LandOwnership, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  
  # randomForest for organism 'i':
  SelectedRF_Model<-RF_TrainingOutput[[i]]
  
  # list for models results :
  myPlots<-list()
  PredictorNames<-c("Degradation", "DistanceResearch","DrySeasonLength", "LandOwnership", "TravelTime")
  # OBS: PredictorNames[4] refere-se à variável categórica 'LandTenure'
  
  # Partial plots for all variables:
  for(v in 1:5){ # v number of variables
    
    # Gráficos para preditores contínuos:
    if(v!=4){
      
      myPlots[[v]]<-partialPlot(x = SelectedRF_Model, # randomForest
                                pred.data = MyData, # training data for the respective randomForest
                                x.var = PredictorNames[v], # variables
                                which.class = NULL, 
                                plot = TRUE, 
                                add = FALSE,
                                rug = TRUE, 
                                ylab="Research Presence",
                                xlab=PredictorNames[v]
      )
      
      myPlots[[v]]$Group<-names(TestingData)[i]
      myPlots[[v]]$Predictor<-PredictorNames[v]
      
    }
    
    # Plots for categorical predictors:
    if(v==4){
      
      myPlots[[v]]<-partialPlot(x = SelectedRF_Model, # object of class randomForest
                                pred.data = MyData, # training data for the respective randomForest
                                x.var = PredictorNames[v], # variable name to be examined
                                which.class = NULL, # for classification data, define which class
                                plot = TRUE, 
                                add = FALSE,
                                rug = TRUE, 
                                ylab="Research Presence",
                                xlab=PredictorNames[v]
      )
      
      myPlots[[v]]$Group<-names(TestingData)[i]
      myPlots[[v]]$Predictor<-PredictorNames[v]
    }
    
  }
  
  Partial_Plots[[i]]<-myPlots
  
}

# Partial_Plots is a set of 4 lists (groups), with 5 results (predictors), one for each group:
# Partial_Plots[[1]] contains five lists, one for each predictor variable
# Each list has the x value and y values of the partial graph for the given variable
Partial_Plots[[1]]



# Export to disk:
save(Partial_Plots, file="RData/PartialPlotsData.RData")



load("RData/PartialPlotsData.RData")



data_frame <- data.frame(
  Coluna1 = unlist(Partial_Plots[[1]]),
  Coluna2 = unlist(Partial_Plots[[2]]),
  Coluna3 = unlist(Partial_Plots[[3]]),
  Coluna4 = unlist(Partial_Plots[[4]])
)

write.csv(data_frame, "Datasets/partial_plots.csv", row.names = FALSE)


# STEP 5 - PREDICT GEOGRAPHICAL PATTERNS OF KNOWLEDGE PROBABILITY----
rm(list=ls()); gc()

# Load model output for training and testing data:
load("RData/TR_and_TS_ModelOutputs.RData")

# Load raster data:
envT <- raster::stack(list.files("Predictors/", full.names=T))
envT[[4]]<-as.factor(envT[[4]]) # convert 'LandOwnership' to factor
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) # reclassify TravelTime

# dataframe:
envT_df<-rasterToPoints(x=envT, spatial=FALSE)
envT_df[is.nan(envT_df)]<-NA # convert NaN values to NA
envT_df<-as.data.table(envT_df) # convert to data.table
summary(envT_df)
envT_df<-envT_df[!is.na(envT_df$Degradation),] # remove NAs
envT_df<-envT_df[!is.na(envT_df$DistanceResearch),] # remove NAs
envT_df<-envT_df[!is.na(envT_df$DrySeasonLength),] # remove NAs
envT_df<-envT_df[!is.na(envT_df$TravelTime ),] # remove NAs
envT_df$LandOwnership<-as.factor(envT_df$LandOwnership) # set 'LandTenure' as factor
summary(envT_df)

### 5.1 - PREDICTING THE PROBABILITY OF KNOWLEDGE FOR EACH FAMILY----
#This sets R to use several cores of your computer (minus 1) to speed up the process
cl<-makePSOCKcluster(detectCores()-1, type="SOCK")
registerDoParallel(cl)
getDoParWorkers()

#Using the 'predict()' function to predict the probability of presence in each cell of the dataframe 'envT_df':
MyPredictions<-foreach(i = 1:length(RF_TrainingOutput), 
                       .export = 'c', 
                       .packages = c("raster", "sp", "randomForest")) %dopar% {
                         predict(object = RF_TrainingOutput[[i]], newdata = envT_df,  type = "response")
                       }
parallel::stopCluster(cl)

group_names<-unique(rbindlist(RF_TestingOutput)$group) #Obtaining the names of the groups/families

#'MyPredictions': list with the predicted probabilities for each cell of the variables and for each family
#'group_names': names of these groups, useful for naming the rasters or maps later.


# convert SpatialPixelsDataFrame:
MyPixelDF<-(data.frame(envT_df[,1:2], 
                       Calliphoridae_flies=MyPredictions[[1]],
                       Mesembrinellidae_flies=MyPredictions[[2]],
                       Model_nulo_nulo=MyPredictions[[3]],
                       Sarcophagidae_flies=MyPredictions[[4]]))

coordinates(MyPixelDF)<- ~ x + y # convert to spatial points
gridded(MyPixelDF)<-TRUE  # griddify your points

plot(MyPixelDF)

# Convert rasters:
MyRasters<-list()
for(i in 1:ncol(MyPixelDF)){
  MyRasters[[i]] <- raster(MyPixelDF, layer=i) 
}
names(MyRasters)<-c("Calliphoridae_flies", "Mesembrinellidae_flies", "Model_nulo_nulo", "Sarcophagidae_flies")
plot(MyRasters)

par(mfrow = c(2, 2)) # 2 linhas, 2 colunas

plot(MyRasters$Calliphoridae_flies)
plot(MyRasters$Mesembrinellidae_flies)
plot(MyRasters$Model_nulo_nulo)
plot(MyRasters$Sarcophagidae_flies)

par(mfrow = c(1, 1)) # reseta para o padrão (opcional)

### 5.2 - CREATING FORECAST RASTER FOR EACH FAMILY ----

## empilhando os rasters
RasterStack<-raster::stack(MyRasters[[1]])
for(i in 2:length(MyRasters)){RasterStack<-raster::stack(RasterStack, MyRasters[[i]])} # faster than using do.call
names(RasterStack)<-names(MyRasters)
raster::crs(RasterStack) <- raster::crs(envT)

plot(RasterStack)


#Remove essas quatro camadas específicas do RasterStack, criando um novo objeto chamado FliesRasters
FliesRasters<-dropLayer(RasterStack, c("Calliphoridae_flies",  
                                       "Mesembrinellidae_flies", 
                                       "Model_nulo_nulo",
                                       "Sarcophagidae_flies"))


plot(FliesRasters)
#Máscara (remoção das áreas fora do habitat)
# raster upland:
UlpandMasck<-raster::raster("rasters/upland.tif")
plot(UlpandMasck)

# use mask layers to remove predictions outside the habitat:
ulpandflies <- raster::mask(RasterStack, UlpandMasck)
#Resultado: um novo raster (ulpandflies) com apenas os pixels dentro da área válida (como um "recorte").

plot(ulpandflies)

# Export the knowledge probability rasters:
dir.create("Projections/", showWarnings = F)
for(i in 1:nlayers(ulpandflies)){terra::writeRaster(terra::rast(ulpandflies[[i]]), filename=paste0("Projections/", names(ulpandflies)[i], ".tif"), overwrite=TRUE)}



# STEP 6 - GET AVERAGE PROBABILITY FROM RESEARCH AND ACROSS ALL ORGANISMS ----

rm(list=ls()); gc()

#listar arquivos
RasterFiles<-list.files(path=paste0("Projections/"), pattern='.tif', full.names=T)
Filenames<-gsub("Projections/", "", RasterFiles)
Filenames<-gsub("\\.tif", "", Filenames)


## para as moascas
# # Separate the raster layers available for each family:
#seleciona os rasters só com flies no nome
FliesFiles<-RasterFiles[RasterFiles %like% "flies"]


# raster:
#Carrega os rasters de moscas
FliesLayers<-lapply(FliesFiles, raster::raster)


# Convert to raster stack and rename layers:
FliesStack<-raster::stack(FliesLayers)
names(FliesStack)<-Filenames[Filenames %like% "flies"]


# Get the average probability of knowledge for each family:
#Calcula a média da probabilidade entre os grupos de moscas
ResearchProb_Flies<-raster::calc(stack(FliesStack), fun=mean, na.rm=T) # it is important to set na.rm=T
names(ResearchProb_Flies)<-"Flies"
plot(ResearchProb_Flies)


# Export the knowledge probability rasters:
dir.create("AvgOutputs/", showWarnings = F)
terra::writeRaster(ResearchProb_Flies, filename=paste0("AvgOutputs/Mean_FLIES.tif"), overwrite=TRUE)



# Obtain the average probability of knowledge for the combination of families:
RasterLayers<-lapply(FliesFiles, raster::raster)
names(FliesStack)<-Filenames[Filenames %like% "flies"]

# Get the average search probability for upland
PixelDataFrame<-list()
for(i in 1:length(RasterLayers)){
  PixelDataFrame[[i]]<-raster::as.data.frame(RasterLayers[[i]], xy=TRUE, na.rm=T)
  PixelDataFrame[[i]]$HabitatOrganism<-names(RasterLayers)[i]
  names(PixelDataFrame[[i]])[3]<-"ResearchProb"
}

# data.table:
PixelDataFrame<-rbindlist(PixelDataFrame)

# Get the average values across all pixels:
PixelDataFrame$x<-round(PixelDataFrame$x, digits = 3)
PixelDataFrame$y<-round(PixelDataFrame$y, digits = 3)
PixelDataFrame$CoordID<-paste0(PixelDataFrame$x, "_", PixelDataFrame$y)
PixelDataFrame<-unique(PixelDataFrame)

# Calculate the average probability of searching among organisms:
AvgOutput<-PixelDataFrame[, .(.N,
                              Avg_x=mean(x, na.rm=T),
                              Avg_y=mean(y, na.rm=T),
                              AvgResearchProb=mean(ResearchProb, na.rm=T)),
                          by=.(CoordID)] 
names(AvgOutput)<-c("CoordID", "N", "x", "y", "ResearchProb")

AvgOutput$HabitatOrganism<-"All"



# Build a raster for all groups:
grid_cells<-data.frame(ResearchProb=AvgOutput$ResearchProb, AvgOutput[,3:4])
grid_cells$x<-round(grid_cells$x, digits = 2)
grid_cells$y<-round(grid_cells$y, digits = 2)
grid_cells<-unique(grid_cells)
sp::coordinates(grid_cells)<- ~ x + y # convert to spatial points
sp::gridded(grid_cells)<-TRUE  # griddify your set of points

plot(grid_cells)

# raster:
grid_cells <- raster(grid_cells)
raster::crs(grid_cells) <- crs(RasterLayers[[1]])
plot(grid_cells)

# Export raster :
terra::writeRaster(terra::rast(grid_cells), filename=paste0("Outputs/AvgResearchProb_Flies.tif"), overwrite=TRUE)



### para o modelo nulo

rm(list=ls()); gc()

#listar arquivos
RasterFiles<-list.files(path=paste0("Projections/"), pattern='.tif', full.names=T)
Filenames<-gsub("Projections/", "", RasterFiles)
Filenames<-gsub("\\.tif", "", Filenames)


NullFile <- RasterFiles[grepl("nulo", RasterFiles, ignore.case = TRUE)]
NullLayer <- raster::raster(NullFile)
names(NullLayer) <- gsub("Projections/|\\.tif", "", NullFile)


terra::writeRaster(NullLayer, filename = "AvgOutputs/Model_Nulo.tif", overwrite = TRUE)

# Converter para data.frame com coordenadas
NullDF <- raster::as.data.frame(NullLayer, xy = TRUE, na.rm = TRUE)
names(NullDF)[3] <- "ResearchProb"
NullDF$HabitatOrganism <- names(NullLayer)
NullDF$x <- round(NullDF$x, 3)
NullDF$y <- round(NullDF$y, 3)
NullDF$CoordID <- paste0(NullDF$x, "_", NullDF$y)

# Média (nesse caso será só ele mesmo, mas mantém a estrutura igual)
NullDF <- data.table::as.data.table(NullDF)
AvgOutputNull <- NullDF[, .(.N,
                            Avg_x = mean(x, na.rm = TRUE),
                            Avg_y = mean(y, na.rm = TRUE),
                            AvgResearchProb = mean(ResearchProb, na.rm = TRUE)),
                        by = .(CoordID)]
names(AvgOutputNull) <- c("CoordID", "N", "x", "y", "ResearchProb")
AvgOutputNull$HabitatOrganism <- "ModelNulo"


# Build a raster for all groups:
grid_null <- data.frame(ResearchProb = AvgOutputNull$ResearchProb, AvgOutputNull[, 3:4])
grid_null$x <- round(grid_null$x, digits = 2)
grid_null$y <- round(grid_null$y, digits = 2)
grid_null <- unique(grid_null)
sp::coordinates(grid_null) <- ~ x + y
sp::gridded(grid_null) <- TRUE
plot(grid_null)


# raster:
grid_null <- raster(grid_null)
raster::crs(grid_null) <- raster::crs(NullLayer)
plot(grid_null)

terra::writeRaster(terra::rast(grid_null), filename = "Outputs/AvgResearchProb_ModelNulo.tif", overwrite = TRUE)



#####################################################################################
############################### SPECIES ########################

# STEP 1 - GET THE TRAINING AND TEST DATA FOLDS ----

rm(list=ls()); gc()

# Read txt file with occurrence data:
occ <- data.table::fread("Datasets_sp/occ_sp.csv", h = T, stringsAsFactors = F)



### 1.1 - SELLING CATEGORICAL DATA: ----
# Get the number of sampling locations per family
occ$group<-paste0(occ$sp,"_")
occ_per_sp<-occ[, .(.N), by=.(group)]
#write.table(occ_per_sp, "occc_per_especie.txt")



### 1.2 - PERFORMING PCA AND TESTING FOR MULTICOLLINEARITY OF VARIABLES ----
# Load environmental variables based on orthogonal PCA axes (PCA variables):
envT <- raster::stack(list.files("Predictors_sp/", full.names=T))
names(envT)
names(envT)<-c("Degradation","DistanceResearch", "DrySeasonLength", "LandOwnership", "TravelTime")


# Reclassify raster values and set 65535 as NA (For traveltime):
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE)

# Convert raster to points:
envT_df<-rasterToPoints(x=envT, spatial=FALSE)
envT_df[is.nan(envT_df)]<-NA # convert NaN values to NA
envT_df<-as.data.table(envT_df)


# Rename columns (predictors) to improve interpretability:
names(envT_df)<-c("x", "y", "Degradation", "DistanceResearch", "DrySeasonLength", "LandOwnership", "TravelTime")
summary(envT_df)

# Remove pixels outside the Legal Amazon:
envT_df<-envT_df[!is.na(envT_df$TravelTime),]


# Check the multicollinearity of the background area:
Vif <- usdm::vif(envT_df[,3:ncol(envT_df)])
data.table::fwrite(Vif, file="Datasets_sp/vif_camadas_nulo_sp.csv")

names(envT_df)


### 1.3 - CREATE SHAPE/LAYERS FOR THE AMAZON BASED ON THE SAMPLE POINTS ----
# Create the shapefile representing the grid cells of your study area:

envT_df$CompleteObs<-rowSums(envT_df[,3:7], na.rm=F)#creating a new column 'CompleteObs', which is the sum columns (3-7) for each row.
length(which(is.na(envT_df$CompleteObs)))/nrow(envT_df) #calculating the proportion of rows with missing data
grid_cells<-envT_df[,c(1:2,8)] #creating a new 'grid_cells' object with coordinates and 'CompleteObs'
coordinates(grid_cells)<- ~ x + y   #convert to spatial points
gridded(grid_cells)<-TRUE  # gridify your point set:
plot(grid_cells)
#save(grid_cells, file="figure_defull/SampleableCells.png")

# Convert the SpatialPixelsDataFrame to a raster object and assign a coordinate reference system (CRS):
sampleable_cells <- raster(grid_cells) #converting 'grid_cells' to raster
raster::crs(sampleable_cells) <- crs(envT) #setting the CRS to that of the variables
sampleable_cells<-raster::reclassify(sampleable_cells, cbind(0, +Inf, 1), right=TRUE) #Reclassifies the raster values: any value greater than 0 will be converted to 1
save(sampleable_cells, file="RData_sp/SampleableCells.RData")


# Mask all predictors to consider only the sampled cells:
cl <- parallel::makePSOCKcluster(parallel::detectCores()-1) #Number of cores on the computer
doParallel::registerDoParallel(cl)
getDoParWorkers()


# rebuild the projection raster to match the exact extent of the sampled cell raster:
# Extract each raster layer in parallel:
load("RData_sp/SampleableCells.RData")
foreach(b = 1:nlayers(envT), .packages = c("raster", "sp", "terra")) %dopar% {
  
  # Load a layer for future projection:
  SingleLayer<-raster::subset(envT, subset=b)
  
  #Crop raster layer to the extent of 'sampleable_cells':
  SingleLayer<-raster::crop(SingleLayer, sampleable_cells)
  
  # Create an empty raster to receive the values from the Single Layer:
  EmptyRaster<-raster(ext=extent(sampleable_cells),
                      crs=crs(sampleable_cells),
                      nrows=dim(sampleable_cells)[1],
                      ncols=dim(sampleable_cells)[2])
  
  # Rewrite the values within the extent of the original predictors:
  RevaluedLayer <- raster::raster(vals=values(SingleLayer),
                                  ext=extent(EmptyRaster),
                                  crs=crs(EmptyRaster),
                                  nrows=dim(EmptyRaster)[1],
                                  ncols=dim(EmptyRaster)[2])
  
  # raster:
  return(writeRaster(RevaluedLayer, filename=paste0("Predictors_sp/", names(envT)[b], ".tif"), format="GTiff", overwrite=TRUE))
  
} # End of foreach:

parallel::stopCluster(cl)


#Replace layers in the predictor folder



### 1.4 - PERFORMING K-FOLD CROSS-VALIDATION FOR THE MODEL ----

envT <- raster::stack(list.files("Predictors_sp/", full.names=T))
names(envT)

# Create a raster of the study area, maintaining the same extent and spatial resolution as the predictor:
StudyArea <- raster::shapefile("Shapefiles/Bullock/Bullock_mask_polygon.shp") # shapefile 
StudyArea_r <- raster::rasterize(StudyArea, y = envT[[1]], background = NA)
plot(StudyArea_r)


occ <- data.table::fread("Datasets_sp/occ_sp.csv", h = T, stringsAsFactors = F)

# Create a new column indicating the presence value:
occ$presence<-1

# Rename columns representing coordinates:
names(occ)[2:3]<-c("x", "y")

# Create an occurrence data list for each family combination:
#occ_xy <- split(occ[,-1], f= paste0(occ$familia))
occ_xy <- split(occ[,-1], f= paste0(occ$sp, "_", occ$group))
spN <- names(occ_xy) # obter nomes de organismos de habitat



### 1.5 - CREATING ABSENCE COUNT (EQUAL TO THE NUMBER OF PRESENCES) ----

Fast_Random_Points <- function(r, n, p) {
  v <- raster::getValues(r) # Get values from the raster:
  v.notNA <- which(!is.na(v)) # remove NA
  v <- v.notNA[!v.notNA%in%raster::cellFromXY(r,p)] # Select pixels not occupied by an occurrence:
  v <- sample(v, n) # Randomly sample n pixels from v:
  pts <- raster::xyFromCell(r, v) # Return as coordinates:
  return(pts)
}


# Get the random point set for each group:
absences<-list()
for(i in 1:length(spN)){
  
  # Set the seed to allow reproducibility:
  set.seed(123)
  
  
  # Create the same number of pseudo-absences as presences:
  new_absences <- Fast_Random_Points(r=StudyArea_r, # máscara pseudo-absences
                                     n=nrow(occ_xy[[i]]), # n pseudo-absences
                                     p=occ_xy[[i]] # list presences
  )
  
  # Create a new column indicating that the presence is false:
  new_absences<-as.data.frame(new_absences)
  new_absences$sp<-unique(occ_xy[[i]]$sp)
  new_absences$presence<-0
  
  # Store in a list:
  absences[[i]]<-new_absences
  rm(new_absences)
  
}


names(absences)<-spN
#Did it create 16 lists, one for each sp, with the same number of presences? If so, okay.



### 1.6 - CREATING TRAINING AND TEST PARTITIONS ----
# For each group, create training and test data partitions for use in cross-validation:
TrainingData<-list()
TestingData<-list()

for(i in 1:length(spN)){
  
  # Set the seed to allow reproducibility:
  set.seed(123)
  
  # Create empty lists to store the training and test data for group 'i':
  
  TR_Data<-list()
  TS_Data<-list()
  
  # Get a vector of row numbers used in each partition for group 'i':
  k_index_1 <- caret::createFolds(y = 1:nrow(occ_xy[[i]]), k = 5, list=TRUE, returnTrain=TRUE) # presences
  k_index_0 <- caret::createFolds(y = 1:nrow(absences[[i]]), k = 5, list=TRUE, returnTrain=TRUE) # absences
  
  for(k in 1:length(k_index_1)){
    
    # Store the training fold 'k' for group 'i':
    TR_presences<-occ_xy[[i]][k_index_1[[k]],]
    TR_absences<-absences[[i]][k_index_0[[k]],]
    
    # test:
    TS_presences<-occ_xy[[i]][-k_index_1[[k]],]
    TS_absences<-absences[[i]][-k_index_0[[k]],]
    
    # store each k partition:
    TR_Data[[k]]<-data.frame(rbind(TR_presences, TR_absences), Kfold=k)
    TS_Data[[k]]<-data.frame(rbind(TS_presences, TS_absences), Kfold=k)
    
    # Remove unnecessary objects before the next iteration:
    rm(TR_presences, TR_absences, TS_presences, TS_absences)
    
  }
  
  # TR_Data e TS_Data:
  names(TR_Data)<-c("Kfold1", "Kfold2", "Kfold3","Kfold4", "Kfold5")
  names(TS_Data)<-c("Kfold1", "Kfold2", "Kfold3","Kfold4", "Kfold5")
  
  
  TrainingData[[i]]<-TR_Data
  TestingData[[i]]<-TS_Data
  
  
  rm(TR_Data, TS_Data, k_index_1, k_index_0)
  
}
names(TrainingData)<-spN
names(TestingData)<-spN

# Export:
save(TestingData, TrainingData, file="RData_sp/TR_and_TS_Data.RData")







# STEP 2 - BUILDING RANDOM FOREST MODELS FOR KNOWLEDGE PROBABILITY ----

rm(list=ls()); gc()

### 2.1 - Load data ----
# Load environmental variables:
envT <- raster::stack(list.files("Predictors_sp/", full.names=T))
names(envT)
envT[[4]]<-as.factor(envT[[4]]) # converte 'LandOwnership' em fator
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) #traveltime


# Read the shapefile with study area boundaries and mask pixels outside it: 
StudyArea <- raster::shapefile("Shapefiles/Bullock/Bullock_mask_polygon.shp") # set diretory
envT <- raster::mask(x = envT, mask = StudyArea)


# Load training and test data:
load("RData_sp/TR_and_TS_Data.RData")

# Create an empty list to store the model results:
RF_TrainingOutput<-list()
RF_TestingOutput<-list()


### 2.2 - RUNNING MODELS ----

for(i in 1:length(TrainingData)){
  
  # Get the environmental data for each training fold k and test fold k of family 'i':
  TR_envData<-lapply(TrainingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  TS_envData<-lapply(TestingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  
  # Create an empty list to store the model results:
  RF_training_k<-list()
  RF_testing_k<-list()
  
  for(k in 1:length(TS_envData)){
    
    # Link the spatial and environmental datasets:
    MyTRData<-cbind(TrainingData[[i]][[k]], TR_envData[[k]])
    MyTSData<-cbind(TestingData[[i]][[k]], TS_envData[[k]])
    MyTRData<-MyTRData[complete.cases(MyTRData),] # remove NA's (if any)
    MyTSData<-MyTSData[complete.cases(MyTSData),] # remove NA's (if any)
    
    
    MyTRData$LandOwnership<-as.factor(MyTRData$LandOwnership)
    MyTSData$LandOwnership<-as.factor(MyTSData$LandOwnership)
    MyTRData$LandOwnership<-factor(MyTRData$LandOwnership, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    MyTSData$LandOwnership<-factor(MyTSData$LandOwnership, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    
    # Build and store model results for partition 'k' and organism 'i':
    RF_training_k[[k]]<-randomForest::tuneRF(x = MyTRData[,c(6:10)], # predictors
                                             y = MyTRData[,4], # response
                                             trace = F,
                                             stepFactor = 1,
                                             ntreeTry = 1000,
                                             doBest = T, # repeat model using tuning parameters
                                             plot = F,
                                             importance = T, # store variable importance
    )
    
    # Add a new list element containing the model outputs for each data fold:
    RF_training_k[[k]][[ (length(RF_training_k[[k]])+1) ]]<-data.frame(tree=1:length(RF_training_k[[k]]$mse),
                                                                       group=names(TrainingData)[i],
                                                                       kfold=k,
                                                                       mse=RF_training_k[[k]]$mse,
                                                                       rsq=RF_training_k[[k]]$rsq)
    names(RF_training_k[[k]])[length(RF_training_k[[k]])]<-"RF_output_k" # name the last element of the list
    
    # the variable importance calculated for each k-fold:
    RF_training_k[[k]][[ (length(RF_training_k[[k]])+1) ]]<-as.data.frame(RF_training_k[[k]]$importance)
    names(RF_training_k[[k]])[length(RF_training_k[[k]])]<-"importance_k" # name the last element of the list
    
    # Predictions for the test data:
    RF_testing_k[[k]]<-data.frame(group=names(TrainingData)[i],
                                  MyTSData, 
                                  Predicted=predict(object=RF_training_k[[k]], newdata=MyTSData[,c(6:10)], type="response"))
    
  } # finish do k for loop
  
  
  
  
  # Store results for organism 'i':
  RF_TrainingOutput[[i]]<-randomForest::combine(RF_training_k[[1]], # get the original lists
                                                RF_training_k[[2]],
                                                RF_training_k[[3]], 
                                                RF_training_k[[4]],
                                                RF_training_k[[5]])
  
  # Add a new list element containing the model outputs for all data folds:
  RF_TrainingOutput[[i]] [[ (length(RF_TrainingOutput[[i]])-1) ]] <-data.frame(rbind(
    RF_training_k[[1]][[(length(RF_training_k[[1]])-1)]],
    RF_training_k[[2]][[(length(RF_training_k[[2]])-1)]],
    RF_training_k[[3]][[(length(RF_training_k[[3]])-1)]],
    RF_training_k[[4]][[(length(RF_training_k[[4]])-1)]],
    RF_training_k[[5]][[(length(RF_training_k[[5]])-1)]]))
  
  # Add a new list element containing variable importance values for all data folds:
  RF_TrainingOutput[[i]] [[ length(RF_TrainingOutput[[i]]) ]] <-data.frame(rbind(
    RF_training_k[[1]][[length(RF_training_k[[1]])]],
    RF_training_k[[2]][[length(RF_training_k[[2]])]],
    RF_training_k[[3]][[length(RF_training_k[[3]])]],
    RF_training_k[[4]][[length(RF_training_k[[4]])]],
    RF_training_k[[5]][[length(RF_training_k[[5]])]]))
  
  RF_TestingOutput[[i]]<-rbindlist(RF_testing_k)
  
  rm(RF_training_k, RF_testing_k, TR_envData, TS_envData, MyTRData, MyTSData)
  
} # finish do loop i for

# Export
names(RF_TrainingOutput)<-c("Chloroprocta_idioidea_Calliphoridae" ,
                            "Chrysomya_albiceps_Calliphoridae" ,
                            "Cochliomyia_macellaria_Calliphoridae" ,
                            "Hemilucilia_semidiaphana_Calliphoridae" ,
                            "Lucilia_eximia_Calliphoridae" ,
                            
                            "Mesembrinella_batesi_Mesembrinellidae" ,
                            "Mesembrinella_bellardiana_Mesembrinellidae"  ,
                            "Mesembrinella_bicolor_Mesembrinellidae" ,   
                            "Mesembrinella_quadrilineata_Mesembrinellidae" , 
                            "Mesembrinella_randa_Mesembrinellidae" ,
                            
                            "Null_model_Amazonia"               ,
                            
                            "Oxysarcodexia_intona_Sarcophagidae"  ,         
                            "Oxysarcodexia_thornax_Sarcophagidae"   ,
                            "Peckia_(Pattonella)_intermutans_Sarcophagidae" ,
                            "Peckia_(Peckia)_chrysostoma_Sarcophagidae"  ,
                            "Peckia_(Sarcodexia)_lambens_Sarcophagidae"   )


names(RF_TestingOutput)<-c("Chloroprocta_idioidea_Calliphoridae" ,
                           "Chrysomya_albiceps_Calliphoridae" ,
                           "Cochliomyia_macellaria_Calliphoridae" ,
                           "Hemilucilia_semidiaphana_Calliphoridae" ,
                           "Lucilia_eximia_Calliphoridae" ,
                           
                           "Mesembrinella_batesi_Mesembrinellidae" ,
                           "Mesembrinella_bellardiana_Mesembrinellidae"  ,
                           "Mesembrinella_bicolor_Mesembrinellidae" ,   
                           "Mesembrinella_quadrilineata_Mesembrinellidae" , 
                           "Mesembrinella_randa_Mesembrinellidae" ,
                           
                           "Null_model_Amazonia"               ,
                           
                           "Oxysarcodexia_intona_Sarcophagidae"  ,         
                           "Oxysarcodexia_thornax_Sarcophagidae"   ,
                           "Peckia_(Pattonella)_intermutans_Sarcophagidae" ,
                           "Peckia_(Peckia)_chrysostoma_Sarcophagidae"  ,
                           "Peckia_(Sarcodexia)_lambens_Sarcophagidae"   )

save(RF_TrainingOutput, RF_TestingOutput, file="RData_sp/TR_and_TS_ModelOutputs.RData")




# STEP 3 - EVALUATING VARIABLE IMPORTANCE AND MODEL PERFORMANCE METRICS ----

rm(list=ls()); gc()

# Load model output for training and test data:
load("RData_sp/TR_and_TS_ModelOutputs.RData")

# Combine into a single dataframe:
TestingOutput<-rbindlist(RF_TestingOutput)

# Results models
ModelPerformance <- lapply(RF_TrainingOutput, purrr::pluck, "RF_output_k")


# Perform model evaluation for the complete set of predicted values using the test data:
### 3.1 - PERFORMING MODEL EVALUATION AND SAVING ----

ModelEvaluation<-list()

TestingOutput$group<-as.factor(TestingOutput$group)
for(i in 1:nlevels(TestingOutput$group)){
  
  ModelEvaluation[[i]]<-ENMTML:::Eval_Jac_Sor_TMLA(
    p = TestingOutput[TestingOutput$group==levels(TestingOutput$group)[i] & TestingOutput$presence==1,]$Predicted,
    a = TestingOutput[TestingOutput$group==levels(TestingOutput$group)[i] & TestingOutput$presence==0,]$Predicted,
    thr = NULL)
  
  names(ModelEvaluation)[i]<-levels(TestingOutput$group)[i]
  
}


##Report the limitations of the models:
ModelEvaluation$Chloroprocta_idioidea_Calliphoridae$SorensenTHR
ModelEvaluation$Mesembrinella_batesi_Mesembrinellidae$SorensenTHR
ModelEvaluation$Null_model_Amazonia$SorensenTHR
ModelEvaluation$Oxysarcodexia_intona_Sarcophagidae$SorensenTHR

###report  sorensenof the models
max(ModelEvaluation$Chloroprocta_idioidea_Calliphoridae$Sorensen, na.rm=T)
max(ModelEvaluation$Mesembrinella_batesi_Mesembrinellidae$Sorensen, na.rm=T)
max(ModelEvaluation$Null_model_Amazonia$Sorensen, na.rm=T)
max(ModelEvaluation$Oxysarcodexia_intona_Sarcophagidae$Sorensen, na.rm=T)


#TestingOutput

######exporta matrix
as.matrix(ModelEvaluation)->ModelEvaluation
write.table (ModelEvaluation, file="Datasets_sp/ModelEvaluation.csv")
save(ModelEvaluation, file="RData_sp/ModelEvaluation.RData")


ModelEvaluation<-rbindlist(ModelEvaluation)
ModelEvaluation
ModelEvaluation<-ModelEvaluation[, .(.N,
                                     avg_Sorensen = mean(Sorensen, na.rm=T),
                                     min_Sorensen = min(Sorensen, na.rm=T),
                                     max_Sorensen = max(Sorensen, na.rm=T),
                                     sd_Sorensen = sd(Sorensen, na.rm=T),
                                     avg_Jaccard= mean(Jaccard, na.rm=T),
                                     min_Jaccard = min(Jaccard, na.rm=T),
                                     max_Jaccard = max(Jaccard, na.rm=T),
                                     sd_Jaccard = sd(Jaccard, na.rm=T)),
                                 by=.(presences,absences)] 
ModelEvaluation

write.table (ModelEvaluation, file="Datasets_sp/ModelEvaluation_medias_final_sp.csv")



### 3.2 - EVALUATING VARIABLE IMPORTANCE ----
# Calculate the average performance metrics of the model for each group in each iteration (tree) and fold (k)::
ModelPerformance<-rbindlist(ModelPerformance)
ModelPerformance
ModelPerformance<-ModelPerformance[, .(.N,
                                       avg_mse = mean(mse, na.rm=T),
                                       min_mse = min(mse, na.rm=T),
                                       max_mse = max(mse, na.rm=T),
                                       sd_mse = sd(mse, na.rm=T),
                                       avg_rsq = mean(rsq, na.rm=T),
                                       min_rsq = min(rsq, na.rm=T),
                                       max_rsq = max(rsq, na.rm=T),
                                       sd_rsq = sd(rsq, na.rm=T)),
                                   by=.(group, tree)] 

ModelPerformance<-ModelPerformance[ModelPerformance$tree==500,]
ModelPerformance


VarImportance <- lapply(RF_TrainingOutput, purrr::pluck, "importance_k")
names(VarImportance)<-ModelPerformance$group 

# Extract variable importance metrics by group:
VarImp_Output<-data.frame()
for(i in 1:length(VarImportance)){
  
  VarImportance[[i]]$Kfold<-c(rep(1,5), rep(2,5), rep(3,5), rep(4,5), rep(5,5))
  VarImportance[[i]]$Variables<-rep(row.names(RF_TrainingOutput[[1]]$importance), 5)
  
  for(k in 1:5){
    
    MyData<-VarImportance[[i]][VarImportance[[i]]$Kfold==k,]
    MyData$group<-ModelPerformance$group[i]
    VarImp_Output<-rbind(VarImp_Output, MyData)
  }
}

write.table (VarImportance, file="Datasets_sp/Importanciadasvariaveis_medias_sp.csv")
rm(VarImportance)

### 3.3 - PREPARING THE DATA FRAME FOR PLOTTING: ----

VarImp_Output<-as.data.table(VarImp_Output)
VarImportance<-VarImp_Output[, .(avg_INP = mean(IncNodePurity, na.rm=T),
                                 min_INP = min(IncNodePurity, na.rm=T),
                                 max_INP = max(IncNodePurity, na.rm=T),
                                 sd_INP = sd(IncNodePurity),
                                 avg_IMSE = mean(X.IncMSE, na.rm=T),
                                 min_IMSE = min(X.IncMSE, na.rm=T),
                                 max_IMSE = max(X.IncMSE, na.rm=T),
                                 sd_IMSE = sd(X.IncMSE, na.rm=T)),
                             by=.(group, Variables)]

# Output:
VarImportance


### 3.4 - VARIABLE IMPORTANCE FIGURE - BARS ----
names(VarImportance$group)
names(VarImportance$group)
# Relabel levels to improve visualization:
VarImportance$group <- factor(VarImportance$group,
                              levels = c("Chloroprocta_idioidea_Calliphoridae" ,
                                         "Chrysomya_albiceps_Calliphoridae" ,
                                         "Cochliomyia_macellaria_Calliphoridae" ,
                                         "Hemilucilia_semidiaphana_Calliphoridae" ,
                                         "Lucilia_eximia_Calliphoridae" ,
                                         
                                         "Mesembrinella_batesi_Mesembrinellidae" ,
                                         "Mesembrinella_bellardiana_Mesembrinellidae"  ,
                                         "Mesembrinella_bicolor_Mesembrinellidae" ,   
                                         "Mesembrinella_quadrilineata_Mesembrinellidae" , 
                                         "Mesembrinella_randa_Mesembrinellidae" ,
                                         
                                         "Null_model_Amazonia"               ,
                                         
                                         "Oxysarcodexia_intona_Sarcophagidae"  ,         
                                         "Oxysarcodexia_thornax_Sarcophagidae"   ,
                                         "Peckia_(Pattonella)_intermutans_Sarcophagidae" ,
                                         "Peckia_(Peckia)_chrysostoma_Sarcophagidae"  ,
                                         "Peckia_(Sarcodexia)_lambens_Sarcophagidae"   ),
                              
                              labels = c("Chloroprocta idioidea - Calliphoridae" ,
                                         "Chrysomya albiceps - Calliphoridae" ,
                                         "Cochliomyia macellaria - Calliphoridae" ,
                                         "Hemilucilia semidiaphana - Calliphoridae" ,
                                         "Lucilia eximia - Calliphoridae" ,
                                         
                                         "Mesembrinella batesi - Mesembrinellidae" ,
                                         "Mesembrinella bellardiana - Mesembrinellidae"  ,
                                         "Mesembrinella bicolor - Mesembrinellidae" ,   
                                         "Mesembrinella quadrilineata - Mesembrinellidae" , 
                                         "Mesembrinella randa - Mesembrinellidae" ,
                                         
                                         "Null model - Amazon"               ,
                                         
                                         "Oxysarcodexia intona - Sarcophagidae"  ,         
                                         "Oxysarcodexia thornax - Sarcophagidae"   ,
                                         "Peckia (Pattonella) intermutans - Sarcophagidae" ,
                                         "Peckia (Peckia) chrysostoma - Sarcophagidae"  ,
                                         "Peckia (Sarcodexia) lambens - Sarcophagidae"   ))

VarImportance$Variables1 <- factor(VarImportance$Variables, 
                                   levels=rev(c("Degradation", "DistanceResearch","DrySeasonLength", "LandOwnership", "TravelTime" )),
                                   labels=rev(c("Degradation", "Distance to\nResearch centers", "Dry season\nLength", "Land\nOwnership", "Travel Time")))


# Plot the normalized IncNodePurity for each variable and group:


cores_especies <- c(
  # Calliphoridae (base: #8dd3c7)
  "#a7e0d4", "#8dd3c7", "#6bc6b9", "#4fbab0", "#3ca99e",
  
  # Mesembrinellidae (base: #bebada)
  "#d9d6e0", "#bebada", "#a59fc5", "#938bb3", "#8177a0",
  
  # Modelo nulo
  "black",
  
  # Sarcophagidae (base: #fdc086)
  "#ffddb2", "#fdc086", "#fdaa5a", "#f48825", "#d86e0b"
)


VarImportance$Family <- case_when(
  grepl("Calliphoridae", VarImportance$group) ~ "Calliphoridae",
  grepl("Mesembrinellidae", VarImportance$group) ~ "Mesembrinellidae",
  grepl("Sarcophagidae", VarImportance$group) ~ "Sarcophagidae",
  grepl("Null", VarImportance$group) ~ "Null"
)

MyPlot1<-
  ggplot(data=VarImportance, aes(x=group, y=avg_IMSE, shape=group)) +
  
  # Add geopoints:
  geom_point(shape=16, aes(colour=group, shape=group, fill=group), size=2) + 
  
  # Add error bars representing min and max values of variable importance:
  geom_errorbar(aes(ymin=min_IMSE, ymax=max_IMSE, col=group), width=0.5, cex=1) +
  
  # Set the colour ramp:
  scale_color_manual(values = c(cores_especies), 
                     labels = c("Chloroprocta idioidea - Calliphoridae" ,
                                "Chrysomya albiceps - Calliphoridae" ,
                                "Cochliomyia macellaria - Calliphoridae" ,
                                "Hemilucilia semidiaphana - Calliphoridae" ,
                                "Lucilia eximia - Calliphoridae" ,
                                
                                "Mesembrinella batesi - Mesembrinellidae" ,
                                "Mesembrinella bellardiana - Mesembrinellidae"  ,
                                "Mesembrinella bicolor - Mesembrinellidae" ,   
                                "Mesembrinella quadrilineata - Mesembrinellidae" , 
                                "Mesembrinella randa - Mesembrinellidae" ,
                                
                                "Null model - Amazon"               ,
                                
                                "Oxysarcodexia intona - Sarcophagidae"  ,         
                                "Oxysarcodexia thornax - Sarcophagidae"   ,
                                "Peckia (Pattonella) intermutans - Sarcophagidae" ,
                                "Peckia (Peckia) chrysostoma - Sarcophagidae"  ,
                                "Peckia (Sarcodexia) lambens - Sarcophagidae"   )) +
  
  # Split geom_points according to levels of 'Variables':
  facet_wrap(~Variables1, strip.position="left", nrow=5, scales = "free_y") +
  
  # Axes titles:
  xlab("") + ylab('% Increase in MSE (Normalized)') +
  
  # Other aesthetics:
  coord_flip() +  # flip coordinates (puts labels on y axis)
  #theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text.y = element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=12, colour="black", face="bold"),
        axis.title.x=element_text(size=12, colour="black", face="bold"),
        legend.position="top",
        #plot.background=element_rect(fill="transparent", colour=NA),
        strip.placement = "outside")

MyPlot1

ggsave("figure_defull_sp/VariableImportance.png", plot=MyPlot1, width=13, height=9)
#ggsave("Figures/VariableImportance2.pdf", plot=MyPlot1, width=8, height=5, units="in", bg="transparent")

# Export:
data.table::fwrite(VarImportance, file="Datasets_sp/VarImportance_sp.csv")



# STEP 4 - INSPECTING THE BEHAVIOR OF RESEARCH PROBABILITY AMONG PREDICTOR VALUES (PARTIAL PLOTS) ----


rm(list=ls()); gc()

# Load training and test data for the model:
load("RData_sp/TR_and_TS_Data.RData")

# Load the model output for training and test data:
load("RData_sp/TR_and_TS_ModelOutputs.RData")

# Load environmental variables:
envT <- raster::stack(list.files("Predictors_sp/", full.names=T))

envT[[4]]<-as.factor(envT[[4]]) # convert 'LandOwnership' to factor
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) # TravelTime - reclassifica os valores raster e define o valor 65535 como NA


# organism 'i':
Partial_Plots<-list()
for(i in 1:length(TestingData)){ 
  
  # Extract environmental data for all presence/pseudo-absence records:  
  All_envData<-lapply(TestingData[[i]], function(x)
    raster::extract(envT, x[,1:2], na.rm=T))
  All_envData<-lapply(All_envData, as.data.frame)
  
  #Combine all data folds used for training:
  MyData<-cbind(rbindlist(TestingData[[i]]), rbindlist(All_envData))
  MyData<-MyData[complete.cases(MyData),] # remove NA's (if any)
  MyData$LandOwnership<-as.factor(MyData$LandOwnership)
  MyData$LandOwnership<-factor(MyData$LandOwnership, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  
  # randomForest for organism 'i':
  SelectedRF_Model<-RF_TrainingOutput[[i]]
  
  # list for models results :
  myPlots<-list()
  PredictorNames<-c("Degradation", "DistanceResearch","DrySeasonLength", "LandOwnership", "TravelTime")
  # OBS: PredictorNames[4] refere-se à variável categórica 'LandOwnership'
  
  # Partial plots for all variables:
  for(v in 1:5){ # v number of variables
    
    # Gráficos para preditores contínuos:
    if(v!=4){
      
      myPlots[[v]]<-partialPlot(x = SelectedRF_Model, # randomForest
                                pred.data = MyData, # training data for the respective randomForest
                                x.var = PredictorNames[v], # variables
                                which.class = NULL, 
                                plot = TRUE, 
                                add = FALSE,
                                rug = TRUE, 
                                ylab="Research Presence",
                                xlab=PredictorNames[v]
      )
      
      myPlots[[v]]$Group<-names(TestingData)[i]
      myPlots[[v]]$Predictor<-PredictorNames[v]
      
    }
    
    # Plots for categorical predictors:
    if(v==4){
      
      myPlots[[v]]<-partialPlot(x = SelectedRF_Model, # object of class randomForest
                                pred.data = MyData, # training data for the respective randomForest
                                x.var = PredictorNames[v], # variable name to be examined
                                which.class = NULL, # for classification data, define which class
                                plot = TRUE, 
                                add = FALSE,
                                rug = TRUE, 
                                ylab="Research Presence",
                                xlab=PredictorNames[v]
      )
      
      myPlots[[v]]$Group<-names(TestingData)[i]
      myPlots[[v]]$Predictor<-PredictorNames[v]
    }
    
  }
  
  Partial_Plots[[i]]<-myPlots
  
}

# Partial_Plots is a set of 16 lists (groups), with 5 results (predictors), one for each group:
# Partial_Plots[[1]] contains five lists, one for each predictor variable
# Each list has the x value and y values of the partial graph for the given variable
Partial_Plots[[1]]



# Export to disk:
save(Partial_Plots, file="RData_sp/PartialPlotsData.RData")



load("RData_sp/PartialPlotsData.RData")



#data_frame <- data.frame(
#  Coluna1 = unlist(Partial_Plots[[1]]),
#  Coluna2 = unlist(Partial_Plots[[2]]),
#  Coluna3 = unlist(Partial_Plots[[3]]),
#  Coluna4 = unlist(Partial_Plots[[4]]),
#  Coluna5 = unlist(Partial_Plots[[5]]),
#  Coluna6 = unlist(Partial_Plots[[6]]),
#  Coluna7 = unlist(Partial_Plots[[7]]),
#  Coluna8 = unlist(Partial_Plots[[8]]),
#  Coluna9 = unlist(Partial_Plots[[9]]),
#  Coluna10 = unlist(Partial_Plots[[10]]),
#  Coluna11 = unlist(Partial_Plots[[11]]),
#  Coluna12 = unlist(Partial_Plots[[12]]),
#  Coluna13 = unlist(Partial_Plots[[13]]),
#  Coluna14 = unlist(Partial_Plots[[14]]),
#  Coluna15 = unlist(Partial_Plots[[15]]),
#  Coluna16 = unlist(Partial_Plots[[16]])
#)

#write.csv(data_frame, "Datasets_sp/partial_plots.csv", row.names = FALSE)



# Função que preenche com NA até o comprimento máximo
pad_with_na <- function(x, max_length) {
  length(x) <- max_length
  return(x)
}

# Define o maior comprimento entre as colunas
max_length <- max(sapply(Partial_Plots[1:16], function(x) length(unlist(x))))

# Cria o data.frame preenchendo com NA onde faltar
data_frame <- data.frame(
  Coluna1 = pad_with_na(unlist(Partial_Plots[[1]]), max_length),
  Coluna2 = pad_with_na(unlist(Partial_Plots[[2]]), max_length),
  Coluna3 = pad_with_na(unlist(Partial_Plots[[3]]), max_length),
  Coluna4 = pad_with_na(unlist(Partial_Plots[[4]]), max_length),
  Coluna5 = pad_with_na(unlist(Partial_Plots[[5]]), max_length),
  Coluna6 = pad_with_na(unlist(Partial_Plots[[6]]), max_length),
  Coluna7 = pad_with_na(unlist(Partial_Plots[[7]]), max_length),
  Coluna8 = pad_with_na(unlist(Partial_Plots[[8]]), max_length),
  Coluna9 = pad_with_na(unlist(Partial_Plots[[9]]), max_length),
  Coluna10 = pad_with_na(unlist(Partial_Plots[[10]]), max_length),
  Coluna11 = pad_with_na(unlist(Partial_Plots[[11]]), max_length),
  Coluna12 = pad_with_na(unlist(Partial_Plots[[12]]), max_length),
  Coluna13 = pad_with_na(unlist(Partial_Plots[[13]]), max_length),
  Coluna14 = pad_with_na(unlist(Partial_Plots[[14]]), max_length),
  Coluna15 = pad_with_na(unlist(Partial_Plots[[15]]), max_length),
  Coluna16 = pad_with_na(unlist(Partial_Plots[[16]]), max_length)
)

write.csv(data_frame, "Datasets_sp/partial_plots.csv", row.names = FALSE)



# STEP 5 - PREDICT GEOGRAPHICAL PATTERNS OF KNOWLEDGE PROBABILITY----
rm(list=ls()); gc()

# Load model output for training and testing data:
load("RData_sp/TR_and_TS_ModelOutputs.RData")

# Load raster data:
envT <- raster::stack(list.files("Predictors_sp/", full.names=T))
envT[[4]]<-as.factor(envT[[4]]) # convert 'LandOwnership' to factor
envT[[5]]<-raster::reclassify(envT[[5]], cbind(65535, +Inf, NA), right=TRUE) # reclassify TravelTime

# dataframe:
envT_df<-rasterToPoints(x=envT, spatial=FALSE)
envT_df[is.nan(envT_df)]<-NA # convert NaN values to NA
envT_df<-as.data.table(envT_df) # convert to data.table
summary(envT_df)
envT_df<-envT_df[!is.na(envT_df$Degradation),] # remove NAs
envT_df<-envT_df[!is.na(envT_df$DistanceResearch),] # remove NAs
envT_df<-envT_df[!is.na(envT_df$DrySeasonLength),] # remove NAs
envT_df<-envT_df[!is.na(envT_df$TravelTime ),] # remove NAs
envT_df$LandOwnership<-as.factor(envT_df$LandOwnership) # set 'LandTenure' as factor
summary(envT_df)

### 5.1 - PREDICTING THE PROBABILITY OF KNOWLEDGE FOR EACH FAMILY----
#This sets R to use several cores of your computer (minus 1) to speed up the process
cl<-makePSOCKcluster(detectCores()-1, type="SOCK")
registerDoParallel(cl)
getDoParWorkers()

#Using the 'predict()' function to predict the probability of presence in each cell of the dataframe 'envT_df':
MyPredictions<-foreach(i = 1:length(RF_TrainingOutput), 
                       .export = 'c', 
                       .packages = c("raster", "sp", "randomForest")) %dopar% {
                         predict(object = RF_TrainingOutput[[i]], newdata = envT_df,  type = "response")
                       }
parallel::stopCluster(cl)

group_names<-unique(rbindlist(RF_TestingOutput)$group) #Obtaining the names of the groups/families

#'MyPredictions': list with the predicted probabilities for each cell of the variables and for each family
#'group_names': names of these groups, useful for naming the rasters or maps later.


species_names <- c(
  "Chloroprocta_idioidea_Calliphoridae",
  "Chrysomya_albiceps_Calliphoridae",
  "Cochliomyia_macellaria_Calliphoridae",
  "Hemilucilia_semidiaphana_Calliphoridae",
  "Lucilia_eximia_Calliphoridae",
  
  "Mesembrinella_batesi_Mesembrinellidae",
  "Mesembrinella_bellardiana_Mesembrinellidae",
  "Mesembrinella_bicolor_Mesembrinellidae",
  "Mesembrinella_quadrilineata_Mesembrinellidae",
  "Mesembrinella_randa_Mesembrinellidae",
  
  "Null_model_Amazonia",
  
  "Oxysarcodexia_intona_Sarcophagidae",
  "Oxysarcodexia_thornax_Sarcophagidae",
  "Peckia_(Pattonella)_intermutans_Sarcophagidae",
  "Peckia_(Peckia)_chrysostoma_Sarcophagidae",
  "Peckia_(Sarcodexia)_lambens_Sarcophagidae"
)


# convert SpatialPixelsDataFrame:
#MyPixelDF<-(data.frame(envT_df[,1:2], 
#                       Calliphoridae_flies=MyPredictions[[1]],
#                       Mesembrinellidae_flies=MyPredictions[[2]],
#                       Model_nulo_nulo=MyPredictions[[3]],
#                       Sarcophagidae_flies=MyPredictions[[4]]))


stopifnot(length(MyPredictions) == length(species_names))

# Suponha que envT_df tem duas colunas: longitude e latitude
MyPixelDF <- data.frame(envT_df[,1:2])  # inicializa só com coords


# Adiciona as colunas com os nomes das espécies
for(i in seq_along(species_names)) {
  MyPixelDF[[species_names[i]]] <- MyPredictions[[i]]
}


coordinates(MyPixelDF)<- ~ x + y # convert to spatial points
gridded(MyPixelDF)<-TRUE  # griddify your points

plot(MyPixelDF)

# Convert rasters:
MyRasters<-list()
for(i in 1:ncol(MyPixelDF)){
  MyRasters[[i]] <- raster(MyPixelDF, layer=i) 
}

names(MyRasters)<-c("Chloroprocta_idioidea_Calliphoridae",
                    "Chrysomya_albiceps_Calliphoridae",
                    "Cochliomyia_macellaria_Calliphoridae",
                    "Hemilucilia_semidiaphana_Calliphoridae",
                    "Lucilia_eximia_Calliphoridae",
                    
                    "Mesembrinella_batesi_Mesembrinellidae",
                    "Mesembrinella_bellardiana_Mesembrinellidae",
                    "Mesembrinella_bicolor_Mesembrinellidae",
                    "Mesembrinella_quadrilineata_Mesembrinellidae",
                    "Mesembrinella_randa_Mesembrinellidae",
                    
                    "Null_model_Amazonia",
                    
                    "Oxysarcodexia_intona_Sarcophagidae",
                    "Oxysarcodexia_thornax_Sarcophagidae",
                    "Peckia_(Pattonella)_intermutans_Sarcophagidae",
                    "Peckia_(Peckia)_chrysostoma_Sarcophagidae",
                    "Peckia_(Sarcodexia)_lambens_Sarcophagidae")

plot(MyRasters)

par(mfrow = c(3, 4)) # 4 linhas, 5 colunas

plot(MyRasters$Chloroprocta_idioidea_Calliphoridae)
plot(MyRasters$Chrysomya_albiceps_Calliphoridae)
plot(MyRasters$Cochliomyia_macellaria_Calliphoridae)
plot(MyRasters$Hemilucilia_semidiaphana_Calliphoridae)
plot(MyRasters$Lucilia_eximia_Calliphoridae)

plot(MyRasters$Mesembrinella_batesi_Mesembrinellidae)
plot(MyRasters$Mesembrinella_bellardiana_Mesembrinellidae)
plot(MyRasters$Mesembrinella_bicolor_Mesembrinellidae)
plot(MyRasters$Mesembrinella_quadrilineata_Mesembrinellidae)
plot(MyRasters$Mesembrinella_randa_Mesembrinellidae)

plot(MyRasters$Null_model_Amazonia)

plot(MyRasters$Oxysarcodexia_intona_Sarcophagidae)
plot(MyRasters$Oxysarcodexia_thornax_Sarcophagidae)
plot(MyRasters$Peckia_(Pattonella)_intermutans_Sarcophagidae)
plot(MyRasters$Peckia_(Peckia)_chrysostoma_Sarcophagidae)
plot(MyRasters$Peckia_(Sarcodexia)_lambens_Sarcophagidae)



par(mfrow = c(1, 1)) # reseta para o padrão (opcional)

### 5.2 - CREATING FORECAST RASTER FOR EACH FAMILY ----

## empilhando os rasters
RasterStack<-raster::stack(MyRasters[[1]])
for(i in 2:length(MyRasters)){RasterStack<-raster::stack(RasterStack, MyRasters[[i]])} # faster than using do.call
names(RasterStack)<-names(MyRasters)
raster::crs(RasterStack) <- raster::crs(envT)

plot(RasterStack)


#Remove essas quatro camadas específicas do RasterStack, criando um novo objeto chamado FliesRasters
FliesRasters<-dropLayer(RasterStack, c("Chloroprocta_idioidea_Calliphoridae",
                                       "Chrysomya_albiceps_Calliphoridae",
                                       "Cochliomyia_macellaria_Calliphoridae",
                                       "Hemilucilia_semidiaphana_Calliphoridae",
                                       "Lucilia_eximia_Calliphoridae",
                                       
                                       "Mesembrinella_batesi_Mesembrinellidae",
                                       "Mesembrinella_bellardiana_Mesembrinellidae",
                                       "Mesembrinella_bicolor_Mesembrinellidae",
                                       "Mesembrinella_quadrilineata_Mesembrinellidae",
                                       "Mesembrinella_randa_Mesembrinellidae",
                                       
                                       "Null_model_Amazonia",
                                       
                                       "Oxysarcodexia_intona_Sarcophagidae",
                                       "Oxysarcodexia_thornax_Sarcophagidae",
                                       "Peckia_(Pattonella)_intermutans_Sarcophagidae",
                                       "Peckia_(Peckia)_chrysostoma_Sarcophagidae",
                                       "Peckia_(Sarcodexia)_lambens_Sarcophagidae"))


plot(FliesRasters)

#Máscara (remoção das áreas fora do habitat)
# raster upland:
UlpandMasck<-raster::raster("rasters/upland.tif")
plot(UlpandMasck)

# use mask layers to remove predictions outside the habitat:
ulpandflies <- raster::mask(RasterStack, UlpandMasck)
#Resultado: um novo raster (ulpandflies) com apenas os pixels dentro da área válida (como um "recorte").

plot(ulpandflies)

# Export the knowledge probability rasters:
dir.create("Projections_sp/", showWarnings = F)
for(i in 1:nlayers(ulpandflies)){terra::writeRaster(terra::rast(ulpandflies[[i]]), filename=paste0("Projections_sp/", names(ulpandflies)[i], ".tif"), overwrite=TRUE)}



# STEP 6 - GET AVERAGE PROBABILITY FROM RESEARCH AND ACROSS ALL ORGANISMS ----

rm(list=ls()); gc()

#listar arquivos
RasterFiles<-list.files(path=paste0("Projections_sp/"), pattern='.tif', full.names=T)
Filenames<-gsub("Projections_sp/", "", RasterFiles)
Filenames<-gsub("\\.tif", "", Filenames)


## para as moscas
# # Separate the raster layers available for each family:
#seleciona os rasters só com flies no nome
FliesFiles<-RasterFiles[RasterFiles %like% "ae"]


# raster:
#Carrega os rasters de moscas
FliesLayers<-lapply(FliesFiles, raster::raster)


# Convert to raster stack and rename layers:
FliesStack<-raster::stack(FliesLayers)
names(FliesStack)<-Filenames[Filenames %like% "ae"]


# Get the average probability of knowledge for each family:
#Calcula a média da probabilidade entre os grupos de moscas
ResearchProb_Flies<-raster::calc(stack(FliesStack), fun=mean, na.rm=T) # it is important to set na.rm=T
names(ResearchProb_Flies)<-"Flies"
plot(ResearchProb_Flies)


# Export the knowledge probability rasters:
dir.create("AvgOutputs_sp/", showWarnings = F)
terra::writeRaster(ResearchProb_Flies, filename=paste0("AvgOutputs_sp/Mean_FLIES_sp.tif"), overwrite=TRUE)



# Obtain the average probability of knowledge for the combination of families:
RasterLayers<-lapply(FliesFiles, raster::raster)
names(FliesStack)<-Filenames[Filenames %like% "ae"]

# Get the average search probability for upland
PixelDataFrame<-list()
for(i in 1:length(RasterLayers)){
  PixelDataFrame[[i]]<-raster::as.data.frame(RasterLayers[[i]], xy=TRUE, na.rm=T)
  PixelDataFrame[[i]]$HabitatOrganism<-names(RasterLayers)[i]
  names(PixelDataFrame[[i]])[3]<-"ResearchProb"
}

# data.table:
PixelDataFrame<-rbindlist(PixelDataFrame)

# Get the average values across all pixels:
PixelDataFrame$x<-round(PixelDataFrame$x, digits = 3)
PixelDataFrame$y<-round(PixelDataFrame$y, digits = 3)
PixelDataFrame$CoordID<-paste0(PixelDataFrame$x, "_", PixelDataFrame$y)
PixelDataFrame<-unique(PixelDataFrame)

# Calculate the average probability of searching among organisms:
AvgOutput<-PixelDataFrame[, .(.N,
                              Avg_x=mean(x, na.rm=T),
                              Avg_y=mean(y, na.rm=T),
                              AvgResearchProb=mean(ResearchProb, na.rm=T)),
                          by=.(CoordID)] 
names(AvgOutput)<-c("CoordID", "N", "x", "y", "ResearchProb")

AvgOutput$HabitatOrganism<-"All"



# Build a raster for all groups:
grid_cells<-data.frame(ResearchProb=AvgOutput$ResearchProb, AvgOutput[,3:4])
grid_cells$x<-round(grid_cells$x, digits = 2)
grid_cells$y<-round(grid_cells$y, digits = 2)
grid_cells<-unique(grid_cells)
sp::coordinates(grid_cells)<- ~ x + y # convert to spatial points
sp::gridded(grid_cells)<-TRUE  # griddify your set of points

plot(grid_cells)

# raster:
grid_cells <- raster(grid_cells)
raster::crs(grid_cells) <- crs(RasterLayers[[1]])
plot(grid_cells)

# Export raster :
terra::writeRaster(terra::rast(grid_cells), filename=paste0("Outputs_sp/AvgResearchProb_Flies_sp.tif"), overwrite=TRUE)



### para o modelo nulo

rm(list=ls()); gc()

#listar arquivos
RasterFiles<-list.files(path=paste0("Projections_sp/"), pattern='.tif', full.names=T)
Filenames<-gsub("Projections_sp/", "", RasterFiles)
Filenames<-gsub("\\.tif", "", Filenames)


NullFile <- RasterFiles[grepl("Amazonia", RasterFiles, ignore.case = TRUE)]
NullLayer <- raster::raster(NullFile)
names(NullLayer) <- gsub("Projections_sp/|\\.tif", "", NullFile)


terra::writeRaster(NullLayer, filename = "AvgOutputs_sp/Model_Nulo.tif", overwrite = TRUE)

# Converter para data.frame com coordenadas
NullDF <- raster::as.data.frame(NullLayer, xy = TRUE, na.rm = TRUE)
names(NullDF)[3] <- "ResearchProb"
NullDF$HabitatOrganism <- names(NullLayer)
NullDF$x <- round(NullDF$x, 3)
NullDF$y <- round(NullDF$y, 3)
NullDF$CoordID <- paste0(NullDF$x, "_", NullDF$y)

# Média (nesse caso será só ele mesmo, mas mantém a estrutura igual)
NullDF <- data.table::as.data.table(NullDF)
AvgOutputNull <- NullDF[, .(.N,
                            Avg_x = mean(x, na.rm = TRUE),
                            Avg_y = mean(y, na.rm = TRUE),
                            AvgResearchProb = mean(ResearchProb, na.rm = TRUE)),
                        by = .(CoordID)]
names(AvgOutputNull) <- c("CoordID", "N", "x", "y", "ResearchProb")
AvgOutputNull$HabitatOrganism <- "ModelNulo"


# Build a raster for all groups:
grid_null <- data.frame(ResearchProb = AvgOutputNull$ResearchProb, AvgOutputNull[, 3:4])
grid_null$x <- round(grid_null$x, digits = 2)
grid_null$y <- round(grid_null$y, digits = 2)
grid_null <- unique(grid_null)
sp::coordinates(grid_null) <- ~ x + y
sp::gridded(grid_null) <- TRUE
plot(grid_null)


# raster:
grid_null <- raster(grid_null)
raster::crs(grid_null) <- raster::crs(NullLayer)
plot(grid_null)

terra::writeRaster(terra::rast(grid_null), filename = "Outputs_sp/AvgResearchProb_ModelNulo.tif", overwrite = TRUE)



### espécies cada

dir.create("Outputs_sp", showWarnings = FALSE)

# Lista de arquivos .tif
RasterFiles <- list.files(path = "Projections_sp/", pattern = "\\.tif$", full.names = TRUE)
Filenames <- gsub("Projections_sp/", "", RasterFiles)
Filenames <- gsub("\\.tif", "", Filenames)

# Lista para guardar os data.frames
SpeciesDFs <- list()

# Loop para processar cada raster
for (i in seq_along(RasterFiles)) {
  
  r <- raster::raster(RasterFiles[i])  # <- aqui você realmente carrega o raster
  sp_name <- Filenames[i]
  
  df <- raster::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "ResearchProb"
  df$HabitatOrganism <- sp_name
  df$x <- round(df$x, 3)
  df$y <- round(df$y, 3)
  df$CoordID <- paste0(df$x, "_", df$y)
  
  df <- data.table::as.data.table(df)
  
  df_avg <- df[, .(.N,
                   Avg_x = mean(x, na.rm = TRUE),
                   Avg_y = mean(y, na.rm = TRUE),
                   AvgResearchProb = mean(ResearchProb, na.rm = TRUE)),
               by = .(CoordID)]
  names(df_avg) <- c("CoordID", "N", "x", "y", "ResearchProb")
  df_avg$HabitatOrganism <- sp_name
  
  SpeciesDFs[[sp_name]] <- df_avg
}


# Criar diretório de saída se não existir
dir.create("Outputs_sp/SpeciesMeanRasters", showWarnings = FALSE)

# Converter para raster e salvar
for (sp_name in names(SpeciesDFs)) {
  
  df <- SpeciesDFs[[sp_name]]
  
  grid <- data.frame(ResearchProb = df$ResearchProb, x = df$x, y = df$y)
  grid$x <- round(grid$x, 2)
  grid$y <- round(grid$y, 2)
  grid <- unique(grid)
  
  sp::coordinates(grid) <- ~ x + y
  sp::gridded(grid) <- TRUE
  
  r <- raster::raster(grid)
  raster::crs(r) <- raster::crs(raster::raster(RasterFiles[1]))  # pega o CRS do primeiro raster
  
  terra::writeRaster(terra::rast(r),
                     filename = paste0("Outputs_sp/SpeciesMeanRasters/", sp_name, "_mean.tif"),
                     overwrite = TRUE)
}

