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
# 1- Random sampling and selection of points, and selection of coordinates to generate the null model
############################################################################################################################



### DIRECTORY AND PACKAGES ----
# First, clear the workspace and set your working directory::
rm(list=ls()); gc()

# Install and load the R packages required to run the analysis:
needed_packages <- c("data.table", "devtools", "dismo", "doParallel", "dplyr", 
                     "foreach", "ggspatial", "ggplot2", "glmnet", "kernlab",
                     "magrittr", "mapview",  "maxnet", "mgcv", "parallel", "plyr", "purrr",
                     "randomForest", "raster",  "rgdal", "sf", "sp", "splitTools", "stringr", "stats", "terra")

#new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(needed_packages, require, character.only = TRUE) 
removeTmpFiles(0.00001) # remove temporary files
rm(needed_packages, new.packages)


# STEP 1 - GET THE TRAINING AND TEST DATA FOLDS ----
# Read txt file with occurrence data:
occ <- data.table::fread("Datasets/occ.csv", h = T, stringsAsFactors = F)

sum(occ$familia=="Sarcophagidae") #2793
sum(occ$familia=="Mesembrinellidae") #1226
sum(occ$familia=="Calliphoridae") #4225
#Como o que tem mais coordenadas e Calliphoridae com 4225 coordenadas, vou usar esse valor como quantidade de pontos sorteados.



r2 <- raster::stack("upland.tif")
envT_df <- data.frame(rasterToPoints(x=r2, spatial=FALSE)) #Pronto, temos os pontos ou coordenadas para sortear, porem precisamos tirar o que e zero (fora da amazonia ou rio).
envT_df2 <- envT_df[envT_df$layer==!0,]
dim(envT_df2) #Tem 2410754 linhas com coordenadas. vamos sortear 4223 delas

num_sort <- sample(1:2410754, 4225) #É preciso salvar esses números para não ficar sorteando sempre!
#write.table(num_sort,"num_sort.txt")
num_sort <- read.table("num_sort.txt", h=T, stringsAsFactors = T)

envT_df_sort <- envT_df2[num_sort$x,] #Coordenadas sorteadas ao acaso. Tambem salvar isso.
#write.table(envT_df_sort,"coordenates_sort.txt")
envT_df_sort <- read.table("coordenates_sort.txt", h=T, stringsAsFactors = T)


r_df <- as.data.frame(r2, xy = TRUE)
library(terra)
ggplot(r_df) +
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data=envT_df_sort, aes(x=x, y=y), color="darkred", alpha=.6) +
  scale_fill_viridis_c(direction = -1, alpha = .4) +
  coord_equal() +
  theme_minimal()

ggsave(filename = "Random_points.pdf", width = 10, height = 7)
