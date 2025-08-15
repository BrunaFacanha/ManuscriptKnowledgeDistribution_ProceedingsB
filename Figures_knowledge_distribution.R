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


# First, clean workspace and set your working directory::
rm(list=ls()); gc()


### FIGURE 2 - RIDGELINE GRAPH OF RESEARCH PROBABILITY----

library(ggridges)

RasterList<-list()
RasterList[[1]]<-raster::raster("AvgOutputs/Mean_UP_flies.tif")
names(RasterList)<-c("Flies")

# For each raster, extract the pixel values as a data table:
PixelDataFrame<-list()
for(i in 1:length(RasterList)){
  PixelDataFrame[[i]]<-raster::as.data.frame(RasterList[[i]], xy=TRUE, na.rm=T)
  PixelDataFrame[[i]]$Habitat<-names(RasterList)[i]
  names(PixelDataFrame[[i]])[3]<-"ResearchProb"
}

# Bind all list in a single data.table:
PixelDataFrame<-rbindlist(PixelDataFrame)
PixelDataFrame$Habitat<-factor(PixelDataFrame$Habitat, levels=c("Flies"))

# Build the plots:
RidgelinePlot<-
  ggplot(PixelDataFrame, aes(x=ResearchProb, y=Habitat, fill=after_stat(x))) + 
  
  # Create the ridges:
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      show.legend=F, scale=2, rel_min_height=0.001) +
  
  # Set the colors:
  scale_fill_gradientn(colours=grDevices::hcl.colors(n = 100, palette="Spectral", rev=T), na.value='transparent') + 
  scale_y_discrete(limits=rev) +
  
  # Add labels for axes and define other aesthetics:
  labs(y="", x="Research probability") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, margin=ggplot2::margin(t=1, r=0, b=0, l=0)),
        axis.title.x = element_text(margin=ggplot2::margin(t=5, r=0, b=0, l=0)),
        axis.title.y = element_text(margin=ggplot2::margin(t=0, r=5, b=0, l=0)))

RidgelinePlot

# Export to disk:
#ggsave("Figures/Figure1_RidglineNoEdits_NOVO.png", plot=RidgelinePlot, width=4, height=4)
#ggsave("Figures/Figure1_RidglineNoEdits_NOVO.pdf", plot=RidgelinePlot, width=4, height=4, units="in", bg = "transparent")


# Note: the maps used to compose figure 2 were created outside of R, in the QGIS environment.


### FIGURE 2.2 - FAMILY RIDGLINES ----

rm(list=ls()); gc()

# Extract the habitat-organism combination from each filename:
setwd("C:/Users/Bruna/OneDrive/Doutorado (2021-2025)/capitulo1/artigo1/Metodo")
RasterFiles<-list.files(path=paste0("Projections/"), pattern='.tif', full.names=T)
Filenames<-gsub("Projections/", "", RasterFiles)
Filenames<-gsub("\\.tif", "", Filenames)

# Load raster for each habitat type:
RasterList<-lapply(RasterFiles, raster::raster)
names(RasterList)<-c("upland_Calliphoridae", "upland_Mesembrinellidae", "upland_Sarcophagidae")

# For each raster, extract the pixel values as a data table:
PixelDataFrame<-list()
for(i in 1:length(RasterList)){
  PixelDataFrame[[i]]<-raster::as.data.frame(RasterList[[i]], xy=TRUE, na.rm=T)
  PixelDataFrame[[i]]$HabitatOrganism<-names(RasterList)[i]
  names(PixelDataFrame[[i]])[4]<-"ResearchProb"
}

# Bind all lists into a single data.table:
PixelDataFrame<-rbindlist(PixelDataFrame)


# Get separate columns to indicate habitat and organism:
PixelDataFrame$Habitat<-gsub("([A-Za-z]+).*", "\\1", PixelDataFrame$HabitatOrganism)
PixelDataFrame$Organism<-gsub(".*_([A-Za-z]+)$", "\\1", PixelDataFrame$HabitatOrganism)

# Build the plots:
#RidgelinePlotAll<-
ggplot(PixelDataFrame, aes(x=upland_Calliphoridae, y=ResearchProb, fill=after_stat(x))) + 
  
  # Create the ridges:
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      show.legend=F, scale=2, rel_min_height=0.001) +
  
  # Set the colors:
  scale_fill_gradientn(colours=grDevices::hcl.colors(n = 100, palette="Spectral", rev=T), na.value='transparent') + 
  scale_y_discrete(limits=rev) +
  
  # Add rectangular boxes in the background:
  # geom_rect(data=PixelDataFrame, xmin=-0.04, xmax=0, ymin=0, ymax=11, fill="#6B9D59") +
  #geom_rect(xmin=-0.04, xmax=0, ymin=4.5, ymax=6.5, fill="#E093C3") +
  #geom_rect(xmin=-0.04, xmax=0, ymin=6.5, ymax=11, fill="#E47352") +
  
  # Add labels for axes and define other aesthetics:
  labs(y="", x="Research probability") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=12),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, margin=ggplot2::margin(t=1, r=0, b=0, l=0)),
        axis.title.x = element_text(margin=ggplot2::margin(t=5, r=0, b=0, l=0)),
        axis.title.y = element_text(margin=ggplot2::margin(t=0, r=5, b=0, l=0)))




RidgelinePlotAll

# Export to disk:
#ggsave("Figures/Figure2_RidglineNoEdits.png", plot=RidgelinePlotAll, width=5, height=5)
#ggsave("Figures/Figure2_RidglinePlot_NoEdits2.pdf", plot=RidgelinePlotAll, width=8, height=8, units="in", bg = "transparent")

# Note: organism silhouettes and vertical colored bands were adjusted outside R, using InkScape software.

### FIGURE 2.3 - RIDGELINE GRAPH OF RESEARCH PROBABILITY - NULL MODEL ----

library(ggridges)

RasterList<-list()
RasterList[[1]]<-raster::raster("AvgOutputs_sp/Model_Nulo.tif")
names(RasterList)<-c("Null model")

# For each raster, extract the pixel values as a data table:
PixelDataFrame<-list()
for(i in 1:length(RasterList)){
  PixelDataFrame[[i]]<-raster::as.data.frame(RasterList[[i]], xy=TRUE, na.rm=T)
  PixelDataFrame[[i]]$Habitat<-names(RasterList)[i]
  names(PixelDataFrame[[i]])[3]<-"ResearchProb"
}



# Bind all list in a single data.table:
PixelDataFrame<-rbindlist(PixelDataFrame)
PixelDataFrame$Habitat<-factor(PixelDataFrame$Habitat, levels=c("Null model"))

# Build the plots:
RidgelinePlot<-
  ggplot(PixelDataFrame, aes(x=ResearchProb, y=Habitat, fill=after_stat(x))) + 
  
  # Create the ridges:
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      show.legend=F, scale=2, rel_min_height=0.001) +
  
  # Set the colors:
  scale_fill_gradientn(colours=grDevices::hcl.colors(n = 100, palette="Spectral", rev=T), na.value='transparent') + 
  scale_y_discrete(limits=rev) +
  
  # Add labels for axes and define other aesthetics:
  labs(y="", x="Knowledge probability") +
  theme(text = element_text(family = "serif"))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=30),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=40),
        axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, margin=ggplot2::margin(t=1, r=0, b=0, l=0)),
        axis.title.x = element_text(margin=ggplot2::margin(t=5, r=0, b=0, l=0)),
        axis.title.y = element_text(margin=ggplot2::margin(t=0, r=5, b=0, l=0)))

RidgelinePlot

# Export to disk:
ggsave("Figures/Figure1_Ridgline_NULL-MODEL_sp.png", plot=RidgelinePlot, width=15, height=10)
#ggsave("Figures/Figure1_Ridgline_NULL-MODEL_sp.pdf", plot=RidgelinePlot, width=15, height=10, units="in", bg = "transparent")


### FIGURE SUPPLEMENTARY 1 - RIDGELINE GRAPH OF RESEARCH PROBABILITY----

library(ggridges)

RasterList<-list()
RasterList[[1]]<-raster::raster("AvgOutputs_sp/Mean_FLIES_sp.tif")
names(RasterList)<-c("Species of flies")

# For each raster, extract the pixel values as a data table:
PixelDataFrame<-list()
for(i in 1:length(RasterList)){
  PixelDataFrame[[i]]<-raster::as.data.frame(RasterList[[i]], xy=TRUE, na.rm=T)
  PixelDataFrame[[i]]$Habitat<-names(RasterList)[i]
  names(PixelDataFrame[[i]])[3]<-"ResearchProb"
}



# Bind all list in a single data.table:
PixelDataFrame<-rbindlist(PixelDataFrame)
PixelDataFrame$Habitat<-factor(PixelDataFrame$Habitat, levels=c("Species of flies"))

# Build the plots:
RidgelinePlot<-
  ggplot(PixelDataFrame, aes(x=ResearchProb, y=Habitat, fill=after_stat(x))) + 
  
  # Create the ridges:
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      show.legend=F, scale=2, rel_min_height=0.001) +
  
  # Set the colors:
  scale_fill_gradientn(colours=grDevices::hcl.colors(n = 100, palette="Spectral", rev=T), na.value='transparent') + 
  scale_y_discrete(limits=rev) +
  
  # Add labels for axes and define other aesthetics:
  labs(y="", x="Knowledge probability") +
  theme(text = element_text(family = "serif"))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=30),
        axis.line = element_line(colour="black"),
        axis.text = element_text(size=40),
        axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5, margin=ggplot2::margin(t=1, r=0, b=0, l=0)),
        axis.title.x = element_text(margin=ggplot2::margin(t=5, r=0, b=0, l=0)),
        axis.title.y = element_text(margin=ggplot2::margin(t=0, r=5, b=0, l=0)))

RidgelinePlot

# Export to disk:
ggsave("Figures/Figure1_Ridgline_FLIES_sp.png", plot=RidgelinePlot, width=15, height=10)
#ggsave("Figures/Figure1_Ridgline_FLIES_sp.pdf", plot=RidgelinePlot, width=15, height=10, units="in", bg = "transparent")

# Note: the maps used to compose figure 2 were created outside of R, in the QGIS environment.



### FIGURE SUPPLEMENTARY 1 - SPECIES RIDGLINES ----

rm(list=ls()); gc()

# Define diretório com rasters das espécies
setwd("C:/Users/blbar/OneDrive/Doutorado (2021-2025)/capitulo1/artigo1/submissao 2 - Proceedings/ressubmissao/data/Projections_sp")
# Lista arquivos .tif
RasterFiles <- list.files(pattern = '.tif$', full.names = TRUE)  # não precisa do path

# Load raster for each habitat type:
RasterList<-lapply(RasterFiles, raster::raster)
names(RasterList)<-c("Chloroprocta_idioidea_Calliphoridae",
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
                     "Peckia_(Sarcodexia)_lambens_Sarcophagidae") # Nomeia por espécie


# Extrai os dados de cada raster
PixelDataFrame <- list()
for(i in seq_along(RasterList)) {
  r <- RasterList[[i]]
  df <- raster::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  
  # Nome da coluna de probabilidade (detecta automaticamente)
  value_col <- setdiff(names(df), c("x", "y"))[1]
  names(df)[names(df) == value_col] <- "ResearchProb"
  
  df$Species <- names(RasterList)[i]
  PixelDataFrame[[i]] <- df
}


names(PixelDataFrame)<-c("Chloroprocta_idioidea_Calliphoridae",
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
                         "Peckia_(Sarcodexia)_lambens_Sarcophagidae") # Nomeia por espécie



# Junta todos os data.frames
PixelDataFrame <- data.table::rbindlist(PixelDataFrame)


library(ggplot2)
library(ggridges)


species_labels <- c(
  "Null_model_Amazonia" = "Null model",
  "Chloroprocta_idioidea_Calliphoridae" = "Chloroprocta idioidea",
  "Chrysomya_albiceps_Calliphoridae" = "Chrysomya albiceps",
  "Cochliomyia_macellaria_Calliphoridae" = "Cochliomyia macellaria",
  "Hemilucilia_semidiaphana_Calliphoridae" = "Hemilucilia semidiaphana",
  "Lucilia_eximia_Calliphoridae" = "Lucilia eximia",
  
  "Mesembrinella_batesi_Mesembrinellidae" = "Mesembrinella batesi",
  "Mesembrinella_bellardiana_Mesembrinellidae" = "Mesembrinella bellardiana",
  "Mesembrinella_bicolor_Mesembrinellidae" = "Mesembrinella bicolor",
  "Mesembrinella_quadrilineata_Mesembrinellidae" = "Mesembrinella quadrilineata",
  "Mesembrinella_randa_Mesembrinellidae" = "Mesembrinella randa",
  
  "Oxysarcodexia_intona_Sarcophagidae" = "Oxysarcodexia intona",
  "Oxysarcodexia_thornax_Sarcophagidae" = "Oxysarcodexia thornax",
  "Peckia_(Pattonella)_intermutans_Sarcophagidae" = "Peckia (Pattonella) intermutans",
  "Peckia_(Peckia)_chrysostoma_Sarcophagidae" = "Peckia (Pattonella) intermutans",
  "Peckia_(Sarcodexia)_lambens_Sarcophagidae" = "Peckia (Sarcodexia) lambens"
  # adicione o resto das espécies aqui
)

PixelDataFrame$SpeciesLabel <- species_labels[as.character(PixelDataFrame$Species)]

# Extraindo a família (última palavra depois do último "_")
PixelDataFrame$Family <- gsub(".*_([A-Za-z]+)$", "\\1", PixelDataFrame$Species)

# Corrige o null model (caso ele vire "Amazonia" ou algo assim)
PixelDataFrame$Family[PixelDataFrame$Species == "Null_model_Amazonia"] <- "Null"

family_colors <- c(
  "Calliphoridae" = "#6B9D59",
  "Mesembrinellidae" = "#E093C3",
  "Sarcophagidae" = "#E47352",
  "Null" = "black"
)

# Cria os rótulos em itálico
# Garante ordem correta com null model primeiro
PixelDataFrame$LabelFormatted <- ifelse(
  PixelDataFrame$SpeciesLabel == "Null model",
  "'Null model'",  # <-- coloco aspas simples para manter como texto plano no eixo
  paste0("italic('", PixelDataFrame$SpeciesLabel, "')")
)

ordered_labels <- c("Null model", sort(setdiff(unique(PixelDataFrame$SpeciesLabel), "Null model")))
PixelDataFrame$LabelFormatted <- factor(PixelDataFrame$LabelFormatted,
                                        levels = ifelse(ordered_labels == "Null model",
                                                        "'Null model'",
                                                        paste0("italic('", ordered_labels, "')")))


# Garantir ordem com o fator do eixo y
PixelDataFrame$LabelFormatted <- factor(PixelDataFrame$LabelFormatted,
                                        levels = rev(levels(PixelDataFrame$LabelFormatted)))

# Criar uma tabela com posições únicas de cada espécie no eixo Y
library(dplyr)

species_positions <- PixelDataFrame %>%
  distinct(LabelFormatted, Family) %>%
  mutate(y_pos = as.numeric(LabelFormatted)) %>%
  arrange(y_pos)

# Identificar os blocos contínuos por família
species_positions <- species_positions %>%
  group_by(Family) %>%
  summarise(ymin = min(y_pos) - 0.5,
            ymax = max(y_pos) + 0.5) %>%
  ungroup()

# Adicionar as cores desejadas por família
species_positions$fill_color <- recode(species_positions$Family,
                                       "Calliphoridae" = "#6B9D59",
                                       "Mesembrinellidae" = "#E093C3",
                                       "Sarcophagidae" = "#E47352",
                                       "Null" = "black"
)


RidgelinePlot_sp <- 
  ggplot(PixelDataFrame, aes(x = ResearchProb, y = LabelFormatted, fill = after_stat(x))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      show.legend = FALSE, scale = 2, rel_min_height = 0.001) +
  scale_fill_gradientn(colours = grDevices::hcl.colors(n = 100, palette = "Spectral", rev = TRUE),
                       na.value = 'transparent') +
  # Add rectangular boxes in the background:
  geom_rect(data = species_positions,
            aes(xmin = -0.04, xmax = -0.03, ymin = ymin, ymax = ymax, fill = NULL),
            inherit.aes = FALSE, fill = species_positions$fill_color) +
  scale_y_discrete(labels = function(x) parse(text = x)) +  # interpreta expression()
  labs(y = "", x = "Knowledge probability") +
  theme(text = element_text(family = "serif"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 30),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                   margin = margin(t = 1)),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)))


RidgelinePlot_sp

setwd("C:/Users/blbar/OneDrive/Doutorado (2021-2025)/capitulo1/artigo1/submissao 2 - Proceedings/ressubmissao/data")

# Export to disk:
ggsave("Figures/Figure1.2_Ridgline_species.png", plot=RidgelinePlot_sp, width=15, height=12)
#ggsave("Figures/Figure1.2_Ridgline_species.pdf", plot=RidgelinePlot_sp, width=15, height=12, units="in", bg = "transparent")

# Note: organism silhouettes and vertical colored bands were adjusted outside R, using InkScape software.



#Fig 3A====
pacman::p_load(openxlsx, ggplot2, ggpubr, scico, OPTS, 
               ggimage)

data1 <- data.table::fread("Datasets/VarImportance_figuras.csv", h = T, stringsAsFactors = F)


# Change the order of levels and labels in the variable group:
#data1$group = factor(data1$group, levels = c("Calliphoridae_flies", "Mesembrinellidae_flies", "Sarcophagidae_flies"),
                     #labels=c("Calliphoridae", "Mesembrinellidae", "Sarcophagidae"))

#groups = levels(as.factor(DM$group))
data1$group = as.factor(data1$group)
data1$Families = as.factor(data1$Families)
data1$avg_IMSE = as.numeric(data1$avg_IMSE)
data1$min_IMSE = as.numeric(data1$min_IMSE)
data1$max_IMSE = as.numeric(data1$max_IMSE)

# Associate each color with a group of organism:
myColors <- c("#6B9D59", #call
              "#E093C3", #mesemb
              "black",
              "#E47352") #sarco

names(myColors)<-c("Calliphoridae", "Mesembrinellidae", "Null model", "Sarcophagidae")
myFillColors<-scale_fill_manual(values = myColors , name="", na.value=NULL)
myColorColors<-scale_color_manual(values = myColors , name="", na.value=NULL)

names(data1)

nome_variaveis <- c(
  "1" = "Degradation",
  "2" = "Dry season\n length",
  "3" = "Land\n ownership",
  "4" = "Travel time",
  "5" = "Distance to\n research centers")


names <- c("Null model",
"Calliphoridae",
"Mesembrinellidae",
"Sarcophagidae")


data1$Families <- factor(data1$Families,
                         levels = c("Null model", "Calliphoridae", "Mesembrinellidae", "Sarcophagidae"))

line_null <- data1[data1$Families == "Null model",]

library(cowplot)
plot_importance <- ggplot(data1, aes(x = avg_IMSE, y = Families, color = Families)) +
  geom_errorbar(aes(xmin = min_IMSE, xmax = max_IMSE), width = 0.25, linewidth = .6, alpha = .5) +
    geom_point(aes(x = avg_IMSE), shape = 16, position=position_dodge(width = 0.52), size = 4, alpha=.9) +
  geom_vline(data = line_null, aes(xintercept = avg_IMSE), linetype = "dashed", linewidth = 0.3, color = "gray30", alpha = .5) +
  scale_y_discrete(labels = names) +
  scale_color_manual(values = myColors) +
  scale_size_continuous(range = c(2.5, 6)) +
  labs(x = "Increase in MSE", y = "", tag = "(a)") +
  xlim(min(data1$min_IMSE) * 1.1, max(data1$max_IMSE) * 1.15) +
  theme_cowplot(font_family = "serif", font_size = 15)+
  theme(legend.position = "top", 
        legend.spacing.x = unit(3, "cm"),
        panel.spacing = unit(1.3, "lines"),
        strip.background = element_rect(fill = NA,color = "black",linewidth = 0.5),
        strip.text = element_text(color = "black"),
        axis.text.y = element_blank())+
  facet_wrap(~ Variables_reoder, scales = "free_x", nrow = 1, strip.position = "top", labeller = labeller(Variables_reoder = nome_variaveis))

plot_importance




#Fig 3B-E====
pacman::p_load(openxlsx, ggplot2, ggpubr, scico, ggimage, rsvg)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(scico)
library(ggimage)
library(rsvg)



DEf = read.xlsx("Datasets/figure5_A-D.xlsx", sheet = 1) #degradation
DMf = read.xlsx("Datasets/figure5_A-D.xlsx", sheet = 2) #drylength
TTf = read.xlsx("Datasets/figure5_A-D.xlsx", sheet = 3) #traveltime
EDf = read.xlsx("Datasets/figure5_A-D.xlsx", sheet = 4) #research
colnames(DMf)

## change values to 0-1
range01 <- function(x){(x-min(x,na.rm=T))/(max(x, na.rm=T)-min(x,na.rm=T))}
DMf$DM = range01(DMf$dry_months)
DEf$DE = range01(DEf$degradation)
EDf$ED = range01(EDf$education)
TTf$TT = range01(TTf$travel_time)

#groups = levels(as.factor(DM$group))
DMf$group = factor(DMf$group, levels = c("Calliphoridae", "Mesembrinellidae","Null model", "Sarcophagidae"))
DEf$group = factor(DEf$group, levels = c("Calliphoridae", "Mesembrinellidae","Null model", "Sarcophagidae"))
EDf$group = factor(EDf$group, levels = c("Calliphoridae", "Mesembrinellidae", "Null model","Sarcophagidae"))
TTf$group = factor(TTf$group, levels = c("Calliphoridae", "Mesembrinellidae","Null model", "Sarcophagidae"))



DMp = 
  ggplot(DMf, aes(DM, Probability.of.Research))+
  geom_line(aes(color = group_env, group = group_env), cex = 1.2, alpha=.8, show.legend = F)+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_color_manual(values = myColors)+
  scale_fill_manual(values = myColors)+
  scale_y_continuous(limits=c(0.10,0.65), breaks=c(0.2, 0.4, 0.6))+
  xlab("Dry season length")+
  ylab("Probability of knowledge")+
  theme_cowplot(font_family = "serif", font_size = 15)
DMp
  

DEp = 
  ggplot(DEf, aes(DE, Probability.of.Research))+
  geom_line(aes(color = group_env, group = group_env), cex = 1.2, alpha=.8, show.legend = F)+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_color_manual(values = myColors)+
  scale_fill_manual(values = myColors)+
  scale_y_continuous(limits=c(0.10,0.65), breaks=c(0.2, 0.4, 0.6))+
  ylab("Probability of knowledge")+
  xlab("Degradation")+
  theme_cowplot(font_family = "serif", font_size = 15)
DEp


EDp = 
  ggplot(EDf, aes(ED, Probability.of.Research))+
  geom_line(aes(color = group_env, group = group_env), cex = 1.2, alpha=.8, show.legend = F)+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_color_manual(values = myColors)+
  scale_fill_manual(values = myColors)+
  scale_y_continuous(limits=c(0.10,0.65), breaks=c(0.2, 0.4, 0.6))+
  ylab("Probability of knowledge")+
  xlab("Distance to research centers")+
  theme_cowplot(font_family = "serif", font_size = 15)
  
EDp


TTp = 
  ggplot(TTf, aes(TT, Probability.of.Research))+
  geom_line(aes(color = group_env, group = group_env), cex=1.2, alpha=.8, show.legend = F)+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_color_manual(values = myColors)+
  scale_fill_manual(values = myColors)+
  scale_y_continuous(limits=c(0.10,0.65), breaks=c(0.2, 0.4, 0.6))+
  ylab("Probability of knowledge")+
  xlab("Travel time")+
  theme_cowplot(font_family = "serif", font_size = 15)
  
TTp

f2 <- plot_grid(TTp, EDp, DEp, DMp, 
          ncol = 2, align = "hv", labels = c("(b)","(c)","(d)","(e)"), label_fontfamily = "serif")

f2



#Fig 3F====
LTf = read.xlsx("Datasets/figure5_E.xlsx")
LTf$group = factor(LTf$group, levels = c("Calliphoridae", "Mesembrinellidae", "Null model", "Sarcophagidae"))
line_null2 <- LTf[LTf$group == "Null model",]

fig_4_land <- ggplot(LTf, aes(y = Probability.of.Research, x = land_tenure_reoder3, group = group, colour=group)) +
  geom_segment(data = line_null2, aes(x = land_tenure_reoder3-0.3, xend = land_tenure_reoder3+0.3, y = Probability.of.Research), linetype = "dashed", linewidth = 0.3, color = "gray30", alpha = .5) +
  geom_point(size=4, alpha=.8, position = position_dodge(width=0.1), show.legend = F) +
  scale_color_manual(values = myColors)+ 
  scale_y_continuous(limits=c(0.17,0.63), breaks=c(0.2, 0.4, 0.6))+
  labs(x=NULL, y="Probability of knowledge", tag = "(f)")+
  scale_x_continuous(limits=c(0.5,9.5), breaks = (9:1),
                     labels=c("Other", "Water \n body","Unassigned \n public", "Private \n land", "Rural \n settlements", "Strict \n Reserve","Protected \n area","Quilombola \n land", "Indigenous \n territories"))+
  theme_cowplot(font_family = "serif", font_size = 15)
fig_4_land 


#Fig 3 prancha====

plot_grid(plot_importance, f2, fig_4_land, 
                ncol = 1, align = "hv", rel_heights = c(1,2.2,1))

#ggsave("Figures/Figure 3_2025-08-08d.png", width = 11, height = 13, dpi = 400, bg = "white")
#ggsave("Figures/Figure 3_2025-08-08d.pdf", width = 11, height = 13)






#Fig 4A====
pacman::p_load(openxlsx, ggplot2, ggpubr, scico, OPTS, 
               ggimage)

data <- data.table::fread("Datasets_sp/VarImportance.csv", h = T, stringsAsFactors = F)

#groups = levels(as.factor(DM$group))
data$group = as.factor(data$group)
data$avg_IMSE = as.numeric(data$avg_IMSE)
data$min_IMSE = as.numeric(data$min_IMSE)
data$max_IMSE = as.numeric(data$max_IMSE)

library(ggplot2)
library(dplyr)
library(forcats)

#family_colors <- c(
#  "Calliphoridae" = "#6B9D59",
#  "Mesembrinellidae" = "#E093C3",
#  "Sarcophagidae" = "#E47352",
#  "Null" = "black")


# 1. Ordem das espécies com null model primeiro
ordem_especies <- c(
  "Null model - Amazon",
  "Chloroprocta idioidea - Calliphoridae",
  "Chrysomya albiceps - Calliphoridae",
  "Cochliomyia macellaria - Calliphoridae",
  "Hemilucilia semidiaphana - Calliphoridae",
  "Lucilia eximia - Calliphoridae",
  "Mesembrinella batesi - Mesembrinellidae",
  "Mesembrinella bellardiana - Mesembrinellidae",
  "Mesembrinella bicolor - Mesembrinellidae",
  "Mesembrinella quadrilineata - Mesembrinellidae",
  "Mesembrinella randa - Mesembrinellidae",
  "Oxysarcodexia intona - Sarcophagidae",
  "Oxysarcodexia thornax - Sarcophagidae",
  "Peckia (Pattonella) intermutans - Sarcophagidae",
  "Peckia (Peckia) chrysostoma - Sarcophagidae",
  "Peckia (Sarcodexia) lambens - Sarcophagidae"
)

# 2. Fator data com a ordem
data$group[data$group == "Null_model"] <- "Null model - Amazon"

data$Species <- factor(data$group, levels = ordem_especies)

# 3. espécies em itálico menos o null model
species_labels <- c(
  "Null model - Amazon" = "Null model",
  "Chloroprocta idioidea - Calliphoridae" = expression(italic("Chloroprocta idioidea")),
  "Chrysomya albiceps - Calliphoridae" = expression(italic("Chrysomya albiceps")),
  "Cochliomyia macellaria - Calliphoridae" = expression(italic("Cochliomyia macellaria")),
  "Hemilucilia semidiaphana - Calliphoridae" = expression(italic("Hemilucilia semidiaphana")),
  "Lucilia eximia - Calliphoridae" = expression(italic("Lucilia eximia")),
  "Mesembrinella batesi - Mesembrinellidae" = expression(italic("Mesembrinella batesi")),
  "Mesembrinella bellardiana - Mesembrinellidae" = expression(italic("Mesembrinella bellardiana")),
  "Mesembrinella bicolor - Mesembrinellidae" = expression(italic("Mesembrinella bicolor")),
  "Mesembrinella quadrilineata - Mesembrinellidae" = expression(italic("Mesembrinella quadrilineata")),
  "Mesembrinella randa - Mesembrinellidae" = expression(italic("Mesembrinella randa")),
  "Oxysarcodexia intona - Sarcophagidae" = expression(italic("Oxysarcodexia intona")),
  "Oxysarcodexia thornax - Sarcophagidae" = expression(italic("Oxysarcodexia thornax")),
  "Peckia (Pattonella) intermutans - Sarcophagidae" = expression(italic("Peckia (Pattonella) intermutans")),
  "Peckia (Peckia) chrysostoma - Sarcophagidae" = expression(italic("Peckia (Peckia) chrysostoma")),
  "Peckia (Sarcodexia) lambens - Sarcophagidae" = expression(italic("Peckia (Sarcodexia) lambens"))
)

nome_variaveis <- c(
  "1" = "Degradation",
  "2" = "Dry season\n length",
  "3" = "Land\n ownership",
  "4" = "Travel time",
  "5" = "Distance to\n research centers"
  # Ajuste os nomes conforme o que você tem em data1$Variables_reoder
)



library(ggplot2)
library(forcats)
library(stringr)
line_null3 <- data[data$Species == "Null model - Amazon",]
data$Species

data2 <- data %>%
  mutate(Species1 = str_trim(str_extract(Species, "^[^-]+")))

ordem_especies2 <- c(
  "Null model",
  "Chloroprocta idioidea",
  "Chrysomya albiceps",
  "Cochliomyia macellaria",
  "Hemilucilia semidiaphana",
  "Lucilia eximia",
  "Mesembrinella batesi",
  "Mesembrinella bellardiana",
  "Mesembrinella bicolor",
  "Mesembrinella quadrilineata",
  "Mesembrinella randa",
  "Oxysarcodexia intona",
  "Oxysarcodexia thornax",
  "Peckia (Pattonella) intermutans",
  "Peckia (Peckia) chrysostoma",
  "Peckia (Sarcodexia) lambens"
)
data2$Species3 <- factor(data2$Species1, levels = ordem_especies2)

species_labels2 <- lapply(ordem_especies2, function(x) parse(text = paste0("italic('", x, "')")))


data2.1 <- data2 %>%
  mutate(
    genero   = str_extract(Species3, "^[^ ]+"),
    epiteto  = str_extract(Species3, "[^ ]+$"),
    abreviado = paste(str_sub(genero, 1, 3), str_sub(epiteto, 1, 4), sep = "."))
data2.1$abreviado2 <- data2.1$abreviado
data2.1[(duplicated(data2.1$abreviado2)),"abreviado2"] <- NA





library(ggrepel)
p1 <- ggplot(data2, aes(x = avg_IMSE, y = Species3, color = Family)) +
  #geom_vline(data = line_null3, aes(xintercept = avg_IMSE), linetype = "dashed", linewidth = 0.3, color = "gray30", alpha = .5) +
  geom_vline(data = data1, aes(xintercept = avg_IMSE, color = Families), linetype = "dashed", linewidth = 0.8, alpha = .3) +
  geom_errorbar(aes(xmin = min_IMSE, xmax = max_IMSE), width = 0.25, linewidth = .6, alpha = 1, color = "gray60") +
  geom_point(position=position_dodge(width = 0.52), size = 3, alpha=.8) +
  geom_text_repel(data = data2.1, aes(x = avg_IMSE, y = Species3, label=abreviado2, color = Family), size=2.5, fontface = "bold.italic", direction = "y", show.legend = F, nudge_y  = 0.3)+
  scale_color_manual("Families", values = myColors) +
  scale_shape_manual("Species", values = c(1, rep(c(16, 15, 17, 18, 6), 3)), labels = species_labels2)+
  labs(x = "Increase in MSE", y = "", tag = "(a)") +
  scale_x_continuous(breaks = c(0, 0.1, 0.2), limits = c(-0.05,0.23))+
  #xlim(min(data1$min_IMSE) * 1, max(data1$max_IMSE) * 1) +
  theme_cowplot(font_family = "serif", font_size = 15)+
  guides(
    shape = guide_legend(ncol = 5),
    color = guide_legend(ncol = 1))+
  theme(legend.position = "top", 
        legend.direction = "horizontal",
        #legend.box = "horizontal",
        #legend.spacing.x = unit(0.1, "mm"),
        #legend.text = element_text(size = 9),
        #legend.title = element_text(size = 9),
        strip.background = element_rect(fill = NA,color = "black",linewidth = 0.5),
        strip.text = element_text(color = "black"),
        axis.text.y = element_blank())+
  facet_wrap(~ Variables_reoder, scales = "free_x", nrow = 1, strip.position = "top", labeller = labeller(Variables_reoder = nome_variaveis))


p1


#ggsave("Figures/Figure 4.1_2025-08-08d.png", width = 15, height = 8, dpi = 400, bg = "white")
#ggsave("Figures/Figure 4.1_2025-08-08d.pdf", width = 15, height = 8)






#Fig 4B-E====
pacman::p_load(openxlsx, ggplot2, ggpubr, scico, ggimage, rsvg)
library(readxl)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(scico)
library(ggimage)
library(rsvg)


dir()
DE = read.xlsx("Datasets_sp/figure5_A-D_sp.xlsx", sheet = 1) #degradation
DM = read.xlsx("Datasets_sp/figure5_A-D_sp.xlsx", sheet = 2) #drylength
ED = read.xlsx("Datasets_sp/figure5_A-D_sp.xlsx", sheet = 3) #research
TT = read.xlsx("Datasets_sp/figure5_A-D_sp.xlsx", sheet = 4) #traveltime
colnames(DM)

## change values to 0-1
range01 <- function(x){(x-min(x,na.rm=T))/(max(x, na.rm=T)-min(x,na.rm=T))}
DM$DM = range01(DM$dry_months)
DE$DE = range01(DE$degradation)
ED$ED = range01(ED$education)
TT$TT = range01(TT$travel_time)

#groups = levels(as.factor(DM$group))
DM$group = factor(DM$group, levels = c( "Chloroprocta_idioidea",
                                        "Chrysomya_albiceps",
                                        "Cochliomyia_macellaria",
                                        "Hemilucilia_semidiaphana",
                                        "Lucilia_eximia",
                                        "Mesembrinella_batesi",
                                        "Mesembrinella_bellardiana",
                                        "Mesembrinella_bicolor",
                                        "Mesembrinella_quadrilineata",
                                        "Mesembrinella_randa",
                                        "Null_model",
                                        "Oxysarcodexia_intona",
                                        "Oxysarcodexia_thornax",
                                        "Peckia_(Pattonella)_intermutans",
                                        "Peckia_(Peckia)_chrysostoma",
                                        "Peckia_(Sarcodexia)_lambens"))
DE$group = factor(DE$group, levels = c("Chloroprocta_idioidea",
                                       "Chrysomya_albiceps",
                                       "Cochliomyia_macellaria",
                                       "Hemilucilia_semidiaphana",
                                       "Lucilia_eximia",
                                       "Mesembrinella_batesi",
                                       "Mesembrinella_bellardiana",
                                       "Mesembrinella_bicolor",
                                       "Mesembrinella_quadrilineata",
                                       "Mesembrinella_randa",
                                       "Null_model",
                                       "Oxysarcodexia_intona",
                                       "Oxysarcodexia_thornax",
                                       "Peckia_(Pattonella)_intermutans",
                                       "Peckia_(Peckia)_chrysostoma",
                                       "Peckia_(Sarcodexia)_lambens"))
ED$group = factor(ED$group, levels = c("Chloroprocta_idioidea",
                                       "Chrysomya_albiceps",
                                       "Cochliomyia_macellaria",
                                       "Hemilucilia_semidiaphana",
                                       "Lucilia_eximia",
                                       "Mesembrinella_batesi",
                                       "Mesembrinella_bellardiana",
                                       "Mesembrinella_bicolor",
                                       "Mesembrinella_quadrilineata",
                                       "Mesembrinella_randa",
                                       "Null_model",
                                       "Oxysarcodexia_intona",
                                       "Oxysarcodexia_thornax",
                                       "Peckia_(Pattonella)_intermutans",
                                       "Peckia_(Peckia)_chrysostoma",
                                       "Peckia_(Sarcodexia)_lambens"))
TT$group = factor(TT$group, levels = c("Chloroprocta_idioidea",
                                       "Chrysomya_albiceps",
                                       "Cochliomyia_macellaria",
                                       "Hemilucilia_semidiaphana",
                                       "Lucilia_eximia",
                                       "Mesembrinella_batesi",
                                       "Mesembrinella_bellardiana",
                                       "Mesembrinella_bicolor",
                                       "Mesembrinella_quadrilineata",
                                       "Mesembrinella_randa",
                                       "Null_model",
                                       "Oxysarcodexia_intona",
                                       "Oxysarcodexia_thornax",
                                       "Peckia_(Pattonella)_intermutans",
                                       "Peckia_(Peckia)_chrysostoma",
                                       "Peckia_(Sarcodexia)_lambens"))




library(dplyr)
library(tidyr)
library(stringr)
DM2 <- DM %>%
  group_by(group_env) %>%
  slice_max(order_by = DM, n = 1, with_ties = FALSE) %>%
  ungroup()
DM3 <- DM2 %>%
  mutate(
    genero   = str_extract(group, "^[^_]+"),
    epiteto  = str_extract(group, "[^_]+$"),
    abreviado = paste(str_sub(genero, 1, 3), str_sub(epiteto, 1, 4), sep = ".")
  )

DMf2 <- DMf %>%
  group_by(group_env) %>%
  slice_max(order_by = DM, n = 1, with_ties = FALSE) %>%
  ungroup()


library(ggrepel)
DMp1 = 
  ggplot(DM, aes(DM, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_line(data=DMf, aes(x=DM, y=Probability.of.Research,,color = group_env, group = group_env), cex = 1.2, alpha=.6, show.legend = F)+ #Familias
  geom_text_repel(data = DM3[-11,], aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  geom_text_repel(data = DMf2, aes(label=group_env, color = group_env), size=3, fontface = "bold.italic", direction = "x", show.legend = F)+ #Familias
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  #scale_y_continuous(limits=c(0.2,0.7), breaks=c(0.2, 0.4, 0.6))+
  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab(NULL)+
  xlab("Dry Season Length")+
  theme_cowplot(font_family = "serif", font_size = 15)
DMp1



DE2 <- DE %>%
  group_by(group_env) %>%
  slice_max(order_by = DE, n = 1, with_ties = FALSE) %>%
  ungroup()
DE3 <- DE2 %>%
  mutate(
    genero   = str_extract(group, "^[^_]+"),
    epiteto  = str_extract(group, "[^_]+$"),
    abreviado = paste(str_sub(genero, 1, 3), str_sub(epiteto, 1, 4), sep = ".")
  )
DEf2 <- DEf %>%
  group_by(group_env) %>%
  slice_max(order_by = DE, n = 1, with_ties = FALSE) %>%
  ungroup()

DEp1 = 
  ggplot(DE, aes(DE, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_line(data=DEf, aes(x=DE, y=Probability.of.Research,,color = group_env, group = group_env), cex = 1.2, alpha=.6, show.legend = F)+ #Familias
  geom_text_repel(data = DE3[-11,], aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  geom_text_repel(data = DEf2, aes(label=group_env, color = group_env), size=3, fontface = "bold.italic", direction = "x", show.legend = F)+ #Familias
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  #scale_y_continuous(limits=c(0.2,0.7), breaks=c(0.2, 0.4, 0.6))+
  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab("Probability of knowledge")+
  xlab("Degradation")+
  theme_cowplot(font_family = "serif", font_size = 15)

DEp1


ED2 <- ED %>%
  group_by(group_env) %>%
  slice_max(order_by = ED, n = 1, with_ties = FALSE) %>%
  ungroup()
ED3 <- ED2 %>%
  mutate(
    genero   = str_extract(group, "^[^_]+"),
    epiteto  = str_extract(group, "[^_]+$"),
    abreviado = paste(str_sub(genero, 1, 3), str_sub(epiteto, 1, 4), sep = ".")
  )
EDf2 <- EDf %>%
  group_by(group_env) %>%
  slice_max(order_by = ED, n = 1, with_ties = FALSE) %>%
  ungroup()

EDp1 = 
  ggplot(ED, aes(ED, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_line(data=EDf, aes(x=ED, y=Probability.of.Research,,color = group_env, group = group_env), cex = 1.2, alpha=.6, show.legend = F)+ #Familias
  geom_text_repel(data = ED3[-11,], aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  geom_text_repel(data = EDf2, aes(label=group_env, color = group_env), size=3, fontface = "bold.italic", direction = "x", show.legend = F)+ #Familias
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  #scale_y_continuous(limits=c(0.2,0.7), breaks=c(0.2, 0.4, 0.6))+
  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab(NULL)+
  xlab("Distance to Research centers")+
  theme_cowplot(font_family = "serif", font_size = 15)
EDp1


TT2 <- TT %>%
  group_by(group_env) %>%
  slice_max(order_by = TT, n = 1, with_ties = FALSE) %>%
  ungroup()
TT3 <- TT2 %>%
  mutate(
    genero   = str_extract(group, "^[^_]+"),
    epiteto  = str_extract(group, "[^_]+$"),
    abreviado = paste(str_sub(genero, 1, 3), str_sub(epiteto, 1, 4), sep = ".")
  )
TTf2 <- TTf %>%
  group_by(group_env) %>%
  slice_max(order_by = TT, n = 1, with_ties = FALSE) %>%
  ungroup()

TTp1 = 
  ggplot(TT, aes(TT, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_line(data=TTf, aes(x=TT, y=Probability.of.Research,,color = group_env, group = group_env), cex = 1.2, alpha=.6, show.legend = F)+ #Familias
  geom_text_repel(data = TT3[-11,], aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  geom_text_repel(data = TTf2, aes(label=group_env, color = group_env), size=3, fontface = "bold.italic", direction = "x", show.legend = F)+ #Familias
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  #scale_y_continuous(limits=c(0.15,0.7), breaks=c(0.2, 0.4, 0.6))+
  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab("Probability of knowledge")+
  xlab("Travel time")+
  theme_cowplot(font_family = "serif", font_size = 15)

TTp1


f2.1 <- plot_grid(TTp1, EDp1, DEp1, DMp1, ncol = 2, align = "hv", labels = c("(a)","(b)","(c)","(d)"), label_fontfamily = "serif")
f2.1






#Fig 4F====
LT = read.xlsx("Datasets_sp/figure5_E_sp.xlsx")
LT$group = factor(LT$group, levels = c("Chloroprocta_idioidea",
                                       "Chrysomya_albiceps",
                                       "Cochliomyia_macellaria",
                                       "Hemilucilia_semidiaphana",
                                       "Lucilia_eximia",
                                       "Mesembrinella_batesi",
                                       "Mesembrinella_bellardiana",
                                       "Mesembrinella_bicolor",
                                       "Mesembrinella_quadrilineata",
                                       "Mesembrinella_randa",
                                       "Null_model",
                                       "Oxysarcodexia_intona",
                                       "Oxysarcodexia_thornax",
                                       "Peckia_(Pattonella)_intermutans",
                                       "Peckia_(Peckia)_chrysostoma",
                                       "Peckia_(Sarcodexia)_lambens"))




LT$group <- factor(LT$group, levels = c("Null_model", setdiff(levels(LT$group), "Null_model")))




legenda_labels <- c(
  "Null_model" = "Null model",
  "Chloroprocta_idioidea" = expression(italic("Chloroprocta idioidea")),
  "Chrysomya_albiceps" = expression(italic("Chrysomya albiceps")),
  "Cochliomyia_macellaria" = expression(italic("Cochliomyia macellaria")),
  "Hemilucilia_semidiaphana" = expression(italic("Hemilucilia semidiaphana")),
  "Lucilia_eximia" = expression(italic("Lucilia eximia")),
  "Mesembrinella_batesi" = expression(italic("Mesembrinella batesi")),
  "Mesembrinella_bellardiana" = expression(italic("Mesembrinella bellardiana")),
  "Mesembrinella_bicolor" = expression(italic("Mesembrinella bicolor")),
  "Mesembrinella_quadrilineata" = expression(italic("Mesembrinella quadrilineata")),
  "Mesembrinella_randa" = expression(italic("Mesembrinella randa")),
  "Oxysarcodexia_intona" = expression(italic("Oxysarcodexia intona")),
  "Oxysarcodexia_thornax" = expression(italic("Oxysarcodexia thornax")),
  "Peckia_(Pattonella)_intermutans" = expression(italic("Peckia (Pattonella) intermutans")),
  "Peckia_(Peckia)_chrysostoma" = expression(italic("Peckia (Peckia) chrysostoma")),
  "Peckia_(Sarcodexia)_lambens" = expression(italic("Peckia (Sarcodexia) lambens"))
)

LT2 <- LT %>%
  left_join(TT3, by = "group")


LT3 <- LT2 %>%
  mutate(
    genero   = str_extract(group, "^[^_]+"),
    epiteto  = str_extract(group, "[^_]+$"),
    abreviado = paste(str_sub(genero, 1, 3), str_sub(epiteto, 1, 4), sep = ".")
  )


line_null4 <- LT[LT$group == "Null_model",]

#Ordenar a variável para organizar a legenda
LT2$environment <- factor(LT2$environment, levels = c("Null model", "Calliphoridae", "Mesembrinellidae", "Sarcophagidae"))

fig_6 <- ggplot(LT2, aes(y = Probability.of.Research.x, x = land_tenure_reoder, colour=environment)) +
  #geom_segment(data = line_null4, aes(x = land_tenure_reoder-0.3, xend = land_tenure_reoder+0.3, y = Probability.of.Research), linetype = "dashed", linewidth = 0.3, color = "gray30", alpha = .5) +
  geom_segment(data = LTf, aes(x = land_tenure_reoder3-0.3, xend = land_tenure_reoder3+0.3, y = Probability.of.Research, colour = group_env), linetype = "dashed", linewidth = 0.7, alpha = .8) +
  geom_point(size=3, alpha=.8, position = position_dodge(width=0.3)) +
  geom_text_repel(data = LT3[-c(102:110),], aes(label=abreviado, color = environment), position = position_dodge(width=0.3), size=2, fontface = "bold.italic", show.legend = F, segment.size=.3, segment.alpha=.5)+
  scale_color_manual("Families", values = myColors)+
  labs(x=NULL, y="Probability of knowledge", tag = "(f)")+
  scale_x_continuous(limits=c(0.5,9.5), breaks = (9:1),
                     labels=c("Other", "Water \n body","Unassigned \n public", "Private \n land", "Rural \n settlements", "Strict \n Reserve","Protected \n area","Quilombola \n land", "Indigenous \n territories"))+
  theme_cowplot(font_family = "serif", font_size = 15)+
  theme(legend.position = "top", 
        legend.spacing.x = unit(3, "cm"),
        panel.spacing = unit(1.3, "lines"))


fig_6
#ggsave("Figures/Figure 4.2_2025-08-08d.png", width = 10, height = 8, dpi = 400, bg = "white")
#ggsave("Figures/Figure 4.2_2025-08-08d.pdf", width = 10, height = 12)



#Fig 4 prancha====
plot_grid(p1, f2.1, fig_6+theme(legend.position = "none"),ncol = 1, align = "hv", rel_heights = c(1.6,2.2,1.5))



#ggsave("Figures/Figure 4_2025-08-08d.png", width = 12, height = 16, dpi = 400, bg = "white")
#ggsave("Figures/Figure 4_2025-08-08d.pdf", width = 12, height = 16)




plot_grid(p1, fig_6+theme(legend.position = "none"),ncol = 1, align = "hv", rel_heights = c(1.3,2.2))




#Fig new pranchas====
p1+labs(x = "Increase in MSE", y = "Species", tag = NULL)
ggsave("Figures/Fig 4_2025-08-14d.pdf", width = 10, height = 7)
ggsave("Figures/Fig 4_2025-08-14d.png", width = 10, height = 7, dpi = 400, bg = "white")




fig_6+labs(tag = NULL)
ggsave("Figures/Fig 6_2025-08-14d.pdf", width = 10, height = 6)
ggsave("Figures/Fig 6_2025-08-14d.png", width = 10, height = 6, dpi = 400, bg = "white")




f2.1
ggsave("Figures/Fig 5(2x2)_2025-08-14d.png", width = 10, height = 7, dpi = 400, bg = "white")
ggsave("Figures/Fig 5(2x2)_2025-08-14d.pdf", width = 10, height = 7)









#Refazendo a 5====
DMp1 = 
  ggplot(DM, aes(DM, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_text_repel(data = DM3, aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab(NULL)+
  xlab("Dry Season Length")+
  theme_cowplot(font_family = "serif", font_size = 15)
DMp1


DEp1 = 
  ggplot(DE, aes(DE, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_text_repel(data = DE3, aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+

  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab(NULL)+
  xlab("Degradation")+
  theme_cowplot(font_family = "serif", font_size = 15)
DEp1


EDp1 = 
  ggplot(ED, aes(ED, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_text_repel(data = ED3, aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab(NULL)+
  xlab("Distance to Research centers")+
  theme_cowplot(font_family = "serif", font_size = 15)
EDp1


TTp1 = 
  ggplot(TT, aes(TT, Probability.of.Research))+
  geom_line(aes(color = environment, linetype = group_env), cex=.5, alpha=.8, show.legend = F)+
  geom_line(aes(color = environment, group = group_env), cex=.5, alpha=.5, show.legend = F)+
  geom_text_repel(data = TT3, aes(label=abreviado, color = environment), size=2, fontface = "bold.italic", direction = "x", show.legend = F)+ #Especies
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_linetype_manual("Species", values = c(1, rep(c(2, 3, 4, 5, 6), 3)))+
  scale_color_manual("Familes", values = myColors) +
  ylab(NULL)+
  xlab("Travel time")+
  theme_cowplot(font_family = "serif", font_size = 15)

TTp1


plot_grid(TTp, TTp1,
          EDp, EDp1,
          DEp, DEp1,
          DMp, DMp1, ncol = 2, align = "hv", labels = c("(a)","(b)","(c)","(d)","(e)","(f)", "(g)", "(h)"), label_fontfamily = "serif")
ggsave("Figures/Fig 5(4x2)_2025-08-11d.png", width = 10, height = 11, dpi = 400, bg = "white")
ggsave("Figures/Fig 5(4x2)_2025-08-11d.pdf", width = 10, height = 11)
