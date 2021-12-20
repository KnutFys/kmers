# AUSTRALIA final R script

# LOAD LIBRARIES
library(vegan) # Needed for distance based rda
library(data.table) # Needed for faster fread funtion
library(ggplot2) # Needed for plotting
library(ggmap) # Needed for creating backgound maps
library(grid) # Needed for inset of grobs
library(gridExtra) # Needed for grid.arange?
library(bit64) # Needed to handle kmer counts in int64
###################################################
# SET PATHS 
###################################################
project_name = 'scripttest'
data.path <- '/media/knut/Storage/Aussie/'
project_path <- paste(c('Projects/', project_name, '/'), sep = '', collapse = '')
####################################################
# SET SEED
####################################################
set.seed(1)
####################################################
# Format file paths
####################################################
bc <- paste(c(
  project_path, project_name, '_', 'braycurtis.csv'), sep = '', collapse = '')
angular <- paste(c(
  project_path, project_name, '_', 'angular.csv'), sep = '', collapse = '')
js <- paste(c(
  project_path, project_name, '_', 'jensenshannon.csv'), sep = '', collapse = '')
metafile <- paste(c(
  project_path, project_name, '_meta.csv'), sep='', collapse='')

plotfile <- paste(c(
  project_path, project_name, '_plot', '.png'), sep='', collapse='')

###################################################
# LOAD DATA
###################################################
data <- as.data.frame(fread(metafile), stringsAsFactors=TRUE)
# Convert kmer_count to billion kmers. R does not play well with 
# int64
data$kmer_count <- data$kmer_count/1000000000
dm.angular <- fread(angular)[,-1]
dm.js <- fread(js)[,-1]
dm.bc <- fread(bc)[,-1]

############################################################
# Subsets
############################################################
# Create Pacific Ocean and Australian subsets
data.australia <- data[data$longitude > 0,]
data.pacific <- data[data$longitude < 0,]
# Create Illumina DNA prep
data.prep <- data[data$library_construction_protocol == 'Illumina DNA Prep',]

#############################################################
# Functions
#############################################################

#########################
# Helper functions
#########################
g_legend <- function(a.gplot){
  # Extract legend grob
  # Function copied from
  # https://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}

split.name <- function(x){
  #Splits geographic name scheme  IE "Australia:Perth" -> "Perth"
  y <- strsplit(x, ':')
  return(y[[1]][2])
}
#######################################################
# Maps
######################################################
create.map <- function(location){
  # Creates maps to plot sample locations on.
  # Downloads from stamen maps. Manually selected area.
  if (location=='Australia') {
    map <- get_stamenmap(
      bbox = c(left=112, bottom=-45, right=160, top= -8), zoom = 4
    )
  }
  if (location=='Pacific') {
    map <- get_stamenmap(
      bbox = c(left=-179, bottom=-67, right=-140, top= 1), zoom = 1
    )
  }
  if (location=='Svalbard') {
    map <- get_stamenmap(
      bbox = c(left=10.8, bottom=78, right=17, top= 79.5), zoom = 4,
      maptype = 'terrain-background'
    )
  }
  map <- ggmap(map)
  return(map)
}

create.mapbox <- function(data, border=0.1, zoom=6) {
  # Helper function to find submaps. Attemts to find 
  # a suitable map for supplied dataset
  lat_min <- round(min(data$latitude) - border, digits=1)
  lat_max <- round(max(data$latitude) + border, digits=1)
  lon_min <- round(min(data$longitude) - border, digits=1)
  lon_max <- round(max(data$longitude) + border, digits=1)
  map.info <- get_stamenmap(
    bbox = c(left=lon_min, bottom=lat_min, right=lon_max, top= lat_max),
    zoom = zoom, maptype = 'terrain-background'
  )
  map <- ggmap(map.info)
  map <- map + theme_void() + theme(legend.position = 'none')
  return(map)
}

fill.mapbox <- function(data, map, aesthetic) {
  # Helper function to plot points on submap
  map <- map + geom_point(data=data, aesthetic, size=3)
  map <- map + theme(panel.border = element_rect(fill = NA, color = 'black'))
  return(map)
}

############################
# Plotters
############################
fourplot <- function(a,b,c,d, out_file='multiplottest.png'){
  # Join four plots with same legend
  a <- a + theme(legend.position = 'bottom', legend.box = "vertical")
  b <- b + theme(legend.position = 'none')
  c <- c + theme(legend.position = 'none')
  d <- d + theme(legend.position = 'none')
  l <- g_legend(a)
  a <- a + theme(legend.position = 'none')
  m <- grid.arrange(arrangeGrob(a ,b, ncol = 2),
                    arrangeGrob(c,d,ncol = 2),
                    top = 'Metagenomic samples from Austalia and the Pacific ',
                    bottom = l)
  ggsave(out_file, plot=m, device = 'png',
         width = 14, height = 14)
}


trend.alpha <- function(data, measure) {
  measure = sym(measure)
  g <- ggplot(data=data, aes(x=kmer_count, y=!!measure))
  g <- g + geom_point(aes(col=library_construction_protocol, shape=sample_type))
  g <- g + geom_smooth(method = 'lm')
  g <- g + labs(col='Library construction protocol') + xlab('Kmer count (Billions)')
  return(g)
}


runNMDS <- function(distanceMatrix, metaData, col, mtitle='Main title',
                    stitle='Sub title') {
  # Creates and plots NMDS representation of Australian
  # metagenomes. Four plots
  col <- sym(col)
  #shape <- sym(shape)
  # Perform multidimentional scaling
  nmds <- metaMDS(distanceMatrix, autotransform = FALSE, try = 10)
  # Add position vectors to data table
  metaData$X <- nmds$points[,1]
  metaData$Y <- nmds$points[,2]
  # Create plot
  g <- ggplot(metaData, aes(x=X, y=Y))
  g <- g + geom_point(aes(col=!!col), size=1.5)
  stress <- paste('Stress: ',
                  round(nmds$stress, digits = 4), sep='',collapse = '')
  g <- g + labs(title = mtitle,
                subtitle = stitle,
                caption = stress,
                col='Sample type')
  g <- g + xlab('NMDS 1') + ylab('NMDS 2')
  return(g)
}

runpcoa <- function(distanceMatrix, metaData, col, mtitle='Main title',
                    stitle='Sub title') {
  # Creates and plots NMDS representation of Australian
  # metagenomes. Four plots
  col <- sym(col)
  #shape <- sym(shape)
  # Perform multidimentional scaling
  pcoa <- cmdscale(distanceMatrix)
  # Add position vectors to data table
  metaData$X <- pcoa[,1]
  metaData$Y <- pcoa[,2]
  # Create plot
  g <- ggplot(metaData, aes(x=X, y=Y))
  g <- g + geom_point(aes(col=!!col), size=1.5)
  g <- g + labs(title = mtitle,
                subtitle = stitle,
                col='Sample type')
  g <- g + scale_colour_manual(values=c('red', 'green', 'black', 'blue'))
  g <- g + xlab('PCOA 1') + ylab('PCOA 2')
  return(g)
}

plotfour<- function() {
  # Create plots for thesis (PCOA)
  set.seed(1) 
  data$sample_type[data$sample_type=='Pelagic' & 
                     data$geo_loc_name=='Pacific Ocean'] <- 'Pelagic (Pacific)'
  a <- runpcoa(dm.angular, data, col='sample_type',mtitle = 'a)',
               stitle = 'Angular distance')
  b <- runpcoa(dm.js, data, col='sample_type',mtitle = 'b)',
               stitle = 'Jensen-Shannon distance')
  c <- runpcoa(dm.bc, data, col='sample_type',mtitle = 'c)',
               stitle = 'Bray-Curtis dissimilarity') 
  ## MAP #################
  aesthetic <- aes(x=longitude, y=latitude, col=sample_type)
  # Create australian map
  main.map <- create.map('Australia')
  #main.map <- fill.mapbox(data.australia, main.map, aesthetic)
  main.map <- main.map + geom_point(data=data, aesthetic, size=2)
  main.map <- main.map + scale_colour_manual(values=c('red', 'green', 'black', 'blue'))
  main.map <- main.map + xlab('Longitude') + ylab('Latitude')
  main.map <- main.map + labs(title = 'Sample sites, Australia',
  caption='Map tiles by Stamen Design, under CC BY 3.0. Data by OpenStreetMap, under ODbL.')
  fourplot(a,b,c,main.map, out_file = 'multiaust.png')
}

runadonis <- function(dm, data) {
  # The function performs permutational anova test on
  # the supplied distance matrix to test if the the grouping
  # of samples by sample_type is significant in explaining the
  # distances. If Pr(>F) < 0.05 performs aditional test 
  # for pairs of sample_types
  result <- adonis(dm ~ sample_type, data=data)
  result$name <- 'main'
  results <- list()
  results[['main']] <- result
  result.sig <- result$aov.tab$`Pr(>F)`[1]
  if (result.sig < 0.05) {
    selection.pc <- data$sample_type == 'Pelagic' | data$sample_type == 'Coastal water'
    selection.ps <- data$sample_type == 'Pelagic' | data$sample_type == 'Sediment'
    selection.sc <- data$sample_type == 'Sediment' | data$sample_type == 'Coastal water'
    
    data.pc <- data[selection.pc,]
    dm.pc <- as.data.frame(dm)[selection.pc, selection.pc]
    result.pc <- adonis(dm.pc ~ sample_type, data=data.pc)
    result.pc$name <- 'PelagicCoastal'
    results[['PelagicCoastal']] <- result.pc
    
    data.ps <- data[selection.ps,]
    dm.ps <- as.data.frame(dm)[selection.ps, selection.ps]
    result.ps <- adonis(dm.ps ~ sample_type, data=data.ps)
    result.ps$name <- 'PelagicSediment'
    results[['PelagicSediment']] <- result.ps
    
    
    data.sc <- data[selection.sc,]
    dm.sc <- as.data.frame(dm)[selection.sc, selection.sc]
    result.sc <- adonis(dm.sc ~ sample_type, data=data.sc)
    result.sc$name <- 'SedimentCoastal'
    results[['SedimentCoastal']] <- result.sc
  }
  return(results)
}

adonisforall <- function(){
  # Runs all oerutational anova tests and prints
  # result to console
  set.seed(1)
  extract <- function(x) {
    for (a in x) {
      cat(
        a$name, ' ',
        a$aov.tab$`Pr(>F)`[1], ' ',
        a$aov.tab$R2, '\n'
      )
    } 
  }
  angular.anova <- runadonis(dm.angular, data)
  js.anova <- runadonis(dm.js, data)
  bc.anova <- runadonis(dm.bc, data)
  cat('Angular distance', '\n')
  extract(angular.anova)
  cat('Jensen-Shannon distance', '\n')
  extract(js.anova)
  cat('Bray-Curtis distance', '\n')
  extract(bc.anova)
}

sample.map <- function() {
  # Create maps for report
  out_file='mapsAusPacauto.png'
  mappointsize = 2.5
  australia <- create.map('Australia')
  aestethic <- aes(x=longitude, y=latitude, col=geo_loc_name,
                   shape = sample_type)
  
  australia <- australia + geom_point(data = data, aestethic,
                                      size=mappointsize)
  australia <- australia + xlab('Longitude') + ylab('Latitude')
  pacific <- create.map('Pacific')
  pacific <- pacific + geom_point(data = data, aestethic,
                                      size=mappointsize)
  pacific <- pacific + xlab('Longitude') + ylab('Latitude') + labs(
    col='Geographic area', shape='Sample type')
  
  #australia <- australia + theme(legend.position = 'bottom', legend.box = "vertical")
  #pacific <- pacific + theme(legend.position = 'none')
  #l <- g_legend(australia)
  australia <- australia + theme(legend.position = 'none')
  m <- grid.arrange(arrangeGrob(australia ,pacific, ncol = 2),
                      top = 'Sample locations Austalia and the Pacific ')
                      #bottom = l)
  ggsave(out_file, plot=m, device = 'png')
           #width = 13, height = 9)
}

############################################
# Run relevant functions
############################################
plotfour()
#adonisforall()
