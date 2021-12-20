# R script used for Svalbard dataset 



dm.path <- './Dm'
comparative.path <- './Comparative'
plot.path <- './Plots'
log.path <- './Logs'
meta.path <- './Metadata'
alpha.path <- './Alpha'


library(vegan)
library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

# Load old metadata used in K-mer analysis
# This data was log transformed in python and lacked some dscriptors
# Used for comparison of RDA based on genus with dbRDA based on kmrs
metaold <- fread(
    paste(meta.path, '/metadata_log.csv', sep='', collapse=''))
alpha.diversity <- fread(paste(
  alpha.path, '/alphadiversity.csv', sep='', collapse=''))[,-1]
# Purtify column names
colnames(metaold)[c(2,3,4,5,6,7)] <- c('TOC','DOC','TN','DN', 'TP', 'DP')
#create a table including alpha diversities
metaalpha <- merge(metaold, alpha.diversity, by.x = 'sample_id', by.y = 'Sample', sort=T)
# Convert JS distance to similarity
metaalpha$Js <- 1 - metaalpha$Js
# Convert Shannons index to pielous index
metaalpha$Pielou <- metaalpha$Shannon/log(4^11/2)
#Remove sample id column
metaold <- metaold[,-1]
#### Load new metadata ###############################
metadata <- fread(paste(meta.path, '/svfi19.csv', sep='', collapse=''))
samplenames <- metadata$sample_id[-10]
#Remove sample 11 (9 allready removed)
metadata <- metadata[-10,]
#Remove unused descriptors
metadata.r <- metadata[, -c(2,3,4,5,7,8,10,11,13,14,15,22,23,24,25,26,28,29,30,31,32)]
#Add age factor def  young < 550m from glacier old >3350m
metadata.r$Age <- with(metadata.r, ifelse(gl_dist < 551, 'Young', 'Med'))
metadata.r$Age <- with(metadata.r, ifelse(gl_dist > 3350, 'Old', Age))
metadata.r$bird <- with(metadata.r, ifelse(bird == 1, 'Yes', 'No'))
setnames(metadata.r, 'bird', 'Birds')
#Set birds and age as factors
metadata.r$Age <- as.factor(metadata.r$Age)
metadata.r$Birds <- as.factor(metadata.r$Birds)
metadata.r <- metadata.r[,-1]

###################################################
# Load Kegg table
###################################################
loadKegg <- function(){
  # Load and format KEGG table 
  kegg <- read.table(paste(
    comparative.path, '/TpmKeggTable.csv', sep='', collapse=''),
    sep = ',', header = TRUE, row.names = 1)
  kegg[is.na(kegg)] = 0
  kegg <- as.data.frame(t(kegg))
  return(kegg)
}

###################################################
# Load distance files
###################################################
load.fulldms <- function(k=11){
    # Loads the 3 main distance matrices and removes a redundant 
    # column containing sample IDs
    if (k==11) {
      dm <- 'Distance matrices'
      dm$js <- fread(paste(
          dm.path, '/full_jensenshannon.csv', sep='', collapse=''))[,-1]
      dm$ang <- fread(paste(
          dm.path, '/full_angular.csv', sep='', collapse=''))[,-1]
      dm$bc <- fread(paste(
          dm.path, '/full_braycurtis.csv', sep='', collapse=''))[,-1]
    }
    if (k==13) {
      dm$js <- fread(paste(
        dm.path, '/JensenShannonFull13.csv', sep='', collapse=''))[,-1]
      dm$ang <- fread(paste(
        dm.path, '/AngularFull13.csv', sep='', collapse=''))[,-1]
      dm$bc <- fread(paste(
        dm.path, '/BrayCurtisFull13.csv', sep='', collapse=''))[,-1]
    }
    return(dm)
}

birdageplot <- function(distmat){
    # Plots bird presence, age and DN conc on
    # Jensen shannon 11mer nmds
    set.seed(1)
    nmds <- metaMDS(distmat, autotransform = F,try=99, trymax = 999)
    tmp <- metadata.r
    tmp$DN <- metadata$DN
    tmp$X <- nmds$points[,1]
    tmp$Y <- nmds$points[,2]
    tmp$Sample <- samplenames
    stress <- paste('Stress: ',
                    round(nmds$stress, digits = 4), sep='',collapse = '')
    g <- ggplot(tmp, aes(x=X, y=Y))
    g <- g + geom_point(aes(col=Age, size=DN, shape=Birds))
    g <- g + geom_text(aes(label=Sample), size=2.5
                       ,position = position_nudge(y=-0.01))
    g <- g + labs(title = ' a) Influence of birds and age on DN',
                  subtitle = 'NMDS based on Jensen-Shannon distance',
                  caption = stress)
    g <- g + xlab('NMDS 1') + ylab('NMDS 2')
    return(g)
}

########################################################
# Load taxonomic abundance data
#######################################################
clean <- function(taxatable) {
  # Removes unclassified reads and viruses from taxonomic tables
  taxatable <- taxatable %>% select(-contains(c(
    'Unclassified', 'Unmapped', 'viruses', 'no genus', 'no order')))
  return(taxatable)
}

loadTaxa <- function(x, tform=FALSE){
    #Load taxonomic abundance data 
    #x='d' selects domain
    #x='p' selects phyla etc down to genus
    
    taxa <- fread(paste(
        comparative.path, '/genus_abund_table.csv', sep='', collapse=''))
    # Set NA to 0
    taxa[is.na(taxa)] <- 0
    # Remove Finse
    taxa <- taxa[,-c(7,8,9,10)]
    
    if (x=='d') {
        # Select only domain
        # note that Domain is labeled as kingdom in dataset
        taxa <- taxa[,-c(2,3,4,5,6)]
        taxa <- aggregate(.~Kingdom, data = taxa, FUN='sum')
    }
    if (x=='p') {
        # Select only phyla
        taxa <- taxa[,-c(1,3,4,5,6)]
        taxa <- aggregate(.~Phyla, data = taxa, FUN='sum')
    }
    if (x=='c') {
        # Select only class
        taxa <- taxa[,-c(1,2,4,5,6)]
        taxa <- aggregate(.~Class, data = taxa, FUN='sum')
    }
    if (x=='o') {
        # Select only order
        taxa <- taxa[,-c(1,2,3,5,6)]
        taxa <- aggregate(.~Order, data = taxa, FUN='sum')
    }
    if (x=='f') {
        # Select only family
        taxa <- taxa[,-c(1,2,3,4,6)]
        taxa <- aggregate(.~Family, data = taxa, FUN='sum')
    }
    if (x=='g') {
        # Select only genus
        taxa <- taxa[,-c(1,2,3,4,5)]
        taxa <- aggregate(.~Genus, data = taxa, FUN='sum')
    }
    taxa <- transpose(taxa)
    colnames(taxa) <- taxa[1,]
    taxa <- taxa[-1,]
    row.names(taxa) = samplenames
    taxa[,1:length(taxa)] <- sapply(taxa[,1:length(taxa)], as.numeric)
    colnames(taxa) <- lapply(
      colnames(taxa), function(x){gsub('^.?_', '',x)})
    # Hellinger transformation of community counts
    if (tform) {
        taxa.t <- decostand(taxa, method = 'hellinger')
        return(taxa.t)
    } else {
        return(taxa)
    }
}
#############################

####################################
# Find best model -- not used
###################################
findbest <- function(community, method){
    if (method=='rda') {
        intercept <- rda(community ~ 1, data=metadata.r.s)
        full <- rda(community ~ ., data=metadata.r.s)
    }
    if (method=='dbrda') {
        intercept <- dbrda(community ~ 1, data=metadata.r.s)
        full <- dbrda(community ~ ., data=metadata.r.s)
    }
    optimal <- ordistep(intercept, scope = full, direction = 'both')
    return(optimal)
}

taxvsk <- function(){
    # Check similarity of distances between
    # various kmer lengths and taxonomic ranks
    dm5 <- fread(paste(
        dm.path, '/5_100000_js.csv', sep='', collapse=''))[,-1]
    dm5 <- as.dist(dm5)
    dm7 <- fread(paste(
        dm.path, '/7_100000_js.csv', sep='', collapse=''))[,-1]
    dm7 <- as.dist(dm7)
    dm9 <- fread(paste(
        dm.path, '/9_100000_js.csv', sep='', collapse=''))[,-1]
    dm9 <- as.dist(dm9)
    dm11 <- fread(paste(
        dm.path, '/11_100000_js.csv', sep='', collapse=''))[,-1]
    dm11 <- as.dist(dm11)
    dmfull <- fread(paste(
        dm.path, '/full_jensenshannon.csv', sep='', collapse=''))[,-1]
    dmfull <- as.dist(dmfull)
    kin <- vegdist(
        decostand(loadTaxa('d', tform = F), method = 'total'), method = 'bray')
    ord <- vegdist(
        decostand(loadTaxa('o', tform = F), method = 'total'), method = 'bray')
    gen <- vegdist(
        decostand(loadTaxa('g', tform = F), method = 'total'), method = 'bray')
    m <- 'Mantel comparisons'
    m$k5kin <- mantel(dm5, kin)
    m$k5ord <- mantel(dm5, ord)
    m$k5gen <- mantel(dm5, gen)
    m$k7kin <- mantel(dm7, kin)
    m$k7ord <- mantel(dm7, ord)
    m$k7gen <- mantel(dm7, gen)
    m$k9kin <- mantel(dm9, kin)
    m$k9ord <- mantel(dm9, ord)
    m$k9gen <- mantel(dm9, gen)
    m$k11kin <- mantel(dm11, kin)
    m$k11ord <- mantel(dm11, ord)
    m$k11gen <- mantel(dm11, gen)
    m$fullkin <- mantel(dmfull, kin)
    m$fullord <- mantel(dmfull, ord)
    m$fullgen <- mantel(dmfull, gen)
    # Collect statistics in a table
    m$statistic <- data.table(
        c(
            m$k5kin$statistic,
            m$k7kin$statistic,
            m$k9kin$statistic,
            m$k11kin$statistic,
            m$fullkin$statistic
        ),
        c(
            m$k5ord$statistic,
            m$k7ord$statistic,
            m$k9ord$statistic,
            m$k11ord$statistic,
            m$fullord$statistic
        ),
        c(
            m$k5gen$statistic,
            m$k7gen$statistic,
            m$k9gen$statistic,
            m$k11gen$statistic,
            m$fullgen$statistic
        )
        
    )
    row.names(m$statistic) = c('K=5','K=7','K=9', 'K=11', 'K=11 Full')
    colnames(m$statistic) = c('Domain','Order','Genus')
    
    m$sig <- data.table(
        c(
            m$k5kin$signif,
            m$k7kin$signif,
            m$k9kin$signif,
            m$k11kin$signif,
            m$fullkin$signif
        ),
        c(
            m$k5ord$signif,
            m$k7ord$signif,
            m$k9ord$signif,
            m$k11ord$signif,
            m$fullord$signif
        ),
        c(
            m$k5gen$signif,
            m$k7gen$signif,
            m$k9gen$signif,
            m$k11gen$signif,
            m$fullgen$signif
        )
        
    )
    row.names(m$sig) = c('K=5','K=7','K=9', 'K=11','K=11 Full')
    colnames(m$sig) = c('Domain','Order','Genus')
    return(m)
}

###################################
# Taxonomic RDA
###################################
genusrda <- function(){
    set.seed(2)
    # Load taxonomic data and remove unclassified genera, 
    # genara with no NCBI match and phages
    tax <- loadTaxa('g')
    tax <- clean(tax)
    # Perform Hellinger transformation
    tax <- decostand(tax, method = 'hellinger')
    rd <- rda(tax ~ . , data=metaold)
    rd$title <- c(
        'd) Alignment based redundancy analysis',
        'Taxonomic rank: Genus'
    )
    rd$anova <- anova(rd)
    return(rd)
}

domainsrda <- function(){
    set.seed(3)
    # Load taxonomic data and remove unc, 
    # viruses and unknown hits
    tax <- loadTaxa('d')
    tax <- clean(tax)
    # Perform Hellinger transformation
    tax <- decostand(tax, method = 'hellinger')
    rd <- rda(tax ~ . , data=metaold)
    rd$title <- c(
        'b) Alignment based redundancy analysis',
        'Taxonomic rank: Domain'
    )
    rd$anova <- anova(rd)
    return(rd)
}

orderrda <- function(){
    set.seed(4)
    # Load taxonomic data and remove unc, 
    # viruses and unknown hits
    tax <- loadTaxa('o')
    tax <- clean(tax)
    # Perform Hellinger transformation
    tax <- decostand(tax, method = 'hellinger')
    rd <- rda(tax ~ . , data=metaold)
    rd$title <- c(
        'c) Alignment based redundancy analysis',
        'Taxonomic rank: Order'
    )
    rd$anova <- anova(rd)
    return(rd)
}

rdaplotter <- function(x,sc=2) {
    # Plot rda object, default scaling  for sites 
    # is weighted average of response variables, 
    # scaling =1 gives sites as a linear comb of
    # explanatory variables
    ordination <- summary(x, scaling=sc)
    plotdata <- data.frame(ordination$sites[,c(1,2)])
    plotdata$names <- samplenames
    plotdata$Age <- metadata.r$Age
    plotdata$Birds <- metadata.r$Birds
    plotdata$DN <- metadata$DN
    # Pick out top species
    species <- as.data.frame(ordination$species[,c(1,2)])
    species$length <- sqrt(rowSums(species^2))
    species.top <- species[species$length > 0.08, ]
    species.top <- species.top[1:min(5,length(species.top))]
    # Get proportion of variance
    rda1_label <- as.character(
        round(ordination$cont$importance[2,1],digits=3)*100)
    rda2_label <- as.character(
        round(ordination$cont$importance[2,2],digits=3)*100)
    rda1_label <- paste(
        'RDA1 ~ ',rda1_label, '%', sep='',collapse = '')
    rda2_label <- paste(
        'RDA2 ~ ',rda2_label, '%', sep='',collapse = '')
    arrowmultiplier = attr(ordination, 'const')
    g <- ggplot(plotdata, aes(x=RDA1, y=RDA2))
    g <- g + xlab(rda1_label) + ylab(rda2_label)
    g <- g + geom_point(aes(col=Age, size=DN, shape=Birds))
    g <- g + geom_text(aes(label=names),
                       size=2.5, position = position_nudge(y=-0.02))
    g <- g + geom_segment(
        data = species.top, aes(x=0,xend=RDA1*attr(ordination,'const'),
                                y=0,yend=RDA2*attr(ordination,'const')), 
        arrow = arrow(length = unit(0.3, 'cm')))
    g <- g + geom_text(
        data = species.top, aes(x=RDA1*arrowmultiplier,y=RDA2*arrowmultiplier,
                                label=row.names(species.top)), size = 2)
    g <- g + labs(title = x$title[1],
                  subtitle = x$title[2])
    return(g)
}

#####################################
# kmer nmds
#####################################
kmernmds <- function(k=11){
  # Run dbrda on kmer 11 o 13 distance matrices. Prints
  # latex formatted results
  dm <- load.fulldms(k)
  set.seed(1)
  js.nmds <- metaMDS(dm$js, try=99, trymax = 999, autotransform = F)
  js.e <- envfit(js.nmds, metaold)
  bc.nmds <- metaMDS(dm$bc, try=99, trymax = 999, autotransform = F)
  bc.e <- envfit(bc.nmds, metaold)
  ang.nmds <- metaMDS(dm$ang, try=99, trymax = 999, autotransform = F)
  ang.e <- envfit(ang.nmds, metaold)
  
  factors <- names(js.e$vectors$pvals)
  js.stars <- symnum(
    js.e$vectors$pvals,
    cutpoints = c(0,0.001,0.01,0.05,0.1,1),
    symbols = c('***','**','*','.',''), na=F, legend = F
  )
  bc.stars <- symnum(
    bc.e$vectors$pvals,
    cutpoints = c(0,0.001,0.01,0.05,0.1,1),
    symbols = c('***','**','*','.',''), na=F, legend = F
  )
  ang.stars <- symnum(
    ang.e$vectors$pvals,
    cutpoints = c(0,0.001,0.01,0.05,0.1,1),
    symbols = c('***','**','*','.',''), na=F, legend = F
  )
  l <- ''
  for (i in 1:length(factors)){
    l <- c(l, paste(
      factors[[i]], 
      '&', round(js.e$vectors$r[[i]], digits = 4),
      '&',round(js.e$vectors$pvals[[i]], digits= 4), '&',
      js.stars[[i]], 
      '&', round(bc.e$vectors$r[[i]], digits = 4),
      '&',round(bc.e$vectors$pvals[[i]], digits= 4), '&',
      bc.stars[[i]], 
      '&', round(ang.e$vectors$r[[i]], digits = 4),
      '&',round(ang.e$vectors$pvals[[i]], digits= 4), '&',
      ang.stars[[i]], 
      '\\\\', sep='', collapse = '') )
  }
  start_table <- c(
    '\\begin{table}[hb]',
    '\\centering',
    '\\caption{Envit results for $k$-mer based NMDS}',
    '\\label{tab:kmernmds}'
  )
  bord <- '\\begin{tabular}{| c | c | c | c | c | c | c | c | c | c |}'
  m <- c(
    '& \\multicolumn{3}{c|}{Jensen-Shannon}&\\multicolumn{3}{c|}{Bray-Curtis} 
     & \\multicolumn{3}{c|}{Angular Distance}\\\\',
    '\\hline')
  desc <- 'Descriptor & $r^2$ & p & Code & $r^2$ & p & Code& $r^2$ & p & Code\\\\'
  start_table <- c(
    start_table,
    bord,
    '\\hline',
    m,
    desc,
    '\\hline'
  )
  end_table <- c(
    '\\hline',
    '\\end{tabular}',
    '\\end{table}'
  )
  cat(start_table, sep = '\n')
  cat(l, sep = '\n')
  cat(end_table, sep = '\n')
}

#####################################
# PLots
#####################################
g_legend <- function(a.gplot){
    # Extract legend grob
    # Function copied from
    # https://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
}

multipl <- function(){
    k <- domainsrda()
    o <- orderrda()
    g <- genusrda()
    k.p <- rdaplotter(k)
    o.p <- rdaplotter(o)
    g.p <- rdaplotter(g)
    g.p <- g.p + theme(legend.position = 'bottom', legend.box = "vertical")
    k.p <- k.p + theme(legend.position = 'none')
    o.p <- o.p + theme(legend.position = 'none')
  #  l <- g_legend(g.p)
   # g.p <- g.p + theme(legend.position = 'none')
    dm = load.fulldms()
    birds <- birdageplot(dm$js)
  #  birds <- birds + theme(legend.position = 'none')
    #m <- grid.arrange(
        #arrangeGrob(k.p ,o.p,
                    #ncol=2), arrangeGrob(g.p,birds,ncol = 2))
#    m <- grid.arrange(k.p ,o.p,
 #                     g.p,birds,
  #                    ncol = 2, nrow=2,
   #                   top = 'Multiplot',
    #                  bottom = l)
    m <- fourplot(birds, k.p, o.p , g.p)
    return(m)
}

plotnmds <- function(){
  dm <- load.fulldms()
  g1 <- birdageplot(dm$js)
  g2 <- birdageplot(dm$bc)
  g3 <- birdageplot(dm$ang)
  g1 <- g1 + theme(legend.position = 'bottom', legend.box = "vertical")
  l <- g_legend(g1)
  g4 <- createKEGGtable()
  g1 <- g1 + theme(legend.position = 'none')
  g2 <- g2 + theme(legend.position = 'none')
  g3 <- g3 + theme(legend.position = 'none')
  g4 <- g4 + theme(legend.position = 'none')
  g1 <- g1 + labs(title = ' a) Jensen-Shannon',
           subtitle = 'NMDS based on K-mers')
  g2 <- g2 + labs(title = ' b) Bray-Curtis',
                  subtitle = 'NMDS based on K-mers')
  g3 <- g3 + labs(title = ' c) Angular Distance',
                  subtitle = 'NMDS based on K-mers')
  g4 <- g4 + labs(title = 'd) KEGG',
           subtitle = 'NMDS based on orthology')
  m <- grid.arrange(arrangeGrob(g1 ,g2, ncol = 2),arrangeGrob(g3 ,g4, ncol = 2),
                    top = 'Influence of birds, age and DN',
                    bottom = l)
  ggsave('nmdskmerkegg.png', plot=m, device = 'png',
         width = 14, height = 14)
}

fourplot <- function(a,b,c,d){
    # Plots domain order genus and JS based NMDS
    # in one plot
    a <- a + theme(legend.position = 'bottom', legend.box = "vertical")
    b <- b + theme(legend.position = 'none')
    c <- c + theme(legend.position = 'none')
    d <- d + theme(legend.position = 'none')
    l <- g_legend(a)
    a <- a + theme(legend.position = 'none')
    m <- grid.arrange(arrangeGrob(a ,b, ncol = 2),
                      arrangeGrob(c,d,ncol = 2),
                      top = 'Multiplot',
                      bottom = l)
    ggsave('birdkingordgen.png', plot=m, device = 'png',
           width = 14, height = 14)
}
######################################
# Convert output to Latex tables
######################################
latexanovaconverter <- function(x, n_rows,row_names, n_objs=1){
    stars <- symnum(
            x$'Pr(>F)',
            cutpoints = c(0.001,0.01,0.05,0.1,1),
            symbols = c('***','**','*','.'), na=F, legend = F
            )
    start_table <- c(
        '\\begin{table}[hb]',
        '\\centering',
        '\\caption{}',
        '\\label{}'
    )
    bord <- '\\begin{tabular}{| c ||'
    desc <- 'Descriptor'
    for (i in 1:n_objs){
        bord <- c(bord,' c | c | c |')
        desc <- c(desc, '& Variance &Pr(>F)& Code')
    }
    bord <- c(bord,'}')
    desc <- paste(c(desc, '\\\\'), sep='', collapse='')
    bord <- paste(bord, sep='',collapse = '')
    start_table <- c(
        start_table,
        bord,
        '\\hline',
        desc,
        '\\hline'
        )
    cat(start_table, sep='\n')
    tot_var = sum(x$Variance)
    for (i in 1:n_rows){
        l <- c(
            row_names[i],
            '&',round(x$Variance[i]/tot_var,digits=3),
            '&',x$'Pr(>F)'[i],
            '&',stars[i],
            '\\\\')
        cat(l,'\n')
    }
    end_table <- c(
        '\\hline',
        '\\end{tabular}',
        '\\end{table}'
    )
    cat(end_table, sep = '\n')
    
}

createKOGtable <- function(){
    # Writes code latex code to console
    # Writes anova results for rda on
    # Genus, Order and domain
    n_objs = 3
    k <- domainsrda()
    k.a <- anova(k, by='terms')
    o <- orderrda()
    o.a <- anova(o, by='terms')
    g <- genusrda()
    g.a <- anova(g, by='terms')
    k.stars <- symnum(
        k.a$'Pr(>F)',
        cutpoints = c(0,0.001,0.01,0.05,0.1,1),
        symbols = c('***','**','*','.',''), na=F, legend = F
    )
    o.stars <- symnum(
        o.a$'Pr(>F)',
        cutpoints = c(0,0.001,0.01,0.05,0.1,1),
        symbols = c('***','**','*','.',''), na=F, legend = F
    )
    g.stars <- symnum(
        g.a$'Pr(>F)',
        cutpoints = c(0,0.001,0.01,0.05,0.1,1),
        symbols = c('***','**','*','.',''), na=F, legend = F
    )
    start_table <- c(
        '\\begin{table}[hb]',
        '\\centering',
        '\\caption{Anova results}',
        '\\label{tab:anovakog}'
    )
    bord <- '\\begin{tabular}{| c ||'
    m <- c(
    '&\\multicolumn{3}{c|}{Domain} & \\multicolumn{3}{c|}{Order}
    &\\multicolumn{3}{c|}{Genus}\\\\',
    '\\hline')
    desc <- 'Descriptor'
    for (i in 1:n_objs){
        bord <- c(bord,' c | c | c |')
        desc <- c(desc, '& Var &Pr(>F)& Code')
    }
    bord <- c(bord,'}')
    desc <- paste(c(desc, '\\\\'), sep='', collapse='')
    bord <- paste(bord, sep='',collapse = '')
    start_table <- c(
        start_table,
        bord,
        '\\hline',
        m,
        desc,
        '\\hline'
    )
    cat(start_table, sep='\n')
    k.tot_var = sum(k.a$Variance)
    o.tot_var = sum(o.a$Variance)
    g.tot_var = sum(g.a$Variance)
    row_names <- row.names(k.a)
    for (i in 1:10){
        l <- c(
            row_names[i],
            '&',round(k.a$Variance[i]/k.tot_var,digits=3),
            '&',k.a$'Pr(>F)'[i],
            '&',k.stars[i],
            '&',round(o.a$Variance[i]/o.tot_var,digits=3),
            '&',o.a$'Pr(>F)'[i],
            '&',o.stars[i],
            '&',round(g.a$Variance[i]/g.tot_var,digits=3),
            '&',g.a$'Pr(>F)'[i],
            '&',g.stars[i],
            '\\\\')
        cat(l,'\n')
    }
    end_table <- c(
        '\\hline',
        '\\end{tabular}',
        '\\end{table}'
    )
    cat(end_table, sep = '\n')
    df <- data.frame(row.names = row.names(k.a)[-10])
    df$k <- k.a$Variance[-10]
    df$o <- o.a$Variance[-10]
    df$g <- g.a$Variance[-10]
    return(df)
}

createKEGGtable <- function(){
  # Writes latex code to console
  # Writes envfit results for nmds on
  # KEGG table and plots birds, DN and age
  set.seed(6)
  kegg <- loadKegg()
  kegg.nmds <- metaMDS(kegg, try=99, trymax = 999)
  kegg.e <- envfit(kegg.nmds, metaold)
  factors <- names(kegg.e$vectors$pvals)
  stars <- symnum(
    kegg.e$vectors$pvals,
    cutpoints = c(0,0.001,0.01,0.05,0.1,1),
    symbols = c('***','**','*','.',''), na=F, legend = F
  )
  l <- ''
  for (i in 1:length(factors)){
    l <- c(l, paste(
      factors[[i]], '&', round(kegg.e$vectors$r[[i]], digits = 4),
      '&',round(kegg.e$vectors$pvals[[i]], digits= 4), '&',
      stars[[i]], '\\\\', sep='', collapse = '') )
  }
  start_table <- c(
    '\\begin{table}[hb]',
    '\\centering',
    '\\caption{Envit results for KEGG orthology fitted to environmental
    variables}',
    '\\label{tab:keggnmds}'
  )
  bord <- '\\begin{tabular}{| c | c | c | c |}'
  desc <- 'Descriptor & $r^2$ & p & Code \\\\'
  start_table <- c(
    start_table,
    bord,
    '\\hline',
    desc,
    '\\hline'
  )
  end_table <- c(
    '\\hline',
    '\\end{tabular}',
    '\\end{table}'
  )
  cat(start_table, sep = '\n')
  cat(l, sep = '\n')
  cat(end_table, sep = '\n')
  tmp <- metadata.r
  tmp$DN <- metadata$DN
  tmp$X <- kegg.nmds$points[,1]
  tmp$Y <- kegg.nmds$points[,2]
  tmp$Sample <- samplenames
  stress <- paste('Stress: ',
                  round(kegg.nmds$stress, digits = 4), sep='',collapse = '')
  g <- ggplot(tmp, aes(x=X, y=Y))
  g <- g + geom_point(aes(col=Age, size=DN, shape=Birds))
  g <- g + geom_text(aes(label=Sample), size=2.5
                     ,position = position_nudge(y=-0.01))
  g <- g + labs(title = ' b) Influence of birds and age on DN',
                subtitle = 'NMDS based on KEGG orthology groups',
                caption = stress)
  g <- g + xlab('NMDS 1') + ylab('NMDS 2')
  return(g)
}

#######################################
# Alpha diversity
######################################

plotalpha <- function(alph, feature) {
  g <- ggplot(data = metaalpha, aes(y=.data[[alph]], x=.data[[feature]]))
  g <- g + geom_point()
  g <- g + geom_smooth(method = 'loess')
  return(g)
}

kogdiversity <- function(){
  # Create dataframe containing evenness/diversity measures for the various
  # samples
  k <- loadTaxa('d', tform = F)
  o <- loadTaxa('o', tform = F)
  g <- loadTaxa('g', tform = F)
  k <- clean(k)
  o <- clean(o)
  g <- clean(g)
  k.d <- diversity(k)
  o.d <- diversity(o)
  g.d <- diversity(g)
  k.e <- k.d/log(specnumber(k))
  o.e <- o.d/log(specnumber(o))
  g.e <- g.d/log(specnumber(g))
  df <- data.frame(row.names = samplenames)
  df[, c('Domain', 'Domain Pielou', 'OrderD','Order Pielou','GenusD','Genus Pielou')] <- c(
    k.d, k.e, o.d, o.e, g.d, g.e)
  return(df)
}

plotvspielou <- function(measure) {
  d <- kogdiversity()
  tmpmeta <- cbind(metaalpha, d)
  g <- ggplot(data = tmpmeta, aes(GenusE))
  g <- g + geom_point(aes(y=.data[[measure]]))
  g <- g + geom_smooth(aes(y=.data[[measure]]), method = 'loess')
  return(g)
}

pairsplotpialou <- function() {
  d <- kogdiversity()
  tmpmeta <- cbind(d,metaalpha)
  pairs.panels(tmpmeta[, c(1,3,5,18,19,20,22,23)])
}

hellbray <- function(x) {
  d <- decostand(x, method = 'hellinger')
  d <- vegdist(d)
  return(d)
}

mantelproctabel <- function(kmer=11) {
  # Creates latex table with mantel procrustes
  # comparisons of distance measures base on taxonomy
  # vs kmers
  dm <- load.fulldms(kmer)
  set.seed(5)
  k <- vegdist(clean(loadTaxa('d')))
  o <- vegdist(clean(loadTaxa('o')))
  g <- vegdist(clean(loadTaxa('g')))
  kegg <- vegdist(loadKegg())
  k.m <- mantel(k, dm$js, permutations = 9999)
  k.p <- protest(k, dm$js, permutations = how(nperm=9999))
  o.m <- mantel(o, dm$js, permutations = 9999)
  o.p <- protest(o, dm$js, permutations = how(nperm=9999))
  g.m <- mantel(g, dm$js, permutations = 9999)
  g.p <- protest(g, dm$js, permutations = how(nperm=9999))
  kegg.m <- mantel(kegg, dm$js, permutations = 9999)
  kegg.p <- protest(kegg, dm$js, permutations = how(nperm=9999))
  
  extractmantel <- function(x) {
    values <- c(round(x$statistic, digits = 4),x$signif)
    return(values)
  }
  extractproc <- function(x) {
    values <- c(round(sqrt(1-x$ss), digits=4), x$signif)
    return(values)
  }
  k.r <- c('Domain',extractmantel(k.m), extractproc(k.p))
  o.r <- c('Order',extractmantel(o.m), extractproc(o.p))
  g.r <- c('Genus',extractmantel(g.m), extractproc(g.p))
  kegg.r <- c('KEGG',extractmantel(kegg.m), extractproc(kegg.p))
  
  start_table <- c(
    '\\begin{table}[hb]',
    '\\centering',
    '\\caption{Results of mantel test and Procrustes test comparing
    dissimilarities between communities based on alignment
    of reads to taxonomic ranks and dissimilarities based on 
    Jensen-Shannon distance for 11-mer counts.}',
    '\\label{tab:taxmantelpro}'
  )
  bord <- '\\begin{tabular}{| c || c | c || c | c |}'
  desc <- 'Rank & Mantel r &  Significance & Procrustes corr. & Significance\\\\'
  start_table <- c(
    start_table,
    bord,
    '\\hline',
    desc,
    '\\hline'
  )
  cat(start_table, sep='\n')
  cat('\\hline', '\\multicolumn{5}{|c|}{Jensen-Shannon}\\\\', sep='\n')
  cat(paste(k.r, sep = '', collapse='&'), '\\\\','\n')
  cat(paste(o.r, sep = '', collapse='&'), '\\\\','\n')
  cat(paste(g.r, sep = '', collapse='&'), '\\\\','\n')
  cat(paste(kegg.r, sep = '', collapse='&'), '\\\\','\n')
  cat('\\hline', '\\multicolumn{5}{|c|}{Bray-Curtis}\\\\', sep='\n')
  
  #Bray Curtis calculations
  k.m <- mantel(k, dm$bc, permutations = 9999)
  k.p <- protest(k, dm$bc, permutations = how(nperm=9999))
  o.m <- mantel(o, dm$bc, permutations = 9999)
  o.p <- protest(o, dm$bc, permutations = how(nperm=9999))
  g.m <- mantel(g, dm$bc, permutations = 9999)
  g.p <- protest(g, dm$bc, permutations = how(nperm=9999))
  kegg.m <- mantel(kegg, dm$bc, permutations = 9999)
  kegg.p <- protest(kegg, dm$bc, permutations = how(nperm=9999))
  k.r <- c('Domain',extractmantel(k.m), extractproc(k.p))
  o.r <- c('Order',extractmantel(o.m), extractproc(o.p))
  g.r <- c('Genus',extractmantel(g.m), extractproc(g.p))
  kegg.r <- c('KEGG',extractmantel(kegg.m), extractproc(kegg.p))
  
  cat(paste(k.r, sep = '', collapse='&'), '\\\\','\n')
  cat(paste(o.r, sep = '', collapse='&'), '\\\\','\n')
  cat(paste(g.r, sep = '', collapse='&'), '\\\\','\n')
  cat(paste(kegg.r, sep = '', collapse='&'), '\\\\','\n')
  
  end_table <- c(
    '\\hline',
    '\\end{tabular}',
    '\\end{table}'
  )
  cat(end_table, sep = '\n')  
  
}

origdbrda <- function(x) {
  # Uses same seed as script used in report. (file redundancy.r)
  # Can thus be used to create tables and plots
  set.seed(1)
  rd <- dbrda(x ~. , data=metaold)
  rd$anova <- anova(rd ,permutations = how(nperm = 10000))
  rd$anovat <- anova(rd, by='terms',permutations = how(nperm = 10000))
  return(rd)
}

origdbrdatable <- function(k=11){
  # Create latex table for dbRDA based on k-mer 
  # length 11 or 13
  dm <- load.fulldms(k)
  bc <- origdbrda(dm$bc)
  js <- origdbrda(dm$js)
  
  n_objs = 2
  bc.stars <- symnum(
    bc$anovat$'Pr(>F)',
    cutpoints = c(0,0.001,0.01,0.05,0.1,1),
    symbols = c('***','**','*','.',''), na=F, legend = F
  )
  js.stars <- symnum(
    js$anovat$'Pr(>F)',
    cutpoints = c(0,0.001,0.01,0.05,0.1,1),
    symbols = c('***','**','*','.',''), na=F, legend = F
  )
  start_table <- c(
    '\\begin{table}[hb]',
    '\\centering',
    '\\caption{Anova results}',
    '\\label{tab:anovadbrda}'
  )
  bord <- '\\begin{tabular}{| c ||'
  m <- c(
    '&\\multicolumn{3}{c|}{Bray-Curtis} & \\multicolumn{3}{c|}{Jensen-Shannon}\\\\',
    '\\hline')
  desc <- 'Descriptor'
  for (i in 1:n_objs){
    bord <- c(bord,' c | c | c |')
    desc <- c(desc, '& SS &$Pr(>F)$& Code')
  }
  bord <- c(bord,'}')
  desc <- paste(c(desc, '\\\\'), sep='', collapse='')
  bord <- paste(bord, sep='',collapse = '')
  start_table <- c(
    start_table,
    bord,
    '\\hline',
    m,
    desc,
    '\\hline'
  )
  cat(start_table, sep='\n')
  bc.tot_ss = sum(bc$anovat$SumOfSqs)
  js.tot_ss = sum(js$anovat$SumOfSqs)
  row_names <- row.names(bc$anovat)
  for (i in 1:10){
    l <- c(
      row_names[i],
      '&',round(bc$anovat$SumOfSqs[i]/bc.tot_ss,digits=3),
      '&',round(bc$anovat$'Pr(>F)'[i], digits=3),
      '&',bc.stars[i],
      '&',round(js$anovat$SumOfSqs[i]/js.tot_ss,digits=3),
      '&',round(js$anovat$'Pr(>F)'[i], digits = 3),
      '&',js.stars[i],
      '\\\\')
    cat(l,'\n')
  }
  end_table <- c(
    '\\hline',
    '\\end{tabular}',
    '\\end{table}'
  )
  cat(end_table, sep = '\n')
  
}

rankings <- function(k=11){
  # Collects scores for all methods used
  # 
  dm <- load.fulldms(k)
  js <- origdbrda(dm$js)
  bc <- origdbrda(dm$bc)
  ang <- origdbrda(dm$ang)
  df <- data.frame(row.names = row.names(js$anovat)[-10])
  df$JS <- js$anovat$SumOfSqs[-10]
  df$BC <- bc$anovat$SumOfSqs[-10]
  df$AD <- ang$anovat$SumOfSqs[-10]
  # Get domain, order and genus RDA results, using same seeds
  kog <- createKOGtable()
  df <- cbind.data.frame(df, kog)
  # Get kegg results with same seed
  set.seed(6)
  kegg <- loadKegg()
  kegg.nmds <- metaMDS(kegg, try=99, trymax = 999)
  kegg.e <- envfit(kegg.nmds, metaold)
  
  set.seed(1)
  js.nmds <- metaMDS(dm$js, try=99, trymax = 999, autotransform = F)
  js.e <- envfit(js.nmds, metaold)
  bc.nmds <- metaMDS(dm$bc, try=99, trymax = 999, autotransform = F)
  bc.e <- envfit(bc.nmds, metaold)
  ang.nmds <- metaMDS(dm$ang, try=99, trymax = 999, autotransform = F)
  ang.e <- envfit(ang.nmds, metaold)
  df$jsnmds <- js.e$vectors$r
  df$bcnmds <- bc.e$vectors$r
  df$angnmds <- ang.e$vectors$r
  df$KEGG <- kegg.e$vectors$r
  
  colnames(df) <- c('JS', 'BC', 'AD', 'Domain', 'Order', 'Genus', 'JS NMDS',
                    'BC NMDS', 'AD NMDS', 'KEGG NMDS')
  rf <- df
  for (i in 1:length(colnames(df))) {
    rf[,i] <- rank(-df[,i])
  }
  start_table <- c(
    '\\begin{table}[hb]',
    '\\centering',
    '\\caption{Ranking of importance of environental variables using
    RDA and NMDS based on kmer counts and alignment based methods}',
    '\\label{tab:rankings}'
  )
  bord <- '\\begin{tabular}{| c | c | c | c | c | c | c | c | c | c | c |}'
  desc <- 'Method & \\multicolumn{6}{c|}{RDA} & \\multicolumn{4}{c|}{NMDS} \\\\'
  desc2 <- 'Source & \\multicolumn{3}{c|}{Kmer} & \\multicolumn{3}{c|}{Taxonomy}
  & \\multicolumn{3}{c|}{Kmer} & Orthology \\\\'
  start_table <- c(
    start_table,
    bord,
    '\\hline',
    desc,
    '\\hline',
    desc2, 
    '\\hline'
  )
  end_table <- c(
    '\\hline',
    '\\end{tabular}',
    '\\end{table}'
  )
  l <- ''
  for (i in 1:length(row.names(df))){
    l <- c(l, paste(row.names(df)[i],'&', rf[i,1],'&',
                    rf[i,2],'&', rf[i,3], '&', rf[i,4], '&', 
                    rf[i,5],'&', rf[i,6], '&', rf[i,7], '&',
                    rf[i,8],'&', rf[i,9], '&', rf[i,10], '\\\\'))
  }
  cat(start_table, sep = '\n')
  cat(l, sep = '\n')
  cat(end_table, sep = '\n')
  return(df)
}

unknown <- function(x) {
  # Find percentage of unknown or unclassified taxa
  k <- loadTaxa('d', tform = F)
  o <- loadTaxa('o', tform = F)
  g <- loadTaxa('g', tform = F)
  colnames(k) <- lapply(
    colnames(k), function(x){gsub('^.?_', '',x)})
  colnames(o) <- lapply(
    colnames(o), function(x){gsub('^.?_', '',x)})
  colnames(g) <- lapply(
    colnames(g), function(x){gsub('^.?_', '',x)})
  df <- data.frame(row.names = samplenames)
  df$'Domain total' <- rowSums(k)
  df$'Domain classified' <- rowSums(clean(k))
  df$'Domain ratio' <- df$'Domain classified'/df$'Domain total'
  df$'Order total' <- rowSums(o)
  df$'Order classified' <- rowSums(clean(o))
  df$'Order ratio' <- df$'Order classified'/df$'Order total'
  df$'Genus total' <- rowSums(g)
  df$'Genus classified' <- rowSums(clean(g))
  df$'Genus ratio' <- df$'Genus classified'/df$'Genus total'
  return(df)
}

compareDiversities <- function(){
  # Calculates correlation of various diversity measures
  # for 
  set.seed(1)
  toCompare <- list('bray', 'gower', 'jaccard', 'kulczynski',
                    'raup', 'cao')
  kmerDist <- load.fulldms()$js
  g <- loadTaxa('g')
  g <- clean(g)
  g <- g[, colSums(g)>0]
  d <- loadTaxa('d')
  d <- clean(d)
  d <- d[, colSums(d)>0]
  dms <- list()
  GenBray <- vegdist(g, method = 'bray')
  DomBray <- vegdist(d, method = 'bray')
  for (i in 2:length(toCompare)){
    mantelResults <- mantel(GenBray, vegdist(g, method = toCompare[i]))
    print(paste('Bray Curtis ', toCompare[i], mantelResults$statistic, sep = ' '))
  }
  mantelResults <- mantel(GenBray, DomBray)
  print(paste('Genus Bray Curtis ', 'Domain Bray Currtis', mantelResults$statistic, sep = ' '))
  return(dms)
}
