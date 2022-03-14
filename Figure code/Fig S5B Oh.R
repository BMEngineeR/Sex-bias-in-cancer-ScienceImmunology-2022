rm( list=ls() )
library(dplyr)
library(cowplot)
library(Seurat)
library(parallel)
library(tools)
library(ggplot2)

# load data

load( file="pooled_CD8.RData" )

# load androgen and interferon gene sets

androgen <- read.table( "HALLMARK_ANDROGEN_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
androgen <- toTitleCase(tolower(androgen))
interferon <- read.table( "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",
  skip=2, header=FALSE, stringsAsFactors=FALSE )[,1]
interferon <- toTitleCase(tolower(interferon))

# annotate sex

sample.ID <- my.integrated@meta.data$Sample_ID

my.integrated@meta.data$sex <- "Male"
my.integrated@meta.data$sex[ 
  sample.ID == "Anti-PD-L1 C tumor CD8" | 
  sample.ID == "Chemo tumor CD8" ] <- "Female"

# extract only anti-PDL1

my.integrated.org <- my.integrated
my.integrated <- subset( my.integrated, 
  cells=grep( "Anti-PD-L1", my.integrated@meta.data$Sample_ID ) )

# extract only the cells corresponding to male

my.integrated.male <- my.integrated
my.integrated.male <- subset( my.integrated.male, 
  cells=which( my.integrated.male@meta.data$sex == "Male" ) )

# correlation analysis

raw.data <- as.matrix(GetAssayData(my.integrated.male, slot = "data"))
loc.androgen <- match( toupper(androgen), rownames(raw.data) )
data.androgen <- unlist(raw.data[ loc.androgen[ !is.na(loc.androgen) ], ])
loc.interferon <- match( toupper(interferon), rownames(raw.data) )
data.interferon <- unlist(raw.data[ loc.interferon[ !is.na(loc.interferon) ], ])

cor.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
gene.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
for ( i in 1:nrow(data.androgen) ) {
  for ( j in 1:nrow(data.interferon) ) {
    x <- data.androgen[i,]
    y <- data.interferon[j,]
    x1 <- x[ -which( x == 0 & y == 0 ) ]
    y1 <- y[ -which( x == 0 & y == 0 ) ]
    cor.mat[i,j] <- cor( x1, y1 )
    gene.mat[i,j] <- paste( rownames(data.androgen)[i], rownames(data.interferon)[j] )
  }
}
boxplot(as.vector(cor.mat))
summary(as.vector(cor.mat))

data.male <- data.frame( as.vector(gene.mat), as.vector(cor.mat) )
data.male <- data.male[ order(data.male[,2]), ]

cors <- list()
cors[["Male"]] <- as.vector(cor.mat)

# extract only the cells corresponding to female

my.integrated.female <- my.integrated
my.integrated.female <- subset( my.integrated.female, 
  cells=which( my.integrated.female@meta.data$sex == "Female" ) )

# correlation analysis

raw.data <- as.matrix(GetAssayData(my.integrated.female, slot = "data"))
loc.androgen <- match( toupper(androgen), rownames(raw.data) )
data.androgen <- unlist(raw.data[ loc.androgen[ !is.na(loc.androgen) ], ])
loc.interferon <- match( toupper(interferon), rownames(raw.data) )
data.interferon <- unlist(raw.data[ loc.interferon[ !is.na(loc.interferon) ], ])

cor.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
gene.mat <- matrix( NA, nrow(data.androgen), nrow(data.interferon) )
for ( i in 1:nrow(data.androgen) ) {
  for ( j in 1:nrow(data.interferon) ) {
    x <- data.androgen[i,]
    y <- data.interferon[j,]
    x1 <- x[ -which( x == 0 & y == 0 ) ]
    y1 <- y[ -which( x == 0 & y == 0 ) ]
    cor.mat[i,j] <- cor( x1, y1 )
    gene.mat[i,j] <- paste( rownames(data.androgen)[i], rownames(data.interferon)[j] )
  }
}
boxplot(as.vector(cor.mat))
summary(as.vector(cor.mat))

data.female <- data.frame( as.vector(gene.mat), as.vector(cor.mat) )
data.female <- data.female[ order(data.female[,2]), ]

cors[["Female"]] <- as.vector(cor.mat)

# boxplot

cors[["Male"]][ cors[["Male"]] > 0.5 ] <- NA
cors[["Female"]][ cors[["Female"]] > 0.5 ] <- NA

pdata <- data.frame(
  c( rep( "Male", length(cors[["Male"]]) ), rep( "Female", length(cors[["Female"]]) ) ),
  c( cors[["Male"]], cors[["Female"]] ) 
)
colnames(pdata) <- c( "Sex", "Corr" )

p <- ggplot( pdata, aes( x=Sex, y=Corr, fill=Sex ) ) + 
  geom_boxplot() + ylab( "Correlation coefficient" ) + ylim( -1, 1 ) +
  geom_hline( yintercept=0, linetype="dashed", color="red", size=2 ) +
  theme(axis.text.y=element_text(size=14),axis.text.x=element_text(size=14),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
        legend.title=element_text(size=14),legend.text=element_text(size=14))
print(p)
