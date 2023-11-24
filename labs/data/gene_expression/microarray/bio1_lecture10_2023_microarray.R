#Basic Microarray Expression analysis
#Dr. T. Ian Simpson, School of Informatics, University of Edinburgh
#ian.simpson@ed.ac.uk
#Bioinformatics 1 - 19th October 2012
#Temporal profiling of Drosophila Peripheral Nervous System Development

# Setting up packages
setRepositories(graphics = F, c(1,2,3,4,5))
install.packages(c('affy','drosophila2cdf','limma','topGO','drosophila2.db'))

# PART ONE - loading the array data ---------------------------------------

#load libraries
library(affy);

#set working directory
setwd('./')

#load Affymetrix CEL files
data <- ReadAffy(); #loading CEL files for wild-type and mutant samples GFP +ve and -ve (16 chips)

#find out a bit about the object
class(data);

#what's an AffyBatch
?AffyBatch

#some access methods
probes(data,which='pm'); #probe perfect match expression
boxplot(data); #box-whisker plots each chip
image(data[,1]); #see a chip image plot

#examine the slots
data@cdfName; #get the chip platform name
data@phenoData; #get the experimental data

#what's an AnnotatedDataFrame
?AnnotatedDataFrame;

#lets annotate on the experimental information
pData(data)$condition <- as.factor(c(rep('wt',8),rep('mut',8)));
pData(data)$gfp <- as.factor(rep(c('p','n'),8));
pData(data)$sample <- as.factor(sort(rep(seq(1,8),2)));

#view the updated phenotype data
pData(data);


# PART TWO - checking array quality and extracting expression valu --------

#library(arrayQualityMetrics);
#arrayQualityMetrics(data,outdir='./AQM_pre/'); #DO NOT RUN THIS COMMAND

#now lets normalise the data
norm <- rma(data); # RMA - the Robust Multi-Chip Average - returns expression in log2

#and extract the expression matrix
exprs <- exprs(norm);

#PART THREE - Principal Component Analysis (PCA)
pca <- prcomp(t(exprs)); # Note, we are doing PCA by chip i.e. trying to see where the variation is between the chips

#look at where the variation is across the Principal Components
summary(pca);

#extract the principal components
pcs <- data.frame(pca$x);

#add in the annotation data by merging (just for the first 2 principal components)
pcs <- merge(pcs[,1:2],pData(data),by=0);

#plot the first two principal components
library(lattice);

#plot by condition using ggplot2 colour mut = red, wt = blue
library(ggplot2);
ggplot(pcs,aes(PC1,PC2,colour=condition))+geom_point()+theme_bw()+scale_colour_manual(values=c('red','blue'));

#plot by gfp status using ggplot2
ggplot(pcs,aes(PC1,PC2,colour=gfp))+geom_point()+theme_bw();

#plot by sample using ggplot2
ggplot(pcs,aes(PC1,PC2,colour=sample))+geom_point()+theme_bw();

#PART FOUR - Differential Expression Calling
#load limma a package for calling differntially expressed genes
library(limma);
#load the Bioconductor annotation file for the Drosophila_2 chips
library(drosophila2.db);

#convert the expresison matrix from log2 scale to base 10
exprs <- data.frame(2^exprs); # we're just extracting the expression slot from the ExpressionSet object, convert to base10

#add group to the annotation data
d <- data.frame(pData(data));
d$group<-interaction(d$gfp,d$sample,sep='');

#rename columns to make easier to identify
colnames(exprs) <- d$group[match(row.names(d),colnames(exprs))];

#these variables can be exported to the environment so that they can be referred to directly
attach(exprs);

#create an enrichment matrix by dividing GFP+ signal by GFP- signal
enrich <- data.frame(s1=p1/n1,s2=p2/n2,s3=p3/n3,s4=p4/n4,s5=p5/n5,s6=p6/n6,s7=p7/n7,s8=p8/n8,row.names=row.names(exprs));

#set up the 'levels' matrix (this shows the groupings for contrasts)
enrich_design <- table(data.frame(row.names=colnames(enrich),unique(d[,1:2])));
#look at it
enrich_design;

#set up the contrast matrix (what we are comparing) - note 'wt-mut'
cont.matrix <- makeContrasts(wt-mut,levels=enrich_design);
#look at it
cont.matrix;

#fit the linear model
fit <- lmFit(enrich, enrich_design);

#calculate the fit coefficients and error for the comparison(s)
fit2 <- contrasts.fit(fit, cont.matrix);

#calculate t-statistics and log-odds of differential expression using emperical Bayes approach
fit2 <- eBayes(fit2);

#have a quick look at the result
topTable(fit2,n=10,sort.by='logFC',p=0.01,adjust='fdr');

#get the differentially expressed genes fdr=1%
diffexp <- topTable(fit2,adjust="fdr",sort.by='logFC',n=1000000000,p=0.01);

#reformat to make the results simpler to view
diffexp_top <- data.frame(diffexp[,c(1,2)],fdr=diffexp[,6],fc=diffexp[,2]);
#bring in the gene names
diffexp_final <- merge(diffexp_top,toTable(drosophila2SYMBOL[rownames(diffexp_top)]),by.x=0,by.y='probe_id',sort=F,all.x=T);
#look at the top10
diffexp_final[1:10,];

#OPTIONAL

#PART FIVE - Functional Annotation
#Using  package to look for functional enrichment in gene lists
library(topGO);

geneList <- fit2$p.value[,1];
names(geneList) <- as.factor(fit2$genes[,1]);

names(geneList) <- as.factor(rownames(fit2));

topDiffGenes <- function(x){
  return(x<=0.01)
}

#now set up the Drosophila GO
library(drosophila2.db);

#create the GOData object that contains all of the gene association information for background and your own gene list
sampleGOdata <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = topDiffGenes,nodeSize = 10,  annot = annFUN.db, affyLib = "drosophila2.db")

#perform a Fisher test to look for statistically significant enrichment of terms
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

#use a Kolmogorov–Smirnov test instead
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")

#apply KS test, with more conservative outcome 'elim'
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

#display all the tests together
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

#PART SIX - Clustering
library(cluster)

#check to see that gthe chips segregate by genotype
cl <- kmeans(t(enrich),2);

#cluster the top 100 differentially expressed genes
cl_group <- diffexp_final[1:100,]

colnames(cl_group)[1] <- 'ID'

#get the ids
ids <- cl_group$ID;

#extract only those expression values
cl_enrich <- enrich[ids,];

#use hierarchical clustering to partition the differentially expressed genes
hc <- hclust(dist(cl_enrich),method='complete');

#show the expression groupings in a heatmap
heatmap(as.matrix(cl_enrich),Rowv=hc$order,)
