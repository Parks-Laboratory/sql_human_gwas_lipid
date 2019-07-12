library(ggplot2)
setwd ("E:/cross_ref/jgrot/masterOutput/output")

#create empty df that will hold all indiv gene dataframes
master = data.frame()

#read input file
genes = read.csv("E:/cross_ref/jgrot/gene.txt", header = FALSE, sep = ",")

#loop to fill master dataframe
for (gene in genes){
  filename = paste(gene, '.txt', sep="")
  geneFrame = read.csv (filename, header = TRUE, sep = "\t")
  master = rbind(master, geneFrame)
}

#variable set up
master_genes<-unique(master$gene_name)
num_genes<-length(master_genes)
num_studies<-numeric(num_genes)

#loop fill num_studies with number of hits per gene
for (i in 1:num_genes){
  gen=subset(master, gene_name==master_genes[i])
  num_studies[i]=length(unique(gen$study_name))
}

#create df that graph will be built with
gene_summary = data.frame(master_genes, num_studies)

#sort
gene_summary<-gene_summary[order(gene_summary$num_studies, decreasing=TRUE), ] 
gene_summary$master_genes<-factor(gene_summary$master_genes, levels=gene_summary$master_genes)

pdf("Summary Plot.pdf") 

#lollipop ggplot 
theme_set(theme_bw())
ggplot(gene_summary, aes(x=master_genes, y=num_studies))+
  geom_point(size=3) +
  geom_segment(aes(x=master_genes,
                   xend=master_genes,
                   y=0,
                   yend=num_studies))+
  labs(title="Gene.txt Summary") + xlab("Genes") + ylab("GWAS Hits")
  
dev.off() 
