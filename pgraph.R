library(ggplot2)
setwd ("E:/cross_ref/jgrot/masterOutput/output")

pgd = read.csv("E:/cross_ref/jgrot/masterOutput/output/pgd.txt", header=TRUE, sep="\t")

#read input file
genes = read.csv("E:/cross_ref/jgrot/gene.txt", header = FALSE, sep = ",")
for (gene in genes){
  filename = paste(gene, '.txt', sep="")
  geneFrame = read.csv (filename, header = TRUE, sep = "\t")
  
  studies = unique(geneFrame$study_name)
  
  num_studies<-length(studies)
  minimums<-numeric(num_studies)
  
  for (i in 1:num_studies){
    tempMin= subset(geneFrame, study_name==studies[i])
    minimums[i] = min(tempMin$p_value)
  }
  
  gene_summary = data.frame(studies, minimums)
  gene_summary<-gene_summary[order(gene_summary$minimums, decreasing=FALSE), ] 
  gene_summary$studies<-factor(gene_summary$studies, levels=gene_summary$studies)
  
  graph_title = paste(toupper(gene), ' P-Val Graph', sep="")
  pdf(graph_title, width= 20)
  
  theme_set(theme_bw())
  print(ggplot(gene_summary, aes(x=studies, y=minimums))+
    geom_point(size=3) +
    geom_segment(aes(x=studies,
                     xend=studies,
                     y=0,
                     yend=minimums))+
    labs(title="pval vs study") + xlab("study") + ylab("pval") +  theme(axis.text.x = element_text(angle=90, vjust=0.6)))
  dev.off()
}

