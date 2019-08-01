library(ggplot2)
setwd ("E:/BRIAN/Cross_Ref_Pipeline/masterOutput/output")

#read input file
genes = scan("E:/BRIAN/Cross_Ref_Pipeline/validGeneFile.txt", what="", sep="\n")

#loop through gene file to produce graph for each gene
for (gene in genes){
  filename = paste(gene, '.txt', sep="")
  geneFrame = read.csv (filename, header = TRUE, sep = "\t")
  
  #create variables
  studies = unique(geneFrame$study_name)
  num_studies<-length(studies)
  minimums<-numeric(num_studies)
  count=numeric(num_studies)
  
  if (num_studies == 0){
    next
  }
  
  #loop to fill minium array with lowest pval of each study
  for (i in 1:num_studies){
    tempMin= subset(geneFrame, study_name==studies[i])
    minimums[i] = tempMin[1, "p_value"]
    count[i] = nrow(tempMin)
  }
  
  #pval -> -log(pval)
  logMin=sapply(minimums, log)
  logMin = -1 * logMin
  
  #create/sort table the graph is formed off
  gene_summary = data.frame(studies, logMin)
  gene_summary<-gene_summary[order(gene_summary$logMin, decreasing=TRUE), ] 
  gene_summary$studies<-factor(gene_summary$studies, levels=gene_summary$studies)
  
  #graph setup
  graph_title = paste(toupper(gene), ' P-Val Graph', sep="")
  pdf(graph_title, width= 20)
  
  theme_set(theme_bw())
  print(ggplot(gene_summary, aes(x=studies, y=logMin))+
    geom_point(size=3) +
    geom_segment(aes(x=studies,
                     xend=studies,
                     y=0,
                     yend=logMin))+
    labs(title="P-Value vs Study", subtitle=toupper(gene)) + xlab("Study") + ylab("-log(pval)") +  theme(axis.text.x = element_text(angle=90, vjust=0.6)))
  dev.off()
}

