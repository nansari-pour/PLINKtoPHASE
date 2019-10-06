setwd("/Users/napour/Documents/Personal/CBRC/raw_data/cbrc_original_dataset/42snps/pedfiles")
files=list.files("/Users/napour/Documents/Personal/CBRC/raw_data/cbrc_original_dataset/42snps/pedfiles",pattern = ".ped")
#######################################
chrs=c(2,7,8,10)
genes=c("UGT1A","CYP3A5","NAT2","CYP2C19-CYP2C9")
#######################################
######### CREATE HAPLO-TABLE from OUT files (not freqs files)
for (i in 1:length(genes)){
  gene_table=data.frame()
  for (j in 1:length(files)){
    name=gsub(".ped","",files[j])
    print(name)
    df=data.frame(rl=readLines(paste(genes[i],name,"out",sep="_")))
    start=NULL
    end=NULL
    for (k in 1:nrow(df)){
      if (df$rl[k]=="BEGIN LIST_SUMMARY"){
        start=k+1}
      if (df$rl[k]=="END LIST_SUMMARY"){
        end=k-1
      }
    }
    haplo=NULL
    freq=NULL
    hap=data.frame(df[start:end,],stringsAsFactors = F)
    for (k in 1:nrow(hap)){
    u=unlist(strsplit(as.character(hap[k,]),split=" "))
    u=u[u!=""]
    haplo=append(haplo,u[2])
    freq=append(freq,as.numeric(u[3]))
    }
    haplo_table=data.frame(Group=name,haplotype=haplo,frequency=freq)
    #sample_size=size=sum(haplo_table$frequency)
    gene_table=rbind(gene_table,haplo_table)
    print(paste("gene_table nrow=", nrow(gene_table)))
  }
  dc=dcast(gene_table,haplotype~Group,value.var="frequency")
  print(dim(dc))
  dc[is.na(dc)]=0
  pivot_table=data.frame(Gene=genes[i],dc[,c("haplotype","w_british","uk_british","GBR","bangladeshi","BEB","b_african","African_in-house_all","African_1000G_all")])
  for (k in 3:ncol(pivot_table)){
    pivot_table[,k]=round(pivot_table[,k]/sum(pivot_table[,k]),digits = 3)
  }
  write.table(format(pivot_table,digits=3),paste0(genes[i],"_haplotype_table.txt"),row.names = F,quote=F)
  print(paste("Gene",genes[i],"done"))
}
