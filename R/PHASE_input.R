setwd("/Users/napour/Documents/Personal/CBRC/raw_data/cbrc_original_dataset/42snps/pedfiles")
files=list.files("/Users/napour/Documents/Personal/CBRC/raw_data/cbrc_original_dataset/42snps/pedfiles",pattern = ".ped")
#######################################
chrs=c(2,7,8,10)
genes=c("UGT1A","CYP3A5","NAT2","CYP2C19-CYP2C9")
#######################################
for (i in 1:length(files)){
  name=gsub(".ped","",files[i])
  snps=read.table("../indian_uk.map",stringsAsFactors = F)
  names(snps)=c("chr","rsid","0","pos")
  SNP=NULL
  for (j in 1:nrow(snps)){
    snp=c(snps$rsid[j],snps$rsid[j])
    SNP=append(SNP,snp)
  }
  #d=read.table(paste("../../ped_files/",files[i],sep=""),stringsAsFactors=F,colClasses="character")
  d=read.table(files[i],stringsAsFactors=F,colClasses="character")
  d=d[,-(3:6)]
  names(d)=c("Group","sample_ID",SNP)
  # remove singleton SNPs (TRVs)#############################################
  d=d[,which(is.na(match(names(d),c("rs72552267","rs55640102"))))]
  # put chrom numbers as names(d)
  snps=snps[which(is.na(match(snps$rsid,c("rs72552267","rs55640102")))),]
  CHROM=NULL
  for (j in 1:nrow(snps)){
    chrom=c(snps$chr[j],snps$chr[j])
    CHROM=append(CHROM,chrom)
  }
  names(d)=c("Group","sample_ID",CHROM)
  # chromosomes with genes with each having at least two SNPs for haplotyping (excluding FMO2)
  for (k in 1:length(chrs)){
    u=d[,which(names(d)==chrs[k])]
    u[u==0]="?"
    u[u=="GA"]="A"
    u[u=="TA"]="A"
    #PHASE input
    if (ncol(u)>2){
      cat(paste0(nrow(u),"\n",ncol(u)/2,"\n",paste0(rep("S",times=ncol(u)/2),collapse=""),"\n"),file=paste0("IN_PHASE/",paste("new",genes[k],name,sep="_")))
      write.table(u,paste0("IN_PHASE/",paste("new",genes[k],name,sep="_")),quote=F,col.names=F,row.names=F,append=T,sep="")
      print(paste(files[i]))
    }
  }
}
