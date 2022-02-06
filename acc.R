title: "ACC on proteoms calculator"
author: "Goulancourt Rebecca"
date: "06/02/2022"

# Installation des librairies

#if (!requireNamespace("BiocManager", quietly   =   TRUE))
#    install.packages("BiocManager")
#BiocManager::install("Biostrings")

#install.packages('protr')
#install.packages('Biostrings')
#install.packages('seqinr')
#install.packages('optparse')
#install.packages('stringr')
#install.packages("Rtsne") 
#install.packages('RColorBrewer')
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("M3C")


# Chargement des librairies

 
library('protr')
library('Biostrings')
library('optparse')
library('stringr')
library('M3C')
library('Rtsne')
library('RColorBrewer')

# Necessary path
path_save <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/img"
path <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/"
files <- list.files(path = path, pattern = (".f")) 
#print(files)

# Data reading & arguments

fastaFile = c()
df_list = NULL
for (f in files) {
  df = data.frame()
  fastaFile <- readDNAStringSet(paste0(path, f))
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  df = assign(paste0("df_", f), df)
  df_list = c(df_list, list(df))
}

for (df in df_list) {
  print(head(df$sequence))
  break
}


#Command line : Rscript --vanilla acc.R -f <fasta-file> -a <column-name> -l <lag> -o <out-file>

option_list  =  list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = "fasta file", metavar = "character"),
  make_option(c("-a", "--accession"), type = "character", default = NULL, 
              help = 'name column sequence', metavar = "character"),
  make_option(c("-l", "--lag"), type = "integer", default = NULL,
              help = "l refers to lag, which is the interval between residues being compared"),
  make_option(c("-o", "--out-file"), type = "character", default = NULL,
              help = "out csv file"),
  make_option(c("-p", "--path"), type = "character", default = NULL,
              help = "path to proteom")
  );

opt_parser  =  OptionParser(option_list = option_list);
opt  =  parse_args(opt_parser);

k = 0
Acc_list = NULL
for (df in df_list) {
  k = k + 1
  mat_vect=c()
  Acc = data.frame()
  #Acc = assign(paste0("Acc_", df), Acc)
for (s in (1:dim(df[1:50,])[1]))
{
  seq = as.character(df[s, 2])
  mat = rbind(AAindex[390,str_sub(seq,1,1)],AAindex[391,str_sub(seq,1,1)],AAindex[392,str_sub(seq,1,1)])
  for (i in (2:nchar(seq)))
  {	
    if (str_sub(seq,i,i)=="Z"|str_sub(seq,i,i)=="U"|str_sub(seq,i,i)=="O"|str_sub(seq,i,i)=="J"|str_sub(seq,i,i)=="B"|str_sub(seq,i,i)=="X" )
    {
      cat('Warning : unrecognized amino acid : ',str_sub(seq,i,i),'\n')
    }
    
    mat = cbind(mat,rbind(AAindex[390,str_sub(seq,i,i)],AAindex[391,str_sub(seq,i,i)],AAindex[392,str_sub(seq,i,i)]))
  }
  mat = t(as.matrix(acc(mat, 4)))
  if (length(mat_vect)==0)
  {
    mat_vect = mat
  }
  else
  {
    mat_vect = rbind(mat_vect,mat)
    Acc = rbind(Acc, mat)
  }
}
  Acc = assign(paste0("Acc_", files[k]), Acc)
  Acc_list = c(Acc_list, list(Acc))
}


setwd(path_save)
pdf("Tsne_ACC.pdf", height = 10,width = 10)
tsne = NULL
col = NULL
col =  palette(rainbow(length(Acc_list))) 
i = 0
for (acc in Acc_list) {
  i = i + 1
  #print(i)
  tsne = Rtsne(acc, labels = as.factor(df$seq_name), perplex = 0.0001, check_duplicates = FALSE)
  if (i == 1) {
  plot(tsne$Y, type = "p", col = col[i])
  }
  #par(new = TRUE)
  else {
  lines(tsne$Y, type = "p", col = col[i])
  }
}
title("Tsne of ACC on proteoms")
legend("topleft", legend = files,
       col = col, lty=1:2, cex=0.8)
dev.off()






