-------------------------------------
#title: "ACC on proteoms calculator"
#author: "Goulancourt Rebecca"
#date: "06/02/2022"
-------------------------------------
# Command line : Rscript --vanilla acc.R 
  
# This script was highly inspired by Clothilde Garrido's work (acc.R) at :  https://github.com/UMR7141/Peptides_Analysis.
  
  
# Librairies

#install.packages('seqinr')
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

library('protr')
library('Biostrings')
library('optparse')
library('stringr')
library('M3C')
library('Rtsne')
library('RColorBrewer')
library('seqinr')


# Necessary paths

path_output <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/output"
path <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/"

# Data reading
list_of_aa = c('M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H')
files <- list.files(path = path, pattern = (".f")) 



reading = function(file) {
  #data = readFASTA(paste0(path, "proteome_diatom.faa"),
  data = readFASTA(file,
                   legacy.mode = TRUE, seqonly = FALSE)
  df = data.frame(nrow = length(data))
  for (i in 1:length(data)) {
    df[i, 1] = cbind(data[i])
  }
  rownames(df) = names(data)
  return(df)
}

fastaFile = c()
df_list = NULL
for (f in files) {
  df = data.frame()
  #fastaFile <- readDNAStringSet(paste0(path, f))
  #seq_name = names(fastaFile)
  #sequence = paste(fastaFile)
  #df <- data.frame(seq_name, sequence)
  print(f)
  df = reading(paste0(path, f))
  print("ok")
  df = assign(paste0("df_", f), df)
  df_list = c(df_list, list(df))
  print("--> ok")
}


Nm = c()
id_seq = c()
k = 0
df2_list = NULL
for (df in df_list) {
  df2_ = data.frame()
  N = 0
  k = k + 1
  for (i in 1:nrow(df)) {
      if (nchar(as.character(df[i,1])) < 100) {
      df2_ = df[-i, 1]
      df[i,] = NA
      N = N + 1
    }
  }
  Nm = c(Nm, N)
  
  df = na.omit(df)
  df2_ = assign(paste0("df2_", files[k]), df)
  df2_list = c(df2_list, list(df2_))
  print(paste0("nb of deleted seq : ", N))
  print(nrow(df))
}
#print(Nm, "total deleted sequences")

for (df in df2_list) {
  print("-------")
  for (i in 1:nrow(df)) {
    if (nchar(as.character(df[i, 1])) < 4) {
      print(i)
    }
  }
}

# Calculation of ACC

k = 0
Acc_list = NULL
#for (df in df_list) {
for (df in df2_list) {
  k = k + 1
  mat_vect=c()
  Acc = data.frame()
  for (s in (1:dim(df)[1]))
  #for (s in (1:50))
{
  #seq = as.character(df[s, 2])
  seq = as.character(df[s, 1])
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

# --> faire un return de l'ACC pour le lire dans python ou les écrire dans un fichier puis le lire via python

setwd(path_output)

for (i in 1:length(files)) {
  write.table(Acc_list[i], file = paste0("Acc_output_", files[i], ".txt"),  
              append = FALSE, sep = "\t", dec = ".", row.names = TRUE,
              col.names = TRUE)
}


for (acc in Acc_list) {
  for (f in files) {
write.table(acc, file = paste0("Acc_output_", f, ".txt"), append = FALSE, sep = "\t",
            dec = ".", row.names = TRUE, col.names = TRUE)
  }
}

# Plot of Acc with Tsne

setwd(path_output)
pdf("Tsne_ACC.pdf", height = 10,width = 10)
tsne = NULL
col = NULL
col =  palette(rainbow(length(Acc_list))) 
i = 0
for (acc in Acc_list) {
  i = i + 1
  tsne = Rtsne(acc, labels = as.factor(df$seq_name), perplex = 0.0001, check_duplicates = FALSE)
  if (i == 1) {
  plot(tsne$Y, type = "p", col = col[i])
  }
  else {
  lines(tsne$Y, type = "p", col = col[i])
  }
}
title("Tsne of ACC on proteoms")
legend("topleft", legend = files,
       col = col, lty=1:2, cex=0.8)
dev.off()






