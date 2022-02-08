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

path_save <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/img"
path <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/"

# Data reading
list_of_aa = c('M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H')

reading = function(file) {
data = readFASTA(paste0(path, "proteome_diatom.faa"),
#data = readFASTA(file,
          legacy.mode = TRUE, seqonly = FALSE)
df = data.frame(nrow = length(data))
for (i in 1:length(data)) {
  df[i, 1] = cbind(data[i])
}
rownames(df) = names(data)
return(df)
}

#hmm = system.file(paste0(path, "proteome_diatom.faa"), package = "seqinr")
#hmm2 = read.fasta((paste0(path, "proteome_diatom.faa")))
files <- list.files(path = path, pattern = (".f")) 


parsing = function(file) {
  
}

id = c()
seq = c()
truc = read.csv(paste0(path, "proteome_diatom.faa"), header = FALSE)
for (elem in truc[1:24, 1]) {
  first = str_sub(elem, 1, 1)
  #print(first)
  if (first == '>') {
    id = c(id, elem)
    s = c()
  }
  else if (first %in% list_of_aa) {
    #while (first %in% list_of_aa) {
      s = paste0(s, elem, sep = "")
      #break
  }
  seq = c(seq, list(s))
  }
  #seq = c(seq, s)
#}

vec2 = c("a", "b", "c")
vec = c()
machin = c()
for (i in 1:vec2) {
machin = paste0(vec, vec2, sep = "")
}
machin


#####################################
#notin = Negate('%in%')
id = c()
seq = c()
truc = read.csv(paste0(path, "proteome_diatom.faa"), header = FALSE)

all = c()
id = c()
s = c()
df = data.frame()
for(seq in truc[1:24, 1]) {
  all = c(all, list(seq))
}

for (item in all) {
  first = str_sub(item, 1, 1)
  if (first == '>') {
    df = rbind(df, as.character(item))
    s = c()
  }
  else if (first %in% list_of_aa) {
    s = c(s, item)
  }
}


for (i in 1:nrow(truc[1:24]-1)) {
  first2 = truc[i+1,1]
  for (line in truc[1:24, 1]) {
  first = str_sub(line, 1, 1)
  
    if (first %in% list_of_aa) {
      #for (i in 1:nrow(truc[1:24, 1])) {
      #  first2 = truc[i,1]
      #}
      #if (nchar(as.character(seq[1])) != 80)
      seq = c(seq, line)
    }
  else if (first == '>') {
    id = c(id, line)
  }
  }
}
  
df = as.data.frame(id, seq)

seq = c()
id = c()
for(i in nrow(truc)[1:24]) {
  first = str_sub(truc[i,1], 1, 1)
  first2 = str_sub(truc[i+1,1], 1, 1)
  #print(first)
  #if ((first %in% list_of_aa) && (first2 %in% list_of_aa)) {
  if ((first %in% list_of_aa) && (first2 != '>')) {
    seq = c(paste0(seq, truc[i,1]))
  }
  #else if (is.na(first)) {
  else if (first == '>') {
    id = c(id, line)
  }
}
####################################################


#df2 = reading(paste0(path, files[2]))

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

# Calculation of ACC

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

# Plot of Acc with Tsne

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






