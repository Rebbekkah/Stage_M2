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

# Necessary path

 
path <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/proteome_diatom.faa"
data2 <-  read.csv(path)

# Data reading & arguments
fastaFile <- readDNAStringSet(path)
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)

#Command line : Rscript --vanilla acc.R -f <fasta-file> -a <column-name> -l <lag> -o <out-file>

option_list  =  list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = "fasta file", metavar = "character"),
  make_option(c("-a", "--accession"), type = "character", default = NULL, 
              help = 'name column sequence', metavar = "character"),
  make_option(c("-l", "--lag"), type = "integer", default = NULL,
              help = "l refers to lag, which is the interval between residues being compared"),
  make_option(c("-o", "--out-file"), type = "character", default = NULL,
              help = "out csv file")
  );

opt_parser  =  OptionParser(option_list = option_list);
opt  =  parse_args(opt_parser);

 
list_of_aa = c('M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H')
print(list_of_aa)


mat_vect=c()
Acc = data.drame()
for (s in (1:dim(df[1:15,])[1]))
{
  #print(s)
  seq = as.character(df[s, 2])
  #print(seq)
  mat = rbind(AAindex[390,str_sub(seq,1,1)],AAindex[391,str_sub(seq,1,1)],AAindex[392,str_sub(seq,1,1)])
  #print(mat)
  #break
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

rownames(Acc) <- df$seq_name[1:15]

#tsne(mat[, 1:ncol(mat)], labels = as.factor(df$seq_name), perplex = 1^-90000000000000000000000000)
#tsne(mat[1,], labels = as.factor(df$seq_name), perplex = )
#RunTSNE(object = mat[1,], perplexity = ...)


Rtsne(mat, labels = as.factor(df$seq_name), perplex = )


Acc = data.frame()
for (seq in df[, 2]) {
  print(seq)
  df[seq, 'sequence'] = acc((seq), lag = 4)
}



######################################################################################
######################################################################################
for (elem in data2[1:100, 1]) {
  if ('>' %in% elem) {
    
  }
}


#data <- read.csv(opt$f, sep = '\t')
#is.data.frame(data2)
`%!in%` = Negate(`%in%`)
df = data.frame()
idt <- c()
seq <- c()
for (elem in data2[1:100, 1]) {
  #vec = c(vec, elem)
  elem2 = str_split(elem, pattern = "")
  #print(elem2)
  #print(typeof(elem2))
  #seq = c(str_split(elem, pattern = ""))
  #print(seq)
  for (letter in elem2) {
    #print(typeof(letter))
    #print(letter) 
    if (letter == '>') {
      idt = c(idt, elem)
    }
    else if (letter %in% list_of_aa) {
      seq = c(seq, elem)
    }
  }
}
idt[1]
seq[1]
#df = data.frame(idt, seq)

#idt_splited = str_split(idt, pattern = " ")
#typeof(idt_splited)
#idt_splited[]
#idt2 = substr(idt_splited, )
idt2











