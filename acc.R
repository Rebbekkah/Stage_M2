# Installation des librairies

 
#if (!requireNamespace("BiocManager", quietly   =   TRUE))
#    install.packages("BiocManager")
#BiocManager::install("Biostrings")

#install.packages('protr')
#install.packages('Biostrings')
#install.packages('seqinr')
#install.packages('optparse')
#install.packages('stringr')
 

# Chargement des librairies

 
library('protr')
library('Biostrings')
library('optparse')
library('stringr')
 

# Necessary path

 
path <- "/Users/rgoulanc/Desktop/Rebecca/FAC/M2BI/Stage/LAFONTAINE/script/proteomes/proteome_diatom.faa"
data2 <-  read.csv(path)

seq_name = names(s)
sequence = paste(s)
df <- data.frame(seq_name, sequence)
# Arguments

 
 


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
 

# Data reading

 
list_of_aa = c('M', 'Q', 'A', 'L', 'S', 'I', 'P', 'K', 'G', 'V', 'R', 'E', 'F', 'D', 'C', 'T', 'N', 'W', 'Y', 'H')
print(list_of_aa)
 
#rbind les lignes pour que 1 ligne = 1 séquence
 
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
df = data.frame(idt, seq)

#idt_splited = str_split(idt, pattern = " ")
#typeof(idt_splited)
#idt_splited[]
#idt2 = substr(idt_splited, )
idt2





