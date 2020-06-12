# Download the required libraries

library("seqinr")
library("R.utils")
library("rBLAST")
library("Biostrings")

# Q1: Download the whole set of E.coli gene DNA sequence.

{
  download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
                destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
}
## Decompress the file.

R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",overwrite=TRUE)

## Create the blast database using makeblast() function, to determine the number of sequence present in the E.coli gene set.

makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl", "-parse_seqids")




# Q2:Download the sample fasta sequences (sample file).
{
  download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
                destfile = "sample.fa")
}

sampleEcoli <- read.fasta("sample.fa") 
str(sampleEcoli) 

# Extracting the allocated gene sequence (43)from the sample fasta sequence.
alEcoli <- sampleEcoli[[43]]
alEcoli
str(alEcoli)

# Determine the sequence length and GC content. 

seqinr::getLength(alEcoli)
seqinr::GC(alEcoli)


# Q3:Creat BLAST database to perform blast searches. 


download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",destfile = "mutblast_function.R")
source("myblastn_tab.R")

# To test the function 


myblastn_tab <- function(myseq=alEcoli)
res <- myblastn_tab (myseq = alEcoli, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")

# To check the blast results.
str(res)
res
head(res)

# Extractdown the sample genesequences, similar to E.coli sequence.

# To determine the first 3 hits.

hits <- as.character(res$sseqid[1:3]) 
hits

# (AAC73878 referes to the Escherichia coli O145:H28 112648 DNA, complete genome)

# Question 04

source("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R")

alEcolimutator<-mutator(alEcoli,30) 
res<- myblastn_tab(myseq=alEcolimutator, db= "tophit.fa")
res


# Question 05


# The sequence for mutation, nmut=nmut as it will recognize the number after the sequence as nmut
myfunc <- function(myseq,nmut) {
  mutseq <- mutator( myseq= alEcolimutator, nmut = nmut) 
  res <- myblastn_tab(myseq= mutseq, db= "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa") #for blast
  if (is.null(res)) {myres= 0} else {myres = 1}
  return(myres)
}
##the sequence for mutation, nmut=nmut as it will recognize the number after the sequence as nmut
myfunc <- function(myseq,nmut)
mutseq <- mutator( myseq= alEcolimutator, nmut = nmut) 

#the sequence for mutation, nmut=nmut as it will recognize the number after the sequence as nmut

# Since the length of gene sequence is 411

myfunc(myseq = allocated_Seq, nmut = 411) 

myfunc(myseq=alEcolimutator,nmut=10)
myfunc(myseq = alEcolimutator,nmut=100)
replicate(n=100,expr = myfunc(myseq=alEcolimutator,nmut=100))
mean(replicate(n=100,myfunc(myseq = 100)))
mean(replicate(n=200,myfunc(myseq = 100)))
mean(replicate(n=300,myfunc(myseq = 100)))
mean(replicate(n=500,myfunc(myseq = 100)))
mean(replicate(n=1000,myfunc(myseq = 100)))