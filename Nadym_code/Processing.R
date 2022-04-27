#Set requires
library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(ape)
library(dplyr)

#Set some incoming data
path_train_set = "/home/alexey/tax_n_refs/silva_nr_v132_train_set.fa.gz"
path_train_set_species = "/home/alexey/tax_n_refs/silva_species_assignment_v132.fa.gz"
map_path = 'nadym_map.csv'
setwd('/home/alexey/Analysis/2020_Nadym/')
path <- 'raw/'


#dada2 raw pipeline
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180), trimLeft=c(19, 20), maxN=0, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, multithread=8)
errF <- learnErrors(filtFs, multithread=15)
errR <- learnErrors(filtRs, multithread=15)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=15, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=15, pool=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
getN <- function(x) sum(getUniques(x))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=15, verbose=TRUE)


##assign taxonomy
taxa.dada2 <- assignTaxonomy(seqtab.nochim, path_train_set , multithread=15)
taxa.dada2 <- addSpecies(taxa.dada2, path_train_set_species)

#Read metadata
mdat <- read.csv(map_path, sep = '\t')
rownames(mdat) <- mdat$SampleID
mdat <- mdat %>% arrange(factor(Filename, levels = rownames(seqtab.nochim)))

if (all(mdat$Filename == rownames(seqtab.nochim))) {
  rownames(seqtab.nochim) <- rownames(mdat)
}


#Make a phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(mdat), 
               tax_table(taxa.dada2))


#Extract dna seqs in taxa_names end export them
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
writeXStringSet(ps@refseq, format = 'fasta', filepath = 'refseqs.fasta')

####################################################################
##Phylogeny? By fasttree or iqtree from QIIME2-plugin

#In QIIME2 run this:

#  qiime tools import --input-path refseqs.fasta --output-path rep-seqs.qza --type 'FeatureData[Sequence]'
#  qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza
#  qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

#  qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree masked-fasttree.qza
#  extract tree file from archives

#read tree and merge
m.fasttree <- read_tree('tree.nwk')
ps <- merge_phyloseq(ps, phy_tree(m.fasttree))
####################################################################
saveRDS(ps, "ps.RData")
