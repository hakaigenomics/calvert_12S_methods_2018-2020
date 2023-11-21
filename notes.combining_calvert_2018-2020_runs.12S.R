#notes for combining calvert 12S runs
#processed by: Evan Morien
#date: Sept 14th, 2021
#working directory: ~/projects/12S_runs/calvert_12S/

####Libraries####
library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(ShortRead)
library(Biostrings)
library(seqinr)

####Environment Setup####
theme_set(theme_bw())
setwd("~/projects/12S_runs/calvert_12S/")


#read in sequence tables from the three runs
#2019-2020
seqtab.2019 <- fread("calvert2019_2020_Run20210809/sequence_table.12S.merged.length_var.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.2019) <- paste0(seqtab.2019[,1], "_2019") #set row names. adding a tag here so that the mergeSequenceTable function won't throw an error. we can merge on row names later, after using the sequence merging capabilities of the dada2 function.
seqtab.2019 <- seqtab.2019[,-1] #remove column with row names in it
seqtab.2019 <- as.matrix(seqtab.2019) #cast the object as a matrix
mode(seqtab.2019) <- "numeric"

#2018
seqtab.2018 <- fread("calvert2018_Re-Run20210903/sequence_table.12S.merged.length_var.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.2018) <- paste0(seqtab.2018[,1], "_2018") #set row names. adding a tag here so that the mergeSequenceTable function won't throw an error. we can merge on row names later, after using the sequence merging capabilities of the dada2 function.
seqtab.2018 <- seqtab.2018[,-1] #remove column with row names in it
seqtab.2018 <- as.matrix(seqtab.2018) #cast the object as a matrix
mode(seqtab.2018) <- "numeric"

#2019-2020 RERUN
seqtab.rerun <- fread("calvert2019_2020_RE-RUN_20210827/sequence_table.12S.merged.length_var.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.rerun) <- paste0(seqtab.rerun[,1], "_redo") #set row names
seqtab.rerun <- seqtab.rerun[,-1] #remove column with row names in it
seqtab.rerun <- as.matrix(seqtab.rerun) #cast the object as a matrix
mode(seqtab.rerun) <- "numeric"


#merge sequence tables
st.all <- mergeSequenceTables(seqtab.2019, seqtab.2018, seqtab.rerun) 

#making unique names for the PCR negatives and blanks, so they don't get merged
row.names(st.all) <- gsub("Negative-PCR_redo", "Negative-PCR_2019-RERUN", row.names(st.all))
row.names(st.all) <- gsub('Negative-PCR_2019$', "Negative-PCR_2019-1stRUN", row.names(st.all))

row.names(st.all) <- gsub("Blank001_redo", "Blank001_RERUN", row.names(st.all))
row.names(st.all) <- gsub("Blank002_redo", "Blank002_RERUN", row.names(st.all))
row.names(st.all) <- gsub("Blank001_2019", "Blank001_1stRUN", row.names(st.all))
row.names(st.all) <- gsub("Blank002_2019", "Blank002_1stRUN", row.names(st.all))

row.names(st.all) <- gsub("NEG1-Rep1_2018", "NEG1-Rep1-2018", row.names(st.all))
row.names(st.all) <- gsub("NEG1-Rep2_2018", "NEG1-Rep2-2018", row.names(st.all))
row.names(st.all) <- gsub("NEG2-Rep1_2018", "NEG2-Rep1-2018", row.names(st.all))
row.names(st.all) <- gsub("NEG3-Rep1_2018", "NEG3-Rep1-2018", row.names(st.all))

#divide into two ASV tables again, by the 2018 run vs new run
a <- which(grepl("2018", row.names(st.all)) == TRUE)
seqtab.new <- as.data.frame(st.all[-a,])
seqtab.2018 <- as.data.frame(st.all[a,])
b <- which(grepl("2019", row.names(seqtab.new)) == TRUE)
seqtab.2019 <- as.data.frame(seqtab.new[b,])
seqtab.redo <- as.data.frame(seqtab.new[-b,])

#remove identifying strings from each separate sequence table
row.names(seqtab.2018) <- gsub("_2018", "", row.names(seqtab.2018))
row.names(seqtab.2019) <- gsub('_2019$', "", row.names(seqtab.2019))
row.names(seqtab.redo) <- gsub("_redo", "", row.names(seqtab.redo))
#make "names" column
seqtab.2018$names <- row.names(seqtab.2018)
seqtab.2019$names <- row.names(seqtab.2019)
seqtab.redo$names <- row.names(seqtab.redo)

#merge data frames on "names" column. this will sum the read counts sample-wise
seqtab.merged <- aggregate(. ~ names, rbind(seqtab.2018,seqtab.2019), sum)
rownames(seqtab.merged) <- seqtab.merged$names

st.all <- aggregate(. ~ names, rbind(seqtab.merged,seqtab.redo), sum)
rownames(st.all) <- st.all$names
st.all$names <- NULL

#### save sequences and do taxonomy assignment with blast ####
##### replace the long ASV names (the actual sequences) with human-readable names####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(st.all)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.length_var.merged_dataset.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "12S_ASV_sequences.length_var.merged_dataset.fasta") #save sequences with new names in fasta format
#IMPORTANT: sanity checks
colnames(st.all) == ASV.seq #only proceed if this tests as true for all elements

#rename your ASVs in the taxonomy table and sequence table objects
colnames(st.all) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(st.all),st.all),"sequence_table.12S.merged.w_ASV_names.length_var.merged_dataset.txt", row.names=FALSE, quote=F, sep="\t")

#run 12S classifier from terrimporter on github (sequences for 12S isolated from mitofish mitochondiral genome repo, classifier trained on these sequences)
java -Xmx248g -jar ~/programs/rdp_classifier_2.13/dist/classifier.jar classify -c 0.8 -t ~/projects/taxonomyDBs/12S_database/terrimporter_12S_fish_classifier/mydata_trained/rRNAClassifier.properties -o taxonomy_table.12S.merged.RDP.txt 12S_ASV_sequences.length_var.fasta

#assign taxonomy with blast NT database at 96% similarity threshold #remember to update blast DB locations when necessary
mkdir blast_96_sim
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2021-11-05/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query 12S_ASV_sequences.length_var.merged_dataset.fasta  -out blast_96_sim/12S_ASV_sequences.length_var.blast.out
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim/12S_ASV_sequences.length_var.blast.out -t ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/rankedlineage.dmp -m ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/merged.dmp  -o blast_96_sim/taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_12S_ASV_sequences.length_var.blast.out > tmp #need header or lca.py breaks
python2 ~/programs/galaxy-tool-lca/lca.species.py -i tmp -o blast_96_sim/taxonomy_table.12S.NCBI_NT.96sim.txt -b 100 -id 96 -cov 50 -t best_hit -tid 98 -tcov 80 -fh environmental,unidentified -flh unclassified

#cleanup
rm blast_96_sim/12S_ASV_sequences.length_var.blast.out #remove blast output without taxonomy
rm taxonomy_12S_ASV_sequences.length_var.blast.out #remove redundant file
mv tmp blast_96_sim/12S_ASV_sequences.length_var.blast.out #replace with taxonomy added blast output


#assign taxonomy with blast NT database at 96% similarity threshold #NO BEST HITS, LCA ONLY
mkdir blast_96_sim_NO_BESTHIT
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2021-11-05/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query 12S_ASV_sequences.length_var.merged_dataset.fasta  -out blast_96_sim_NO_BESTHIT/12S_ASV_sequences.length_var.blast.out
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim_NO_BESTHIT/12S_ASV_sequences.length_var.blast.out -t ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/rankedlineage.dmp -m ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/merged.dmp  -o blast_96_sim_NO_BESTHIT/taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_12S_ASV_sequences.length_var.blast.out > tmp #need header or lca.py breaks
python2 ~/programs/galaxy-tool-lca/lca.species.py -i tmp -o blast_96_sim_NO_BESTHIT/taxonomy_table.12S.NCBI_NT.96sim.txt -b 100 -id 96 -cov 50 -t only_lca -fh environmental,unidentified -flh unclassified

#cleanup
rm blast_96_sim_NO_BESTHIT/12S_ASV_sequences.length_var.blast.out #remove blast output without taxonomy
rm taxonomy_12S_ASV_sequences.length_var.blast.out #remove redundant file
mv tmp blast_96_sim_NO_BESTHIT/12S_ASV_sequences.length_var.blast.out #replace with taxonomy added blast output
