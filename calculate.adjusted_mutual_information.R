#!/usr/bin/Rscript
################################################################################
#Calculate pairwise AMI, NMI and ARI between pairs of OTU sets
#
#AMI		=> Adjusted Mutual Information
#NMI		=> Normalized Mutual Information
#ARI		=> Adjusted Rand Index
#
#All measures are calculated following the formulas given in
#Vinh et al, 2009, Information Theoretic Measures for Clustering Comparison [...]
#Schmidt et al, 2014, Limits to Robustness & Reproducibility [...] (-> Supplementary Text S1)
#
#2016-04-10
#sebastian.schmidt@imls.uzh.ch
################################################################################


################################################################################
################################################################################
# Load Packages
library("parallel", warn.conflicts=F, quietly=T);

#
#MOVE TO WORKING DIRECTORY (WHERE THE INPUT FILES ARE)
#

#Set global parameters
use.cores = 40;				#how many CPUs to use during parallelization
seq.count = 887870;		#total number of sequences in the dataset; a known variable

#Define helper function
#=> needed later for more convenient data parsing
split.feature.vector <- function (X) {eval(parse(text=paste("c(", substr(X, 1, nchar(X)-1), ")")))}
################################################################################
################################################################################


################################################################################
################################################################################
#Precalculate sum(n.ij) matrix for expected mutual information (EMI)
#This is convenient, because the EMI depends only on OTU sizes and can later
#be reutilized when calculating AMI values. A lookup is much faster than recalculating
#EMI values at each step all over again.
#EMI is calculated for all combinations of OTU sizes from 1 to 1000, providing a
#handy 1000x1000 (triangular) lookup matrix (or list, in this case).
#EMI is calculated as described in formula 24a in Vinh et al, 2009.
n <- seq.count;
#Timing
writeLines(paste(date(), "=> Calculating EMI lookup table."), sep = "\n");
#Nested listwise computation
emi.lookup <- mclapply(
  seq(1, 1000),
  function (a.i) {
    unlist(lapply(
      seq(1, 1000),
      function (b.j) {
        n.ij <- seq(max(a.i + b.j - n, 1), min(a.i, b.j));
        sum((n.ij / n) * log((n.ij * n) / (a.i * b.j)) * exp(lfactorial(a.i) + lfactorial(b.j) + lfactorial(n - a.i) + lfactorial(n - b.j) - lfactorial(n) - lfactorial(n.ij) - lfactorial(a.i - n.ij) - lfactorial(b.j - n.ij) - lfactorial(n - a.i - b.j + n.ij)))
      }
    ));
  },
  mc.cores = use.cores
);
#Timer stop
writeLines(paste(date(), "=> Done."), sep = "\n");
################################################################################
################################################################################


################################################################################
################################################################################
#Load data
#Importantly, OTU datasets are loaded in two different formats. In this case, we are
#comparing the average linkage (al) clustering to complete linkage (cl) clustering
#at a similarity threshold of 97%. Let's call the al-dataset "set.1" and the cl-dataset "set.2".
#Both sets will be loaded in different formats, for R convenience.
#
#set.1 is in an OTU-wise format, like so:
#[otu_idx]	[sequence_indices, comma-separated]
#=> listing for every OTU all the sequences it contains.
#For example:
#1	1,10,12,15,20,
#2	2,7,8,9,11,25,27,
#...
#For an OTU set where the first OTU contains sequences 1,10,12,15 and 20;
#the second OTU contains seqs 2,7,8,etc. This data format is called "OTU mapping",
#and for parsing it, the convenience function "split.feature.vector" is used (see above).
#
#set.2 is in a sequence-wise format, like so:
#[seq_idx]	[otu_idx]
#=> listing for every sequence the OTU into which it clustered.
#For example (using the above example):
#1	1
#2	2
#...
#7	2
#8	2
#9	2
#10	1
#...

#Load OTU mapping for set.1 ("mapping.1")
#=> in list format
tmp.data <- scan(file = "al.97.otu_mapping.otu", what = list(integer(), character()), sep = "\t");
mapping.1 <- list(); length(mapping.1) <- length(tmp.data[[1]]);
mapping.1[tmp.data[[1]]] <- mclapply(tmp.data[[2]], split.feature.vector, mc.cores = 8);			#This step may take a moment for large files...
#Get OTU sizes for set.1
sizes.1 <- unlist(lapply(mapping.1, length));
#Get total number of OTUs for set.1
otu_count.1 <- length(mapping.1);

#Load sequence mapping for set.2 ("seq.mapping.2")
tmp.data <- scan(file = "cl.97.seq_mapping.otu", what = list(integer(), integer()), sep = "\t");
seq.mapping.2 <- numeric(max(tmp.data[[2]])); seq.mapping.2[tmp.data[[1]]] <- tmp.data[[2]];
#Get OTU sizes for set.2
sizes.2 <- tabulate(seq.mapping.2);
#Get total number of OTUs for set.2
otu_count.2 <- length(sizes.2);
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate indices

#Compute shared mapping (n.ij)
#=> the sparse "contingency table"
#The output list "current.shared" contains a sublist for each OTU in set.1
#(length(current.shared) == otu_count.1)
#This sublist has two entries, "$n.ij" and "sizes.2". With the example data:
#
#> current.shared[[1]]
# $n.ij
# [1] 573 118 188  53  52   3  12   1   4   8   1   1   1   1   1   1
# 
# $sizes.2
# [1] 573 118 188  53  52   3  12   1   4   8   1   1   1   1   1   1
#
#...this means that the 1018 sequences in OTU 1 of set.1 (check sizes.1[1]) clustered
#into 16 different OTUs in set.2 (try length(current.shared[[1]]$n.ij)). For each of these
#OTUS, "$n.ij" holds the number of concordant sequences (those in OTU 1 of set.1 and the respective OTU in set.2).
#At the same time, "$sizes.2" has the size of these OTUs in set.2. In the current case, all 16 OTUs in set.2
#contain *only* sequences of OTU 1 in set.1. In other words, cl clustering split the al OTU 1 into 16 smaller OTUs
#which are completely contained in that one OTU from set.1.
#
#Note: this is usually the most time-consuming step.
writeLines(paste(date(), "Calculating current shared mapping..."));
current.shared <- mclapply(mapping.1, function (map.1) {
  tmp <- tabulate(seq.mapping.2[map.1]);
  out <- list();
  out$n.ij <- as.numeric(tmp[tmp != 0]);
  out$sizes.2 <- sizes.2[which(tmp != 0)];
  out
}, mc.cores = use.cores);

#Compute ARI
#=> equation (3) in Vinh et al., 2009
writeLines(paste(date(), "Calculating Adjusted Rand Index..."));
#Calculate size-dependent terms ("sums") for the expected index
sum.1 <- sum(unlist(lapply(sizes.1, function(X) {choose(X, 2)})));
sum.2 <- sum(unlist(lapply(sizes.2, function(X) {choose(X, 2)})));
#Calculate expected Rand Index for these OTU size distributions
expected.index <- (sum.1 * sum.2) / choose(n, 2);
#Calculate Rand Index based on shared mapping (contingency table)
rand.index <- sum(unlist(mclapply(current.shared, function (shared) {choose(shared$n.ij, 2)}, mc.cores = use.cores)));
#Calculate ARI from all of the above
ari <- (rand.index - expected.index) / (0.5*(sum.1 + sum.2) - expected.index);
#With the example data, we find an ARI of 0.585

#Compute Mutual Information & Entropies (for NMI & AMI calculation)
#=> equations (4) & (5) in Vinh et al., 2009
writeLines(paste(date(), "Calculating Mutual Information and Entropies..."));
mutual.information <- sum(unlist(mclapply(current.shared, function (shared) {sum((shared$n.ij / n) * log((shared$n.ij * n) / (sum(shared$n.ij) * shared$sizes.2)))}, mc.cores = use.cores)));
entropy.1 <- sum((sizes.1 / n) * log(sizes.1 / n));
entropy.2 <- sum((sizes.2 / n) * log(sizes.2 / n));

#Compute NMI from mutual information and entropies
nmi <- (-2 * mutual.information) / (entropy.1 + entropy.2);
#Should be 0.944

#Compute Expected Mutual Information (EMI)
#=> for the AMI formula
#Now, this is a little involved. Above, we calculated a lookup table for EMI values for
#certain combinations of OTU sizes (1,..,1000). However, some OTUs in the set are larger,
#and for these, we have to calculate the EMI on the way.
writeLines(paste(date(), "Calculating Expected Mutual Information..."));
tab.sizes.1 <- as.numeric(tabulate(sizes.1)); idx.sizes.1 <- which(tab.sizes.1 != 0);
tab.sizes.2 <- as.numeric(tabulate(sizes.2)); idx.sizes.2 <- which(tab.sizes.2 != 0);
idx.small.sizes.2 <- idx.sizes.2[idx.sizes.2 <= 1000];
idx.large.sizes.2 <- idx.sizes.2[idx.sizes.2 >  1000];
expected.mutual.information <- sum(unlist(mclapply(
  idx.sizes.1,
  function (a.i) {
    #Treat a.i <= 1000
    if (a.i <= 1000) {
      #Sum up n.ij (precalculated) for all b.j <= 1000
      small.sum <- sum(tab.sizes.1[a.i] * tab.sizes.2[idx.small.sizes.2] * emi.lookup[[a.i]][idx.small.sizes.2]);
      #Additionally, calculate n.ij for b.j > 1000
      large.sum <- sum(tab.sizes.1[a.i] * tab.sizes.2[idx.large.sizes.2] * unlist(lapply(
        idx.large.sizes.2,
        function (b.j) {
          n.ij <- as.numeric(seq(max(a.i + b.j - n, 1), min(a.i, b.j)));
          num.a.i <- as.numeric(a.i); num.b.j <- as.numeric(b.j);
          sum((n.ij / n) * log((n.ij * n) / (num.a.i * num.b.j)) * exp(lfactorial(num.a.i) + lfactorial(num.b.j) + lfactorial(n - num.a.i) + lfactorial(n - num.b.j) - lfactorial(n) - lfactorial(n.ij) - lfactorial(num.a.i - n.ij) - lfactorial(num.b.j - n.ij) - lfactorial(n - num.a.i - num.b.j + n.ij)))
        }
      )))
      #Pass output
      small.sum + large.sum
    }
    else {
      #Sum over n.ij for all a.i >= 1000 and all b.j
      sum(tab.sizes.1[a.i] * tab.sizes.2[idx.sizes.2] * unlist(lapply(
        idx.sizes.2,
        function (b.j) {
          n.ij <- as.numeric(seq(max(a.i + b.j - n, 1), min(a.i, b.j)));
          num.a.i <- as.numeric(a.i); num.b.j <- as.numeric(b.j);
          sum((n.ij / n) * log((n.ij * n) / (num.a.i * num.b.j)) * exp(lfactorial(num.a.i) + lfactorial(num.b.j) + lfactorial(n - num.a.i) + lfactorial(n - num.b.j) - lfactorial(n) - lfactorial(n.ij) - lfactorial(num.a.i - n.ij) - lfactorial(num.b.j - n.ij) - lfactorial(n - num.a.i - num.b.j + n.ij)))
        }
      )))
    }
  },
  mc.cores = use.cores
)));
#For the example data, this gives an EMI of 2.935. This is the mutual information
#that can be expected by chance for two OTU sets with the OTU size distributions of set.1 and set.2.

#Compute AMI
ami <- (mutual.information - expected.mutual.information) / (sqrt(entropy.1 * entropy.2) - expected.mutual.information);
#Should be 0.913

#Stop the timer
writeLines(paste(date(), "=> Done."));
################################################################################
################################################################################


q()

