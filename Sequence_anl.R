library(tidyr)
library(dplyr)
library(TraMineR)
library(openxlsx)

#0. setting up ----------------------
data1 <- read.csv("R_sequence_list_0424.csv")

#View the data 
head(data1)

#Remove rows where Surface == "End"
data_clean<-data1 %>% filter (Fomite !="End")


#Reshape data: wide format where each row is a participant (ID), and columns are the touch sequence
data_wide <- data_clean %>%
  arrange(ID, Sequence) %>%
  group_by(ID) %>%
  #mutate(seq_order = row_number()) %>%
  pivot_wider(names_from = Sequence, values_from = Fomite)

# View reshaped data
print(data_wide)

#Remove ID column and define the sequence object------------------
seq_data<-seqdef(data_wide[,-1])

#View the sequence
seq_data

#1. Most frequent complete sequences (including all sequences)--------------
seqf<-seqtab(seq_data,idxs=1:20)
print(seqf)

##Excel file save
seqf.df<-as.data.frame(seqf)
write.xlsx(seqf.df, "top_sequences.xlsx")


#2. Seqence length >=2 analysis-----------------------------

#Count the length of each sequence
seq_lengths<-seqlength(seq_data)

#Select sequence having more than 1 sequence length
seq_data_long<- seq_data[seq_lengths >1,]

#Most frequent complete sequences (only more than 1 sequence)
seqf_long<-seqtab(seq_data_long,idxs=1:25)
print(seqf_long)

##Excel file save
seqf_long.df<-as.data.frame(seqf_long)
write.xlsx(seqf_long.df, "top_sequences_long.xlsx")


#3. Subsequence Pattern Analysis-----------------------
eseq<-seqecreate(data1, id=data1$ID, timestamp=data1$Sequence, event= data1$Fomite)
# Frequent sub-sequence patterns
subseq_freq <- seqefsub(eseq, min.support = 10 , constraint = seqeconstraint(window.size = 1))  # Adjust threshold if needed (e.g. more than 10% supported pattern)

# (window size in here restricts how close the events in a pattern must be to each other.) 
# (e.g. window.size=2 means that: The events in the subsequence can occur within a maximum distance of 2 positions from each other, so A-B-C can be counted (A)-(C) subsquence (distance from A to C is below 2) but A-B-D-C is not (A)-(C) distance is 3 (over windowsize=2))
# View results
print(subseq_freq)

##Excel file save
subseqf.df<-as.data.frame(subseq_freq)
write.xlsx(subseqf.df, "top_Subsequences.xlsx")
###-----------------making plot-----------------###
# Plot most frequent sequences
#windows(width = 9, height = 6)

#tiff("seqfplot.tiff", width = 1200, height = 800, res = 300)
#seqfplot(seq_data, with.legend="right", border = NA) #legend.prop=0.4), cex.legend=0.8)

#dev.off()

# Plot most frequent sequences (long sequence)
#windows(width = 9, height = 6)

#tiff("seqfplot_long.tiff", width = 1200, height = 800, res = 300, po)
#seqfplot(seq_data_long, with.legend="right", border = NA) #legend.prop=0.4), cex.legend=0.8)

#dev.off()





