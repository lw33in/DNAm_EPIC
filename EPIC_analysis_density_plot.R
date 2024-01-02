# DNA methylation EPIC analysis & color density plots

# Load data --------------------------------------------------------------------
# Load beta values
beta <- read.csv("/EPICarray/beta_values-noob.csv")
head(beta)
# Save each sample for GEO 
condition1_beta_geo <- beta [c(1,2) , ]
# Load manifest
manifest <- read.csv("/EPICarray/MethylationEPIC_v-1-0_B4.csv")
head(manifest)

# Data pre-processing ----------------------------------------------------------
# Assign "IlmnID" (cgxxxxx) to "UCSC_RefGene_Name" and "UCSC_RefGene_Group"
cg_w_gene_location <- manifest [ , c("IlmnID" , "UCSC_RefGene_Name" , "UCSC_RefGene_Group")  ]
# Overlap "cg_w_gene_location" with our "beta output" 
beta_cleaned <- beta %>%
  colnames(c( "IlmnID" , "condition1" , "condition2" , "condition3" , "wt"))
beta_cleaned <- beta_cleaned [ -c(1,2,3) ,]

# Add cg/gene/location information to beta 
cg_w_gene_location # from manifest 
beta_w_gene_location <- merge(cg_w_gene_location,
                             beta_cleaned,
                             by=c("IlmnID"))
beta_w_gene_location[c(1:60),]
beta_w_gene_location$UCSC_RefGene_Group

# UCSC_RefGene_Group: Gene region feature category describing the CpG position, from UCSC.
# Features listed in the same order as the target gene transcripts.
  # TSS200 = 0–200 bases upstream of the transcriptional start site (TSS).
  # TSS1500 = 200–1500 bases upstream of the TSS.
  # 5'UTR = Within the 5' untranslated region, between the TSS and the ATG start site.
  # Body = Between the ATG and stop codon; irrespective of the presence of introns, exons, TSS, or promoters.
  # 3'UTR = Between the stop codon and poly A signal.

# All probes for the gene (5'-3')
beta_w_gene_location_all <- beta_w_gene_location
beta_w_gene_location_all$UCSC_RefGene_Group <- gsub(";.*","",beta_w_gene_location_all$UCSC_RefGene_Group)

# Remove probes at gene body and 3'UTR
beta_w_gene_location  <- beta_w_gene_location[!grepl("Body", beta_w_gene_location$UCSC_RefGene_Group),]  %>%
  beta_w_gene_location[!grepl("3'UTR", beta_w_gene_location$UCSC_RefGene_Group),]  
# beta_w_gene_location  <- beta_w_gene_location[!grepl("1stExon", beta_w_gene_location$UCSC_RefGene_Group),]   
beta_w_gene_location[1:50,]
beta_w_gene_location$UCSC_RefGene_Group = gsub(";.*","",beta_w_gene_location$UCSC_RefGene_Group)
dim(beta_w_gene_location) 

# Clean up gene names 
beta_w_gene_location [1:70,]
beta_w_gene_location$UCSC_RefGene_Name = gsub(";.*","",beta_w_gene_location$UCSC_RefGene_Name)

write.table ( beta_w_gene_location, 
             "/Users/WeiyaNi/Desktop/FarnhamLab/ZFX_project/EPICarray/promoter_probe_gene_location_beta.txt",
              quote = FALSE,
              sep="\t",
              col.names=T, row.names=F)

# Improt promoter probes : TSS1500 , TSS200, 5'UTR, 1st exon (excluding body and 3'UTR)
beta_w_gene_location <- read.delim("/EPICarray/promoter_probe_gene_location_beta.txt", header = T)

# Overlap with TF bound promoters ----------------------------------------------
down_bound_846 <- read.delim("/Knockout/down_bound_bothTF.txt",header = F)
noFC_bound_850 <- read.delim("/Knockout/noFC_dound_bothTF.txt",header=F)

down_bound_846_EPIC <- beta_w_gene_location [ beta_w_gene_location$UCSC_RefGene_Name %in% down_bound_846$V1 , ]
# Check the number of genes has EPIC data
length(unique(down_bound_846_EPIC$UCSC_RefGene_Name))
# Extract upstream probes
down_bound_846_EPIC_upstream  <- down_bound_846_EPIC[!grepl("1stExon", down_bound_846_EPIC$UCSC_RefGene_Group),] %>%
  down_bound_846_EPIC_upstream[!grepl("5'UTR", down_bound_846_EPIC_upstream$UCSC_RefGene_Group),] 
# Extract downstream probes
down_bound_846_EPIC_downstream  <- down_bound_846_EPIC[!grepl("TSS1500", down_bound_846_EPIC$UCSC_RefGene_Group),] %>% 
  down_bound_846_EPIC_downstream[!grepl("TSS200", down_bound_846_EPIC_downstream$UCSC_RefGene_Group),] 
# Update `down_bound_846_EPIC` with genes with matched EPIC probes
down_bound_846_EPIC <- rbind(down_bound_846_EPIC_upstream,
                             down_bound_846_EPIC_downstream)

# All promoter probes in condition1
down_bound_846_EPIC_condition1 <- down_bound_846_EPIC[, c("condition1" , "wt")]
down_bound_846_EPIC_condition1$change <- abs(down_bound_846_EPIC_condition1$Seven - down_bound_846_EPIC_condition1$WT)

# Color density plot -----------------------------------------------------------
# Plot out beta change OVER 0.2 
down_bound_846_EPIC_condition1 <- subset(down_bound_846_EPIC_condition1, down_bound_846_EPIC_condition1$change >0.2 )
condition1_kennel <- matrix(down_bound_846_EPIC_seven[,"Condition1"])
wt_kennel <- matrix(down_bound_846_EPIC_seven[,"wt"])

Lab.palette <- colorRampPalette( c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), space = "Lab")
pdf("/plots/condition1.846.allProbe.pdf",width=3.7,height=4,paper='special')
i.s <- smoothScatter(wt_kennel, seven_kennel,
                     colramp = Lab.palette,
                     pch=NA, 
                     nrpoints = Inf,
                     xlab = "wt",ylab = "Condition1")
dev.off()

sessionInfo()
