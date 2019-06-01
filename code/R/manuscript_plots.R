###############################
#### Generate PCA on HBCC #####
#### and MSSM #################
###############################


#MSSM
downloadFile <- function(id){
  fread(synGet(id)$path, data.table = F)
}
downloadFile_version <- function(id , version){
  fread(synGet(id, version = version)$path, data.table = F)
}
# Download reprocessed counts (DLPFC)
COUNT_ID = 'syn17346208'
ALL_USED_IDs = COUNT_ID
COUNT = downloadFile_version(COUNT_ID, version = 2) %>% data.frame()
# rownames(COUNT) = COUNT$gene_id
# COUNT$gene_id = NULL
COUNT$transcript_id.s. = NULL

# Download gene lengths (DLPFC)
GENE.LEN = downloadFile_version('syn17346397', version = 2) %>%
  tidyr::gather(sampleID, Length, -gene_id, -`transcript_id(s)`) %>%
  group_by(gene_id) %>%
  summarise(Length = median(Length, na.rm = T)) %>%
  ungroup() %>% data.frame()
ALL_USED_IDs = c(ALL_USED_IDs, 'syn16783530')

# Get ancestry vector calculated using gemtools
ANCESTRY_ID = 'syn2511399'
ALL_USED_IDs = c(ALL_USED_IDs, ANCESTRY_ID)
ANCESTRY = downloadFile_version(ANCESTRY_ID, version = 3) %>%
  plyr::rename(c('DNA_report..Genotyping.Sample_ID' = 'SNP_report:Genotyping_Sample_ID'))

# Get genotype ids from synapse
GENOTYPE_ID = 'syn16816490'
ALL_USED_IDs = c(ALL_USED_IDs, GENOTYPE_ID)
GENOTYPE = downloadFile_version(GENOTYPE_ID, version = 3) %>%
  dplyr::select(Individual_ID, `SNP_report:Genotyping_Sample_ID`) %>%
  dplyr::inner_join(ANCESTRY)

#HBCC
# # Download reprocessed counts (DLPFC)
COUNT_ID = 'syn17894685'
ALL_USED_IDs = c(ALL_USED_IDs, COUNT_ID)
COUNT_HBCC = downloadFile_version(COUNT_ID, version = 4) %>% data.frame()
# rownames(COUNT_HBCC) = COUNT_HBCC$gene_id
# COUNT_HBCC$gene_id = NULL
COUNT_HBCC$transcript_id.s. = NULL

# Download gene lengths (DLPFC)
GENE.LEN_HBCC = downloadFile_version('syn18324060', version = 3) %>% 
  tidyr::gather(sampleID, Length, -gene_id, -`transcript_id(s)`) %>%
  group_by(gene_id) %>%
  summarise(Length = median(Length, na.rm = T)) %>%
  ungroup() %>% data.frame()
ALL_USED_IDs = c(ALL_USED_IDs, 'syn16783074')

# Get ancestry vector calculated using gemtools
ANCESTRY_ID = 'syn9922992'
ALL_USED_IDs = c(ALL_USED_IDs, ANCESTRY_ID)
ANCESTRY_HBCC = downloadFile_version(ANCESTRY_ID, version = 2 ) %>% 
  plyr::rename(c(ID = 'SNP_report:Genotyping_Sample_ID'))

# Get genotype ids from synapse
GENOTYPE_ID = 'syn16816490'
ALL_USED_IDs = c(ALL_USED_IDs, GENOTYPE_ID)
GENOTYPE_HBCC = downloadFile_version(GENOTYPE_ID, version = 4) %>% 
  dplyr::select(Individual_ID, `SNP_report:Genotyping_Sample_ID`, `SNP_report:Exclude?`) %>% 
  dplyr::inner_join(ANCESTRY_HBCC)

#Metadata
# Get clinical metadata 
CLINICAL_ID = 'syn3354385'
ALL_USED_IDs = c(ALL_USED_IDs, CLINICAL_ID)
CLINICAL = downloadFile_version(CLINICAL_ID, version = 4)
# Get RNASeq QCmetadata
METADATA_QC_DLPFC_ID = 'syn18358379' 
ALL_USED_IDs[length(ALL_USED_IDs)+1] = METADATA_QC_DLPFC_ID
METADATA_QC_DLPFC = downloadFile_version(METADATA_QC_DLPFC_ID, version = 2)
# Get unreleased metadata to map rRNA rate 
md <- downloadFile_version("syn16816488", version = 13) %>% 
  dplyr::select("Individual_ID","Sample_RNA_ID", "rnaSeq_report:rRNA_Rate") %>% 
  rename(rRNA_Rate = `rnaSeq_report:rRNA_Rate`)

METADATA = right_join(CLINICAL, METADATA_QC_DLPFC, by = c("Individual ID" = "Individual_ID")) %>% 
  left_join(., md, by = c("Individual ID" = "Individual_ID","Sample_RNA_ID")) %>% 
  left_join(., md, by = c("Individual ID" = "Individual_ID")) %>% 
  mutate(rRNA_Rate = coalesce(`rRNA_Rate.x`,`rRNA_Rate.y`)) %>% 
  dplyr::select(-one_of("Sample_RNA_ID.y", "rRNA_Rate.x", "rRNA_Rate.y")) %>% 
  rename(Sample_RNA_ID = `Sample_RNA_ID.x`) %>% 
  distinct()

METADATA = rename(METADATA, SampleID = Sample_RNA_ID)

GENO = bind_rows(GENOTYPE, GENOTYPE_HBCC)

ANCESTRY = bind_rows(ANCESTRY, ANCESTRY_HBCC)

MD = left_join(METADATA, GENO, by = c("Individual ID" = "Individual_ID"))

COUNT = full_join(COUNT, COUNT_HBCC, by = c("gene_id"))

# #Filter on missing metadata variables
# ind = METADATA$SampleID [is.na(METADATA$Ethnicity) | is.na(METADATA$Institution) | is.na(METADATA$Dx)]
# writeLines(paste('Following', length(ind), 'counts are missing any metadata'))
# writeLines(paste(ind, collapse = ', '))
# METADATA <- METADATA  %>% dplyr::filter(!(SampleID %in% ind)) 
# 
# ind = METADATA$SampleID [is.na(METADATA$PMI)]
# writeLines(paste('Following', length(ind), 'counts are missing PMI'))
# writeLines(paste(ind, collapse = ', '))
# METADATA <- METADATA  %>% dplyr::filter(!(SampleID %in% ind)) 
# 
# ind = METADATA$SampleID [is.na(METADATA$Reported_Gender)]
# writeLines(paste('Following', length(ind), 'counts are missing gender'))
# writeLines(paste(ind, collapse = ', '))
# METADATA <- METADATA  %>% dplyr::filter(!(SampleID %in% ind)) 
# 
# ind = METADATA$SampleID [is.na(METADATA$Age_of_Death)]
# writeLines(paste('Following', length(ind), 'counts are missing age of death'))
# writeLines(paste(ind, collapse = ', '))
# METADATA <- METADATA  %>% dplyr::filter(!(SampleID %in% ind))
# 
# ind = METADATA$SampleID [is.na(METADATA$EV.1)]
# writeLines(paste('Following', length(ind), 'counts are missing ancestry information'))
# writeLines(paste(ind, collapse = ', '))
# METADATA <- METADATA  %>% dplyr::filter(!(SampleID %in% ind))


## Get GC content from biomart
backgroundGenes = data.frame(gene_id = COUNT$gene_id) %>%
  dplyr::mutate(id = gene_id) %>%
  tidyr::separate(id, c('ensembl_gene_id','position'), sep = '\\.')

# Define biomart object
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "oct2018.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")

# Query biomart
Ensemble2HGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "percentage_gene_gc_content", "gene_biotype", "chromosome_name"),
                       filters = "ensembl_gene_id", values = backgroundGenes$ensembl_gene_id,
                       mart = mart)

GENE.GC.CONT = Ensemble2HGNC %>%
  dplyr::left_join(backgroundGenes) %>% 
  dplyr::select(gene_id, percentage_gene_gc_content, chromosome_name) %>%
  unique
rownames(GENE.GC.CONT) = GENE.GC.CONT$gene_id

SEX.COUNTS = GENE.GC.CONT %>% 
  left_join(COUNT) %>%
  dplyr::select(-one_of("percentage_gene_gc_content")) %>%
  filter(chromosome_name == "X" |chromosome_name == "Y") %>% 
  tidyr::gather(key = item, value = value, -c(gene_id, chromosome_name)) %>%
  mutate(value = log(value)) %>%
  rename(`counts(log)`= value) %>% 
  rename(SampleID = item) %>%
  left_join(METADATA[,c("SampleID", "Sex", "Institution")])


my.theme <- theme_bw() %+replace% theme(legend.position = 'top', axis.text.x = element_text(angle = 90, hjust = 1), plot.title=element_text(hjust=0.5))


##XIST and UTY expression 
#ENSG00000229807.10 and ENSG00000183878.15 
FILT <- SEX.COUNTS %>% 
  filter(gene_id == "ENSG00000229807.10" | gene_id == "ENSG00000183878.15") %>% 
  dplyr::select(-one_of("chromosome_name")) %>% 
  tidyr::spread(key = gene_id, value = `counts(log)`) %>% 
  mutate(XIST = as.numeric(`ENSG00000229807.10`)) %>% 
  mutate(UTY = as.numeric(`ENSG00000183878.15`))

p = ggplot(FILT, aes (x= XIST, y = UTY)) 

p = p + geom_point(aes(color=`Sex`, shape=Institution)) 
p


FactorCovariates <- c('Individual ID', "Institution", "Reported Gender", "Library_Batch", "Dx", "Flowcell_Batch")
ContCovariates <- c("Age of Death", "PMI (in hours)", "RIN", "Mapped_Reads", "Intragenic_Rate", "Intronic_Rate", "Intergenic_Rate",
                    "Genes_Detected", "Expression_Profiling_Efficiency", "rRNA_Rate", "Total_Reads")
                    
                    #("EV.1", "EV.2", "EV.3", "EV.4", "EV.5")

# Find inter relation between factor covariates
COVARIATES = METADATA[,c(FactorCovariates,ContCovariates),drop=F]
COVARIATES[,FactorCovariates] <- data.frame(lapply(COVARIATES[,FactorCovariates],function(x){
  x <- sapply(x,function(y){str_replace_all(as.character(y),'[^[:alnum:]]','_')})}))
rownames(COVARIATES) <- METADATA$SampleID

# Convert factor covariates to factors
COVARIATES[,FactorCovariates] = lapply(COVARIATES[,FactorCovariates], factor)
COVARIATES[,ContCovariates] = lapply(COVARIATES[,ContCovariates], function(x){
  x = as.numeric(as.character(gsub('[\\,\\%]','',x)))
})

# Add in RIN^2 values
# COVARIATES$RIN2 = COVARIATES$RIN^2
# ContCovariates = c(ContCovariates, 'RIN2')

my.theme <- theme_bw() %+replace% theme(legend.position = 'top', axis.text.x = element_text(angle = 90, hjust = 1), plot.title=element_text(hjust=0.5))

# RIN
p = list()
p[[1]] = ggplot(COVARIATES, aes(x = Dx, y = RIN)) + geom_boxplot()
p[[1]] = p[[1]] + ggtitle('RIN') + my.theme

# Age of Death
p[[2]] = ggplot(COVARIATES, aes(x = Dx, y = Age_of_Death)) + geom_boxplot()
p[[2]] = p[[2]] + ggtitle('AgeOfDeath') + my.theme

# PMI
p[[3]] = ggplot(COVARIATES, aes(x = Dx, y = PMI)) + geom_boxplot()
p[[3]] = p[[3]] + ggtitle('PMI (in hours)') + my.theme

# Intronic Rate
p[[4]] = ggplot(COVARIATES, aes(x = Dx, y = IntronicRate)) + geom_boxplot()
p[[4]] = p[[4]] + ggtitle('Intronic Rate') + my.theme

# IntergenicRate
p[[5]] = ggplot(COVARIATES, aes(x = Dx, y = IntergenicRate)) + geom_boxplot()
p[[5]] = p[[5]] + ggtitle('Intergenic Rate') + my.theme

# Transcripts Detected
p[[6]] = ggplot(COVARIATES, aes(x = Dx, y = IntragenicRate)) + geom_boxplot()
p[[6]] = p[[6]] + ggtitle('Intragenic Rate') + my.theme

# Mapped Reads
p[[7]] = ggplot(COVARIATES, aes(x = Dx, y = MappedReads)) + geom_boxplot()
p[[7]] = p[[7]] + ggtitle('Mapped Reads') + my.theme

# PercentAligned
p[[8]] = ggplot(COVARIATES, aes(x = Dx, y = TotalReads)) + geom_boxplot()
p[[8]] = p[[8]] + ggtitle('Total Reads') + my.theme

# rRNARate
p[[9]] = ggplot(COVARIATES, aes(x = Dx, y = rRNARate)) + geom_boxplot()
p[[9]] = p[[9]] + ggtitle('rRNARate') + my.theme

multiplot(plotlist = p, cols = 3)

