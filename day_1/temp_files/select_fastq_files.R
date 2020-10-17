library(tidyverse)

dat <- read_tsv("/Users/nt246/Box/Lab\ work/Genomic\ Libraries/Nextera/Silversides_Master_SampleSheet_incl_ReSeq.txt")

View(dat)
table(dat$Population)
geo <- c("Magdalen_Is", "Minas_Basin", "Patchogue", "Jekyll_Is")

pop <- dat %>% 
  filter(Population %in% geo, !is.na(`Nextera1-2_Run`)) %>% 
  select(ID, Population, `Nextera1-2_GNomEx_ID`, Nextera3_Lane, Nextera3_GNomEx_ID) %>% 
  arrange(ID)

View(pop)
pop

write_tsv(pop, path = "/Users/nt246/Dropbox/To\ transfer/Physalia-course/temp/samples_for_fastq_exercise.txt")


qc <- read_tsv("/Users/nt246/Box/Lab\ work/Genomic\ Libraries/Nextera/AllSamples/SummaryTables/QC_Stats_With_Sample_And_Mapping_Info_AllNextera.txt")

View(qc)

qc %>% 
  filter(Population %in% geo, !is.na(`Nextera1.2_Run`)) %>% 
  select(ID, Population, `Nextera1.2_Run`, Nextera3_Lane, MeanFragLength) %>% 
  head(20)

