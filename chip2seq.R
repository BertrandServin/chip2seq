library(tidyverse)
library(seqinr)
library(Rsamtools)
library(BSgenome)

WD = 'capra_hircus/'
Manifile = 'Goat_IGGC_65K_v2_15069617X365016_A2.csv'
BWADB = 'Capra_hircus'
manifest = read_csv(paste0(WD,'manifests/',Manifile),skip=7) %>% drop_na(AlleleA_ProbeSeq)

##testm=head(manifest)

## align probeseq Sequances
## Maybe try using blaster package at some point. This can take a while.
write.fasta(sequences=as.list(manifest$AlleleA_ProbeSeq), names=manifest$Name, paste0(WD,"fasta/alleleA_probeseq.fa"),open='w')
cmd = paste0("bwa mem -o ",WD,"sam/alleleA_probseq.sam ",WD,"bwadb/",BWADB," ",WD,"fasta/alleleA_probeseq.fa")
system(cmd)

##  align context sequence
context.seq = manifest %>% mutate(contextSeq=str_replace(SourceSeq,"\\[[ATGC]\\/[ATGC]\\]","N")) %>% select(Name,contextSeq)
write.fasta(sequences=as.list(context.seq$contextSeq), names=context.seq$Name, paste0(WD,"fasta/context.seq.fa"),open="w")
cmd = paste0("bwa mem -o ",WD,"sam/context.seq.sam ",WD,"bwadb/",BWADB," ",WD,"fasta/context.seq.fa")
system(cmd)

Aprobe.sam = read_table("capra_hircus/sam/alleleA_probseq.sam",
    comment="@",
    col_names=c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","TLEN","SEQ","QUAL")) %>%
        mutate(Name=QNAME,strand = ifelse(FLAG==0,"plus","minus")) %>% select(Name,RNAME,POS,FLAG,strand)
                        
manifest.annot = manifest %>% left_join(Aprobe.sam) %>% filter(IlmnStrand %in% c('TOP','BOT'))
manifest.annot %>% select(Name,IlmnStrand,strand,Chr,RNAME,MapInfo,POS,SNP)
## Get Chromosome numbers from NCBI REF
assembly  = read_table('capra_hircus/fasta/GCF_001704415.1_ARS1_assembly_report.txt',comment='#',col_names=F)
assembly = assembly %>% select(X3,X7) 
colnames(assembly) = c('CHRNAME', 'RNAME')

## Get Reference Alleles
annot = manifest.annot %>% select(Name,IlmnStrand,strand,Chr,MapInfo,SNP,RNAME,POS) %>% left_join(assembly)
annot = annot %>% mutate(CHRNAME = replace_na(CHRNAME,'0')) %>% mutate(refpos = ifelse(strand=='minus',POS-1,POS+50))
annot.gd = annot %>% filter(CHRNAME!='0')
fasta_file=FaFile("capra_hircus/fasta/Capra_hircus")
gr1 = GRanges(annot.gd$RNAME,IRanges(start=annot.gd$refpos,end=annot.gd$refpos))
refbase=getSeq(fasta_file,gr1)
annot.gd$REF = as.data.frame(refbase)$x
write_csv(annot.gd, file=paste0(WD,'manifests/',Manifile,'.annot'))


