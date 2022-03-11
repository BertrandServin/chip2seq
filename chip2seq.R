library(tidyverse)
WD = 'capra_hircus/'
Manifile = 'Goat_IGGC_65K_v2_15069617X365016_A2.csv'
BWADB = 'Capra_hircus'
manifest = read_csv(paste0(WD,'manifests/',Manifile),skip=7)

##testm=head(manifest)
library(seqinr)

## align probeseq Sequances
## Maybe try using blaster package at some point. This can take a while.
write.fasta(sequences=as.list(manifest$AlleleA_ProbeSeq), names=manifest$Name, paste0(WD,"fasta/alleleA_probeseq.fa"),open='w')
cmd = paste0("bwa mem -o ",WD,"sam/alleleA_probseq.sam ",WD,"bwadb/",BWADB," ",WD,"fasta/alleleA_probeseq.fa")
system(cmd)
##system("blastn -query data/alleleA_probeseq.fa -db data/my_db -max_target_seqs 1 -max_hsps 5 -outfmt \"10 qseqid qstart qend saccver sseqid sstart send evalue length mismatch sstrand\"  > data/alleleA_probseq.blast)

##  align context sequence
context.seq = manifest %>% mutate(contextSeq=str_replace(SourceSeq,"\\[[ATGC]\\/[ATGC]\\]","N")) %>% select(Name,contextSeq)
write.fasta(sequences=as.list(context.seq$contextSeq), names=context.seq$Name, paste0(WD,"fasta/context.seq.fa"),open="w")
cmd = paste0("bwa mem -o ",WD,"sam/context.seq.sam ",WD,"bwadb/",BWADB," ",WD,"fasta/context.seq.fa")
system(cmd)
##system("blastn -query data/context.seq.fa -db data/my_db  -qcov_hsp_perc 99 -outfmt \"10 qseqid saccver sseqid sstart send \" | grep  \",1$\"  > data/context.seq.blast")

Aprobe.sam = read_table("capra_hircus/sam/alleleA_probseq.sam",
    comment="@",
    col_names=c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","TLEN","SEQ","QUAL")) %>%
        mutate(Name=QNAME,strand = ifelse(FLAG==0,"plus","minus")) %>% select(Name,RNAME,POS,FLAG,strand)
                        
manifest.annot = manifest %>% left_join(Aprobe.sam) %>% filter(IlmnStrand %in% c('TOP','BOT'))



## Read in blast data
## probeseq.blast=read_csv('data/alleleA_probseq.blast',col_names=c('qseqid','qstart','qend','saccver','sseqid','sstart','send','evalue','length','mismatch','sstrand'))
 
## context.blast=read_csv('data/context.seq.blast', col_names=c('qseqid','c.saccver','c.sseqid','c.sstart','c.send','c.evalue', 'c.length', 'c.mismatch'))

## blastres = left_join(probeseq.blast,context.blast) %>% mutate(matchpos = (sstart>=c.sstart)&(sstart<=c.send)&(send>=c.sstart)&(send<=c.send))

##  %>% filter(mismatch<3) %>% group_by(qseqid) %>% filter(row_number()==1) %>%
##     mutate(strand=if_else((send-sstart)>0,'+','-') )

