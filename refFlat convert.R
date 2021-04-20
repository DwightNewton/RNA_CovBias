#refFlat file rearranging
#Input is refFlat file derived from "gtfToGenePred" function, source: UCSC
options(stringsAsFactors = FALSE)

#load refFlat and ensembl metadata
refFlat <- read.table("GRCh38-r98-refFlat.txt")
ensemblMeta <- read.csv("martquery_1104155809_616.txt")

refFlat$Gene <- ensemblMeta$Gene.stable.ID[match(refFlat$V1, ensemblMeta$Transcript.stable.ID)]

#Check it merged correctly - it did
ensemblMeta$Transcript.stable.ID[ensemblMeta$Gene.stable.ID == "ENSG00000223972"]
ensemblMeta$Gene.stable.ID[ensemblMeta$Transcript.stable.ID == "ENST00000456328"]

#Rearrange as required: GENE_NAME, TRANSCRIPT_NAME, CHROMOSOME, STRAND, TX_START, TX_END, CDS_START, CDS_END, EXON_COUNT, EXON_STARTS, EXON_ENDS
#If additional columns are generated, they need to be removed, check manually 
refFlat_rearranged <- refFlat[,c(11,1:10)]

#Save as .txt without row/colnames and as tab-seperated 
write.table(refFlat_rearranged, file="Grch38-r98-refFlat-rearranged.txt", col.names = FALSE, row.names = FALSE, sep="\t", quote = FALSE)