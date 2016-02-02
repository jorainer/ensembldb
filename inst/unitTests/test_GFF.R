notrun_test_builds <- function(){
    input <- "/Users/jo/Projects/EnsDbs/83/Homo_sapiens.GRCh38.83.gtf.gz"
    fromGtf <- ensDbFromGtf(input, outfile=tempfile())
    ## provide wrong ensembl version
    fromGtf <- ensDbFromGtf(input, outfile=tempfile(), version="75")
    ## provide wrong genome version
    fromGtf <- ensDbFromGtf(input, outfile=tempfile(), genomeVersion="75")
    EnsDb(fromGtf)
    ## provide wrong organism
    fromGtf <- ensDbFromGtf(input, outfile=tempfile(), organism="blalba")
    EnsDb(fromGtf)
    ## GFF
    input <- "/Users/jo/Projects/EnsDbs/83/Homo_sapiens.GRCh38.83.chr.gff3.gz"
    fromGff <- ensDbFromGff(input, outfile=tempfile())
    EnsDb(fromGff)
    fromGff <- ensDbFromGff(input, outfile=tempfile(), version="75")
    EnsDb(fromGff)
    fromGff <- ensDbFromGff(input, outfile=tempfile(), genomeVersion="bla")
    EnsDb(fromGff)
    fromGff <- ensDbFromGff(input, outfile=tempfile(), organism="blabla")
    EnsDb(fromGff)

    ## AH
    library(AnnotationHub)
    ah <- AnnotationHub()
    fromAh <- ensDbFromAH(ah["AH47963"], outfile=tempfile())
    EnsDb(fromAH)
    fromAh <- ensDbFromAH(ah["AH47963"], outfile=tempfile(), version="75")
    EnsDb(fromAH)
    fromAh <- ensDbFromAH(ah["AH47963"], outfile=tempfile(), genomeVersion="bla")
    EnsDb(fromAH)
    fromAh <- ensDbFromAH(ah["AH47963"], outfile=tempfile(), organism="blabla")
    EnsDb(fromAH)
}



notrun_test_ensdbFromGFF <- function(){
    library(ensembldb)
    ##library(rtracklayer)
    ## VERSION 83
    gtf <- "/Users/jo/Projects/EnsDbs/83/Homo_sapiens.GRCh38.83.gtf.gz"
    fromGtf <- ensDbFromGtf(gtf, outfile=tempfile())
    egtf <- EnsDb(fromGtf)

    gff <- "/Users/jo/Projects/EnsDbs/83/Homo_sapiens.GRCh38.83.gff3.gz"
    fromGff <- ensDbFromGff(gff, outfile=tempfile())
    egff <- EnsDb(fromGff)

    ## Compare EnsDbs
    ensembldb:::compareEnsDbs(egtf, egff)
    ## OK, only Entrezgene ID "problems"

    ## Compare with the one built with the Perl API
    library(EnsDb.Hsapiens.v83)
    edb <- EnsDb.Hsapiens.v83

    ensembldb:::compareEnsDbs(egtf, edb)

    ensembldb:::compareEnsDbs(egff, edb)
    ## OK, I get different genes...
    genes1 <- genes(egtf)
    genes2 <- genes(edb)

    only2 <- genes2[!(genes2$gene_id %in% genes1$gene_id)]

    ## That below was before the fix to include feature type start_codon and stop_codon
    ## to the CDS type.
    ## Identify which are the different transcripts:
    txGtf <- transcripts(egtf)
    txGff <- transcripts(egff)
    commonIds <- intersect(names(txGtf), names(txGff))
    haveCds <- commonIds[!is.na(txGtf[commonIds]$tx_cds_seq_start) & !is.na(txGff[commonIds]$tx_cds_seq_start)]
    diffs <- haveCds[txGtf[haveCds]$tx_cds_seq_start != txGff[haveCds]$tx_cds_seq_start]
    head(diffs)

    ## What could be reasons?
    ## 1) alternative CDS?
    ## Checking the GTF:
    ## tx ENST00000623834: start_codon: 195409 195411.
    ##                     first CDS: 195259 195411.
    ##                     last CDS: 185220 185350.
    ##                     stop_codon: 185217 185219.
    ## So, why the heck is the stop codon OUTSIDE the CDS???
    ## library(rtracklayer)
    ## theGtf <- import(gtf, format="gtf")
    ## ## Apparently, the GTF contains the additional elements start_codon/stop_codon.
    ## theGff <- import(gff, format="gff3")


    ## transcripts(egtf, filter=TxidFilter(diffs[1]))
    ## transcripts(egff, filter=TxidFilter(diffs[1]))


    ## VERSION 81
    ## Try to get the same via AnnotationHub
    gff <- "/Users/jo/Projects/EnsDbs/81/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz"
    fromGff <- ensDbFromGff(gff, outfile=tempfile())
    egff <- EnsDb(fromGff)

    gtf <- "/Users/jo/Projects/EnsDbs/81/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz"
    fromGtf <- ensDbFromGtf(gtf, outfile=tempfile())
    egtf <- EnsDb(fromGtf)

    ## Compare those two:
    ensembldb:::compareEnsDbs(egff, egtf)
    ## Why are there some differences in the transcripts???
    trans1 <- transcripts(egff)
    trans2 <- transcripts(egtf)
    onlyInGtf <- trans2[!(trans2$tx_id %in% trans1$tx_id)]

    ##gtfGRanges <- ah["AH47963"]

    library(AnnotationHub)
    ah <- AnnotationHub()
    fromAh <- ensDbFromAH(ah["AH47963"], outfile=tempfile())  ## That's human...
    eah <- EnsDb(fromAh)

    ## Compare it to gtf:
    ensembldb:::compareEnsDbs(eah, egtf)
    ## OK. Same cds starts and cds ends.

    ## Compare it to gff:
    ensembldb:::compareEnsDbs(eah, egff)
    ## hm.

    ## Compare to EnsDb
    library(EnsDb.Hsapiens.v81)
    edb <- EnsDb.Hsapiens.v81
    ensembldb:::compareEnsDbs(edb, egtf)
    ## Problem with CDS
    ensembldb:::compareEnsDbs(edb, egff)
    ## That's fine.

    ## Summary:
    ## GTF and AH are the same.
    ## GFF and Perl API are the same.

    ## OLD STUFF BELOW.

    ##fromAh <- EnsDbFromAH(ah["AH47963"], outfile=tempfile(), organism="Homo sapiens", version=81)

    ## Try with a fancy species:
    gff <- "/Users/jo/Projects/EnsDbs/83/gadus_morhua/Gadus_morhua.gadMor1.83.gff3.gz"
    fromGtf <- ensDbFromGff(gff, outfile=tempfile())

    gff <- "/Users/jo/Projects/EnsDbs/83/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.83.gff3.gz"
    fromGff <- ensDbFromGff(gff, outfile=tempfile())
    ## That works.

    ## Try with a file from AnnotationHub: Gorilla gorilla.
    library(AnnotationHub)
    ah <- AnnotationHub()
    ah <- ah["AH47962"]

    res <- ensDbFromAH(ah, outfile=tempfile())
    edb <- EnsDb(res)
    genes(edb)


    ## ensRel <- query(ah, c("GTF", "ensembl"))

    ## gtf <- "/Users/jo/Projects/EnsDbs/83/Homo_sapiens.GRCh38.83.gtf.gz"
    ## ## GTF
    ## dir.create("/tmp/fromGtf")
    ## fromGtf <- ensDbFromGtf(gtf, path="/tmp/fromGtf", verbose=TRUE)
    ## ## GFF
    ## dir.create("/tmp/fromGff")
    ## fromGff <- ensembldb:::ensDbFromGff(gff, path="/tmp/fromGff", verbose=TRUE)

    ## ## ZBTB16:
    ## ## exon: ENSE00003606532 is 3rd exon of tx: ENST00000335953
    ## ## exon: ENSE00003606532 is 3rd exon of tx: ENST00000392996
    ## ## the Ensembl GFF has 2 entries for this exon.

}



