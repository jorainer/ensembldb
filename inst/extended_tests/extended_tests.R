## This script comprises extended tests.
##*****************************************************************
## Gviz stuff
notrun_test_genetrack_df <- function(){
    do.plot <- FALSE
    if(do.plot){
        ##library(Gviz)
        options(ucscChromosomeNames=FALSE)
        data(geneModels)
        geneModels$chromosome <- 7
        chr <- 7
        start <- min(geneModels$start)
        end <- max(geneModels$end)
        myGeneModels <- getGeneRegionTrackForGviz(edb, chromosome=chr,
                                                  start=start,
                                                  end=end)
        ## chromosome has to be the same....
        gtrack <- GenomeAxisTrack()
        gvizTrack <- GeneRegionTrack(geneModels, name="Gviz")
        ensdbTrack <- GeneRegionTrack(myGeneModels, name="ensdb")
        plotTracks(list(gtrack, gvizTrack, ensdbTrack))
        plotTracks(list(gtrack, gvizTrack, ensdbTrack), from=26700000,
                   to=26780000)
        ## Looks very nice...
    }
    ## Put the stuff below into the vignette:
    ## Next we get all lincRNAs on chromosome Y
    Lncs <- getGeneRegionTrackForGviz(edb,
                                      filter=list(SeqNameFilter("Y"),
                                                  GeneBiotypeFilter("lincRNA")))
    Prots <- getGeneRegionTrackForGviz(edb,
                                       filter=list(SeqNameFilter("Y"),
                                                   GeneBiotypeFilter("protein_coding")))
    if(do.plot){
        plotTracks(list(gtrack, GeneRegionTrack(Lncs, name="lincRNAs"),
                        GeneRegionTrack(Prots, name="proteins")))
        plotTracks(list(gtrack, GeneRegionTrack(Lncs, name="lincRNAs"),
                        GeneRegionTrack(Prots, name="proteins")),
                   from=5000000, to=7000000, transcriptAnnotation="symbol")
    }
    ## is that the same than:
    TestL <- getGeneRegionTrackForGviz(edb,
                                       filter=list(GeneBiotypeFilter("lincRNA")),
                                       chromosome="Y", start=5000000, end=7000000)
    TestP <- getGeneRegionTrackForGviz(edb,
                                       filter=list(GeneBiotypeFilter("protein_coding")),
                                       chromosome="Y", start=5000000, end=7000000)
    if(do.plot){
        plotTracks(list(gtrack, GeneRegionTrack(Lncs, name="lincRNAs"),
                        GeneRegionTrack(Prots, name="proteins"),
                        GeneRegionTrack(TestL, name="compareL"),
                        GeneRegionTrack(TestP, name="compareP")),
                   from=5000000, to=7000000, transcriptAnnotation="symbol")
    }
    expect_true(all(TestL$exon %in% Lncs$exon))
    expect_true(all(TestP$exon %in% Prots$exon))
    ## Crazy amazing stuff
    ## system.time(
    ##     All <- getGeneRegionTrackForGviz(edb)
    ## )
}



notrun_test_getSeqlengthsFromMysqlFolder <- function() {
    ## Test this for some more seqlengths.
    library(EnsDb.Rnorvegicus.v79)
    db <- EnsDb.Rnorvegicus.v79
    seq_info <- seqinfo(db)
    seq_lengths <- ensembldb:::.getSeqlengthsFromMysqlFolder(
        organism = "Rattus norvegicus", ensembl = 79,
        seqnames = seqlevels(seq_info))
    sl <- seqlengths(seq_info)
    sl_2 <- seq_lengths$length
    names(sl_2) <- rownames(seq_lengths)
    checkEquals(sl, sl_2)
    ## Mus musculus
}

notrun_test_ensDbFromGtf_Gff_AH <- function() {
    gtf <- paste0("/Users/jo/Projects/EnsDbs/80/caenorhabditis_elegans/",
                  "Caenorhabditis_elegans.WBcel235.80.gtf.gz")
    outf <- tempfile()
    db <- ensDbFromGtf(gtf = gtf, outfile = outf)
    ## use Gff
    gff <- paste0("/Users/jo/Projects/EnsDbs/84/canis_familiaris/gff3/",
                  "Canis_familiaris.CanFam3.1.84.gff3.gz")
    outf <- tempfile()
    db <- ensDbFromGff(gff, outfile = outf)

    ## Checking one from ensemblgenomes:
    gtf <- paste0("/Users/jo/Projects/EnsDbs/ensemblgenomes/30/",
                  "solanum_lycopersicum/",
                  "Solanum_lycopersicum.GCA_000188115.2.30.chr.gtf.gz"
                  )
    outf <- tempfile()
    db <- ensDbFromGtf(gtf = gtf, outfile = outf)
    gtf <- paste0("/Users/jo/Projects/EnsDbs/ensemblgenomes/30/",
                  "solanum_lycopersicum/",
                  "Solanum_lycopersicum.GCA_000188115.2.30.gtf.gz"
                  )
    outf <- tempfile()
    db <- ensDbFromGtf(gtf = gtf, outfile = outf)

    ## AH
    library(AnnotationHub)
    ah <- AnnotationHub()
    query(ah, c("release-83", "gtf"))
    ah_1 <- ah["AH50418"]
    db <- ensDbFromAH(ah_1, outfile = outf)
    ah_2 <- ah["AH50352"]
    db <- ensDbFromAH(ah_2, outfile = outf)
}

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
    db <- EnsDb.Hsapiens.v83

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


    ## transcripts(egtf, filter=TxIdFilter(diffs[1]))
    ## transcripts(egff, filter=TxIdFilter(diffs[1]))


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



############################################################
## Can not perform these tests right away, as they require a
## working MySQL connection.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

dontrun_test_useMySQL <- function() {
    edb_mysql <- useMySQL(edb, user = "anonuser", host = "localhost", pass = "")
}

dontrun_test_connect_EnsDb <- function() {
    library(RMySQL)
    con <- dbConnect(MySQL(), user = "anonuser", host = "localhost", pass = "")

    ensembldb:::listEnsDbs(dbcon = con)
    ## just with user.
    ensembldb:::listEnsDbs(user = "anonuser", host = "localhost", pass = "",
                           port = 3306)

    ## Connecting directly to a EnsDb MySQL database.
    con <- dbConnect(MySQL(), user = "anonuser", host = "localhost", pass = "",
                     dbname = "ensdb_hsapiens_v75")
    edb_mysql <- EnsDb(con)
}

notrun_compareEnsDbs <- function() {
    res <- ensembldb:::compareEnsDbs(edb, edb)
}

############################################################
## Massive test validating the cds:
## compare the length of the CDS with the length of the encoded protein.
## Get the CDS sequence, translate that and compare to protein sequence.
notrun_massive_cds_test <- function() {
    ## Get all CDS:
    tx_cds <- cdsBy(edb, by = "tx", filter = SeqNameFilter(c(1:22, "X", "Y")))
    prots <- proteins(edb, filter = TxIdFilter(names(tx_cds)),
                      return.type = "AAStringSet")
    checkTrue(all(names(tx_cds) %in% mcols(prots)$tx_id))
    tx_cds <- tx_cds[mcols(prots)$tx_id]
    ## Check that the length of the protein sequence is length of CDS/3
    diff_width <- sum(width(tx_cds)) != width(prots) * 3
    ## Why??? I've got some many differences here???
    sum(diff_width)
    ## Check some of them manually in Ensembl

    ## 1st: - strand.
    tx_1 <- tx_cds[diff_width][1]
    ## Protein: 245aa
    prots[diff_width][1]
    ## OK.
    ## Tx 2206bp:
    exns <- exonsBy(edb, filter = TxIdFilter(names(tx_1)))
    sum(width(exns))
    ## OK.
    ## Now to the CDS:
    cds_ex1 <- "ATGGCGTCCCCGTCTCGGAGACTGCAGACTAAACCAGTCATTACTTGTTTCAAGAGCGTTCTGCTAATCTACACTTTTATTTTCTGG"
    cds_ex2 <- "ATCACTGGCGTTATCCTTCTTGCAGTTGGCATTTGGGGCAAGGTGAGCCTGGAGAATTACTTTTCTCTTTTAAATGAGAAGGCCACCAATGTCCCCTTCGTGCTCATTGCTACTGGTACCGTCATTATTCTTTTGGGCACCTTTGGTTGTTTTGCTACCTGCCGAGCTTCTGCATGGATGCTAAAACTG"
    cds_ex3 <- "TATGCAATGTTTCTGACTCTCGTTTTTTTGGTCGAACTGGTCGCTGCCATCGTAGGATTTGTTTTCAGACATGAG"
    cds_ex4 <- "ATTAAGAACAGCTTTAAGAATAATTATGAGAAGGCTTTGAAGCAGTATAACTCTACAGGAGATTATAGAAGCCATGCAGTAGACAAGATCCAAAATACG"
    cds_ex5 <- "TTGCATTGTTGTGGTGTCACCGATTATAGAGATTGGACAGATACTAATTATTACTCAGAAAAAGGATTTCCTAAGAGTTGCTGTAAACTTGAAGATTGTACTCCACAGAGAGATGCAGACAAAGTAAACAATGAA"
    cds_ex6 <- "GGTTGTTTTATAAAGGTGATGACCATTATAGAGTCAGAAATGGGAGTCGTTGCAGGAATTTCCTTTGGAGTTGCTTGCTTCCAA"
    cds_ex7 <- "CTGATTGGAATCTTTCTCGCCTACTGCCTCTCTCGTGCCATAACAAATAACCAGTATGAGATAGTGTAA"
    cds_seq <- c(cds_ex1, cds_ex2, cds_ex3, cds_ex4, cds_ex5, cds_ex6, cds_ex7)
    nchar(cds_seq)
    width(tx_1)
    ## The length should be identical:
    checkEquals(sum(nchar(cds_seq)), sum(width(tx_1)), checkNames = FALSE)
    ## OK; so WHAT???
    sum(width(tx_1)) / 3
    ## So, start codon is encoded into a methionine.
    ## Stop codon is either a TAA, TGA or TAG. UAG can be encoded into Sec (U), UAG into Pyl (O)
    library(Biostrings)
    dna_s <- DNAString(paste0(cds_ex1, cds_ex2, cds_ex3, cds_ex4, cds_ex5, cds_ex6, cds_ex7))
    translate(dna_s)
    ## Look at that!
    translate(DNAString("TAA")) ## -> translates into *
    translate(DNAString("TGA")) ## -> translates into *
    translate(DNAString("TAG")) ## -> translates into *
    translate(DNAString("ATG")) ## -> translates into M

    ## Assumption:
    ## If the mRNA ends with a TAA, the protein sequence will be 1aa shorter than
    ## length(CDS)/3.
    ## If the mRNA ends with a TAG, UAG the AA length is length(CDS)/3

    ## Check one of the mRNA where it fits:
    tx_2 <- tx_cds[!diff_width][1]
    prots[!diff_width][1]
    ## AA is 137 long, ends with I.
    sum(width(tx_2)) / 3  ## OK
    ## Check Ensembl:
    tx_2_1 <- "ATGCTAAAACTG"
    tx_2_2 <- "TATGCAATGTTTCTGACTCTCGTTTTTTTGGTCGAACTGGTCGCTGCCATCGTAGGATTTGTTTTCAGACATGAG"
    tx_2_3 <- "ATTAAGAACAGCTTTAAGAATAATTATGAGAAGGCTTTGAAGCAGTATAACTCTACAGGAGATTATAGAAGCCATGCAGTAGACAAGATCCAAAATACG"
    tx_2_4 <- "TTGCATTGTTGTGGTGTCACCGATTATAGAGATTGGACAGATACTAATTATTACTCAGAAAAAGGATTTCCTAAGAGTTGCTGTAAACTTGAAGATTGTACTCCACAGAGAGATGCAGACAAAGTAAACAATGAA"
    tx_2_5 <- "GGTTGTTTTATAAAGGTGATGACCATTATAGAGTCAGAAATGGGAGTCGTTGCAGGAATTTCCTTTGGAGTTGCTTGCTTCCAA"
    tx_2_6 <- "GACATT"
    tx_2_cds <- paste0(tx_2_1, tx_2_2, tx_2_3, tx_2_4, tx_2_5, tx_2_6)
    nchar(tx_2_cds)
    sum(width(tx_2))
    ## OK.
    translate(DNAString(tx_2_cds))

    ## Next assumption:
    ## If we don't have a 3' UTR the AA sequence corresponds to length(CDS)/3
    tx_cds <- cdsBy(edb, by = "tx", filter = SeqNameFilter(c(1:22, "X", "Y")),
                    columns = c("tx_seq_start", "tx_seq_end", "tx_cds_seq_start",
                                "tx_cds_seq_end"))
    prots <- proteins(edb, filter = TxIdFilter(names(tx_cds)),
                      return.type = "AAStringSet")
    checkTrue(all(names(tx_cds) %in% mcols(prots)$tx_id))
    tx_cds <- tx_cds[mcols(prots)$tx_id]
    ## Calculate the CDS width.
    tx_cds_width <- sum(width(tx_cds))
    txs <- transcripts(edb, filter = TxIdFilter(names(tx_cds)))
    txs <- txs[names(tx_cds)]
    ## Subtract 3 from the width if we've got an 3'UTR.
    to_subtract <- rep(3, length(tx_cds_width))
    to_subtract[((end(txs) == txs$tx_cds_seq_end) &
                 as.logical(strand(txs) == "+"))
                | ((start(txs) == txs$tx_cds_seq_start)
                    & as.logical(strand(txs) == "-"))] <- 0
    tx_cds_width <- tx_cds_width - to_subtract
    ## Check that the length of the protein sequence is length of CDS/3
    diff_width <- tx_cds_width != width(prots) * 3
    ## Why??? I've got some many differences here???
    sum(diff_width)
    length(diff_width)
    ## AAAA, still have some that don't fit!!!
    tx_3 <- tx_cds[diff_width][1]
    prots[diff_width][1]
    ## AA is 259aa long, ends with T., TX is: ENST00000371584
    ## WTF, we've got no START CODON!!!
    sum(width(tx_3)) / 3 ## OMG!!!

    ## Now, exclude those without a 5' UTR:
    no_five <- ((start(txs) == txs$tx_cds_seq_start) &
                as.logical(strand(txs) == "+")) |
        ((end(txs) == txs$tx_cds_seq_end) &
         as.logical(strand(txs) == "-"))
    still_prot <- (diff_width & !no_five)
}

notrun_test_getGenomeFaFile <- function(){
    library(EnsDb.Hsapiens.v82)
    edb <- EnsDb.Hsapiens.v82

    ## We know that there is no Fasta file for that Ensembl release available.
    Fa <- getGenomeFaFile(edb)
    ## Got the one from Ensembl 81.
    genes <- genes(edb, filter=SeqNameFilter("Y"))
    geneSeqsFa <- getSeq(Fa, genes)
    ## Get the transcript sequences...
    txSeqsFa <- extractTranscriptSeqs(Fa, edb, filter=SeqNameFilter("Y"))

    ## Get the TwoBitFile.
    twob <- ensembldb:::getGenomeTwoBitFile(edb)
    ## Get thegene sequences.
    ## ERROR FIX BELOW WITH UPDATED VERSIONS!!!
    geneSeqs2b <- getSeq(twob, genes)

    ## Have to fix the seqnames.
    si <- seqinfo(twob)
    sn <- unlist(lapply(strsplit(seqnames(si), split=" ", fixed=TRUE), function(z){
        return(z[1])
    }))
    seqnames(si) <- sn
    seqinfo(twob) <- si

    ## Do the same with the TwoBitFile
    geneSeqsTB <- getSeq(twob, genes)

    ## Subset to all genes that are encoded on chromosomes for which
    ## we do have DNA sequence available.
    genes <- genes[seqnames(genes) %in% seqnames(seqinfo(Dna))]

    ## Get the gene sequences, i.e. the sequence including the sequence of
    ## all of the gene's exons and introns.
    geneSeqs <- getSeq(Dna, genes)

    library(AnnotationHub)
    ah <- AnnotationHub()
    quer <- query(ah, c("release-", "Homo sapiens"))
    ## So, I get 2bit files and toplevel stuff.
    Test <- ah[["AH50068"]]

}



notrun_test_extractTranscriptSeqs <- function(){
    ## Note: we can't run that by default as we can not assume everybody has
    ## AnnotationHub and the required ressource installed.
    ## That's how we want to test the transcript seqs.
    genome <- getGenomeFaFile(edb)
    ZBTB <- extractTranscriptSeqs(genome, edb, filter=GenenameFilter("ZBTB16"))
    ## Load the sequences for one ZBTB16 transcript from FA.
    faf <- system.file("txt/ENST00000335953.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000335953"
    ## cDNA
    checkEquals(unname(as.character(ZBTB[tx])),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxIdFilter(tx))
    CDS <- extractTranscriptSeqs(genome, cBy)
    checkEquals(unname(as.character(CDS)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(genome, fBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(genome, tBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))


    ## Another gene on the reverse strand:
    faf <- system.file("txt/ENST00000200135.fa.gz", package="ensembldb")
    Seqs <- readDNAStringSet(faf)
    tx <- "ENST00000200135"
    ## cDNA
    cDNA <- extractTranscriptSeqs(genome, edb, filter=TxIdFilter(tx))
    checkEquals(unname(as.character(cDNA)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## do the same, but from other strand
    exns <- exonsBy(edb, "tx", filter=TxIdFilter(tx))
    cDNA <- extractTranscriptSeqs(genome, exns)
    checkEquals(unname(as.character(cDNA)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    strand(exns) <- "+"
    cDNA <- extractTranscriptSeqs(genome, exns)
    checkTrue(unname(as.character(cDNA)) !=
              unname(as.character(Seqs[grep(names(Seqs), pattern="cdna")])))
    ## CDS
    cBy <- cdsBy(edb, "tx", filter=TxIdFilter(tx))
    CDS <- extractTranscriptSeqs(genome, cBy)
    checkEquals(unname(as.character(CDS)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="cds")])))
    ## 5' UTR
    fBy <- fiveUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(genome, fBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr5")])))
    ## 3' UTR
    tBy <- threeUTRsByTranscript(edb, filter=TxIdFilter(tx))
    UTR <- extractTranscriptSeqs(genome, tBy)
    checkEquals(unname(as.character(UTR)),
                unname(as.character(Seqs[grep(names(Seqs), pattern="utr3")])))
}

notrun_test_getCdsSequence <- function(){
    ## That's when we like to get the sequence from the coding region.
    genome <- getGenomeFaFile(edb)
    tx <- extractTranscriptSeqs(genome, edb, filter=SeqNameFilter("Y"))
    cdsSeq <- extractTranscriptSeqs(genome, cdsBy(edb, filter=SeqNameFilter("Y")))
    ## that's basically to get the CDS sequence.
    ## UTR sequence:
    tutr <- extractTranscriptSeqs(genome, threeUTRsByTranscript(edb, filter=SeqNameFilter("Y")))
    futr <- extractTranscriptSeqs(genome, fiveUTRsByTranscript(edb, filter=SeqNameFilter("Y")))
    theTx <- "ENST00000602770"
    fullSeq <- as.character(tx[theTx])
    ## build the one from 5', cds and 3'
    compSeq <- ""
    if(any(names(futr) == theTx))
        compSeq <- paste0(compSeq, as.character(futr[theTx]))
    if(any(names(cdsSeq) == theTx))
        compSeq <- paste0(compSeq, as.character(cdsSeq[theTx]))
    if(any(names(tutr) == theTx))
        compSeq <- paste(compSeq, as.character(tutr[theTx]))
    checkEquals(unname(fullSeq), compSeq)
}

notrun_test_cds <- function(){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    cds <- cds(txdb)
    cby <- cdsBy(txdb, by="tx")

    gr <- cby[[7]][1]
    seqlevels(gr) <- sub(seqlevels(gr), pattern="chr", replacement="")
    tx <- transcripts(edb, filter=GRangesFilter(gr, condition="overlapping"))
    cby[[7]]

    ## Note: so that fits! And we've to include the stop_codon feature for GTF import!
    ## Make an TxDb from GTF:
    gtf <- "/Users/jo/Projects/EnsDbs/75/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
    library(GenomicFeatures)
    Test <- makeTxDbFromGFF(gtf, format="gtf", organism="Homo sapiens")
    scds <- cdsBy(Test, by="tx")
    gr <- scds[[7]][1]
    tx <- transcripts(edb, filter=GRangesFilter(gr, condition="overlapping"))
    scds[[7]]
    ## Compare:
    ## TxDb form GTF has: 865692 879533
    ## EnsDb: 865692 879533

    ## Next test:
    gr <- scds[[2]][1]
    tx <- transcripts(edb, filter=GRangesFilter(gr, condition="overlapping"))
    tx
    scds[[2]]
    ## start_codon: 367659 367661, stop_codon: 368595 368597 CDS: 367659 368594.
    ## TxDb from GTF includes the stop_codon!
}


dontrun_benchmark_ordering_genes <- function() {
    .withR <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- TRUE
        genes(x, ...)
    }
    .withSQL <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- FALSE
        genes(x, ...)
    }
    library(microbenchmark)
    microbenchmark(.withR(edb), .withSQL(edb), times = 10)  ## same
    microbenchmark(.withR(edb, columns = c("gene_id", "tx_id")),
                   .withSQL(edb, columns = c("gene_id", "tx_id")),
                   times = 10)  ## R slightly faster.
    microbenchmark(.withR(edb, columns = c("gene_id", "tx_id"),
                          SeqNameFilter("Y")),
                   .withSQL(edb, columns = c("gene_id", "tx_id"),
                            SeqNameFilter("Y")),
                   times = 10)  ## same.
}

## We aim to fix issue #11 by performing the ordering in R instead
## of SQL. Thus, we don't want to run this as a "regular" test
## case.
dontrun_test_ordering_cdsBy <- function() {
    doBench <- FALSE
    if (doBench)
        library(microbenchmark)
    .withR <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- TRUE
        cdsBy(x, ...)
    }
    .withSQL <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- FALSE
        cdsBy(x, ...)
    }
    res_sql <- .withSQL(edb)
    res_r <- .withR(edb)
    checkEquals(res_sql, res_r)
    if (dobench)
        microbenchmark(.withSQL(edb), .withR(edb),
                       times = 3)  ## R slightly faster.
    res_sql <- .withSQL(edb, filter = SeqNameFilter("Y"))
    res_r <- .withR(edb, filter = SeqNameFilter("Y"))
    checkEquals(res_sql, res_r)
    if (dobench)
        microbenchmark(.withSQL(edb, filter = SeqNameFilter("Y")),
                       .withR(edb, filter = SeqNameFilter("Y")),
                       times = 10)  ## R 6x faster.
}

dontrun_test_ordering_exonsBy <- function() {
    doBench <- FALSE
    if (doBench)
        library(microbenchmark)
    .withR <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- TRUE
        exonsBy(x, ...)
    }
    .withSQL <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- FALSE
        exonsBy(x, ...)
    }
    res_sql <- .withSQL(edb)
    res_r <- .withR(edb)
    checkEquals(res_sql, res_r)
    if (doBench)
        microbenchmark(.withSQL(edb), .withR(edb),
                       times = 3)  ## about the same; R slightly faster.
    ## with using a SeqNameFilter in addition.
    res_sql <- .withSQL(edb, filter = SeqNameFilter("Y"))
    res_r <- .withR(edb, filter = SeqNameFilter("Y")) ## query takes longer.
    checkEquals(res_sql, res_r)
    if (doBench)
        microbenchmark(.withSQL(edb, filter = SeqNameFilter("Y")),
                       .withR(edb, filter = SeqNameFilter("Y")),
                       times = 3)  ## SQL twice as fast.
    ## Now getting stuff by gene
    res_sql <- .withSQL(edb, by = "gene")
    res_r <- .withR(edb, by = "gene")
    ## checkEquals(res_sql, res_r) ## Differences due to ties
    if (doBench)
        microbenchmark(.withSQL(edb, by = "gene"),
                       .withR(edb, by = "gene"),
                       times = 3)  ## SQL faster; ???
    ## Along with a SeqNameFilter
    res_sql <- .withSQL(edb, by = "gene", filter = SeqNameFilter("Y"))
    res_r <- .withR(edb, by = "gene", filter = SeqNameFilter("Y"))
    ## Why does the query take longer for R???
    ## checkEquals(res_sql, res_r) ## Differences due to ties
    if (doBench)
        microbenchmark(.withSQL(edb, by = "gene", filter = SeqNameFilter("Y")),
                       .withR(edb, by = "gene", filter = SeqNameFilter("Y")),
                       times = 3)  ## SQL faster.
    ## Along with a GeneBiotypeFilter
    if (doBench)
        microbenchmark(.withSQL(edb, by = "gene", filter = GeneBiotypeFilter("protein_coding"))
                     , .withR(edb, by = "gene", filter = GeneBiotypeFilter("protein_coding"))
                     , times = 3)
}

dontrun_test_ordering_transcriptsBy <- function() {
    .withR <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- TRUE
        transcriptsBy(x, ...)
    }
    .withSQL <- function(x, ...) {
        ensembldb:::orderResultsInR(x) <- FALSE
        transcriptsBy(x, ...)
    }
    res_sql <- .withSQL(edb)
    res_r <- .withR(edb)
    checkEquals(res_sql, res_r)
    microbenchmark(.withSQL(edb), .withR(edb), times = 3) ## same speed

    res_sql <- .withSQL(edb, filter = SeqNameFilter("Y"))
    res_r <- .withR(edb, filter = SeqNameFilter("Y"))
    checkEquals(res_sql, res_r)
    microbenchmark(.withSQL(edb, filter = SeqNameFilter("Y")),
                   .withR(edb, filter = SeqNameFilter("Y")),
                   times = 3) ## SQL slighly faster.
}

dontrun_query_tune <- function() {
    ## Query tuning:
    library(RSQLite)
    con <- dbconn(edb)

    Q <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from gene join tx on (gene.gene_id=tx.gene_id) join tx2exon on (tx.tx_id=tx2exon.tx_id) join exon on (tx2exon.exon_id=exon.exon_id) where gene.seq_name = 'Y'"
    system.time(dbGetQuery(con, Q))

    Q2 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from exon join tx2exon on (tx2exon.exon_id = exon.exon_id) join tx on (tx2exon.tx_id = tx.tx_id) join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y'"
    system.time(dbGetQuery(con, Q2))

    Q3 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from tx2exon join exon on (tx2exon.exon_id = exon.exon_id) join tx on (tx2exon.tx_id = tx.tx_id) join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y'"
    system.time(dbGetQuery(con, Q3))

    Q4 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from tx2exon join exon on (tx2exon.exon_id = exon.exon_id) join tx on (tx2exon.tx_id = tx.tx_id) join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y' order by tx.tx_id"
    system.time(dbGetQuery(con, Q4))

    Q5 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from tx2exon inner join exon on (tx2exon.exon_id = exon.exon_id) inner join tx on (tx2exon.tx_id = tx.tx_id) inner join gene on (gene.gene_id=tx.gene_id) where gene.seq_name = 'Y' order by tx.tx_id"
    system.time(dbGetQuery(con, Q5))

    Q6 <- "select distinct tx2exon.exon_id,exon.exon_seq_start,exon.exon_seq_end,gene.seq_name,tx2exon.tx_id,gene.seq_strand,tx2exon.exon_idx from gene inner join tx on (gene.gene_id=tx.gene_id) inner join tx2exon on (tx.tx_id=tx2exon.tx_id) inner join exon on (tx2exon.exon_id=exon.exon_id) where gene.seq_name = 'Y' order by tx.tx_id asc"
    system.time(dbGetQuery(con, Q6))
}

## implement:
## .checkOrderBy: checks order.by argument removing columns that are
## not present in the database
## orderBy columns are added to the columns.
## .orderDataFrameBy: orders the dataframe by the specified columns.

notrun_test_protein_domains <- function() {
    res <- ensembldb:::getWhat(edb, columns = c("protein_id", "tx_id", "gene_id",
                                                "gene_name"),
                               filter = list(ProtDomIdFilter("PF00096")))
}

notrun_compare_full <- function(){
    ## That's on the full thing.
    ## Test if the result has the same ordering than the transcripts call.
    allTx <- transcripts(edb)
    txLen <- transcriptLengths(edb, with.cds_len=TRUE, with.utr5_len=TRUE,
                               with.utr3_len=TRUE)
    checkEquals(names(allTx), rownames(txLen))
    system.time(
        futr <- fiveUTRsByTranscript(edb)
    )
    ## 23 secs.
    futrLen <- sum(width(futr))  ## do I need reduce???
    checkEquals(unname(futrLen), txLen[names(futrLen), "utr5_len"])
    ## 3'
    system.time(
        tutr <- threeUTRsByTranscript(edb)
    )
    system.time(
        tutrLen <- sum(width(tutr))
    )
    checkEquals(unname(tutrLen), txLen[names(tutrLen), "utr3_len"])
}

notrun_compare_to_genfeat <- function(){
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    system.time(
        Len <- transcriptLengths(edb)
    )
    ## Woa, 52 sec
    system.time(
        txLen <- lengthOf(edb, "tx")
    )
    ## Faster, 31 sec
    checkEquals(Len$tx_len, unname(txLen[rownames(Len)]))
    system.time(
        Len2 <- transcriptLengths(txdb)
    )
    ## :) 2.5 sec.
    ## Next.
    system.time(
        Len <- transcriptLengths(edb, with.cds_len = TRUE)
    )
    ## 56 sec
    system.time(
        Len2 <- transcriptLengths(txdb, with.cds_len=TRUE)
    )
    ## 4 sec.

    ## Calling the transcriptLengths of GenomicFeatures on the EnsDb.
    system.time(
        Def <- GenomicFeatures::transcriptLengths(edb)
    ) ## 26.5 sec

    system.time(
        WithCds <- GenomicFeatures::transcriptLengths(edb, with.cds_len=TRUE)
    ) ## 55 sec

    system.time(
        WithAll <- GenomicFeatures::transcriptLengths(edb, with.cds_len=TRUE,
                                                      with.utr5_len=TRUE,
                                                      with.utr3_len=TRUE)
    ) ## 99 secs

    ## Get my versions...
    system.time(
        MyDef <- ensembldb:::.transcriptLengths(edb)
    ) ## 31 sec
    system.time(
        MyWithCds <- ensembldb:::.transcriptLengths(edb, with.cds_len=TRUE)
    ) ## 44 sec
    system.time(
        MyWithAll <- ensembldb:::.transcriptLengths(edb, with.cds_len=TRUE,
                                                    with.utr5_len=TRUE,
                                                    with.utr3_len=TRUE)
    ) ## 63 sec

    ## Should be all the same!!!
    rownames(MyDef) <- NULL
    checkEquals(Def, MyDef)
    ##
    rownames(MyWithCds) <- NULL
    MyWithCds[is.na(MyWithCds$cds_len), "cds_len"] <- 0
    checkEquals(WithCds, MyWithCds)
    ##
    rownames(MyWithAll) <- NULL
    MyWithAll[is.na(MyWithAll$cds_len), "cds_len"] <- 0
    MyWithAll[is.na(MyWithAll$utr3_len), "utr3_len"] <- 0
    MyWithAll[is.na(MyWithAll$utr5_len), "utr5_len"] <- 0
    checkEquals(WithAll, MyWithAll)
}
