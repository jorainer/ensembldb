library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

notrun_test_getGenomeFaFile <- function(){
    library(EnsDb.Hsapiens.v82)
    edb <- EnsDb.Hsapiens.v82

    ## We know that there is no Fasta file for that Ensembl release available.
    Fa <- getGenomeFaFile(edb)
    ## Got the one from Ensembl 81.
    genes <- genes(edb, filter=SeqnameFilter("Y"))
    geneSeqsFa <- getSeq(Fa, genes)
    ## Get the transcript sequences...
    txSeqsFa <- extractTranscriptSeqs(Fa, edb, filter=SeqnameFilter("Y"))

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


