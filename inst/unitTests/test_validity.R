library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_validity_functions <- function() {
    OK <- ensembldb:::dbHasRequiredTables(dbconn(edb))
    checkTrue(OK)
    ## Check the tables
    OK <- ensembldb:::dbHasValidTables(dbconn(edb))
    checkTrue(OK)
}

test_validateEnsDb <- function() {
    checkTrue(ensembldb:::validateEnsDb(edb))
}

test_compareProteins <- function() {
    if (hasProteinData(edb)) {
        res <- ensembldb:::compareProteins(edb, edb)
        checkEquals(res, "OK")
    }
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

