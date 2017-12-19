#' @include proteinToX.R

#' @title Map transcript-relative coordinates to amino acid residues of the
#'     encoded protein
#'
#' @description
#'
#' `transcriptToProtein` maps within-transcript coordinates to the corresponding
#' coordinates within the encoded protein sequence. The provided coordinates
#' have to be within the coding region of the transcript (excluding the stop
#' codon) but are supposed to be relative to the first nucleotide of the
#' transcript (including the 5' UTR).
#'
#' @details
#'
#' Transcript-relative coordinates are mapped to the amino acid residues they
#' encode. As an example, positions within the transcript that corrspond to
#' nucleotides 1 to 3 in the CDS are mapped to the first position in the
#' protein sequence (see examples for more details).
#'
#' @param x `IRanges` with the coordinates within the transcript. Coordinates
#'     are counted from the start of the transcript (including the 5' UTR). The
#'     Ensembl IDs of the corresponding transcripts have to be provided either
#'     as `names` of the `IRanges`, or in one of its metadata columns.
#'
#' @param db `EnsDb` object.
#'
#' @param id `character(1)` specifying where the transcript identifier can be
#'     found. Has to be either `"name"` or one of `colnames(mcols(prng))`.
#' 
#' @return
#'
#' `IRanges` with the same length (and order) than the input `IRanges`
#' `x`. Each element in `IRanges` provides the coordinates within the
#' protein sequence, names being the (Ensembl) IDs of the protein. The
#' original transcript ID and the transcript-relative coordinates are provided
#' as metadata columns. Metadata columns `"cds_ok"` indicates whether the
#' length of the transcript's CDS matches the length of the encoded protein.
#' `IRanges` with a start coordinate of `-1` is returned for transcript
#' coordinates that can not be mapped to protein-relative coordinates
#' (either no transcript was found for the provided ID, the transcript
#' does not encode a protein or the provided coordinates are not within
#' the coding region of the transcript).
#'
#' @family coordinate mapping functions
#'
#' @author Johannes Rainer
#' 
#' @md
#'
#' @examples
#' 
#' library(EnsDb.Hsapiens.v75)
#' ## Restrict all further queries to chromosome x to speed up the examples
#' edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "X")
#'
#' ## Define an IRanges with the positions of the first 2 nucleotides of the
#' ## coding region for the transcript ENST00000381578
#' txpos <- IRanges(start = 692, width = 2, names = "ENST00000381578")
#'
#' ## Map these to the corresponding residues in the protein sequence
#' ## The protein-relative coordinates are returned as an `IRanges` object,
#' ## with the original, transcript-relative coordinates provided in metadata
#' ## columns tx_start and tx_end
#' transcriptToProtein(txpos, edbx)
#'
#' ## We can also map multiple ranges. Note that for any of the 3 nucleotides
#' ## encoding the same amino acid the position of this residue in the
#' ## protein sequence is returned. To illustrate this we map below each of the
#' ## first 4 nucleotides of the CDS to the corresponding position within the
#' ## protein.
#' txpos <- IRanges(start = c(692, 693, 694, 695),
#'     width = rep(1, 4), names = rep("ENST00000381578", 4))
#' transcriptToProtein(txpos, edbx)
#'
#' ## If the mapping fails, an IRanges with negative start position is returned.
#' ## Mapping can fail (as below) because the ID is not known.
#' transcriptToProtein(IRanges(1, 1, names = "unknown"), edbx)
#'
#' ## Or because the provided coordinates are not within the CDS
#' transcriptToProtein(IRanges(1, 1, names = "ENST00000381578"), edbx)
transcriptToProtein <- function(x, db, id = "name") {
    if (missing(x) || !is(x, "IRanges"))
        stop("Argument 'x' is required and has to be an 'IRanges' object")
    if (missing(db) || !is(db, "EnsDb"))
        stop("Argument 'db' is required and has to be an 'EnsDb' object")
    .txs_to_proteins(x, db = db, id = id)
}

#' This version performs the mapping for each IRange separately! Results are
#' the same as in .txs_to_protins. While the code is easier to read and to
#' debug, it is also considerably slower.
#'
#' @noRd
.tx_to_protein <- function(x, db, id = "name") {
    if (length(x) > 1)
        return(unlist(IRangesList(
            lapply(split(x, 1:length(x)), FUN = .tx_to_protein,
                   db = db, id = id)), use.names = FALSE))
    id <- match.arg(id, c("name", colnames(mcols(x))))
    if (id == "name")
        ids <- names(x)
    else ids <- mcols(x)[, id]
    if (any(is.null(ids)))
        stop("All IDs are NULL")
    empty_ranges <- IRanges(start = -1, width = 1)
    mcols(empty_ranges) <- DataFrame(tx_id = ids,
                                     tx_start = start(x),
                                     tx_end = end(x),
                                     cds_ok = NA)
    tx_lens <- .transcriptLengths(db,
                                  filter = TxIdFilter(ids),
                                  with.cds_len = TRUE,
                                  with.utr5_len = TRUE,
                                  with.utr3_len = TRUE)
    ## Does the tx exist?
    if (nrow(tx_lens) == 0) {
        warning("Transcript '", ids, "' can not be found", call. = FALSE)
        return(empty_ranges)
    }
    ## Is it a protein coding tx?
    if (is.na(tx_lens$cds_len)) {
        warning("Transcript '", ids, "' does not encode a protein",
                call. = FALSE)
        return(empty_ranges)
    }
    if (is.na(tx_lens$utr5_len))
        tx_lens$utr5_len <- 0
    if (is.na(tx_lens$utr3_len))
        tx_lens$utr3_len <- 0
    ## Are the coordinates within the CDS? Exclude the stop codon!
    if (start(x) <= tx_lens$utr5_len |
        end(x) > (tx_lens$utr5_len + tx_lens$cds_len - 3)) {
        warning("The provided coordinates for '", ids,
                "' are not within the coding region", call. = FALSE)
        return(empty_ranges)
    }
    ## Define the coordinates
    end(empty_ranges) <- ceiling((end(x) - tx_lens$utr5_len) / 3)
    start(empty_ranges) <- ceiling((start(x) - tx_lens$utr5_len) / 3)
    ## Now check if the CDS length matches the protein sequence length.
    prt <- proteins(db, filter = TxIdFilter(ids), columns = "protein_sequence")
    names(empty_ranges) <- prt$protein_id
    cds_ok <- nchar(prt$protein_sequence) * 3 == (tx_lens$cds_len - 3)
    if (!cds_ok)
        warning("The CDS of '", ids, "' does not match the length of the ",
                "encoded protein. Returned protein coordinates might not be",
                " correct", call. = FALSE)
    mcols(empty_ranges)$cds_ok <- cds_ok
    empty_ranges
}

#' @details
#'
#' Transcript-relative coordinates are mapped to the amino acid residues they
#' encode. As an example, positions within the transcript that corrspond to
#' nucleotides 1 to 3 in the CDS are mapped to the first position in the
#' protein sequence (see examples for more details).
#' (e.g. a transcript relative position corresponding to the first nucleotide
#' of the transcript's CDS is mapped to the first position in the protein
#' sequence, same as if the transcri(no matter which transcript relative
#' position corresponding to the first,
#' second or third nucleotide of the CDS are provided, all (any) of  
#'
#' @author Johannes Rainer
#' 
#' @md
#'
#' @noRd
.txs_to_proteins <- function(x, db, id = "name") {
    ## Fetch all in one.
    id <- match.arg(id, c("name", colnames(mcols(x))))
    if (id == "name")
        ids <- names(x)
    else ids <- mcols(x)[, id]
    if (any(is.null(ids)))
        stop("One or more of the provided IDs are NULL", call. = FALSE)
    names(x) <- ids
    ## Define internal IDs - we'll need them to return the result in the
    ## correct order.
    internal_ids <- paste0(names(x), ":", start(x), ":", end(x))
    tx_lens <- ensembldb:::.transcriptLengths(db,
                                  filter = TxIdFilter(unique(ids)),
                                  with.cds_len = TRUE,
                                  with.utr5_len = TRUE,
                                  with.utr3_len = TRUE)
    ## Process transcripts not found
    not_found <- !(ids %in% tx_lens$tx_id)
    fail_res <- IRanges()
    if (any(not_found)) {
        ## issue #69 does not apply here, even if redundant ids are used,
        ## all of them are considered.
        warning("Transcript(s) ", paste0("'", ids[not_found], "'",
                                         collapse = ", "),
                " could not be found", call. = FALSE)
        fail_res <- IRanges(start = rep(-1, sum(not_found)),
                            width = rep(1, sum(not_found)))
        ## Add metadata columnds.
        mcols(fail_res) <- DataFrame(tx_id = ids[not_found],
                                     tx_start = start(x[not_found]),
                                     tx_end = end(x[not_found]),
                                     cds_ok = NA)
        if (all(not_found))
            return(fail_res[match(internal_ids,
                                  paste0(mcols(fail_res)$tx_id, ":",
                                         mcols(fail_res)$tx_start, ":",
                                         mcols(fail_res)$tx_end))])
    }
    ## Remove non-coding transcripts
    not_coding <- is.na(tx_lens$cds_len)
    if (any(not_coding)) {
        not_coding_id <- tx_lens$tx_id[not_coding]
        warning("Transcript(s) ", paste0("'", not_coding_id, "'",
                                         collapse = ", "),
                " do/does not encode a protein", call. = FALSE)
        ## Subset x to the non-coding ones (issue #69)
        x_sub <- x[ids %in% not_coding_id]
        not_coding_res <- IRanges(start = rep(-1, length(x_sub)),
                                  width = rep(1, length(x_sub)))
        ## Add metadata columnds.
        mcols(not_coding_res) <- DataFrame(tx_id = names(x_sub),
                                           tx_start = start(x_sub),
                                           tx_end = end(x_sub),
                                           cds_ok = NA)
        fail_res <- c(fail_res, not_coding_res)
        if (all(not_coding))
            return(fail_res[match(internal_ids,
                                  paste0(mcols(fail_res)$tx_id, ":",
                                         mcols(fail_res)$tx_start, ":",
                                         mcols(fail_res)$tx_end))])
        tx_lens <- tx_lens[!not_coding, , drop = FALSE]
    }
    ## Ensure we have no missing 3' or 5' UTRs
    tx_lens$utr3_len[is.na(tx_lens$utr3_len)] <- 0
    tx_lens$utr5_len[is.na(tx_lens$utr5_len)] <- 0
    x <- x[ids %in% tx_lens$tx_id]
    id_idx <- match(names(x), tx_lens$tx_id)
    ## Are the coordinates within the CDS? Exclude the stop codon!
    outside_cds <- start(x) <= tx_lens$utr5_len[id_idx] |
        end(x) > (tx_lens$utr5_len[id_idx] + tx_lens$cds_len[id_idx] - 3)
    if (any(outside_cds)) {
        ## Also save against issue #69.
        outside_res <- IRanges(start = rep(-1, sum(outside_cds)),
                               width = rep(1, sum(outside_cds)))
        ## Add metadata columnds.
        mcols(outside_res) <- DataFrame(tx_id = names(x)[outside_cds],
                                           tx_start = start(x[outside_cds]),
                                           tx_end = end(x[outside_cds]),
                                           cds_ok = NA)
        fail_res <- c(fail_res, outside_res)
        warning("Provided coordinates for ",
                paste0("'", names(x)[outside_cds], "'", collapse = ", "),
                " are not within the coding region", call. = FALSE)
        if (all(outside_cds))
            return(fail_res[match(internal_ids,
                                  paste0(mcols(fail_res)$tx_id, ":",
                                         mcols(fail_res)$tx_start, ":",
                                         mcols(fail_res)$tx_end))])
        x <- x[!outside_cds]
    }
    ## Now check if the CDS length matches the protein sequence length.
    prt <- proteins(db, filter = TxIdFilter(names(x)),
                    columns = "protein_sequence")
    id_idx_prt <- match(names(x), prt$tx_id)
    id_idx <- match(names(x), tx_lens$tx_id)
    res <- IRanges(start = ceiling((start(x) - tx_lens$utr5_len[id_idx]) / 3),
                   end = ceiling((end(x) - tx_lens$utr5_len[id_idx]) / 3),
                   names = prt$protein_id[id_idx_prt])
    cds_ok <- nchar(prt$protein_sequence[id_idx_prt]) * 3 ==
        (tx_lens$cds_len[id_idx] - 3)
    if (any(!cds_ok))
        warning("The CDS of ", paste0("'", unique(names(x)[!cds_ok]), "'",
                                      collapse = ", "),
                " does not match the length of the encoded protein. Returned",
                " protein coordinates for this/these transcript(s) might not",
                " be correct", call. = FALSE)
    mcols(res) <- DataFrame(tx_id = names(x), tx_start = start(x),
                            tx_end = end(x), cds_ok = cds_ok)
    res <- c(fail_res, res)
    res[match(internal_ids, paste0(mcols(res)$tx_id, ":",
                                   mcols(res)$tx_start, ":",
                                   mcols(res)$tx_end))]
}

#' @title Map transcript-relative coordinates to genomic coordinates
#'
#' @description
#'
#' `transcriptToGenome` maps transcript-relative coordinates to genomic
#' coordinates.
#'
#' @inheritParams transcriptToProtein
#' 
#' @return
#'
#' `GRangesList` with the same length (and order) than the input `IRanges`
#' `x`. Each `GRanges` in the `GRangesList` provides the genomic coordinates
#' corresponding to the provided within-transcript coordinates. The
#' original transcript ID and the transcript-relative coordinates are provided
#' as metadata columns as well as the ID of the individual exon(s). An empty
#' `GRanges` is returned for transcripts that can not be found in the database.
#' 
#' @md
#' 
#' @author Johannes Rainer
#'
#' @family coordinate mapping functions
#'
#' @examples
#' 
#' library(EnsDb.Hsapiens.v75)
#' ## Restrict all further queries to chromosome x to speed up the examples
#' edbx <- filter(EnsDb.Hsapiens.v75, filter = ~ seq_name == "X")
#'
#' ## Below we map positions 1 to 5 within the transcript ENST00000381578 to
#' ## the genome. The ID of the transcript has to be provided either as names
#' ## or in one of the IRanges' metadata columns
#' txpos <- IRanges(start = 1, end = 5, names = "ENST00000381578")
#'
#' transcriptToGenome(txpos, edbx)
#' ## The object returns a GRangesList with the genomic coordinates, in this
#' ## example the coordinates are within the same exon and map to a single
#' ## genomic region.
#'
#' ## Next we map nucleotides 501 to 505 of ENST00000486554 to the genome
#' txpos <- IRanges(start = 501, end = 505, names = "ENST00000486554")
#'
#' transcriptToGenome(txpos, edbx)
#' ## The positions within the transcript are located within two of the
#' ## transcripts exons and thus a `GRanges` of length 2 is returned.
#'
#' ## Next we map multiple regions, two within the same transcript and one
#' ## in a transcript that does not exist.
#' txpos <- IRanges(start = c(501, 1, 5), end = c(505, 10, 6),
#'     names = c("ENST00000486554", "ENST00000486554", "some"))
#'
#' res <- transcriptToGenome(txpos, edbx)
#'
#' ## The length of the result GRangesList has the same length than the
#' ## input IRanges
#' length(res)
#'
#' ## The result for the last region is an empty GRanges, because the
#' ## transcript could not be found in the database
#' res[[3]]
#'
#' res
transcriptToGenome <- function(x, db, id = "name") {
    if (missing(x) || !is(x, "IRanges"))
        stop("Argument 'x' is required and has to be an 'IRanges' object")
    if (missing(db) || !is(db, "EnsDb"))
        stop("Argument 'db' is required and has to be an 'EnsDb' object")
    res <- .tx_to_genome(x, db, id = id)
    not_found <- sum(lengths(res) == 0)
    if (not_found > 0)
        warning(not_found, " transcript(s) could not be found in the database")
    res
}

.tx_to_genome <- function(x, db, id = "name") {
    id <- match.arg(id, c("name", colnames(mcols(x))))
    if (id == "name")
        ids <- names(x)
    else ids <- mcols(x)[, id]
    if (any(is.null(ids)))
        stop("One or more of the provided IDs are NULL", call. = FALSE)
    names(x) <- ids
    exns <- exonsBy(db, filter = TxIdFilter(unique(ids)))
    ## Loop over x
    res <- lapply(split(x, f = 1:length(x)), function(z) {
        if (any(names(exns) == names(z))) {
            gen_map <- .to_genome(exns[[names(z)]], z)
            mcols(gen_map) <- DataFrame(mcols(gen_map), tx_start = start(z),
                                        tx_end = end(z))
            gen_map[order(gen_map$exon_rank)]
        } else GRanges()
    })
    names(res) <- names(x)
    as(res, "GRangesList")
}
