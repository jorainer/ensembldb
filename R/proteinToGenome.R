#' @include Methods.R

## Code related to mapping of coordinates within proteins to genomic
## coordinates.

#' @description Convert within-protein coordinates to transcript (CDS) relative
#'    coordinates.
#'
#' @param x `IRanges` object with coordinates within a protein sequence.
#'
#' @return `IRanges` with the coordinates within the coding region (!) of the
#'     encoding transcript.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' prt <- IRanges(start = 23, end = 27)
#' .proteinCoordsToTx(prt)
.proteinCoordsToTx <- function(x) {
    end(x) <- end(x) * 3
    start(x) <- 1 + (start(x) - 1) * 3
    x
}

#' @description
#'
#' Fetch the CDS for all transcripts encoding the specified protein.
#'
#' @param x `EnsDb` object.
#'
#' @param id `character` with the protein ID(s).
#'
#' @param idType `character(1)` defining what type of IDs are provided. Has to
#'     be one of `"protein_id"` (default), `"uniprot_id"` or `"tx_id"`.
#'
#' @return a `list` with the CDS of the encoding transcript(s) for each provided
#'     id (as a `GRangesList`). Names of the `list` are the ids, if no
#'     transcript was found `NULL` is returned.
#'
#' @author Johannes Rainer
#' 
#' @md
#'
#' @noRd
.cds_for_id <- function(x, id, idType = "protein_id") {
    cds <- lapply(id,
                  FUN = function(x_id) {
                      suppressWarnings(
                          res <-
                              cdsBy(x, by = "tx",
                                    filter = .filter_for_idType(x_id, idType),
                                    columns = unique(c(idType,
                                                       "tx_id",
                                                       "protein_id")))
                      )
                      if (length(res) == 0)
                          warning("No CDS found for '", x_id, "'",
                                  call. = FALSE)
                      res
                  })
    names(cds) <- id
    cds
}

#' @description
#'
#' Uses .cds_for_id to fetch CDS for protein identifiers and in addition
#' checks that the CDS has a length that fits the range.
#'
#' @param x `EndDb` object.
#'
#' @param range The `IRange` with the position within protein.
#'
#' @param id `character(1)` specifying where the protein identifier can be found.
#'     Has to be either `"name"` or one of `colnames(mcols(prng))`.
#'
#' @param idType `character(1)` defining what type of IDs are provided. Has to
#'     be one of `"protein_id"` (default), `"uniprot_id"` or `"tx_id"`.
#'
#' @return `list` of `GRangesList`.
#' 
#' @author Johannes Rainer
#' 
#' @md
#'
#' @noRd
.cds_for_id_range <- function(x, range, id = "name", idType = "protein_id") {
    idType <- match.arg(idType, c("protein_id", "uniprot_id", "tx_id"))
    id <- match.arg(id, c("name", colnames(mcols(range))))
    if (id == "name")
        ids <- names(range)
    else ids <- mcols(range)[, id]
    cds <- .cds_for_id(x, ids, idType = idType)
    ## check returned CDS if their width matches the provided protein ranges
    mapply(cds, end(range) * 3,
                  FUN = function(x, x_end) {
                      if (!is.null(x))
                          x[sum(width(x)) >= x_end]
                  })
}

#' @description
#'
#' Fetches the protein sequence using the provided `"protein_id"` in `cds` and
#' selects checks whether the cds lenghts match the protein sequence. The
#' function returns a single CDS for each protein/transcript (the first of
#' those that have the correct length) and throws a warning if none of the
#' specified CDS matches the protein sequence.
#'
#' @details
#'
#' Mismatch between the CDS and the AA sequence are in most instances caused by
#' incomplete (3' or 5') CDS sequences.
#'
#' @param x `EnsDb` object.
#'
#' @param cds `GRangesList`
#'
#' @return `list` of `GRangesList` with one CDS per element/protein.
#' 
#' @author Johannes Rainer
#' 
#' @md
#'
#' @noRd
.cds_matching_protein <- function(x, cds) {
    ## Fetch the protein sequences, all in one go.
    ## Loop through the cds
    prot_ids <- unique(unlist(lapply(cds, function(z) {
        lapply(z, function(y) y$protein_id)
    }), use.names = FALSE))
    prot_seqs <- proteins(x, filter = ProteinIdFilter(prot_ids),
                          columns = c("protein_id", "protein_sequence"))
    ## Calculate the expected CDS length (add +3 to add the stop codon).
    exp_cds_len <- nchar(prot_seqs$protein_sequence) * 3 + 3
    names(exp_cds_len) <- prot_seqs$protein_id
    mapply(cds, names(cds), FUN = function(z, nm) {
        if (!is.null(z)) {
            cds_lens <- sum(width(z))
            prt_ids <- unlist(lapply(z, function(y) unique(y$protein_id)))
            diffs <- cds_lens - exp_cds_len[prt_ids]
            ## Select one for which the
            idx <- order(abs(diffs))[1]
            if (diffs[idx] != 0)
                warning("Could not find a CDS whith the expected length for ",
                        "protein: '", nm, "'. The returned genomic ",
                        "coordinates might thus not be correct for this ",
                        "protein.", call. = FALSE)
            ## Alternatively we could align the RNA and AA sequences and trim
            ## the protein sequence...
            z[idx]
        }
    })
}

#' @description
#'
#' Map protein-relative coordinates to coordinates within the encoding
#' transcript.
#'
#' @noRd
proteinToTranscript <- function(protein, x, id = "name",
                                idType = "protein_id") {
}


#' @description
#'
#' Map within-protein coordinates to genome coordinates.
#' 
#' @param protein `IRanges` with the coordinates within the protein(s). The
#'     object has also to provide some means to identify the protein (see
#'     details).
#'
#' @param genome `EnsDb` object.
#'
#' @param id `character(1)` specifying where the protein identifier can be
#'     found. Has to be either `"name"` or one of `colnames(mcols(prng))`.
#'
#' @param idType `character(1)` defining what type of IDs are provided. Has to
#'     be one of `"protein_id"` (default), `"uniprot_id"` or `"tx_id"`.
#'
#' @author Johannes Rainer
#' 
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' irs <- IRanges(start = c(14, 17, 13), end = c(14, 23, 13))
#' names(irs) <- c("ALAT2_HUMAN", "H0YIP2_HUMAN", "what")
proteinToGenome <- function(protein, genome, id = "name",
                            idType = "protein_id") {
    coords_cds <- .proteinCoordsToTx(protein)
    ## 1) retrieve CDS for each protein
    message("Fetching CDS for ", length(protein), " proteins ... ",
            appendLF = FALSE)
    cds_genome <- .cds_for_id_range(genome, protein, id = id, idType = idType)
    message("OK")
    ## Convert exons to tx relative, i.e. 1...10,11...15 etc.
    ## identify the exons with the start end of the coords_cds
    ## calculate the coords.
    res <- mapply(cds_genome, as(coords_cds, "IRangesList"),
                  function(x, y) {
                      cds_rel <- .splice(x)
                      start_exon <- which()
                  })
}

#' @description
#'
#' Function to map coordinates within a CDS to genomic coordinates, based on the
#' genomic coordinates of the exons of the CDS. See examples below.
#'
#' @note
#'
#' While designed for cds-relative coordinates, the function should also
#' work for `g_coords` that represent the exons of a transcript and
#' `cds_coords` being transcript-relative.
#' 
#' @param g_coords `GRanges` with the exons encoding the CDS. These have to be
#'     provided in the correct order (i.e. first element being the first exon
#'     etc, also for genes on the reverse strand!). If `g_coords` is a
#'     `GRangesList`, `cds_coords` is mapped to genomic coordinates using each
#'     `GRanges` in the `GRangesList`.
#'
#' @param cds_coords `IRanges` with CDS-relative coordinates.
#'
#' @return A `GRanges` with the result from the within-cds coordinate mapping
#'     to the genome.
#' 
#' @md
#'
#' @noRd
#'
#' @author Johannes Rainer based on code from Laurent Gatto and Sebastian Gibb
#'
#' @examples
#'
#' g_coords <- GRanges("1",
#'     ranges = IRanges(start = c(3, 8, 15, 19), end = c(5, 12, 16, 21)),
#'     strand = "+")
#'
#' cds_coords <- IRanges(start = 5, end = 12)
#' ## Expect: 9:20
#' 
#' .to_genome(g_coords, cds_coords)
#'
#' ## Reverse strand example
#' g_c <- GRanges("1", IRanges(start = c(28, 21, 16, 10, 3),
#'     end = c(30, 25, 17, 12, 6)), strand = "-")
#'
#' c_c <- IRanges(start = 2, end = 2)
#' ## Expect: 29:29
#' .to_genome(g_c, c_c)
#'
#' c_c <- IRanges(start = 8, end = 16)
#' ## Expect: 4:6, 10:12, 16:17, 21:21
#' .to_genome(g_c, c_c)
#'
#' ## What if we've got a GRangesList?
#' g_l <- GRangesList(g_c, g_c[2:5])
#' c_c <- IRanges(start = 4, end = 9)
#' ## Expect: 17:17, 21:25
#' ##         11:12, 16:17, 21:22
#' .to_genome(g_l, c_c)
.to_genome <- function(g_coords, cds_coords) {
    if (is(g_coords, "GRangesList"))
        return(lapply(g_coords, .to_genome, cds_coords = cds_coords))
    if (!is(g_coords, "GRanges"))
        stop("'g_coords' is supposed to be a 'GRanges' object")
    cds_rel <- .splice(g_coords)
    strnd <- unique(as.character(strand(g_coords)))
    seqlvl <- unique(seqlevels(g_coords))
    if (length(strnd) > 1)
        stop("All exons are expected to be located on the same strand")
    if (length(seqlvl) > 1)
        stop("All exons are expected to be located on the same chromosome")
    if (strnd == "-")
        strnd_num <- -1L
    else
        strnd_num <- 1L
    cums <- cumsum(width(cds_rel))
    ## Find the exons in which the start and end is located.
    start_exon <- which(start(cds_rel) <= start(cds_coords) &
                        end(cds_rel) >= start(cds_coords))
    end_exon <- which(start(cds_rel) <= end(cds_coords) &
                      end(cds_rel) >= end(cds_coords))
    if (length(start_exon) == 0 | length(end_exon) == 0)
        stop("The within transcript/CDS coordinates are outside the region ",
             "defined by the provided exons")
    ## Convert within CDS coordinates into within exon coordinates.
    exon_rel_start <- start(cds_coords) + width(cds_rel)[start_exon] -
        cums[start_exon]
    exon_rel_end <- end(cds_coords) + width(cds_rel)[end_exon] -
        cums[end_exon]
    ## Convert into genomic coordinates.
    if (strnd_num > 0) {
         genome_start <- start(g_coords)[start_exon] + exon_rel_start - 1
         genome_end <- start(g_coords)[end_exon] + exon_rel_end - 1
    } else {
        genome_end <- end(g_coords)[start_exon] - exon_rel_start + 1
        genome_start <- end(g_coords)[end_exon] - exon_rel_end + 1
    }
    genome_coords <- GRanges(seqnames = seqlvl,
                             IRanges(start = genome_start, end = genome_end),
                             strand = strnd)
    seqinfo(genome_coords) <- seqinfo(g_coords)
    intersect(genome_coords, g_coords)
}

.filter_for_idType <- function(x, idType) {
    if (idType == "protein_id")
        return(ProteinIdFilter(x))
    if (idType == "uniprot_id")
        return(UniprotFilter(x))
    if (idType == "tx_id")
        return(TxIdFilter(x))
    NULL
}

#' @description
#'
#' *Splices* the provided `GRanges`/`IRanges` by removing all intronic ranges.
#'
#' @param x `IRanges`
#' 
#' @author Johannes Rainer
#' 
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' ir <- IRanges(start = c(5, 10, 15), end = c(7, 13, 16))
#' .splice(ir)
.splice <- function(x) {
    if (is(x, "IRangesList") | is(x, "list") | is(x, "GRangesList"))
        lapply(x, .splice)
    else
        IRanges(start = c(1, cumsum(width(x)[-length(x)]) + 1),
                width = width(x))
}
