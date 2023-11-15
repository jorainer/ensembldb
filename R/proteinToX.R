#' @include Methods.R

#' @title Map protein-relative coordinates to positions within the transcript
#'
#' @description
#'
#' `proteinToTranscript` maps protein-relative coordinates to positions within
#' the encoding transcript. Note that the returned positions are relative to
#' the complete transcript length, which includes the 5' UTR.
#'
#' Similar to the [proteinToGenome()] function, `proteinToTranscript` compares
#' for each protein whether the length of its sequence matches the length of
#' the encoding CDS and throws a warning if that is not the case. Incomplete
#' 3' or 5' CDS of the encoding transcript are the most common reasons for a
#' mismatch between protein and transcript sequences.
#'
#' @details
#'
#' Protein identifiers (supported are Ensembl protein IDs or Uniprot IDs) can
#' be passed to the function as `names` of the `x` `IRanges` object, or
#' alternatively in any one of the metadata columns (`mcols`) of `x`.
#'
#' @note
#'
#' While mapping of Ensembl protein IDs to Ensembl transcript IDs is 1:1, a
#' single Uniprot identifier can be annotated to several Ensembl protein IDs.
#' `proteinToTranscript` calculates in such cases transcript-relative
#' coordinates for each annotated Ensembl protein.
#'
#' Mapping using Uniprot identifiers needs also additional internal checks that
#' can have a significant impact on the performance of the function. It is thus
#' strongly suggested to first identify the Ensembl protein identifiers for the
#' list of input Uniprot identifiers (e.g. using the [proteins()] function and
#' use these as input for the mapping function.
#'
#' @inheritParams proteinToGenome
#' 
#' @return
#'
#' `IRangesList`, each element being the mapping results for one of the input
#' ranges in `x`. Each element is a `IRanges` object with the positions within
#' the encoding transcript (relative to the start of the transcript, which
#' includes the 5' UTR). The transcript ID is reported as the name of each
#' `IRanges`. The `IRanges` can be of length > 1 if the provided
#' protein identifier is annotated to more than one Ensembl protein ID (which
#' can be the case if Uniprot IDs are provided). If the coordinates can not be
#' mapped (because the protein identifier is unknown to the database) an
#' `IRanges` with negative coordinates is returned.
#'
#' The following metadata columns are available in each `IRanges` in the result:
#' + `"protein_id"`: the ID of the Ensembl protein for which the within-protein
#'   coordinates were mapped to the genome.
#' + `"tx_id"`: the Ensembl transcript ID of the encoding transcript.
#' + `"cds_ok"`: contains `TRUE` if the length of the CDS matches the length
#'    of the amino acid sequence and `FALSE` otherwise.
#' + `"protein_start"`: the within-protein sequence start coordinate of the
#'   mapping.
#' + `"protein_end"`: the within-protein sequence end coordinate of the mapping.
#'
#' @family coordinate mapping functions
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' library(EnsDb.Hsapiens.v86)
#' ## Restrict all further queries to chromosome x to speed up the examples
#' edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "X")
#'
#' ## Define an IRange with protein-relative coordinates within a protein for
#' ## the gene SYP
#' syp <- IRanges(start = 4, end = 17)
#' names(syp) <- "ENSP00000418169"
#' res <- proteinToTranscript(syp, edbx)
#' res
#' ## Positions 4 to 17 within the protein span are encoded by the region
#' ## from nt 23 to 64.
#'
#' ## Perform the mapping for multiple proteins identified by their Uniprot
#' ## IDs.
#' ids <- c("O15266", "Q9HBJ8", "unexistant")
#' prngs <- IRanges(start = c(13, 43, 100), end = c(21, 80, 100))
#' names(prngs) <- ids
#'
#' res <- proteinToTranscript(prngs, edbx, idType = "uniprot_id")
#'
#' ## The result is a list, same length as the input object
#' length(res)
#' names(res)
#'
#' ## No protein/encoding transcript could be found for the last one
#' res[[3]]
#'
#' ## The first protein could be mapped to multiple Ensembl proteins. The
#' ## region within all transcripts encoding the region in the protein are
#' ## returned
#' res[[1]]
#'
#' ## The result for the region within the second protein
#' res[[2]]
proteinToTranscript <- function(x, db, id = "name",
                                idType = "protein_id") {
    if (missing(x) || !is(x, "IRanges"))
        stop("Argument 'x' is required and has to be an 'IRanges' object")
    if (missing(db) || !is(db, "EnsDb"))
        stop("Argument 'db' is required and has to be an 'EnsDb' object")
    coords_cds <- .proteinCoordsToTx(x)
    ## 1) retrieve CDS for each protein
    message("Fetching CDS for ", length(x), " proteins ... ",
            appendLF = FALSE)
    cds_genome <- .cds_for_id_range(db, x, id = id, idType = idType)
    miss <- lengths(cds_genome) == 0
    if (any(miss))
        warning("No CDS found for: ", paste0(names(cds_genome)[miss],
                                             collapse = ", "))
    message(sum(!miss), " found")
    ## 2) ensure that the CDS matches the AA sequence length
    message("Checking CDS and protein sequence lengths ... ", appendLF = FALSE)
    cds_genome <- .cds_matching_protein(db, cds_genome)
    are_ok <- vapply(cds_genome, function(z) {
        if (is(z, "GRangesList"))
            all(z[[1]]$cds_ok)
        else NA
    }, FUN.VALUE = logical(1))
    are_ok <- are_ok[!is.na(are_ok)]
    ## We've got now a list of GRanges
    message(sum(are_ok), "/", length(are_ok), " OK")
    ## Get for each transcript it's 5' UTR and add its width to the coords_cds
    tx_ids <- unique(unlist(lapply(cds_genome, names), use.names = FALSE))
    if (length(tx_ids)) {
        five_utr <- fiveUTRsByTranscript(db, filter = TxIdFilter(tx_ids))
        ## Calculate 5' widths for these
        five_width <- sum(width(five_utr))
    }
    as(mapply(
        cds_genome, as(coords_cds, "IRangesList"), split(x, 1:length(x)),
        FUN = function(gnm, cds, prt) {
            if (is.null(gnm)) {
                ## Define the metadata columns
                mc <- DataFrame(protein_id = NA_character_,
                                tx_id = NA_character_,
                                cds_ok = NA,
                                protein_start = start(prt),
                                protein_end = end(prt))
                if (idType == "uniprot_id")
                    mc$uniprot_id <- names(prt)
                else mc$protein_id <- names(prt)
                ir <- IRanges(start = -1, width = 1)
                mcols(ir) <- mc
                if(!is.null(mcols(prt))) 
                    mcols(ir) <- cbind(mcols(ir),mcols(prt))
                ir
            } else {
                ids <- names(gnm)
                res <- IRanges(start = start(cds) + five_width[ids],
                               end = end(cds) + five_width[ids],
                               names = ids)
                ## Populate mcols
                mc <- DataFrame(
                    protein_id = unlist(
                        lapply(gnm, function(z) z$protein_id[1])),
                    tx_id = ids,
                    cds_ok = gnm[[1]]$cds_ok[1],
                    protein_start = start(prt),
                    protein_end = end(prt)
                )
                if (idType == "uniprot_id")
                    mc$uniprot_id <- names(prt)
                mcols(res) <- mc
                if(!is.null(mcols(prt))) 
                    mcols(res) <- cbind(mcols(res),mcols(prt))
                res
            }
        }), "IRangesList")
}


#' @title Map within-protein coordinates to genomic coordinates
#'
#' @name proteinToGenome
#'
#' @aliases proteinToGenome proteinToGenome,EnsDb-method
#'
#' @description
#'
#' `proteinToGenome` maps protein-relative coordinates to genomic coordinates
#' based on the genomic coordinates of the CDS of the encoding transcript. The
#' encoding transcript is identified using protein-to-transcript annotations
#' (and eventually Uniprot to Ensembl protein identifier mappings) from the
#' submitted `EnsDb` object (and thus based on annotations from Ensembl).
#'
#' Not all coding regions for protein coding transcripts are complete, and the
#' function thus checks also if the length of the coding region matches the
#' length of the protein sequence and throws a warning if that is not the case.
#'
#' The genomic coordinates for the within-protein coordinates, the Ensembl
#' protein ID, the ID of the encoding transcript and the within protein start
#' and end coordinates are reported for each input range.
#'
#' @details
#'
#' Protein identifiers (supported are Ensembl protein IDs or Uniprot IDs) can
#' be passed to the function as `names` of the `x` `IRanges` object, or
#' alternatively in any one of the metadata columns (`mcols`) of `x`.
#'
#' @note
#'
#' While the mapping for Ensembl protein IDs to encoding transcripts (and
#' thus CDS) is 1:1, the mapping between Uniprot identifiers and encoding
#' transcripts (which is based on Ensembl annotations) can be one to many. In
#' such cases `proteinToGenome` calculates genomic coordinates for
#' within-protein coordinates for all of the annotated Ensembl proteins and
#' returns all of them. See below for examples.
#'
#' Mapping using Uniprot identifiers needs also additional internal checks that
#' have a significant impact on the performance of the function. It is thus
#' strongly suggested to first identify the Ensembl protein identifiers for the
#' list of input Uniprot identifiers (e.g. using the [proteins()] function and
#' use these as input for the mapping function.
#'
#' A warning is thrown for proteins which sequence does not match the coding
#' sequence length of any encoding transcripts. For such proteins/transcripts
#' a `FALSE` is reported in the respective `"cds_ok"` metadata column.
#' The most common reason for such discrepancies are incomplete 3' or 5' ends
#' of the CDS. The positions within the protein might not be correclty
#' mapped to the genome in such cases and it might be required to check
#' the mapping manually in the Ensembl genome browser.
#'
#' @param x `IRanges` with the coordinates within the protein(s). The
#'     object has also to provide some means to identify the protein (see
#'     details).
#'
#' @param db `EnsDb` object to be used to retrieve genomic coordinates of
#'     encoding transcripts.
#'
#' @param id `character(1)` specifying where the protein identifier can be
#'     found. Has to be either `"name"` or one of `colnames(mcols(prng))`.
#'
#' @param idType `character(1)` defining what type of IDs are provided. Has to
#'     be one of `"protein_id"` (default), `"uniprot_id"` or `"tx_id"`.
#'
#' @return
#'
#' `list`, each element being the mapping results for one of the input
#' ranges in `x` and names being the IDs used for the mapping. Each
#' element can be either a:
#' + `GRanges` object with the genomic coordinates calculated on the
#'   protein-relative coordinates for the respective Ensembl protein (stored in
#'   the `"protein_id"` metadata column.
#' + `GRangesList` object, if the provided protein identifier in `x` was
#'   mapped to several Ensembl protein IDs (e.g. if Uniprot identifiers were
#'   used). Each element in this `GRangesList` is a `GRanges` with the genomic
#'   coordinates calculated for the protein-relative coordinates from the
#'   respective Ensembl protein ID.
#'
#' The following metadata columns are available in each `GRanges` in the result:
#' + `"protein_id"`: the ID of the Ensembl protein for which the within-protein
#'   coordinates were mapped to the genome.
#' + `"tx_id"`: the Ensembl transcript ID of the encoding transcript.
#' + `"exon_id"`: ID of the exons that have overlapping genomic coordinates.
#' + `"exon_rank"`: the rank/index of the exon within the encoding transcript.
#' + `"cds_ok"`: contains `TRUE` if the length of the CDS matches the length
#'    of the amino acid sequence and `FALSE` otherwise.
#' + `"protein_start"`: the within-protein sequence start coordinate of the
#'   mapping.
#' + `"protein_end"`: the within-protein sequence end coordinate of the mapping.
#'
#' Genomic coordinates are returned ordered by the exon index within the
#' transcript.
#'
#' @family coordinate mapping functions
#'
#' @seealso \code{\link[GenomicFeatures]{proteinToGenome}} in the
#'     \pkg{GenomicFeatures} package for methods that operate on a
#'     TxDb or GRangesList object.
#'
#' @author Johannes Rainer based on initial code from Laurent Gatto and
#'     Sebastian Gibb
#'
#' @md
#'
#' @examples
#'
#' library(EnsDb.Hsapiens.v86)
#' ## Restrict all further queries to chromosome x to speed up the examples
#' edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "X")
#'
#' ## Define an IRange with protein-relative coordinates within a protein for
#' ## the gene SYP
#' syp <- IRanges(start = 4, end = 17)
#' names(syp) <- "ENSP00000418169"
#' res <- proteinToGenome(syp, edbx)
#' res
#' ## Positions 4 to 17 within the protein span two exons of the encoding
#' ## transcript.
#'
#' ## Perform the mapping for multiple proteins identified by their Uniprot
#' ## IDs.
#' ids <- c("O15266", "Q9HBJ8", "unexistant")
#' prngs <- IRanges(start = c(13, 43, 100), end = c(21, 80, 100))
#' names(prngs) <- ids
#'
#' res <- proteinToGenome(prngs, edbx, idType = "uniprot_id")
#'
#' ## The result is a list, same length as the input object
#' length(res)
#' names(res)
#'
#' ## No protein/encoding transcript could be found for the last one
#' res[[3]]
#'
#' ## The first protein could be mapped to multiple Ensembl proteins. The
#' ## mapping result using all of their encoding transcripts are returned
#' res[[1]]
#'
#' ## The coordinates within the second protein span two exons
#' res[[2]]
setMethod("proteinToGenome", "EnsDb",
  function(x, db, id = "name", idType = "protein_id") {
    if (missing(x) || !is(x, "IRanges"))
        stop("Argument 'x' is required and has to be an 'IRanges' object")
    if (missing(db) || !is(db, "EnsDb"))
        stop("Argument 'db' is required and has to be an 'EnsDb' object")
    coords_cds <- .proteinCoordsToTx(x)
    ## 1) retrieve CDS for each protein
    message("Fetching CDS for ", length(x), " proteins ... ",
            appendLF = FALSE)
    cds_genome <- .cds_for_id_range(db, x, id = id, idType = idType)
    miss <- lengths(cds_genome) == 0
    if (any(miss))
        warning("No CDS found for: ", paste0(names(cds_genome)[miss],
                                             collapse = ", "))
    message(sum(!miss), " found")
    ## 2) ensure that the CDS matches the AA sequence length
    message("Checking CDS and protein sequence lengths ... ", appendLF = FALSE)
    cds_genome <- .cds_matching_protein(db, cds_genome)
    are_ok <- unlist(lapply(cds_genome, function(z) {
        if (is(z, "GRangesList") & length(z) > 0)
            all(z[[1]]$cds_ok)
        else NA
    }))
    message(sum(are_ok, na.rm = TRUE), "/", length(are_ok), " OK")
    are_ok <- are_ok[!is.na(are_ok)]
    ## We've got now a list of GRanges
    ## Perform the mapping for each input range with each mapped cds
    x_IRangesList <- as(x, "IRangesList") ## preserve the metadata
    x_IRangesList@unlistData@elementMetadata <- x@elementMetadata
    res <- mapply(
        cds_genome, as(coords_cds, "IRangesList"), x_IRangesList,
        FUN = function(gnm, cds, prt) {
            if (!length(gnm)) { ## addresses NULL and empty elements
                GRanges()
            } else {
                ## Unlist because we'd like to have a GRanges here. Will split
                ## again later.
                maps <- unlist(.to_genome(gnm, cds))
                ## Don't want to have GRanges names!
                names(maps) <- NULL
                mcols(maps) <- cbind(mcols(maps),
                                     protein_start = start(prt),
                                     protein_end = end(prt))
                if(!is.null(mcols(prt))) 
                    mcols(maps) <- cbind(mcols(maps),mcols(prt))
                maps[order(maps$exon_rank)]
            }
        })
    ## Split each element again, if there are more than one protein_id. Names
    ## of the elements are then the protein_id.
    lapply(res, function(z) {
        if (length(unique(z$protein_id)) > 1)
            split(z, f = z$protein_id)
        else z
    })
  })

#' @rdname proteinToGenome
#'
#' @aliases proteinToGenome proteinToGenome,Preloaded-method
#'
#' @param db For the method for `EnsDb` objects: An `EnsDb` object to be used to
#'     retrieve genomic coordinates of encoding transcripts.<br><br>
#'     For the method for `CompressedGRangesList` objects: A `CompressedGRangesList` object 
#'     generated by [cdsBy()] where by = 'tx' and columns = c(listColumns(edb,'tx')
#'     ,'protein_id','uniprot_id','protein_sequence').
#'
#' @md
#'
#' @examples
#'
#' ## Meanwhile, this function can be called in parallel processes if you preload
#' ## the CDS data with desired data columns
#' cds <- cdsBy(edb,columns = c(listColumns(edb,'tx'),'protein_id','uniprot_id','protein_sequence'))
#' cds <- cdsBy(edb,columns = c(listColumns(edb,'tx'),'protein_id','protein_sequence'))
#' cds <- cdsBy(edb,columns = c('tx_id','protein_id','protein_sequence'))
#' ## Define an IRange with protein-relative coordinates within a protein for
#' ## the gene SYP
#' syp <- IRanges(start = 4, end = 17)
#' names(syp) <- "ENSP00000418169"
#' res <- proteinToGenome(syp, cds)
#' res
#' ## Positions 4 to 17 within the protein span two exons of the encoding
#' ## transcript.
#'
#' ## Perform the mapping for multiple proteins identified by their Uniprot
#' ## IDs.
#' ids <- c("O15266", "Q9HBJ8", "unexistant")
#' prngs <- IRanges(start = c(13, 43, 100), end = c(21, 80, 100))
#' names(prngs) <- ids
#'
#' res <- proteinToGenome(prngs, cds, idType = "uniprot_id")
setMethod("proteinToGenome", "CompressedGRangesList",
  function(x, db, id = "name", idType = "protein_id") {
    if (missing(x) || !is(x, "IRanges"))
        stop("Argument 'x' is required and has to be an 'IRanges' object")
    if (missing(db) || !is(db, "CompressedGRangesList"))
        stop("Argument 'db' is required and has to be an 'CompressedGRangesList' object")
    coords_cds <- .proteinCoordsToTx(x)
    ## 1) retrieve CDS for each protein
    message("Fetching CDS for ", length(x), " proteins ... ",
            appendLF = FALSE)
    cds_genome <- .cds_for_id_range(db, x, id = id, idType = idType)
    miss <- lengths(cds_genome) == 0
    if (any(miss))
        warning("No CDS found for: ", paste0(names(cds_genome)[miss],
                                             collapse = ", "))
    message(sum(!miss), " found")
    ## 2) ensure that the CDS matches the AA sequence length
    message("Checking CDS and protein sequence lengths ... ", appendLF = FALSE)
    cds_genome <- .cds_matching_protein(cds_genome)
    are_ok <- unlist(lapply(cds_genome, function(z) {
        if (is(z, "GRangesList") & length(z) > 0)
            all(z[[1]]$cds_ok)
        else NA
    }))
    message(sum(are_ok, na.rm = TRUE), "/", length(are_ok), " OK")
    are_ok <- are_ok[!is.na(are_ok)]
    ## We've got now a list of GRanges
    ## Perform the mapping for each input range with each mapped cds
    x_IRangesList <- as(x, "IRangesList") ## preserve the metadata
    x_IRangesList@unlistData@elementMetadata <- x@elementMetadata
    res <- mapply(
        cds_genome, as(coords_cds, "IRangesList"), x_IRangesList,
        FUN = function(gnm, cds, prt) {
            if (!length(gnm)) { ## addresses NULL and empty elements
                GRanges()
            } else {
                ## Unlist because we'd like to have a GRanges here. Will split
                ## again later.
                maps <- unlist(.to_genome(gnm, cds))
                ## Don't want to have GRanges names!
                names(maps) <- NULL
                mcols(maps) <- cbind(mcols(maps),
                                     protein_start = start(prt),
                                     protein_end = end(prt))
                if(!is.null(mcols(prt))) 
                    mcols(maps) <- cbind(mcols(maps),mcols(prt))
                maps[order(maps$exon_rank)]
            }
        })
    ## Split each element again, if there are more than one protein_id. Names
    ## of the elements are then the protein_id.
    lapply(res, function(z) {
        if (length(unique(z$protein_id)) > 1)
            split(z, f = z$protein_id)
        else z
    })
  })

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

#' @description Fetch the CDS for all transcripts encoding the specified protein.
#'
#' @note
#'
#' Use one query to fetch CDS for all (unique) input IDs. If input IDs are
#' Uniprot identifiers we have to perform additional checks and data
#' re-organizations because one transcript (and thus CDS) can be associated
#' with multiple Uniprot identifiers.
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
.cds_for_id2 <- function(x, id, idType = "protein_id") {
    if (idType != "tx_id") {
        map <- transcripts(x, filter = .filter_for_idType(unique(id), idType),
                           columns = c("tx_id", idType), order.by = "tx_id",
                           return.type = "data.frame")
        tx_id <- split(map$tx_id, map[, idType])
    } else {
        tx_id <- id
        tx_id <- split(tx_id, id)
    }
    if (length(tx_id)) {
        suppressWarnings(
            cds <- cdsBy(
                x, columns = unique(c(idType, "tx_id", "protein_id")), by = "tx",
                filter = TxIdFilter(unique(unlist(tx_id, use.names = FALSE))))
        )
        if (length(cds)) {
            if (idType == "uniprot_id") {
                ## Additional (in)sanity checks for Uniprot identifiers
                ## This has a significant impact on performance, but otherwise
                ## we end up with duplicated transcript entries!
                tmp <- unlist(cds)
                tmp <- as.list(split(unname(tmp), tmp$uniprot_id))
                cds <- lapply(tmp, function(z) split(z, z$tx_id))
            } else
                cds <- lapply(tx_id, function(z) cds[z])
        }
        else cds <- list()
    } else cds <- list()
    cds <- cds[id]
    names(cds) <- id # to add also names for elements not found
    cds
}

#' @description Fetch the CDS for all transcripts encoding the specified protein using 
#'     preloaded data.
#'
#' @note
#'
#' Use one query to fetch CDS for all (unique) input IDs using preloaded data. If input 
#' IDs are Uniprot identifiers we have to perform additional checks and data
#' re-organizations because one transcript (and thus CDS) can be associated
#' with multiple Uniprot identifiers.
#' 
#' @param x `CompressedGRangesList` object generated by [cdsBy()] where 
#'     by = 'tx' and columns = c(listColumns(edb,'tx'),'protein_id','uniprot_id','protein_sequence').
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
#' @author Johannes Rainer, Boyu Yu
#'
#' @md
#'
#' @noRd
setMethod(".cds_for_id2", "CompressedGRangesList",
  function(x, id, idType = "protein_id") {
    ## No need for get the transcripts first, just use the metadata of CDS
    ## if you find this modification well, you may want to change the 
    ## generic .cds_for_id2 function above
    if (!all(c("tx_id", "protein_id",idType,"protein_sequence") %in% colnames(mcols(x[[1]]))))
        stop(
            paste(
                c(
                    "CDS must have", 
                    unique(c("'tx_id'","'protein_sequence'", "'protein_id'",paste("'",idType,"'",sep=''))), 
                    "in the metadata."
                ), 
            collapse = ' ')
        )
    map <- x@unlistData[as.data.frame(x@partitioning)[,1]]@elementMetadata[,unique(c('tx_id', 'protein_id', idType))]
    map <- map[map[,idType] %in% id,]
    if (idType != "tx_id") {
        tx_id <- split(map$tx_id, map[, idType])
    } else {
        tx_id <- id
        tx_id <- split(tx_id, id)
    }
    if (length(tx_id)) {
        cds <- x[names(x) %in% unique(unlist(tx_id, use.names = FALSE))]
        if (length(cds)) {
            if (idType == "uniprot_id") {
                ## Additional (in)sanity checks for Uniprot identifiers
                ## This has a significant impact on performance, but otherwise
                ## we end up with duplicated transcript entries!
                tmp <- unlist(cds)
                tmp <- as.list(split(unname(tmp), tmp$uniprot_id))
                cds <- lapply(tmp, function(z) split(z, z$tx_id))
            } else if ("uniprot_id" %in% colnames(mcols(x[[1]]))) {
                warning(
                    "Multiple uniprot ids were preserved for the same transcript entry. ",
                    "Please don't include uniprot_id column in cdsBy() if unneeded."
                )
                cds <- lapply(tx_id, function(z) cds[z])
            } else {
                cds <- lapply(tx_id, function(z) cds[z])
            }
        }
        else cds <- list()
    } else cds <- list()
    cds <- cds[id]
    names(cds) <- id # to add also names for elements not found
    cds
})

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
    cds <- .cds_for_id2(x, ids, idType = idType)
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
#' checks whether the cds lenghts match the protein sequence. The
#' function returns all CDS for each protein/transcript that have the correct
#' length) and throws a warning if none of the specified CDS matches the
#' protein sequence. In the latter case a single CDS (the one with the CDS
#' length closest to the expected length). Whether a CDS has the expected length
#' or not is also reported in the metadata column `"cds_ok"`.
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
#' @return `list` of `GRangesList`s with one list element per input protein,
#'    and the `GRangesList` containing the CDS of all matching (or one not
#'    matching) transcripts.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
setGeneric(".cds_matching_protein", signature="x",
    function(x, ...) standardGeneric(".cds_matching_protein")
)
setMethod(".cds_matching_protein", "EnsDb",
  function(x, cds) {
    ## Fetch the protein sequences, all in one go.
    ## Loop through the cds
    prot_ids <- unique(unlist(lapply(cds, function(z) {
        lapply(z, function(y) y$protein_id)
    }), use.names = FALSE))
    if (length(prot_ids) == 0)
        return(lapply(cds, function(z) NULL))
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
            ## Return all for which diff is 0, or otherwise the one with the
            ## smallest diff.
            if (any(diffs == 0)) {
                z <- lapply(z, function(grng) {
                    mcols(grng)$cds_ok <- TRUE
                    grng
                })
                GRangesList(z[diffs == 0])
            } else {
                warning("Could not find a CDS whith the expected length for ",
                        "protein: '", nm, "'. The returned genomic ",
                        "coordinates might thus not be correct for this ",
                        "protein.", call. = FALSE)
                ## Alternatively we could align the RNA and AA sequences and
                ## trim the protein sequence...
                z <- lapply(z, function(grng) {
                    mcols(grng)$cds_ok <- FALSE
                    grng
                })
                GRangesList(z[which.min(abs(diffs))])
            }
        }
    })
})

#' @description
#'
#' Fetches the protein sequence using the provided `"protein_id"` in `cds` and
#' checks whether the cds lenghts match the protein sequence. The
#' function returns all CDS for each protein/transcript that have the correct
#' length) and throws a warning if none of the specified CDS matches the
#' protein sequence. In the latter case a single CDS (the one with the CDS
#' length closest to the expected length). Whether a CDS has the expected length
#' or not is also reported in the metadata column `"cds_ok"`.
#'
#' @details
#'
#' Mismatch between the CDS and the AA sequence are in most instances caused by
#' incomplete (3' or 5') CDS sequences.
#'
#' @param x `list` object generated by .cds_for_id2.
#'
#'
#' @return `list` of `GRangesList`s with one list element per input protein,
#'    and the `GRangesList` containing the CDS of all matching (or one not
#'    matching) transcripts.
#'
#' @author Johannes Rainer, Boyu Yu
#'
#' @md
#' 
#' @noRd
setMethod(".cds_matching_protein", "list",
  function(x) {
    mapply(x, names(x), FUN = function(z, nm) {
        if (!is.null(z)) {
            cds_lens <- sum(width(z))
            exp_cds_len <- unlist(lapply(z, function(y)  nchar(unique(y$protein_sequence)) * 3 + 3))
            ## Calculate the expected CDS length (add +3 to add the stop codon).
            diffs <- cds_lens - exp_cds_len
            ## Return all for which diff is 0, or otherwise the one with the
            ## smallest diff.
            if (any(diffs == 0)) {
                z <- lapply(z, function(grng) {
                    mcols(grng)$cds_ok <- TRUE
                    grng
                })
                GRangesList(z[diffs == 0])
            } else {
                warning("Could not find a CDS whith the expected length for ",
                        "protein: '", nm, "'. The returned genomic ",
                        "coordinates might thus not be correct for this ",
                        "protein.", call. = FALSE)
                ## Alternatively we could align the RNA and AA sequences and
                ## trim the protein sequence...
                z <- lapply(z, function(grng) {
                    mcols(grng)$cds_ok <- FALSE
                    grng
                })
                GRangesList(z[which.min(abs(diffs))])
            }
        }
    })
})

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
    ## TODO: preserve mcols.
    if (is(g_coords, "GRangesList"))
        return(GRangesList(
            lapply(g_coords, .to_genome, cds_coords = cds_coords)))
    if (!is(g_coords, "GRanges"))
        stop("'g_coords' is supposed to be a 'GRanges' object")
    cds_rel <- .splice(g_coords)
    strnd <- unique(as.character(strand(g_coords)))
    seqlvl <- unique(seqnames(g_coords))
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
    if (length(start_exon) == 0 | length(end_exon) == 0) {
        warning("The within transcript/CDS coordinates are outside the region ",
                "defined by the provided exons", call. = FALSE)
        return(GRanges())
    }
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
    genome_coords <- intersect(genome_coords, g_coords)
    ## Now grab all of the metadata columns from the second...
    mcols(genome_coords) <- mcols(g_coords[findOverlaps(genome_coords,
                                                        g_coords,
                                                        select = "first")])
    genome_coords
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
