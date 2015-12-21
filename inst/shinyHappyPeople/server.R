## list all packages...
packs <- installed.packages()
epacks <- packs[grep(packs, pattern="^Ens")]

## library(EnsDb.Hsapiens.v75)
## edb <- EnsDb.Hsapiens.v75

TheFilter <- function(input){
    Cond <- input$condition
    ## check if we've got something to split...
    Vals <- input$geneName
    ## check if we've got ,
    if(length(grep(Vals, pattern=",")) > 0){
        ## don't want whitespaces here...
        Vals <- gsub(Vals, pattern=" ", replacement="", fixed=TRUE)
        Vals <- unlist(strsplit(Vals, split=","))
    }
    if(length(grep(Vals, pattern=" ", fixed=TRUE)) > 0){
        Vals <- unlist(strsplit(Vals, split=" ", fixed=TRUE))
    }
    if(input$type=="Gene name"){
        return(GenenameFilter(Vals, condition=Cond))
    }
    if(input$type=="Chrom name"){
        return(SeqnameFilter(Vals, condition=Cond))
    }
    if(input$type=="Gene biotype"){
        return(GenebiotypeFilter(Vals, condition=Cond))
    }
    if(input$type=="Tx biotype"){
        return(TxbiotypeFilter(Vals, condition=Cond))
    }
}

## checkSelectedPackage <- function(input){
##     if(is.null(input$package)){
##         return(FALSE)
##     }else{
##         require(input$package, character.only=TRUE)
##         message("Assigning ", input$package, " to variable edb.")
##         assign("edb", get(input$package), envir=globalenv())
##         return(TRUE)
##     }
## }

## Based on the given EnsDb package name it loads the library and returns
## the object.
getEdb <- function(x){
    require(x, character.only=TRUE)
    return(get(x))
}

## Define server logic required to draw a histogram
shinyServer(function(input, output) {

    ## Generate the select field for the package...
    output$packages <- renderUI(
        selectInput("package", "Select installed EnsDb package", as.list(epacks))
    )

    selectedPackage <- reactive({
        names(epacks) <- epacks
        ## epacks <- sapply(epacks, as.symbol)
        ## load the package.
        if(length(input$package) > 0){
            require(input$package, character.only=TRUE)
            ## Actually, should be enough to just return input$package...
            ##return(switch(input$package, epacks))
            return(getEdb(input$package))
        }else{
            return(NULL)
        }
    })

    ## Metadata infos
    output$metadata_organism <- renderText({
        edb <- selectedPackage()
        if(!is.null(edb)){
            ## db <- getEdb(edb)
            paste0("Organism: ", organism(edb))
        }

    })
    output$metadata_ensembl <- renderText({
        edb <- selectedPackage()
        if(!is.null(edb)){
            ## db <- getEdb(edb)
            md <- metadata(edb)
            rownames(md) <- md$name
            paste0("Ensembl version: ", md["ensembl_version", "value"])
        }
    })
    output$metadata_genome <- renderText({
        edb <- selectedPackage()
        if(!is.null(edb)){
            ## db <- getEdb(edb)
            md <- metadata(edb)
            rownames(md) <- md$name
            paste0("Genome build: ", md["genome_build", "value"])
        }
    })

    output$genename <- renderText({
        if(length(input$geneName) > 0){
            input$geneName
        }else{
            return()
        }
    })
    ## That's the actual queries for genes, transcripts and exons...
    output$Genes <- renderDataTable({
        ## if(!checkSelectedPackage(input))
        ##     return()
        if(length(input$package) == 0)
            return
        if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
            edb <- selectedPackage()
            res <- genes(edb, filter=TheFilter(input),
                         return.type="data.frame")
            assign(".ENS_TMP_RES", res, envir=globalenv())
            return(res)
        }
    })
    output$Transcripts <- renderDataTable({
        if(length(input$package) == 0)
            return
        if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
            edb <- selectedPackage()
            res <- transcripts(edb, filter=TheFilter(input),
                         return.type="data.frame")
            assign(".ENS_TMP_RES", res, envir=globalenv())
            return(res)
        }
    })
    output$Exons <- renderDataTable({
        if(length(input$package) == 0)
            return
        if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
            edb <- selectedPackage()
            res <- exons(edb, filter=TheFilter(input),
                         return.type="data.frame")
            assign(".ENS_TMP_RES", res, envir=globalenv())
            return(res)
        }
    })
    observe({
        if(input$closeButton > 0){
            ## OK, now, gather all the data and return it in the selected format.
            edb <- selectedPackage()
            resType <- input$returnType
            resTab <- input$resultTab
            res <- NULL
            ## If result type is data.frame we just return what we've got.
            if(resType == "data.frame"){
                res <- get(".ENS_TMP_RES")
            }else{
                ## Otherwise we have to fetch a little bit more data, thus, we perform the
                ## query again and return it as GRanges.
                if(resTab == "Genes")
                    res <- genes(edb, filter=TheFilter(input), return.type="GRanges")
                if(resTab == "Transcripts")
                    res <- transcripts(edb,filter=TheFilter(input), return.type="GRanges")
                if(resTab == "Exons")
                    res <- exons(edb,filter=TheFilter(input), return.type="GRanges")
            }
            rm(".ENS_TMP_RES", envir=globalenv())
            stopApp(res)
        }
    })
})



## ## Define server logic required to draw a histogram
## shinyServer(function(input, output) {

##     ## generate the select field for the package...
##     output$packages <- renderUI(
##         selectInput("package", "Select EnsDb package", as.list(epacks))
##     )

##     ## generating metadata info.
##     output$metadata_organism <- renderText({
##         if(!checkSelectedPackage(input))
##             return()
##         paste0("Organism: ", organism(edb))
##     })
##     output$metadata_ensembl <- renderText({
##         if(!checkSelectedPackage(input))
##             return()
##         md <- metadata(edb)
##         rownames(md) <- md$name
##         paste0("Ensembl version: ", md["ensembl_version", "value"])
##     })
##     output$metadata_genome <- renderText({
##         if(!checkSelectedPackage(input))
##             return()
##         md <- metadata(edb)
##         rownames(md) <- md$name
##         paste0("Genome build: ", md["genome_build", "value"])
##     })
##     ## output$genename <- renderText({
##     ##     if(!checkSelectedPackage(input))
##     ##         return()
##     ##     input$geneName
##     ## })
##     ## That's the actual queries for genes, transcripts and exons...
##     output$Genes <- renderDataTable({
##         if(!checkSelectedPackage(input))
##             return()
##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
##             res <- genes(edb, filter=TheFilter(input),
##                          return.type="data.frame")
##             return(res)
##         }
##     })
##     output$Transcripts <- renderDataTable({
##         if(!checkSelectedPackage(input))
##             return()
##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
##             res <- transcripts(edb, filter=TheFilter(input),
##                          return.type="data.frame")
##             return(res)
##         }
##     })
##     output$Exons <- renderDataTable({
##         if(!checkSelectedPackage(input))
##             return()
##         if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
##             res <- exons(edb, filter=TheFilter(input),
##                          return.type="data.frame")
##             return(res)
##         }
##     })
##     ## observe({
##     ##     if(input$Close > 0){
##     ##         stopApp("AAARGHHH")
##     ##     }
##     ## })
## })


