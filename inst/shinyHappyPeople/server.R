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

checkSelectedPackage <- function(input){
    if(is.null(input$package)){
        return(FALSE)
    }else{
        require(input$package, character.only=TRUE)
        assign("edb", get(input$package), envir=globalenv())
        return(TRUE)
    }
}

## Define server logic required to draw a histogram
shinyServer(function(input, output) {

    ## generate the select field for the package...
    output$packages <- renderUI(
        selectInput("package", "Select EnsDb package", as.list(epacks))
    )

    ## generating metadata info.
    output$metadata_organism <- renderText({
        if(!checkSelectedPackage(input))
            return()
        paste0("Organism: ", organism(edb))
    })
    output$metadata_ensembl <- renderText({
        if(!checkSelectedPackage(input))
            return()
        md <- metadata(edb)
        rownames(md) <- md$name
        paste0("Ensembl version: ", md["ensembl_version", "value"])
    })
    output$metadata_genome <- renderText({
        if(!checkSelectedPackage(input))
            return()
        md <- metadata(edb)
        rownames(md) <- md$name
        paste0("Genome build: ", md["genome_build", "value"])
    })
    ## output$genename <- renderText({
    ##     if(!checkSelectedPackage(input))
    ##         return()
    ##     input$geneName
    ## })
    ## That's the actual queries for genes, transcripts and exons...
    output$Genes <- renderDataTable({
        if(!checkSelectedPackage(input))
            return()
        if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
            res <- genes(edb, filter=TheFilter(input),
                         return.type="data.frame")
            return(res)
        }
    })
    output$Transcripts <- renderDataTable({
        if(!checkSelectedPackage(input))
            return()
        if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
            res <- transcripts(edb, filter=TheFilter(input),
                         return.type="data.frame")
            return(res)
        }
    })
    output$Exons <- renderDataTable({
        if(!checkSelectedPackage(input))
            return()
        if(!is.na(input$geneName) & length(input$geneName) > 0 & input$geneName!=""){
            res <- exons(edb, filter=TheFilter(input),
                         return.type="data.frame")
            return(res)
        }
    })
})


