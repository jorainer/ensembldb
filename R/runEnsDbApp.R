## running the shiny web app.
runEnsDbApp <- function(...){
    if(requireNamespace("shiny", quietly=TRUE)){
        message("Starting the EnsDb shiny web app. Use Ctrl-C to stop.")
        shiny::runApp(appDir=system.file("shinyHappyPeople", package="ensembldb"), ...)
    }else{
        stop("Package shiny not installed!")
    }
}

