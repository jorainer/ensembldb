loadEnsDb <- function( x ){
    ## con <- ensDb( x )
    ## EDB <- new( "EnsDb", ensdb=con )
    return( EnsDb( x ) )
}
