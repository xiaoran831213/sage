## check the extracted genes
source('src/dsg_tgz.R')
main <- function(fs = NULL)
{
    if(is.null(fs))
       fs <- dir('raw/GT/4', '.*tgz', full.names = T)
    N <- length(fs)
    ex <- list()
    for(i in 1L:N)
    {
        f <- fs[i]
        g <- sub('.tgz', '', basename(f))
        r <- try(readDTZ(f), silent = T)
        if(inherits(r, 'try-error'))
            ex[[g]] <- r
        cat('\r', sprintf("%05X/%05X, NE=%05X", i, N, length(ex)))
    }
    cat('\n')
    invisible(ex)
}
