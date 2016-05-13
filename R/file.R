# File compatibility version. MUST be updated whenever a package change might
# render the previous saved files unusable.
.COMPAT = "0.0.1-alpha1"

#' File persistence of simulations
#' 
#' Functions to save and load simulated data from the disk
#'
#' The function \code{GillesComToFile} saves all of the relevant
#' data about a simulated community to a file, including the random
#' number generator seed. This allows for retaking the simulation later
#' or even in other computers with reproducible results.
#' 
#' The function \code{GillesComFromFile} loads all of the relevant
#' data about a simulated community from a file. It returns the saved random
#' number generator seed as an invisible vector, so if you want to continue the
#' simulation with the saved key, use ".Random.seed <- GillesComFromFile()".
#' @export
#' @rdname file
#' @param file filename to which save / from which load the data
GillesComToFile <- function(file="GillesCom.rda") {
    save_int = save_int()
    # Checking save_int (as it cannot be set zero)
    if (save_int == 0) stop("You must first initialize the community!")
    abundance = abundance()
    compat = .COMPAT
    seed = .Random.seed # Must be GLOBALLY assigned to restore the seed generator
    date = Sys.time()
    history = history()
    interaction = interaction()
    K = K()
    d0 = d0()
    birth = birth()
    migration = migration()
    time = elapsed_time()
    save(compat, seed, date, abundance, history, interaction, K, d0, birth, migration, time, save_int, file=file)
    cat("File saved. Size:", format.h(file.info(file)$size), "\n")
}

#' @export
#' @rdname file
GillesComFromFile <- function(file="GillesCom.rda") {
    # to avoid NOTEs at R check:
    compat <- NULL; seed <- NULL;
    load(file=file)
    cat("File loaded. Size:", format.h(file.info(file)$size), "\nSimulation date/time:", format.Date(date),"\n")
    if (compat != .COMPAT) warning("NOTE: Incompatible file type!\nExpected ", .COMPAT, ", got ", compat)
    load_community(abundance, history, interaction, K, d0, birth, migration, time, save_int)
    return(invisible(seed));
}


format.h <- function (x) {
    if (x < 1024) return(paste(x, "bytes"))
    if (x < 1024*1024) return (paste(floor(x/1024), "Kbytes"))
    if (x < 1024*1024*1024) return (paste(floor(x/1024/1024), "Mbytes"))
    #else 
    return (paste(floor(x/1024/1024/1024), "Gbytes"))
}

