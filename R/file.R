# File compatibility version. MUST be updated whenever a package change might
# render the previous saved files unusable.
.COMPAT = "0.1.2"

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
#' data about a simulated community from a file. If you also want to
#' restore the number generator seed, use ".Random.seed <- SeedFromFile()".
#' @param file filename to which save / from which load the data
#' @param community An object of class \code{\link{Community}}
#' @examples
#' filename = tempfile()
#' # Generate a new community
#' a = Init_Community(50); bdm(a, 10000, progress="none")
#' # Saves the simulated data
#' GillesComToFile(a, filename)
#' # Loads the simulated data
#' b = GillesComFromFile(filename)
#' # Notice that a and b are different objects...
#' identical (a,b)
#' # But their content is the same
#' identical (a$trajectories, b$trajectories)
#' # To load the random seed as well...
#' .Random.seed <- SeedFromFile(filename)
#' @export
#' @rdname file
GillesComToFile <- function(community, file="GillesCom.rda") {
    save_int = tryCatch( community$save_int, 
                        error = function(e) stop("You must initialize the community! See ?Init_Community")
                        )
    abundance = community$abundance

    compat = .COMPAT
    seed = .Random.seed # Must be GLOBALLY assigned to restore the seed generator
    date = Sys.time()
    trajectories = community$trajectories
    interaction = community$interaction
    stochastic = community$stochastic
    K = community$K
    d0 = community$d0
    b = community$b
    m = community$m
    time = community$time
    cycles = community$cycles
    save(compat, seed, date, abundance, trajectories, interaction, K, d0, b, m, time, save_int, cycles, stochastic, file=file)
    cat("File saved. Size:", format.h(file.info(file)$size), "\n")
}

#' @export
#' @rdname file
GillesComFromFile <- function(file="GillesCom.rda") {
    # to avoid NOTEs at R check:
    abundance <- NULL; trajectories <- NULL; interaction <- NULL; K <- NULL; d0 <- NULL
    b <- NULL; m <- NULL; time <- NULL; save_int <- NULL; 
    compat <- NULL; seed <- NULL; cycles <- NULL; stochastic <- NULL;
    # Does the actual loading
    load(file=file)
    cat("File loaded. Size:", format.h(file.info(file)$size), "\nSimulation date/time:", format.Date(date),"\n")
    if (compat != .COMPAT) warning("NOTE: Incompatible file type!\nExpected ", .COMPAT, ", got ", compat)
    a = new (Community, abundance, interaction, save_int)
    a$b = b; a$cycles = cycles; a$d0 = d0; a$K = K; a$m = m; a$stochastic = stochastic; a$time = time
    a$trajectories = trajectories
    return(a);
}

#' @export
#' @rdname file
SeedFromFile <- function(file="GillesCom.rda") {
    # to avoid NOTEs at R check:
    seed <- NULL
    # Does the actual loading
    load(file=file)
    return(seed);
}

format.h <- function (x) {
    if (x < 1024) return(paste(x, "bytes"))
    if (x < 1024*1024) return (paste(floor(x/1024), "Kbytes"))
    if (x < 1024*1024*1024) return (paste(floor(x/1024/1024), "Mbytes"))
    #else 
    return (paste(floor(x/1024/1024/1024), "Gbytes"))
}

