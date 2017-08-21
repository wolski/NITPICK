#' calculate accurate isotopic masses and abundances
#' 
#' \code{ams.mercury.mercury} uses the mercury7 algortihm for
#' the quick and accurate calculation of isotopic masses and
#' abundances
#' 
#' 
#' @usage ams.mercury.mercury(composition, charge, limit)
#' @param composition a vector of integers; needs to have \code{length==MAXLENGTH} of the libmercury++ library.
#'  This is a current shortcoming that is going to be resolvedin the near future.
#' @param charge  a single integer giving the desired charge
#' @param limit a pruning limit, recommended values are between 1e-25 and 1e-30
#' @details See the C++ source
#' @return 	Returns a two-column matrix. The first columns holds the m/z values, the second the corrsponding abundances.
#' @references Rockwood AL, Haimi P, Efficient calculation of accurate masses of isotopic peaks, J Am Soc Mass Spectrom, 17, 415-419, (2006).
#' @author  Marc Kirchner
#' @useDynLib Mercury
#' @export
#' @examples
#' 
#' # CH4
#' x <- ams.mercury.mercury(c(4,1,0,0,0), 1, 0);
#' plot(x[,1], x[,2], type="h")
#' x <- ams.averagine(1000)
#' 
#' x <- ams.mercury.mercury(unlist(x$model), 1, 0);
#' x <- ams.averagine(3000)
#' x <- ams.mercury.mercury(unlist(x$model), 1, 0);
#' plot(x[,1], x[,2], type="h")
#' 
ams.mercury.mercury <- function(composition, charge, limit)
{
	# currently. we expect H-C-N-O-S
	stopifnot(length(composition)==5);
	tmp <- .Call("Rmercury", 
		as.numeric(composition), 
		as.integer(charge), 
		as.numeric(limit),
		PACKAGE="Mercury"
	);
	return(cbind(tmp$mz, tmp$abundance));		
}

