"ams.mercury.mercury" <-
function(composition, charge, limit)
{
	# currently. we expect H-C-N-O-S
	stopifnot(length(composition)==5);
	tmp <- .Call("Rmercury", 
		as.numeric(composition), 
		as.integer(charge), 
		as.numeric(limit),
		PACKAGE="amsmercury"
	);
	return(cbind(tmp$mz, tmp$abundance));		
}

