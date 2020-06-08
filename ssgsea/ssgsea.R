
#' single-sample gene set enrichment analysis
 # from Barbie et al., 2009
ssgsea <- function(E, S, alpha = .25, compute.p = 0L) {
	# get dimensions
	m <- NROW(E)
	n <- NCOL(E)

	# create output
	stat <- data.frame( es = NaN, p = NaN, E.es = NaN, V.es = NaN )
	stat <- stat[ rep(1L, n), ]
	rownames(stat) <- colnames(E)

	# compute the ssGSEA ES score
	get.es <- function(r, w.ind, pos) {
		lr <- length(r)
		lw <- length(w.ind)
		# NB. rolled out cumsums here..
		es <- (( sum( (lr-(w.ind-1L)) * ( r[w.ind]^alpha + pos/(lr-lw) ) ) + 
			-lr*(lr+1.)*.5 * pos/(lr-lw) )) / pos
	}

	# loop over each sample
	for (j in seq_len(n)) {
		# get ranks
		r <- sort( rank( E[ !is.na(E[, j]), j ], ties = 'average' ),
			decreasing = T )
		w.ind <- which( names(r) %in% S )

		# precompute constants
		lr <- length(r)
		lw <- length(w.ind)
		pos <- sum( r[w.ind] ^ alpha )

		# get observed ES
		stat[j, ]$es <- get.es( r, w.ind, pos )

		# get expected ES
		if (compute.p > 0L) {
			# Monte Carlo sample
			es.ref <- rep(NaN, compute.p)
			for (k in seq_len(compute.p))
				es.ref[k] <- get.es(r, sample(lr, lw), pos)

			# get statistics
			stat[j, ]$p <- .5*mean(sign( stat[j, ]$es - es.ref )) + .5
			stat[j, ]$E.es <- mean(es.ref)
			stat[j, ]$V.es <- var(es.ref)
		}
	}

	return (stat)
}
