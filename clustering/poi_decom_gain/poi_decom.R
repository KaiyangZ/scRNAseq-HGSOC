
#
# fits a model:
#   Y ~ Poi( (( X %*% D + Z %*% C )) %*% G )
# where Y is m-by-n, X is m-by-r, D is r-by-n, Z is m-by-k, C is k-by-n, G = diag(g) is n-by-n
#
# Y are output profiles: total expression
# X are nuisance: 
# D is the design for nuisances: weights matrix
# Z are latent profiles: 
# C are latent classes (binary)
# G are sample gains (e.g. total RNA)
#
# m is features (e.g. genes)
# n is samples
# r is categories in nuisance
# k is categories in latent
#

poi_decom_control <- function(Y, D, k, g = 1., w = 1., max.iter = 100L, min.delta = 1.48e-9) {
	# get dimensions
	m <- NROW(Y)
	n <- NCOL(Y)
	r <- NROW(D)

#cat(sprintf('m (genes) = %d, n (samps) = %d, r (design rank) = %d, k (clusters) = %d', m, n, r, k), sep = '\n', file = stderr())
	
	# check inputs
	stopifnot(identical(dim(Y), c(m, n)))
	stopifnot(identical(dim(D), c(r, n)))
	stopifnot(is.vector(g) && length(g) == n)

	# if we have a weighted problem, convert to unweighted
	if (!all(w == 1.)) {
		Y <- .diag.mul.right(Y, w)
		g <- g * w
	}

	# make initial guess
	t_X_Z <- matrix(1., r+k, m)
	L <- sample(k, n, replace = T)

	# precompute stuff
	Y <- Y + 0.   # NB. cast Y to double
	t_Y <- t(Y)
	D_C <- rbind( D, matrix(NaN, k, n) )

	# loop
	sum.res <- Inf
	for (iter in seq_len(max.iter)) {
		# solve (X, Z)
		D_C[seq_len(k)+r, ] <- outer( 1:k, L, `==` )
		t_X_Z <- poi_decom_control_solve(t_X_Z, .diag.mul.right(D_C, g), t_Y, 1L, min.delta)
	
		# relabel
		res <- poi_decom_control_label(D_C, t_X_Z, Y, g = g, s = r)
		L <- res$L

		# update cost
		last.sum.res <- sum.res
		sum.res <- res$cost
#cat(sprintf('iter = %d, sum.res = %g (last = %g), delta = %-g', iter, sum.res, last.sum.res, sum.res - last.sum.res), sep = '\n', file = stderr())
		if (!(last.sum.res - sum.res > min.delta))
			# early exit
			break;
	}

	# solve final parameters
	t_X_Z <- poi_decom_control_solve(t_X_Z, .diag.mul.right(D_C, g), t_Y)

	return (list(X = t(t_X_Z[seq_len(r), , drop = F]), Z = t(t_X_Z[seq_len(k)+r, , drop = F]), L = L, res = sum.res))
}

.diag.mul.right <- function(A, d)
	return (t(t(A) * d))

dyn.load('poi_decom.so')

# solve B ~ Poi( A %*% X ) for X 
poi_decom_control_solve <- function(X, t_A, B, max.iter = 100L, min.delta = 1.48e-9)
	return (.Call('poi_decom_control_solve_R', X, t_A, B, max.iter, min.delta))

# solve B ~ Poi( (( A[,1:s] %*% X[1:s,] + X[s+L,] )) %*% G ) for L
poi_decom_control_label <- function(X, t_A, B, g = 1., s = 0L) {
	L <- rep(0L, NCOL(B))
	cost <- .Call('poi_decom_control_label_R', L, X, t_A, B, g * rep(1., NCOL(B)), s)
	return (list(L = L, cost = cost))
}

# this is analogous to sqared distance in normal
# not symmetric !!! must put the data & rates correct way
# @param A rate from the model (like X %*% D %*% G + Z %*% C %*% G)
# @param B readouts (data)
poi_dist <- function(A, B) {
	lhs <- B * -log(A / B)
	lhs[!(B > 0.)] <- 0.
	return (lhs - (B - A))
}

# this computes a centroid for any subset of the data
# this is expensive-- just use res$Z if you want the leave centroids after estimation
# @param A effective rate for the control part (like X %*% D)
# @param B readouts (data)
poi_decom_control_centroid <- function(A, B, g = 1., w = 1.) {
	if (!all(w == 1.)) {
		B <- .diag.mul.right(B, w)
		g <- g * w
	}
	return (.Call('poi_decom_control_centroid_R', t(A), t(B), g * rep(1., NCOL(B))))
}

# select a model using BIC
#BIC = log(n) * k - 2*log()
bic_select <- function(Y, g, w = NULL, fits = NULL) {
    # get number of samples
    samps <- length(Y)
    # extract parameter count & RSSs from the fits
    if (!is.null(fits)) {
        # count the parameters
        dofs <- sapply(fits, function(f) {
            dofs <- 0L
            if (!is.null(f$X))
                dofs <- dofs + length(f$X)
            if (!is.null(f$Z))
                dofs <- dofs + length(f$Z)

	    if (!is.null(w))
		dofs <- dofs + length(w)
	    dofs <- dofs + length(g)
            return (dofs)
        })
    }
    #bics <- log(samps)* (dofs + 2) + 2*sapply(fits, function(x) x$res)
    bics <- log(samps)* (dofs) + 2*sapply(fits, function(x) x$res)
    return(bics)
}


# compute pairwise distances between each sample
# @param A effective rate for the control part (like X %*% D)
# @param B readouts (data)
poi_cellDist <- function(A, B, g = 1., w = 1., gw = 1.) {
	if (!all(w == 1.)) {
		B <- .diag.mul.right(B, w)
		g <- g * w
	}
	return (.Call('poi_cellDist_R', t(A), t(B), g * rep(1., NCOL(B)), gw))
}




