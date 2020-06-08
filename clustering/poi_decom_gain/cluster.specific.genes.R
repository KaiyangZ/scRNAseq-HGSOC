#!/usr/bin/env Rscript

cluster.specific.genes <- function(X, D, g, Y, inds1, inds2 ) {

    centroid.1 <- poi_decom_control_centroid(X %*% D[, inds1], Y[, inds1, drop = F], g[inds1])
    centroid.2 <- poi_decom_control_centroid(X %*% D[, inds2], Y[, inds2, drop = F], g[inds2])

    centroid <- poi_decom_control_centroid(X %*% D[, c(inds1, inds2)], Y[, c(inds1, inds2), drop = F], g[c(inds1, inds2)])

    cost <- function(X, Y, D, g, inds, centroid){
      return(rowSums(poi_dist( .diag.mul.right(( X %*% D[, inds] + centroid), g[inds]), Y[, inds, drop = F])))
    }

    cost1 <- cost( X, Y, D, g, inds1, centroid.1) + cost( X, Y, D, g, inds2, centroid.2)
    cost0 <- cost( X, Y, D, g, c(inds1, inds2), centroid)

    centroidDiff <- centroid.1 - centroid.2

    n <- length(inds1) + length(inds2)
    p1 <- 2L
    p0 <- 1L

    rss2 <- cost1
    rss1 <- cost0
    df2 <- n-p1
    df1 <- p1-p0

    # get F-statistic
    f <- ( (rss1-rss2)/df1 / ( rss2/df2 ) )
    #f <- ((( length(inds1) * (centroid.1 -centroid)^2)  + ( length(inds2) * (centroid.2 -centroid)^2))/df1) / ( rss2/df2 ) 

    f[!(df1 > 0 & df2 > 0)] <- NaN

    pf <- pf(f, df1 = df1, df2 = df2, lower.tail = F, log.p = F)

    lik <- (rss1-rss2)/df1
    lik[!(df1 > 0 & df2 > 0)] <- NaN

    p.val <- pchisq( 2*lik , df1 , lower.tail = F, log.p = F ) # cumulative probability
    p.adj <- p.adjust(p.val, method = "BH")

		chisq.compare <- function(x1, x2) {
			# this uses a chi-squared test to compare if two discrete samples x1, x2
			 # follow an equal distribution determined by the epmf of c(x1, x2)
			make.table <- function(x, p) {
				y <- rep(0L, p)
				for (j in seq_along(x))
					y[x[j]] <- y[x[j]] + 1L
				return (y)
			}
			stat <- function(t1, t0) {
				O <- t1
				E <- t0 / sum(t0) * sum(t1)
				v <- O-E; v <- v*v   # NB. v <- (O-E)^2/E
				v[v>0] <- v[v>0] / E[v>0]
				return (v)
			}
			p <- max( x1+1L, x2+1L )
			t1 <- make.table( x1+1L, p )
			t2 <- make.table( x2+1L, p )
			t12 <- t1 + t2
			return (pchisq( sum(stat(t1, t12)) + sum(stat(t2, t12)), sum(t1 > 0L)-1L + sum(t2 > 0L)-1L, lower.tail = F))
		}
                                
		# p-value for data being equal-- WARN: does not take model into account
		Y.eq.p.val <- rep(NaN, nrow(Y))
		for (i in seq_len(nrow(Y)))
			Y.eq.p.val[i] <- chisq.compare(Y[i, inds1], Y[i, inds2])

		# p-value for centroid being significant
    costZ <- cost( X, Y, D, g, c(inds1, inds2), 0.)
    Z.p.val <- pchisq( 2*(costZ-cost0) , 1L , lower.tail = F, log.p = F )

		 # p-value for centroid being significant in *either* (or *both*)
		Z.p.val.e <- pchisq( 2*(costZ-cost1) , 2L, lower.tail = F, log.p = F)

		 # p-value for centroid being significant in *both* 
		 Z.p.val.b <- pmax( pchisq( 2*( cost(X, Y, D, g, inds1, 0.) - cost(X, Y, D, g, inds1, centroid.1) ), 1L, lower.tail = F, log.p = F ),
		 	pchisq( 2*( cost(X, Y, D, g, inds2, 0.) - cost(X, Y, D, g, inds2, centroid.2) ), 1L, lower.tail = F, log.p = F) )

    return (data.frame(fstats = f, lik = lik, centroidDiff = centroidDiff, centroid.1 = centroid.1, centroid.2 = centroid.2, centroid = centroid, p.val = p.val, p.adj = p.adj, Y.eq.p.val = Y.eq.p.val, Z.p.val = Z.p.val, Z.p.val.e = Z.p.val.e, Z.p.val.b = Z.p.val.b, stringsAsFactors = F))
}



