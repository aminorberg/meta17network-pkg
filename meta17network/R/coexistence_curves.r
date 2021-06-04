#' Coexistence and species-area curves
#'
#' Calculate coexistence and species-area curves with different starting points
#' @param Y Virus community data
#' @param xy Matrix or data.frame of x and y coordinates
#' @param sample_ids Vector of row identifiers that matches with Y and xy
#' @param dirs Directories
#' @param nsim Number of simulations
#' @return List: 1) (True) coexistences always within the same plant, 2) pooled coexistences, 3) cumulative species richness 
#' @export

coexistence_curves <- function(Y, xy, sample_ids, dirs, nsim = 100) {

	# order plants according to their spatial proximity (first plant selected randomly)
	distsXY <- as.matrix(dist(xy, diag = FALSE))
	rownames(distsXY) <- sample_ids
	colnames(distsXY) <- rownames(distsXY)
	diag(distsXY) <- NA

	Y_dist_ord_ls <- list()
	startPlant <- list()
	for (n in 1:nsim) {
		startPlant[[n]] <- sample(1:nrow(Y), 1)
		Y_dist_ord <- Y[startPlant[[n]],]
		distsXY_tmp <- distsXY
		for (i in 1:(nrow(distsXY)-1)) {
			focalPlant <- which(distsXY_tmp[startPlant[[n]], ] == min(distsXY_tmp[startPlant[[n]], ], 
																  na.rm = TRUE))
			if (length(focalPlant) > 1) {
				focalPlant <- focalPlant[1]
			}
			distsXY_tmp[startPlant[[n]],focalPlant] <- NA
			Y_dist_ord <- rbind(Y_dist_ord, Y[focalPlant,])        
		}
	Y_dist_ord_ls[[n]] <- Y_dist_ord
	}
	Y_dist_ord_ls_mat <- lapply(Y_dist_ord_ls, as.matrix)

	insidSumS2 <- insidSumS <- NA
	sp_area <- insids2 <- insids <- list()
	for (m in 1:length(Y_dist_ord_ls)) {        
		tmp <- insids[[1]] <- as.matrix(Y_dist_ord_ls_mat[[m]][1,]) %*% t(Y_dist_ord_ls_mat[[m]][1,])
		tmp[upper.tri(tmp, diag = TRUE)] <- 0
		insidsSum <- sum(tmp > 0)
		tmp_comm_mat <- 1 * (apply(Y_dist_ord_ls_mat[[m]], 2, cumsum) > 0)
		insids2[[1]] <- as.matrix(tmp_comm_mat[1,]) %*% t(tmp_comm_mat[1,])
		tmp2 <- insids2[[1]]
		tmp2[upper.tri(tmp2, diag = TRUE)] <- 0
		insidsSum2 <- sum(tmp2 > 0)    
		for (n in 2:nrow(Y_dist_ord_ls_mat[[1]])) {
			tmp <- insids[[n]] <- t(Y_dist_ord_ls_mat[[m]][1:n,]) %*% Y_dist_ord_ls_mat[[m]][1:n,]        
			tmp[upper.tri(tmp, diag = TRUE)] <- 0
			tmp <- sum(tmp > 0)
			insidsSum <- c(insidsSum, tmp)        
			insids2[[n]] <- t(tmp_comm_mat[1:n,]) %*% tmp_comm_mat[1:n,]
			tmp2 <- insids2[[n]]
			tmp2[upper.tri(tmp2, diag = TRUE)] <- 0
			tmp2 <- sum(tmp2 > 0)
			insidsSum2 <- c(insidsSum2, tmp2)
		}
	sp_area[[m]] <- apply(Y_dist_ord_ls_mat[[m]], 2, cumsum)
	insidSumS <- cbind(insidSumS, insidsSum)
	insidSumS2 <- cbind(insidSumS2, insidsSum2)
	}
	# insidSumS: only co-occurrence inside the same plant is counted as coexistence
	# insidSumS2: also co-occurrence between plants is counted as coexistence

	tmp <- lapply(lapply(sp_area, function(x){(x > 0) * 1}), rowSums)
	sp_area <- simplify2array(tmp)

	res <- list("coex_within_hosts" = insidSumS[,-1], 
			 	"coex_pooled" = insidSumS2[,-1], 
			 	"sps_area" = sp_area)
	return(res)

}
