# Load required libraries
library(sf)
library(units)
library(betapart)
library(mgcv)

# Get spatial distances between sample locations -------------------
get_dists <- function(IBA_locs){
  # Select distinct trap IDs from the IBA_locs dataset, ensuring no missing trapID values
  IBA_samples <- IBA_locs |> 
    drop_na(trapID) |>          # Remove rows with NA in trapID
    select(trapID) |>           # Select only the trapID column
    distinct()                  # Keep only distinct trapIDs
  
  # Calculate pairwise distances between all traps
  pDist <- st_distance(IBA_samples$geometry , IBA_samples$geometry)
  pDist <- set_units(pDist , "km")  # Convert distances to kilometers
  
  # Create a data frame of pairwise distances for all unique trap pairs
  ind       <- which( lower.tri(pDist,diag=FALSE) , arr.ind = TRUE )
  trap_dist <- data.frame( trap_1 = IBA_samples$trapID[ind[,1]],
                           trap_2 = IBA_samples$trapID[ind[,2]],
                           distance = pDist[ind]) |>
    mutate(trap_ID = interaction(trap_1 , trap_2 , sep = "_")) |> 
    droplevels()
  
  return(trap_dist)  # Return the data frame with distances
}

# Partition beta diversity between samples --------------
partition_beta_diversity <- function(sp_matrix_total , trap_dist) {
  # Calculate beta diversity components
  bpart_core <- betapart.core(sp_matrix_total)  # Core calculation of beta diversity
  beta_dist <- beta.pair(bpart_core , index.family = "jaccard")  # Pairwise beta diversity using Jaccard index
  
  # Convert Jaccard turnover component to matrix and extract lower triangle
  IBA_jtu <- beta_dist$beta.jac |>as.matrix()  
  ind_jtu <- which( lower.tri(IBA_jtu,diag=FALSE ), arr.ind = TRUE )
  
  # Create a data frame with Jaccard turnover values for each trap pair
  IBA_turnover <- data.frame( 
    
    trap_1 = rownames(sp_matrix_total)[ind_jtu[,1]],
    trap_2 = rownames(sp_matrix_total)[ind_jtu[,2]],
    jaccard = IBA_jtu[ind_jtu]) |>
    mutate(trap_ID = interaction(trap_1 , trap_2 , sep = "_")) |> droplevels()
  
  # Merge the trap distance data with the Jaccard turnover data
  IBA_betapart <- inner_join(trap_dist , IBA_turnover , by = "trap_ID")
  IBA_betapart <- drop_units(IBA_betapart)  # Remove units from the distance column for simplicity
  
  return(IBA_betapart)  # Return the merged data frame
}

# Fit spline with monotonic assumption 
monotonic_gam <- function(beta_dist,nK=5){
  # Fits a single cubic regression spline to the distance / turnover data
  
  # predict from unconstrained GAM fit
  fitUC <- gam(jaccard~s(distance,k=nK,bs="cr"),data=beta_dist)
  
  # Construct smooth term's design matrix
  # Find linear constraints sufficient for monotonicity of a cubic regression spline
  sm <- smoothCon(s(distance,k=nK,bs="cr"),beta_dist,knots=NULL)[[1]]
  DM <- mono.con(sm$xp) # Monotonic constraint
  
  # Construct necessary list of info for pcls() function
  G <- list(
    X=sm$X,                               # design matrix of smooth term
    C=matrix(0,0,0),                      # [0 x 0] matrix (no equality constraints)
    sp=fitUC$sp,                          # smoothing parameter estimates (taken from unconstrained model)
    p=sm$xp,                              # array of feasible initial parameter estimates
    y = beta_dist$jaccard,                # response vector         
    w =rep(1 , times = nrow(beta_dist)),  # weights for data  (all set to 1)  
    Ain = DM$A,                           # matrix for the inequality constraints
    bin = DM$b,                           # vector for the inequality constraints
    off = 0,                              # Offset values   
    S = sm$S                              # list of penalty matrices
  )
  
  # Use penalised constrained least squares fitting to get the monotonic fit.
  p <- pcls(G)    
  
  # predict 
  outData <- data.frame(distance = seq(1 , max(beta_dist$distance) , l = 1e3))
  outData$pred_fit <- Predict.matrix(sm, data.frame(distance = outData$distance)) %*% p
  return(outData)
  
}
 
