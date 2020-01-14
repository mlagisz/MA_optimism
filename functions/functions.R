
pacman::p_load(purrr, dplyr, tidyverse, magrittr) #fmsb

### Support functions for# "MA of cognitive bias (optimism) in animal studies"


####################### LATENCY
#' @title Calculating Cohens d and Hedges g for latency data
#' @description Function for calculating Cohens d and Hedges g for latency data using log-transformed means and SD, for between and within-subject study designs
#' @param data Dataframe object containing effect sizes, their variance, and sample sizes
#' @param Worse Name of the variable (vector) containing means of the "Worse" condition group
#' @param WorseSD Name of the variable (vector) containing standard deviations of the "Worse" condition group
#' @param WorseN Name of the variable (vector) containing sample sizes of the "Worse" condition group
#' @param Better Name of the variable (vector) containing means of the "Better" condition group
#' @param BetterSD Name of the variable (vector) containing standard deviations of the "Better" condition group
#' @param BetterN Name of the variable (vector) containing sample sizes of the "Better" condition group
#' @param type Optional parameter indicating whether we assume the data is log-normally distributed (default or "lnorm") or a delta method should be used ("delta").
#' @param WithinBetween Optional parameter indicating whether formula for between-subject study design should be used (default or "between") or a within-subject formula should be used ("within").
#' @param adjusted Optional logical parameter indicating whether d should be adjusted for small sample sizes, i.e. calculating Hedges d. Default value is "TRUE"
#' @export
#Value: adds columns with calculated d and Vd to the dataframe

calc_ES_latency <- function(data, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type=c("lnorm", "delta")){
  
  if (missing(data)) 
    stop("Must specify dataframe via 'data' argument.")
  if (missing(type)) 
    type <- "lnorm" 
  
  
  # turning parameters that are column names into strings
  Worse <- data[[deparse(substitute(Worse))]]
  WorseSD <- data[[ deparse(substitute(WorseSD))]]
  WorseN <- data[[deparse(substitute(WorseN))]]
  Better <- data[[deparse(substitute(Better))]]
  BetterSD <- data[[deparse(substitute(BetterSD))]]
  BetterN <- data[[deparse(substitute(BetterN))]]
  WithinBetween <- data[[deparse(substitute(WithinBetween))]]
  
  if(type == "delta"){  
    M1 <- log(Better) - (BetterSD^2)/(2*Better^2)
    M2 <- log(Worse) - (WorseSD^2)/(2*Worse^2)
    V1 <- (BetterSD^2)/(Better^2)
    V2 <- (WorseSD^2)/(Worse^2)    
  }
  
  if(type == "lnorm"){  
    M1 <- log(Better) - log( sqrt(1 + ( (BetterSD^2)/(Better^2) ) ) )
    M2 <- log(Worse) - log( sqrt(1 + ( (WorseSD^2)/(Worse^2) ) ) )
    V1 <- log( 1 + ( (BetterSD^2)/(Better^2) ) )
    V2 <- log( 1 + ( (WorseSD^2)/(Worse^2) ) )
  }
  
  N1 <- BetterN
  N2 <- WorseN  
  #calculate sPooled
  sPooled <- sqrt( ((N1-1) * V1 + (N2-1) * V2) / (N1 + N2 - 2) )
  
  #calculate Cohen's d and Hedges g (d adjusted for small sample sizes)
  d <- (M1-M2)/sPooled
  if(adjusted==TRUE){
    d <- d * (1-(3/ (4 * (N1 + N2 -2) -1))) #calculate Hedges d (adjusted for small N)
  }
  
  # use vectorization to switch between calcualting standard error for d with within or between-subject study design option
  se_choose <- function(N2, N1, WithinBetween, d) {
    if (WithinBetween == "between") {
      SEd <- sqrt(((N2 + N1) / (N2 * N1) ) + ( (d^2) / (2 * (N2 + N1 - 2))))
    } #for between-subject study design
    else {
      SEd <- sqrt(1/N2 + ((d^2) / (2 * (N2  - 1))))
    } #for within-subject study design
  }
  
  SEd <- pmap_dbl(list(N2, N1, WithinBetween, d), se_choose)
  
  #calculate Variance for d 
  Vd <- SEd^2
  
  #Reversing sign on latency data to match proportion data. 
  #Now positive values means optimistic for both lat and prop data.
  d <- (-d)   
  
  # putting all together in a tibble and combine
  es <- tibble(M1, V1, N1, M2, V2, N2, d, Vd)
  data <- bind_cols(data, es)
  return(data)  
}



###TESTS
# test <- tibble(Worse = c(1,1),
#                WorseSD = c(2,2),
#                WorseN = c(20, 20),
#                Better = c(1.1,1.1),
#                BetterSD = c(2, 2),
#                BetterN = c(20, 20),
#                WithinBetween = c("within", "between"))
# 
# calc_ES_latency(test, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="lnorm")
# test %<>% calc_ES_latency(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="lnorm")
# test2 <- test %>% calc_ES_latency(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="lnorm") 
# 



####################### PERCENT / LOGIT
#' @title Calculating Cohens d and Hedges g for proportion (or logit-transformed proportion) data
#' @description Function for calculating Cohens d and Hedges g for proportion data using logit-transformed means and SD, for between and within-subject study designs
#' @param data Dataframe object containing effect sizes, their variance and sample sizes
#' @param Worse Name of the variable (vector) containing means of the "Worse" condition group (proportions of positive responses)
#' @param WorseSD Name of the variable (vector) containing standard deviations of the "Worse" condition group
#' @param WorseN Name of the variable (vector) containing sample sizes of the "Worse" condition group
#' @param Better Name of the variable (vector) containing means of the "Better" condition group (proportions of positive responses)
#' @param BetterSD Name of the variable (vector) containing standard deviations of the "Better" condition group
#' @param BetterN Name of the variable (vector) containing sample sizes of the "Better" condition group
#' @param type Optional parameter indicating whether we assume the data is percent or logit.
#' @param WithinBetween Optional parameter indicating whether formula for between-subject study design should be used (default or "between") or a within-subject formula should be used ("within").
#' @param adjusted Optional logical parameter indicating whether d should be adjusted for small sample sizes, i.e. calculating Hedges d. Default value is "TRUE"
#' @export
#Value: adds columns with calculated d and Vd to the dataframe

calc_ES_proportion <- function(data, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type=c("proportion", "logit")){
  
  if (missing(data)) 
    stop("Must specify dataframe via 'data' argument.")
  if (missing(type)) 
    type <- "logit" 
  
  # turning parameters that are column names into strings
  Worse <- data[[deparse(substitute(Worse))]]
  WorseSD <- data[[ deparse(substitute(WorseSD))]]
  WorseN <- data[[deparse(substitute(WorseN))]]
  Better <- data[[deparse(substitute(Better))]]
  BetterSD <- data[[deparse(substitute(BetterSD))]]
  BetterN <- data[[deparse(substitute(BetterN))]]
  WithinBetween <- data[[deparse(substitute(WithinBetween))]]
  
  if(type == "proportion"){  
      M1 <- qlogis(Better) + ( (BetterSD^2/2) * ( (1/(1-Better)^2 ) - (1/Better^2) ) )
      M2 <- qlogis(Worse)  + ( (WorseSD^2)/2) * ( (1/(1-Worse)^2 ) - (1/(Worse^2) ) )
      V1 <- (BetterSD^2) * ( 1/Better + 1/(1-Better) )^2 
      V2 <- (WorseSD^2)  * ( 1/Worse +  1/(1-Worse) )^2 
  }

  if(type == "logit"){  
    M1 <- Better
    V1 <- BetterSD^2
    M2 <- Worse
    V2 <- WorseSD^2
  }
  
  N1 <- BetterN
  N2 <- WorseN  

  #calculate sPooled
  sPooled <- sqrt( ((N1-1) * V1 + (N2-1) * V2) / (N1 + N2 - 2) )
  
  #calculate Cohen's d and Hedges g (d adjusted for small sample sizes)
  d <- (M1-M2)/sPooled   # positive values means optimistic for both lat and prop data
  if(adjusted==TRUE){
    d <- d * (1-(3/ (4 * (N1 + N2 - 2) -1))) #calculate Hedges d (adjusted for small N)
  }
  
  # use vectorization to switch between calcualting standard error for d with within or between-subject study design option
  se_choose <- function(N2, N1, WithinBetween, d) {
    if (WithinBetween == "between") {
      SEd <- sqrt(((N2 + N1) / (N2 * N1) ) + ( (d^2) / (2 * (N2 + N1 - 2))))
    } #for between-subject study design
    else {
      SEd <- sqrt(1/N2 + ((d^2) / (2 * (N2  - 1))))
    } #for within-subject study design
  }
  
  SEd <- pmap_dbl(list(N2, N1, WithinBetween, d), se_choose)
  
  #calculate Variance for d 
  Vd <- SEd^2
  
  # putting all together in a tibble and combine
  es <- tibble(M1, V1, N1, M2, V2, N2, d, Vd)
  data <- bind_cols(data, es)
  return(data)  
}



###TESTS
# test <- tibble(Worse = c(0.5,0.5), 
#                WorseSD = c(0.5,0.5), 
#                WorseN = c(20, 20), 
#                Better = c(0.6,0.6), 
#                BetterSD = c(0.5, 0.5), 
#                BetterN = c(20, 20), 
#                WithinBetween = c("within", "between"))
# 
# calc_ES_proportion(test, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="proportion") 
# test %<>% calc_ES_proportion(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="proportion")
# test2 <- test %>% prepare_df(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="proportion") 



############################ VCV MATRIX

#' @title Variance-covariance and correlation matrix function basing on shared level ID
#' @description Function for generating simple variance-covariance and correlation matrices 
#' @param data Dataframe object containing effect sizes, their variance, unique IDs and clustering variable
#' @param V Name of the variable (vector) containing effect size variances variances
#' @param cluster Name of the variable (vector) indicating which effects belong to the same cluster. Same value of 'cluster' are assumed to be nonindependent (correlated).
#' @param obs Name of the variable (vector) containing individual IDs for each value in the V (Vector of variances). If this parameter is missing, label will be labelled with consecutive integers starting from 1.
#' @param rho Known or assumed correlation value among effect sizes sharing same 'cluster' value. Default value is 0.5.
#' @param type Optional logical parameter indicating whether a full variance-covariance matrix (default or "vcv") is needed or a correlation matrix ("cor") for the non-independent blocks of variance values.
#' @export
#Value: Labelled full variance-covariance or correlation matrice of the size and labels matching initial dataframe will be returned 

make_VCV_matrix <- function(data, V, cluster, obs, type=c("vcv", "cor"), rho=0.5){
  
  if (missing(data)) 
    stop("Must specify dataframe via 'data' argument.")
  if (missing(V)) 
    stop("Must specify name of the variance variable via 'V' argument.")
  if (missing(cluster)) 
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  if (missing(obs)) 
    obs <- 1:length(V)   
  if (missing(type)) 
    type <- "vcv" 
  
  new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
  rownames(new_matrix) <- data[ ,obs]
  colnames(new_matrix) <- data[ ,obs]
  # find start and end coordinates for the subsets
  shared_coord <- which(data[ ,cluster] %in% data[duplicated(data[ ,cluster]), cluster]==TRUE)
  # matrix of combinations of coordinates for each experiment with shared control
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]   #add the diagonal
  }
  
  if(type == "cor"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1   #add the diagonal of 1
  }
  
  return(new_matrix)
}




# General modeling functions 
# Functions for I2

#' Title Function to obtain total and separate I2 from multilevel-meta-analytic model
#'
#' @param model 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
I2 <- function(model, method = c("Wolfgang", "Shinichi")){
  
  ## evaluate choices
  method <- match.arg(method)
  
  # Wolfgang's method
  if(method == "Wolfgang"){
    W <- solve(model$V) 
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    I2_each  <- model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
    
    # or my way
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1/model$vi) * (model$k-1) / (sum(1/model$vi)^2 - sum((1/model$vi)^2)) 
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + sigma2_v) #s^2_t = total variance
    I2_each  <- model$sigma2 / (sum(model$sigma2) + sigma2_v)
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
  }
  return(I2s)
}

# test <- dataset$fit4.1[[3]]
# I2(test, method = "Wolfgang")
# I2(test, method = "Shinichi")


#' Title: R2 based on Nakagawa & Schielzeth 2013
#'
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
R2 <- function(model){
  #warning("Conditional R2 is not meaningful and the same as marginal R2\n")
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  R2
  #Rm <- round(100*R2m, 3)
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) / 
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  return(R2s)
}


#' Title: the function to get estimates from rma objects (metafor)
#'
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_est <- function (model, mod = " ") {
  
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub 
  
  table <- tibble(name = name, estimate = estimate, lowerCL = lowerCL, upperCL = upperCL)
}


#' Title: the function to get prediction intervals (crediblity intervals) from rma objects (metafor)
#'
#' @param model: rma.mv object 
#' @param mod: the name of a moderator 
get_pred <- function (model, mod = " ") {
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  len <- length(name)
  
  if(len != 1){
  newdata <- matrix(NA, ncol = len, nrow = len)
  for(i in 1:len) {
    # getting the position of unique case from X (design matrix)
    pos <- which(model$X[,i] == 1)[[1]]
    newdata[, i] <- model$X[pos,]
    }
  pred <- predict.rma(model, newmods = newdata)
  }
  else {
    pred <- predict.rma(model)
    }
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub 
  
  table <- tibble(name = name, lowerPR = lowerPR, upperPR = upperPR)
}

#Here are links for how to do confidence regions for rma.mv regression lines
#https://www.rdocumentation.org/packages/metafor/versions/1.9-9/topics/predict.rma
#https://stackoverflow.com/questions/50804464/out-of-sample-prediction-for-rma-object-in-metafor


#' Title: Contrast name geneator
#'
#' @param name: a vector of character strings
cont_gen <- function (name) {
  combination <- combn(name,2)
  name_dat <- t(combination)
  names <- paste(name_dat[ ,1], name_dat[, 2], sep = "-")
  return(names)
}



# ############################ SEQUENTIAL VIF
# 
# ## function for sequentially building the model that accounts for collinearity among the explanatory variables
# # followng https://www.r-bloggers.com/collinearity-and-stepwise-vif-selection/
# 
# vif_func<-function(in_frame,thresh=10,trace=T,...){
#   
#   library(fmsb)
#   
#   if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
#   
#   #get initial vif value for all comparisons of variables
#   vif_init<-NULL
#   var_names <- names(in_frame)
#   for(val in var_names){
#     regressors <- var_names[-which(var_names == val)]
#     form <- paste(regressors, collapse = '+')
#     form_in <- formula(paste(val, '~', form))
#     vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
#   }
#   vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
#   
#   if(vif_max < thresh){
#     if(trace==T){ #print output of each iteration
#       prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
#       cat('\n')
#       cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
#     }
#     return(var_names)
#   }
#   else{
#     
#     in_dat<-in_frame
#     
#     #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
#     while(vif_max >= thresh){
#       
#       vif_vals<-NULL
#       var_names <- names(in_dat)
#       
#       for(val in var_names){
#         regressors <- var_names[-which(var_names == val)]
#         form <- paste(regressors, collapse = '+')
#         form_in <- formula(paste(val, '~', form))
#         vif_add<-VIF(lm(form_in, data = in_dat, ...))
#         vif_vals<-rbind(vif_vals,c(val,vif_add))
#       }
#       max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
#       
#       vif_max<-as.numeric(vif_vals[max_row,2])
#       
#       if(vif_max<thresh) break
#       
#       if(trace==T){ #print output of each iteration
#         prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
#         cat('\n')
#         cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
#         flush.console()
#       }
#       
#       in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
#       
#     }
#     
#     return(names(in_dat))
#     
#   }
#   
# }
# 
# #TEST
# #rand.vars <- dplyr::select(dat, Captive_Wild.caught, Sex, Age, WithinBetween, Blind, Automated, FoodDeprived, TaskType, CueTypeCat, ReinforcementCat, AffectManipCat, AffectManipTiming, AmbigReinforced, MeasureType, ScalePoint)
# #vif_func(in_frame=rand.vars,thresh=5,trace=T) #does not work for factors!



