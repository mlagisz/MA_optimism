
pacman::p_load(purrr, dplyr, tidyverse, magrittr)

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
  
  if(type == "lnorm"){  
    M1 <- log(Better) - (BetterSD^2)/(2*Better^2)
    V1 <- (BetterSD^2)/(Better^2)
    M2 <- log(Worse) - (WorseSD^2)/(2*Worse^2)
    V2 <- (WorseSD^2)/(Worse^2)    
  }
  
  if(type == "delta"){  
    M1 <- log(Better) - log( sqrt(1 + (BetterSD^2)/(Better^2) ) )
    V1 <- log( 1 + (BetterSD^2)/(Better^2) )
    M2 <- log(Worse) - log( sqrt(1 + (WorseSD^2)/(Worse^2) ) )
    V2 <- log( 1 + (WorseSD^2)/(Worse^2) )
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
test <- tibble(Worse = c(1,1), 
               WorseSD = c(2,2), 
               WorseN = c(20, 20), 
               Better = c(1.1,1.1), 
               BetterSD = c(2, 2), 
               BetterN = c(20, 20), 
               WithinBetween = c("within", "between"))

calc_ES_latency(test, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="lnorm") 
#test %<>% calc_ES_latency(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="lnorm")
#test2 <- test %>% calc_ES_latency(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="lnorm") 




####################### PERCENT / LOGIT
#' @title Calculating Cohens d and Hedges g for proportion (percent or logit) data
#' @description Function for calculating Cohens d and Hedges g for proportion data using logit-transformed means and SD, for between and within-subject study designs
#' @param data Dataframe object containing effect sizes, their variance and sample sizes
#' @param Worse Name of the variable (vector) containing means of the "Worse" condition group
#' @param WorseSD Name of the variable (vector) containing standard deviations of the "Worse" condition group
#' @param WorseN Name of the variable (vector) containing sample sizes of the "Worse" condition group
#' @param Better Name of the variable (vector) containing means of the "Better" condition group
#' @param BetterSD Name of the variable (vector) containing standard deviations of the "Better" condition group
#' @param BetterN Name of the variable (vector) containing sample sizes of the "Better" condition group
#' @param type Optional parameter indicating whether we assume the data is percent or logit.
#' @param WithinBetween Optional parameter indicating whether formula for between-subject study design should be used (default or "between") or a within-subject formula should be used ("within").
#' @param adjusted Optional logical parameter indicating whether d should be adjusted for small sample sizes, i.e. calculating Hedges d. Default value is "TRUE"
#' @export
#Value: adds columns with calculated d and Vd to the dataframe

calc_ES_proportion <- function(data, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type=c("percent", "logit")){
  
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
  
  if(type == "percent"){  
      Better <- Better/100 #convert form percentage to proportion
      BetterSD <- BetterSD/100 #convert form percentage to proportion
      Worse <- Worse/100 #convert form percentage to proportion
      WorseSD <- WorseSD/100 #convert form percentage to proportion
      M1 <- qlogis(Better) + (BetterSD^2)/2 * (1/((1-Better)^2) -1/(Better^2))
      V1 <- (BetterSD^2) * ( 1/Better + 1/(1-Better)^2 )
      M2 <- qlogis(Worse) + (WorseSD^2)/2 * (1/((1-Worse)^2) -1/(Worse^2))
      V2 <- (WorseSD^2) * ( 1/Worse + 1/(1-Worse)^2 )
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
  
  # putting all together in a tibble and combine
  es <- tibble(M1, V1, N1, M2, V2, N2, d, Vd)
  data <- bind_cols(data, es)
  return(data)  
}



###TESTS
test <- tibble(Worse = c(50,50), 
               WorseSD = c(50,50), 
               WorseN = c(20, 20), 
               Better = c(60,60), 
               BetterSD = c(50, 50), 
               BetterN = c(20, 20), 
               WithinBetween = c("within", "between"))

calc_ES_proportion(test, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="percent") 
#test %<>% calc_ES_proportion(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="percent")
#test2 <- test %>% prepare_df(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type="percent") 







#' ########################################################################################## SUNDAY MESS
#' ##########################################################################################
#' 
#' ### PREPARE DATA
#' prepare_df <- function(data, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, MeasureType, DataScale, type=c("lnorm", "delta")){
#'   # turning parameters that are column names into strings
#'   Worse <- data[[deparse(substitute(Worse))]]
#'   WorseSD <- data[[ deparse(substitute(WorseSD))]]
#'   WorseN <- data[[deparse(substitute(WorseN))]]
#'   Better <- data[[deparse(substitute(Better))]]
#'   BetterSD <- data[[deparse(substitute(BetterSD))]]
#'   BetterN <- data[[deparse(substitute(BetterN))]]
#'   MeasureType <- data[[deparse(substitute(MeasureType))]]
#'   DataScale <- data[[deparse(substitute(DataScale))]]
#'   
#'   type <- match.arg(type)
#'   
#'   
#'   ##for latency
#'   if(MeasureType == "latency" && type == "lnorm") {  
#'     M1 <- log(Better) - (BetterSD^2)/(2*Better^2)
#'     V1 <- (BetterSD^2)/(Better^2)
#'     M2 <- log(Worse) - (WorseSD^2)/(2*Worse^2)
#'     V2 <- (WorseSD^2)/(Worse^2)    
#'   } else {
#'     if(MeasureType == "latency" && type == "delta") {  
#'     M1 <- log(Better) - log( sqrt(1 + (BetterSD^2)/(Better^2) ) )
#'     V1 <- log( 1 + (BetterSD^2)/(Better^2) )
#'     M2 <- log(Worse) - log( sqrt(1 + (WorseSD^2)/(Worse^2) ) )
#'     V2 <- log( 1 + (WorseSD^2)/(Worse^2) )
#'   } else {
#' 
#'     #for proportion
#'   if(MeasureType == "proportion" && DataScale == "natural"){
#'     Better <- Better/100 #convert form percentage to proportion
#'     BetterSD <- BetterSD/100 #convert form percentage to proportion
#'     Worse <- Worse/100 #convert form percentage to proportion
#'     WorseSD <- WorseSD/100 #convert form percentage to proportion
#'     M1 <- qlogis(Better) + (BetterSD^2)/2 * (1/((1-Better)^2) -1/(Better^2))
#'     V1 <- (BetterSD^2) * ( 1/Better + 1/(1-Better)^2 )
#'     M2 <- qlogis(Worse) + (WorseSD^2)/2 * (1/((1-Worse)^2) -1/(Worse^2))
#'     V2 <- (WorseSD^2) * ( 1/Worse + 1/(1-Worse)^2 )
#'   } else {
#' 
#'     M1 <- Better
#'     V1 <- BetterSD^2
#'     M2 <- Worse
#'     V2 <- WorseSD^2}
#'   }
#' }
#' 
#'   
#'   N1 <- BetterN
#'   N2 <- WorseN
#'   # putting all together in a tibble and combine
#'   MVN <- as_tibble(cbind(data, M1, V1, N1, M2, V2, N2))
#'   return(MVN)
#' } 
#' 
#' 
#' ###TESTS
#' #dat2 <- dat %>% prepare_df(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, MeasureType, DataScale, type="lnorm") 
#' 
#' test <- tibble(Worse = c(1,2), 
#'                WorseSD = c(2,2), 
#'                WorseN = c(20, 20), 
#'                Better = c(2, 4), 
#'                BetterSD = c(2, 2), 
#'                BetterN = c(20, 20), 
#'                MeasureType = c("latency", "proportion"), 
#'                DataScale = c("natural", "logit")
#'   
#' )
#' 
#' prepare_df(test, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, MeasureType, DataScale, type="lnorm") 
#' 
#' 
#' # test <- tibble(
#' #   Wores = 1,
#' #   ....
#' # )
#' 
#' 
#' ##########################################################################################
#' ##########################################################################################
#' ###
#' 
#' 
#'   
#' 
#' ##########################################################################################
#' 
#'   #' @title Calculating Cohens d and Hedges g for proportion data
#'   #' @description Function for calculating Cohens d and Hedges g for proportion-like data using logit-transformed means and SD, for between and within-subject study designs
#'   #' @param data Dataframe object containing effect sizes, their variance, unique IDs and clustering variable
#'   #' @param Worse Name of the variable (vector) containing means of the "Worse" condition group
#'   #' @param WorseSD Name of the variable (vector) containing standard deviations of the "Worse" condition group
#'   #' @param WorseN Name of the variable (vector) containing sample sizes of the "Worse" condition group
#'   #' @param Better Name of the variable (vector) containing means of the "Better" condition group
#'   #' @param BetterSD Name of the variable (vector) containing standard deviations of the "Better" condition group
#'   #' @param BetterN Name of the variable (vector) containing sample sizes of the "Better" condition group
#'   #' @param type Optional parameter indicating whether the data is "natural" (default, percentage, gets converted to proportion) or "logit"
#'   #' @param WithinBetween Optional parameter indicating whether formula for between-subject study design should be used (default or "between") or a within-subject formula should be used ("within").
#'   #' @param adjusted Optional logical parameter indicating whether d should be adjusted for small sample sizes, i.e. calculating Hedges d. Default value is "TRUE"
#'   #' @export
#'   #Value: adds columns with calculated d and Vd to the dataframe
#'   
#'   calc_d_proportion <- function(data, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, DataScale, adjusted=TRUE){
#'     
#'     if (missing(data)) 
#'       stop("Must specify dataframe via 'data' argument.")
#'     
#'     # turning parameters that are column names into strings
#'     Worse <- data[[deparse(substitute(Worse))]]
#'     WorseSD <- data[[ deparse(substitute(WorseSD))]]
#'     WorseN <- data[[deparse(substitute(WorseN))]]
#'     Better <- data[[deparse(substitute(Better))]]
#'     BetterSD <- data[[deparse(substitute(BetterSD))]]
#'     BetterN <- data[[deparse(substitute(BetterN))]]
#'     DataScale <- data[[deparse(substitute(DataScale))]]
#'     WithinBetween <- data[[deparse(substitute(WithinBetween))]]
#'     
#'     # if(DataScale == "natural"){  
#'     #   Better <- Better/100 #convert form percentage to proportion
#'     #   BetterSD <- BetterSD/100 #convert form percentage to proportion
#'     #   Worse <- Worse/100 #convert form percentage to proportion
#'     #   WorseSD <- WorseSD/100 #convert form percentage to proportion
#'     #   M1 <- qlogis(Better) + (BetterSD^2)/2 * (1/((1-Better)^2) -1/(Better^2))
#'     #   V1 <- (BetterSD^2) * ( 1/Better + 1/(1-Better)^2 )
#'     #   M2 <- qlogis(Worse) + (WorseSD^2)/2 * (1/((1-Worse)^2) -1/(Worse^2))
#'     #   V2 <- (WorseSD^2) * ( 1/Worse + 1/(1-Worse)^2 )
#'     # }
#'     # 
#'     # if(DataScale == "logit"){  
#'     #   M1 <- Better
#'     #   V1 <- BetterSD^2
#'     #   M2 <- Worse
#'     #   V2 <- WorseSD^2
#'     # }
#'     
#'     # use vectorization to switch between calculalting standard error for d with within or between-subject study design option
#'     se_choose <- function(data, Worse, WorseSD, Better, BetterSD, DataScale) {
#'       if(DataScale == "logit"){  
#'         M1 <- Better
#'         V1 <- BetterSD^2
#'         M2 <- Worse
#'         V2 <- WorseSD^2
#'       }
#'       else {
#'         Better <- Better/100 #convert form percentage to proportion
#'         BetterSD <- BetterSD/100 #convert form percentage to proportion
#'         Worse <- Worse/100 #convert form percentage to proportion
#'         WorseSD <- WorseSD/100 #convert form percentage to proportion
#'         M1 <- qlogis(Better) + (BetterSD^2)/2 * (1/((1-Better)^2) -1/(Better^2))
#'         V1 <- (BetterSD^2) * ( 1/Better + 1/(1-Better)^2 )
#'         M2 <- qlogis(Worse) + (WorseSD^2)/2 * (1/((1-Worse)^2) -1/(Worse^2))
#'         V2 <- (WorseSD^2) * ( 1/Worse + 1/(1-Worse)^2 )
#'       }
#'       # # putting all together in a tibble and combine
#'       MV <- tibble(M1, V1, M2, V2)
#'       data <- bind_cols(data, MV)
#'       return(data)
#'     }
#'     ok <- pmap_dfr(list(data, Worse, WorseSD, Better, BetterSD, DataScale), se_choose)
#'     return(ok)
#'     
#'     # #calculate sPooled
#'     # sPooled <- sqrt( ((BetterN-1) * V1 + (WorseN-1) * V2) / (BetterN + WorseN - 2) )
#'     # 
#'     # #calculate Cohen's d and Hedges g (d adjusted for small sample sizes)
#'     # d <- (M1-M2)/sPooled
#'     # if(adjusted==TRUE){
#'     #   d <- d * (1-(3/ (4 * (BetterN + WorseN -2) -1))) #calculate Hedges d (adjusted for small N)
#'     # }
#'     # 
#'     # # use vectorization to switch between calculalting standard error for d with within or between-subject study design option
#'     # se_choose <- function(WorseN, BetterN, WithinBetween, d) {
#'     #   if (WithinBetween == "between") {
#'     #     SEd <- sqrt(((WorseN + BetterN) / (WorseN * BetterN) ) + ( (d^2) / (2 * (WorseN + BetterN - 2))))
#'     #   } #for between-subject study design
#'     #   else {
#'     #     SEd <- sqrt(1/WorseN + ((d^2) / (2 * (WorseN  - 1))))
#'     #   } #for within-subject study design
#'     # }
#'     # SEd <- pmap_dbl(list(WorseN, BetterN, WithinBetween, d), se_choose)
#'     # 
#'     # #calculate Variance for d 
#'     # Vd <- SEd^2
#'     # 
#'     # # putting all together in a tibble and combine
#'     # es <- tibble(d, Vd)
#'     # data <- bind_cols(data, es)
#'     # return(data)  
#'   }
#'   
#'   
#'   # TEST USE
#'   str(dat_prop) 
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, DataScale)
#' 
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=TRUE, type=c("lnorm"))
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=FALSE, type=c("lnorm"))
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, type=c("natural"))
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, type=c("logit"))
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, adjusted=FALSE)
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN)
#'   calc_d_proportion(dat_prop, Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween, type)
#'   
#'   #pacman::p_load(tidyverse, magrittr)
#'   testdf <- dat_lat
#'   str(testdf) 
#'   testdf %<>% calc_d_proportion(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, WithinBetween)
#'   str(testdf) 
#'   
#'   
  
# ##########################################################################################
# ##########################################################################################
# ###OLD
# 
# ## Function for calculating effect sizes for meta-analysis and standard error for d. 
# # Equations are taken from Nakagawa and Cuthill 2007 Biol. Reviews.
# # Calc.d calculates Hedges'd (sometimes g), which is 'unbiased for small sample size'. 
# # Cohen's d can be returned using adjusted = F.
# Calc.d <- function(CMean, CSD, CN, EMean, ESD, EN, adjusted=T){
#   varianceoftreatment <- (ESD)^2
#   varianceofcontrols <- (CSD)^2
#   sPooled<-sqrt(((EN - 1)*varianceoftreatment + (CN - 1)*varianceofcontrols) / (EN + CN - 2))	
#   d <- (EMean - CMean) / sPooled 
#   H.d <- d * (1 - (3 / (4 * (EN + CN - 2) - 1)))
#   if(adjusted==T){return(H.d)}
#   if(adjusted==F){return(d)}
# }
# 
# # Calc.SE.d calculates for d values given the necessery n's - IT NEEDS TO BE ^2 to be V !
# Calc.SE.d <- function(CN, EN, d){
#   SE <- sqrt(( (CN + EN) / (as.numeric(CN) * as.numeric(EN)) ) + ( (d^2) / (2 * (CN + EN - 2) ) ) )
#   return(SE)
# }


## function to claculate total heterogeneity I-squared, 
# followng http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
# r - ; res - 
calc.I2 <- function (vi, res) 
{
  W <- diag(1/vi)
  X <- model.matrix(res)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  return(100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P))))
}

#you can also estimate how much of the total variance can be attributed to between- and within-cluster heterogeneity separately:
#100 * res$sigma2 / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P)))


 
# ##########################################################################################
# ##########################################################################################
# 
# ##SN other project example function
# es_cal <- function(data, m1, m2, sd1, sd2, n1, n2) {
#   
#   # required packages
#   #require(dplyr)
#   
#   # turning them into strings
#   m1 <- data[[deparse(substitute(m1))]]
#   m2 <- data[[ deparse(substitute(m1))]]
#   sd1 <- data[[deparse(substitute(sd1))]]
#   sd2 <- data[[deparse(substitute(sd2))]]
#   n1 <- data[[deparse(substitute(n1))]]
#   n2 <- data[[deparse(substitute(n2))]]
#   
#   # SMD
#   # the pooled standard devision
#   sdp <- sqrt( ( (n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2 ) /
#                  (n1 + n2 - 2) )
#   # small N correction
#   J <- 1 - 3 / (4 * (n1 + n2) - 9)
#   smd <- ((m2 - m1) / sdp) * J
#   
#   #Var(SMD)
#   vsmd <- (n1 + n2) / (n1*n2) + smd ^ 2 /(2 * (n1 + n2))
#   
#   # SSDD
#   ssdd <-  (sd2 - sd1) / sdp
#   # Var(SSDD)
#   vssdd <- (1 / sdp ^ 2) * ( (sd1 ^ 2 / ( 2 * (n1 - 1) ) ) +
#                                (sd2 ^ 2 / (2 * (n2 - 1) ) ) +
#                                (sd2 - sd1)^2/(2*(n1 + n2 - 2) ) )
#   
#   # lnRR (modify later using the correction formula)
#   lnrr <- log(m2 / m1)
#   
#   # Var(lnRR)
#   vlnrr <- sd1^2/(n1 * m1^2) + sd2^2/(n2 * m2^2)
#   
#   # lnVR
#   lnvr <- log(sd2 / sd1) + 1 / (2 * (n1 - 1) ) - 1 / (2 * (n1 - 1) )
#   
#   # Var(lnVR)
#   vlnvr <- 1 / (2 * (n2 - 1) ) + 1 / (2 * (n1 - 1) )
#   
#   # lnVR
#   lncvr <- log(sd2 / m2) - log(sd1 / m1) + 1 / (2 * (n1 - 1) ) - 1 / (2 * (n1 - 1) )
#   
#   # Var(lnVR)
#   vlncvr <- 1 / ( 2 * (n2 - 1) ) + 1 / ( 2 * (n1 - 1) ) + sd1^2/(n1 * m1^2) + sd2^2/(n2 * m2^2)
#   
#   # putting all together in a tibble and combine
#   es <- tibble(smd, ssdd, lnrr, lnvr, lncvr, vsmd, vssdd, vlnrr, vlnvr, vlncvr)
#   data <- bind_cols(data, es)
#   return(data)
# }
# 
# # USE
# pacman::p_load(tidyverse, gridExtra, purrr, magrittr, metafor, cowplot, moments)
# dat %<>% es_cal(mean1_r, mean2_r, sd1_r, sd2_r, n1_r, n2_r)
# 


# #### ORIGINAL custom functions
# #### Functions for calculating effect size for proportion data modified by ML&SN from the function (originally by A M Senior): we use qlogis(meanproportion) and SD=sqrt(pi^2/3)
# #### Calc.d calculates Hedges'd (sometimes g), which is 'unbiased for small sample size'. Cohen's d can be returned using adjusted = F.
# #### Calc.SE.d calculates d values given the necessery n's
# 
# ## Calc.d.percA calculates the effect size d for proportion data with fixed SD:
# # Calc.d.propA <- function(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, adjusted=T){
# #   #logit SD, divide by 100 to get proportions from procentage data
# #   logitBetter <- qlogis(Better/100)
# #   logitWorse <- qlogis(Worse/100)
# #   #variance from SD
# #   VarianceBetter <- pi^2/3
# #   VarianceWorse <- pi^2/3
# #   #sPooled 
# #   sPooled <- sqrt(((BetterN-1)*VarianceBetter + (WorseN-1)*VarianceWorse)/(BetterN + WorseN - 2))
# #   #d and Hedges d
# #   d <- (logitBetter-logitWorse)/sPooled
# #   H.d <- d * (1-(3/ (4 * (BetterN + WorseN -2) -1)))
# #   if(adjusted==T){return(H.d)}
# #   if(adjusted==F){return(d)}
# # }
# 
# 
# ## Calc.d.propB calculates the effect size d for proportion data with non-fixed SD:
# Calc.d.propB <- function(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, adjusted=T, DataType="natural"){
#   #if natural data then logit means and SDs, divide by 100 to get proportions from procentage data
#   if (DataType == "natural")
#   {
#     logitBetter <- qlogis(Better/100)
#     logitWorse <- qlogis(Worse/100)
#     SDlogitBetter <- (BetterSD/100) * ( 1/(Better/100) + 1/(1-(Better/100)) )
#     SDlogitWorse <- (WorseSD/100) * ( 1/(Worse/100) + 1/(1-(Worse/100)) )
#   }
#   #if logit data then no need to transform
#   if (DataType == "logit")
#   {
#     logitBetter <- Better
#     logitWorse <- Worse
#     SDlogitBetter <- BetterSD
#     SDlogitWorse <- WorseSD
#   }
#   #variance
#   VarianceBetter <- SDlogitBetter^2
#   VarianceWorse <- SDlogitWorse^2
#   #sPooled 
#   sPooled <- sqrt(((BetterN-1)*VarianceBetter + (WorseN-1)*VarianceWorse)/(BetterN + WorseN - 2))
#   #d and Hedges d
#   d <- (logitBetter-logitWorse)/sPooled
#   H.d <- d * (1- (3/ (4 * (BetterN + WorseN -2) -1)))
#   if(adjusted==T){return(H.d)}  
#   if(adjusted==F){return(d)}
# }
# 
# 
# ## Calc.d.lat calculates the effect size d for latency data:
# Calc.d.lat_raw <- function(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, adjusted=T){
#   VarianceBetter <- BetterSD^2
#   VarianceWorse <- WorseSD^2
#   #sPooled is calculated with the logged variables
#   sPooled <- sqrt(((BetterN-1)*VarianceBetter + (WorseN-1)*VarianceWorse)/(BetterN + WorseN -2))
#   #d and Hedges d
#   d <- (Better-Worse)/sPooled
#   H.d <- d * (1-(3/ (4 * (BetterN + WorseN -2) -1)))
#   if(adjusted==T){return(-H.d)} #Reversing sign on latency data to match proportion data. #Now positive values means optimistic for both lat and prop data.
#   if(adjusted==F){return(-d)} #Reversing sign on latency data to match proportion data. #Now positive values means optimistic for both lat and prop data.
# }
# 
# ## Calc.d.lat calculates the effect size d for latency data (after log-transformation):
# Calc.d.lat_log <- function(Worse, WorseSD, WorseN, Better, BetterSD, BetterN, adjusted=T){
#   #log mean and SDn
#   logBetter <- log(Better)
#   logWorse <- log(Worse)
#   #variance
#   VarianceBetter <- (BetterSD/Better)^2
#   VarianceWorse <- (WorseSD/Worse)^2
#   #sPooled is calculated with the logged variables
#   sPooled <- sqrt(((BetterN-1)*VarianceBetter + (WorseN-1)*VarianceWorse)/(BetterN + WorseN -2))
#   #d and Hedges d
#   d <- (logBetter-logWorse)/sPooled
#   H.d <- d * (1-(3/ (4 * (BetterN + WorseN -2) -1)))
#   if(adjusted==T){return(-H.d)} #Reversing sign on latency data to match proportion data. #Now positive values means optimistic for both lat and prop data.
#   if(adjusted==F){return(-d)} #Reversing sign on latency data to match proportion data. #Now positive values means optimistic for both lat and prop data.
# }
# 
# # #Standard error for d
# # Calc.SE.d <- function(WorseN, BetterN, d){
# #   SE <- sqrt(((WorseN + BetterN) / (WorseN * BetterN) ) + ( (d^2) / (2 * (WorseN + BetterN - 2)))) #between study design
# #   SE <- sqrt(( (d^2) / (2 * (WorseN  - 1)))) #within study design - assuming the perfect correlation (conservative estimate)
# #   return(SE)
# # }
# 
# 
# #Standard error for d with within - between subject study design option
# Calc.SE.d <- function(WorseN, BetterN, d, WithinBetween="between"){
#   if (WithinBetween == "between") {SE <- sqrt(((WorseN + BetterN) / (WorseN * BetterN) ) + ( (d^2) / (2 * (WorseN + BetterN - 2))))} #for between-subject study design
#   if (WithinBetween == "within") {SE <- sqrt(1/WorseN + ((d^2) / (2 * (WorseN  - 1))))} #for within-subject study design - assuming 0.5 correlation (conservative estimate)
#   return(SE)
# }
# 
