
#EXACT


#' @title Enumerate the Permutations of the Elements of a Vector When Some of those Elements are Identical
#' @description This function enumerates the possible combinations of \code{n} elements where the first element is repeated \code{n1} times, the sencond element is repeated \code{n2} times, the third \code{n3} times, ...
#' @param Sample_Sizes Numeric vector (\code{n1},...,\code{nk}) that indicates the number of times each element is repeated.
#' @return Returns a matrix where each row contains a permutation.
#' @section Warning:
#' The number of permutations increases rapidly with lenght (Sets). The computational time increases exponetianly with the number of elements and with the number of sets.
#' @examples
#' Sample_Sizes <- c(2,2,2)
#' Permutations_With_Repetition(Sample_Sizes)
#'
#' Sample_Sizes <- c(2,3)
#' Permutations_With_Repetition(Sample_Sizes)
#'
#' Permutations_With_Repetition(c(3,2,3))
#'
#' Permutations_With_Repetition(c(10,5))
#' @export
Permutations_With_Repetition <- function(Sample_Sizes){

  if(!is.numeric(Sample_Sizes))
    stop("Some elements of 'Sample_Sizes' are not numeric")

  N_Samples <- length(Sample_Sizes)
  node <- NULL
  node[[1]] <- list(N_Samples=Sample_Sizes, Sample_Sizes=Sample_Sizes, padre=0, n_hijos=0, hijo=rep(-1,N_Samples), seq=rep(0,sum(Sample_Sizes)))
  continue <- TRUE
  cur_node <- 1
  last_node <- 1
  N_permutations <- factorial(sum(Sample_Sizes))

  for(i in 1:N_Samples) N_permutations <- N_permutations/factorial(Sample_Sizes[i])

  x <- 'y'
  if (N_permutations > 500000) {
    x <- readline(prompt="Computational time will be very high with these sample sizes. Type 'y' to continue or anything else to exit: ")
  }

  permutations <- NULL

  if ((x == 'y') || (x == 'Y')){

    permutations <- matrix(rep(0,sum(Sample_Sizes)*N_permutations), nrow=N_permutations, ncol=sum(Sample_Sizes))

    j <- 1
    while(continue){

      for(i in 1:N_Samples){

        if((node[[cur_node]]$N_Samples[i]!=0)&&(node[[cur_node]]$N_Samples[i]!=sum(node[[cur_node]]$N_Samples))){
          node[[last_node+1]] <- list(N_Samples=Sample_Sizes, Sample_Sizes=Sample_Sizes, padre=cur_node, n_hijos=0, hijo=rep(-1,N_Samples))
          node[[last_node+1]]$N_Samples <- node[[cur_node]]$N_Samples
          node[[last_node+1]]$N_Samples[i] <- node[[cur_node]]$N_Samples[i]-1
          node[[cur_node]]$n_hijos <- node[[cur_node]]$n_hijos + 1
          node[[last_node+1]]$padre <- cur_node
          last_node <- last_node+1
          node[[cur_node]]$hijo[i] <- last_node
          node[[last_node]]$seq <- node[[cur_node]]$seq
          node[[last_node]]$seq[sum(Sample_Sizes)-sum(node[[cur_node]]$N_Samples)+1] <- i
        }else{
          if(node[[cur_node]]$N_Samples[i]>0){
            for(k in 1: node[[cur_node]]$N_Samples[i])
              node[[cur_node]]$seq[sum(Sample_Sizes)-sum(node[[cur_node]]$N_Samples)+k] <- i
          }
        }
      }

      if(node[[cur_node]]$n_hijos == 0 ){
        permutations[j,] <- node[[cur_node]]$seq
        j <- j+1
      }

      if(cur_node==last_node) continue <- FALSE

      cur_node<- cur_node + 1
    }
  }
  return(permutations)
}


#' @import Rglpk
#' @title Linear Ordering Problem (LOP)
#' @description This function computes the solution of the Linear Ordering Problem.
#' @param mat_LOP Preference matrix defining the Linear Ordering Problem. A numeric square matrix for which we want to obtain the permutation of rows/columns that maximizes the sum of the elements above the main diagonal.
#' @return The function returns a list with the following elements:
#' \enumerate{
#'  \item{ \code{obj_val}: Optimal value of the solution of the Linear Ordering Problem, i.e., the sum of the elements above the main diagonal under the permutation rows/cols solution.}
#'  \item{ \code{permutation}: Solution of the Linear Ordering Problem, i.e., the rows/cols permutation.}
#'  \item{ \code{permutation_matrix}: Optimal permutation matrix of the Linear Ordering Problem. }
#' }
#' @section References:
#' Martí, R. and Reinelt, G. The Linear Ordering Problem: Exact and Heuristic Methods in Combinatorial Optimization. Springer, first edition 2011.
#'@examples
#'## Square matrix
#'##
#'##  |  1    2    2  |
#'##  |  2    3    3  |
#'##  |  3    2    2  |
#'##
#'## The optimal permutation of rows/cols is (2,3,1),
#'## and the solution of the Linear Ordering Problem is 8.
#'## Te permutation matrix of the solution is
#'##  |  0    0    0  |
#'##  |  1    0    1  |
#'##  |  1    0    0  |
#'
#' mat_LOP <- matrix(c(1,2,3,2,3,2,2,3,2), nrow=3)
#' LOP(mat_LOP)
#' @export
LOP <- function(mat_LOP){

  if(!is.matrix(mat_LOP) || !is.numeric(mat_LOP))
    stop("'mat_LOP' must be a numeric square matrix")

  if (nrow(mat_LOP) != ncol(mat_LOP))
    stop("'mat_LOP' must be a square matrix")

  n_var <- length(mat_LOP[1,])*length(mat_LOP[1,])-length(mat_LOP[1,])

  if(length(mat_LOP[1,])>2)
    n_con <- n_var +  factorial(length(mat_LOP[1,]))/factorial(length(mat_LOP[1,])-3)

  if(length(mat_LOP[1,])!= length(mat_LOP[,1])) stop("matLOP is not a square matrix ")

  if(length(mat_LOP[1,]) == 2)
    n_con <- n_var

  obj <- c(rep(0,n_var))
  count <- 0
  mat <- matrix(rep(0,n_var*n_con),nrow=(n_con))
  dir <- c(rep("",n_con))
  rhs <- c(rep(0,n_con))
  types <- c(rep("B",n_var))
  index <- matrix(rep(0,length(mat_LOP[1,])*length(mat_LOP[1,])), nrow=(length(mat_LOP[1,])))

  for(i in 1:length(mat_LOP[1,])) {
    for(j in 1:length(mat_LOP[1,])) {
      if((i!=j)){
        count <- count + 1
        obj[count] <- mat_LOP[i,j]
        index[i,j] <- count
      }
    }
  }

  count <- 0
  for(i in 1:length(mat_LOP[1,])) {
    for(j in 1:length(mat_LOP[1,])){
      if((i!=j)){
        count<-count + 1
        mat[count,index[i,j]] <- 1
        mat[count,index[j,i]] <- 1
        dir[count] <- "=="
        rhs[count] <- 1
      }
    }
  }

  if(length(mat_LOP[1,])>2) {
    for(i in 1:length(mat_LOP[1,])){
      for(j in 1:length(mat_LOP[1,])){
        for(k in 1:length(mat_LOP[1,])){
          if((k!=i)&&(k!=j)&&(j!=i)){
            count<-count + 1
            mat[count,index[i,j]] <- 1
            mat[count,index[j,k]] <- 1
            mat[count,index[k,i]] <- 1
            dir[count] <-"<="
            rhs[count] <- 2
          }
        }
      }
    }
  }

  max <- TRUE
  sol <- Rglpk_solve_LP(obj, mat, dir, rhs, types=types, max=TRUE)
  obj_val <- sol$optimum
  solution <- sol$solution
  permutation_matrix <- matrix(rep(0,length(mat_LOP[1,])*length(mat_LOP[1,])), nrow=length(mat_LOP[1,]))

  permutation <- c(1:length(mat_LOP[1,]))
  sum_row <- c(rep(0,length(mat_LOP[1,])))

  for(i in 1:length(mat_LOP[1,])){
    for(j in 1:length(mat_LOP[1,])){
      if((i!=j)){
        if(solution[index[i,j]]==1) permutation_matrix[i,j] <- 1
      }
    }
  }

  for(i in 1:length(mat_LOP[1,]))
    sum_row[i] <-   sum(permutation_matrix[i,])

  z <- cbind(sum_row,permutation)
  permutation <- z[order(-sum_row,permutation),][,2]

  return(list(obj_val=obj_val, permutation=permutation, permutation_matrix=permutation_matrix))

}


#' @title Probability Distribution of Concordance Coefficient and Kruskal Wallis Statistic
#' @description This function computes the probability distribution tables of Concordance Coefficient and Kruskal Wallis Statistic.
#' @param Sample_Sizes Numeric vector (\code{n1},...,\code{nk}) containing the number of repetitions of each element, i.e., the size of each sample in the experiment.
#' @param H 0 by default. If set to 1, the probability distribution table of Kruskal-Wallis Statistic is also calculated and returned.
#' @return The function returns a list with the following elements:
#' \enumerate{
#'  \item{ \code{C_freq}: Matrix with the probability distribution of Concordance Coefficient. Each row in the matrix contains the disorder, the value of the statistic, the frequency and its probability.}
#'  \item{ \code{H_freq}: Matrix with the probability distribution of Kruskal Wallis Statistic. Each row in the matrix contains the value of the statistic, the frequency and its probability (only if H = 1).}
#' }
#' @section Warning:
#' The number of permutations increases rapidly with lenght (Sets). The computational time increases exponetianly with the number of elements and with the number of sets.
#'@examples
#' Sample_Sizes <- c(5,4)
#' CT_Distribution(Sample_Sizes)
#' CT_Distribution(Sample_Sizes, H = 1)
#'
#' CT_Distribution(c(3,3,3), H = 1)
#' @export
CT_Distribution<-function(Sample_Sizes, H = 0){

  if(!is.numeric(Sample_Sizes))
    stop("Some elements of 'Sample_Sizes' are not numeric")

  N_g <- length(Sample_Sizes)
  if (N_g < 2)
    stop("'Sample_Sizes' must be a numeric vector with at least 2 elements")

  if ((H != 0) && (H != 1))
    warning("The probability distribution table of Kruskal Wallis Statistic is not calculated. Set H = 1 for it")

  N <- sum(Sample_Sizes)
  b <- 0
  GP <- 0

  N_T <-0
  for(i in 1:(N_g-1)){
    for (j in (i+1):N_g){
      N_T<- N_T + Sample_Sizes[i]*Sample_Sizes[j]
    }
  }

  for(i in 1:N_g) if((Sample_Sizes[i])%%2!=0) b <- b+1

  N_min <- 0
  for(i in 1:(N_g-1)){
    for (j in (i+1):N_g){
      N_min <- N_min + trunc(Sample_Sizes[i]*Sample_Sizes[j]/2)
    }
  }

  GB <- 0
  if(b>1){
    if(b%%2==0){
      l <- b/2
      GB <- l*(3*l-1)/2
    }else{
      l<-(b-1)/2
      GB<-l*(3*l+1)/2
    }
  }
  GB
  if(b==1){GB <- 0}
  N_min <- N_min + GB

  p_p <- Permutations_With_Repetition(Sample_Sizes)

  x <- 'y'
  if (length(p_p[,1]) > 15000) {
    x <- readline(prompt=" Computational time will be very high with these sample sizes. Type 'y' to continue or anything else to exit: ")
  }

  C_freq <- NULL
  H_freq <- NULL

  if ((x == 'y') || (x == 'Y')){

    if(!is.null(p_p)) {

      pref <- list()
      for(k in 1:length(p_p[,1])){
        pref[[k]] <- matrix(rep(0,N_g*N_g), nrow=N_g, ncol=N_g)
        for(i in 1:(N-1)){
          for(j in (i+1):N){
            pref[[k]][p_p[k,i],p_p[k,j]] <-  pref[[k]][p_p[k,i],p_p[k,j]] + 1
          }
        }
      }

      opt<-c(rep(0,length(p_p[,1])))
      rho<-c(rep(0,length(p_p[,1])))

      for(l in 1:length(p_p[,1])){
        opt[l] <- LOP(pref[[l]])$obj_val
        opt[l] <- N_T-opt[l]
        rho[l] <- 1 - opt[l]/(N_T-N_min)
      }

      if (any(sapply(rho, is.na))) {
        C_freq <-  cbind(NA,NA,NA,NA)
        colnames(C_freq)<-c("disorder", "Concordance Coefficient", "Frequency", "Probability")

        H_freq <- cbind(NA,NA,NA)
        colnames(H_freq)<-c("H Statistic", "Frequency", "Probability")
      }
      else {

        if (H == 1) {

          KW <- c(rep(0,length(p_p[,1])))

          for(l in 1:length(p_p[,1])){
            KW[l] <- sum(tapply(1:N, p_p[l,], "sum")^2 / tapply(1:N, p_p[l,], "length"))
            KW[l] <- 12*KW[l]/(N*(N+1)) -3*(N+1)
          }
        }

        rho <- round(sort(rho),2)
        count <- 1
        count2 <- 1
        rho2 <- 0
        rho2[1] <- rho[1]
        C_freq <- 0
        C_freq[1] <- count2
        for(i in 2:length(rho)){
          if (rho[i]==rho[i-1]){ count <- count + 1 }
          if (rho[i]>rho[i-1]) {  C_freq[count2] <- count; count2 <- count2+1; rho2[count2] <- rho[i]; count <- 1 }
        }
        C_freq[count2] <- count

        if (H == 1) {
          KW <- round(sort(KW),2)
          count <- 1
          count2 <- 1
          KW2 <- 0
          KW2[1] <- rho[1]
          H_freq <- 0
          H_freq[1] <- count2

          for(i in 2:length(KW)){
            if(KW[i]==KW[i-1]){ count <- count +1}
            if(KW[i]>KW[i-1]){  H_freq[count2] <- count ;count2 <- count2+1; KW2[count2] <- KW[i]; count <-  1 }
          }
          H_freq[count2] <- count

          H_freq <- cbind(KW2, H_freq, round((H_freq/sum(H_freq)),4))
          colnames(H_freq)<-c("H Statistic", "Frequency", "Probability")
        }

        C_freq <- cbind(round(((N_T-N_min)*(1-rho2)),0), rho2, C_freq, round((C_freq/sum(C_freq)),4))
        colnames(C_freq)<-c("disorder", "Concordance Coefficient", "Frequency", "Probability")
      }
    }
  }

  if (H == 1) return(list(C_freq=C_freq,H_freq=H_freq))
  else return(list(C_freq=C_freq))
}


#' @title Critical Values of Concordance and Kruskal Wallis Tests. Exact Method
#' @description This function computes the critical values of the Concordance and Kruskal Wallis Tests. Exact p-value of desired significance levels of .10, .05 and .01.
#' @param Sample_Sizes Numeric vector (\code{n1},...,\code{nk}) containing the number of repetitions of each element, i.e., the size of each sample in the experiment.
#' @param H 0 by default. If set to 1, the critical values of the Kruskal-Wallis Test are also calculated and returned.
#' @return The function returns a list with the following elements:
#' \enumerate{
#'  \item{ \code{C_results}: Concordance Coefficient results. Critical values and exact p-values for a desired significance levels of 0.1, .05 and .01. }
#'  \item{ \code{H_results}: Kruskal Wallis results. Critical values and exact p-values for a desired significance levels of 0.1, .05 and .01 (only if H = 1). }
#' }
#' @section Warning:
#' The number of permutations increases rapidly with lenght (Sets). The computational time increases exponetianly with the number of elements and with the number of sets.
#'@examples
#' Sample_Sizes <- c(4,4,2)
#' CT_Critical_Values(Sample_Sizes, H = 1)
#' @export
CT_Critical_Values<-function(Sample_Sizes, H = 0){

  if(!is.numeric(Sample_Sizes))
    stop("Some elements of 'Sample_Sizes' are not numeric")

  if (length(Sample_Sizes) < 2)
    stop("'Sample_Sizes' must be a numeric vector with at least 2 elements")

  if ((H != 0) && (H != 1)){
    warning("The critical values of the Kruskal Wallis Test are not calculated. Set H = 1 for it")
    H <- 0
  }

  distributions <- CT_Distribution(Sample_Sizes, H)

  C_freq <- distributions$C_freq

  C_results <- NULL
  H_results <- NULL

  if(!is.null(C_freq)) {

    if (is.na(C_freq[1,1])) {
      C_results <- matrix(c(rep(NA,9)), nrow=3, byrow=TRUE)
      colnames(C_results)<-c("|  disorder", "|  Concordance Coefficient ","|  p-value")
      rownames(C_results)<-c("Sig level .10", "Sig level .05","Sig level .01")
      H_results <- matrix(c(rep(NA,6)), nrow=3, byrow=TRUE)
      colnames(H_results)<-c("|  H Statistic", "| p-value")
      rownames(H_results)<-c("Sig level .10", "Sig level .05", "Sig level .01")
    }
    else {

      C_result <- c(rep(0,9))
      for(k in 1:3){
        if(k==1)  alpha <- 0.1
        if(k==2)  alpha <- 0.05
        if(k==3)  alpha <- 0.01

        sum <- 0
        for(i in 1:(length(C_freq[,4])-1)){
          sum <- sum+C_freq[i,4]
          if( sum > (1-alpha) ){
            C_result[1+(k-1)*3] <- C_freq[i+1,1]
            C_result[2+(k-1)*3] <- C_freq[i+1,2]
            C_result[3+(k-1)*3] <- 1-sum
            alpha<- -1000
          }
        }
        i <- length(C_freq[,3])

        if(C_result[1+(k-1)*3]==0){
          if(k==1)  alpha <- 0.1
          if(k==2)  alpha <- 0.05
          if(k==3)  alpha <- 0.01
          if(C_freq[i,4]< alpha){
            C_result[1+(k-1)*3] <- C_freq[i,1]
            C_result[2+(k-1)*3] <- C_freq[i,2]
            C_result[3+(k-1)*3] <- C_freq[i,4]
          }
          else {
            C_result[1+(k-1)*3] <- NA
            C_result[2+(k-1)*3] <- NA
            C_result[3+(k-1)*3] <- NA

          }
        }
      }

      if (H == 1) {

        H_freq  <- distributions$H_freq

        H_result <- c(rep(0,6))
        for(k in 1:3){
          if(k==1)  alpha <- 0.1
          if(k==2)  alpha <- 0.05
          if(k==3)  alpha <- 0.01

          sum <- 0
          for(i in 1:(length(H_freq[,3])-1)){
            sum <- sum+H_freq[i,3]
            if( sum > (1-alpha) ){
              H_result[1+(k-1)*2] <- H_freq[i+1,1]
              H_result[2+(k-1)*2] <- 1-sum
              alpha<- -1000
            }
          }

          i <- length(H_freq[,3])

          if(H_result[1+(k-1)*2] == 0){
            if(k==1)  alpha <- 0.1
            if(k==2)  alpha <- 0.05
            if(k==3)  alpha <- 0.01
            if(H_freq[i,3]< alpha){
              H_result[1+(k-1)*2] <- H_freq[i,1]
              H_result[2+(k-1)*2] <- H_freq[i,3]
            }
            else {
              H_result[1+(k-1)*2] <- NA
              H_result[2+(k-1)*2] <- NA
            }
          }
        }
        H_results <- matrix(H_result, nrow=3, byrow=TRUE)
        colnames(H_results)<-c("|  H Statistic", "| p-value")
        rownames(H_results)<-c("Sig level .10", "Sig level .05", "Sig level .01")
      }

      C_results <- matrix(C_result, nrow=3, byrow=TRUE)
      colnames(C_results)<-c("|  disorder", "|  Concordance Coefficient ","|  p-value")
      rownames(C_results)<-c("Sig level .10", "Sig level .05","Sig level .01")
    }
  }

  if (H == 1) return(list(C_results=C_results, H_results=H_results))
  else return(list(C_results=C_results))

}

#' @import graphics
#' @title Probability Plot for Concordance Coefficient and Kruskal Wallis Statistic
#' @description This function performs the graphical visualization of the probability distribution of Concordance Coefficient and Kruskal Wallis Statistic.
#' @param C_freq Frecuency and probability distribution of Concordance Coefficient obtained from the functions: \code{\link{CT_Distribution}}, \code{\link{CT_Distribution_Sim}}.
#' @param H_freq Frecuency and probability distribution of Kruskal Wallis statistic obtained from the functions: \code{\link{CT_Distribution}}, \code{\link{CT_Distribution_Sim}}.
#' @examples
#' Sample_Sizes <- c(5,4)
#' Distributions <-  CT_Distribution(Sample_Sizes, H = 1)
#' C_freq <- Distributions$C_freq
#' H_freq <- Distributions$H_freq
#' CT_Probability_Plot(C_freq)
#' CT_Probability_Plot(C_freq, H_freq)
#'
#' Sample_Sizes <- c(5,5,5)
#' Distributions <-  CT_Distribution_Sim(Sample_Sizes, Num_Sim = 1000, H = 1)
#' C_freq <- Distributions$C_freq
#' H_freq <- Distributions$H_freq
#' CT_Probability_Plot(C_freq, H_freq)
#' @export
CT_Probability_Plot<-function(C_freq = NULL, H_freq = NULL){

  if ((!is.null(C_freq)) && (!is.matrix(C_freq) || (ncol(C_freq) != 4)))
    stop("'C_freq' must be a matrix obtained from the function CT_Distribution or CT_Distribution_Sim")

  if ((!is.null(H_freq)) && (!is.matrix(H_freq) || (ncol(H_freq) != 3)))
    stop("'H_freq' must be a matrix obtained from the function CT_Distribution or CT_Distribution_Sim")

  if(length(C_freq)>1){
    x <- C_freq[,2]
    y <- C_freq[,4]
    plot(x,y,type="h", xlab=list("1-relative disorder",cex=1.5), ylim = c(0, max(y)), xlim = c(0,1 ), ylab=list("Prob",cex=1.5), frame="false", axes="false", main="Concordance Coefficient")
    axis(1, at = c(0,1), labels = c(0,1), tick = TRUE,lwd=1)
    axis(2, at = c(0, max(y)), labels = c(0,max(y)), tick = TRUE, lwd = 2)
  }

  if(length(H_freq)>1){
    x <- H_freq[,1]
    y <- H_freq[,3]
    plot(x,y,type="h",xlab=list("H statistic",cex=1.5), ylim = c(0, max(y)), xlim = c(0,max(x) ), ylab=list("Prob",cex=1.5), frame="false", axes="false", main="Kruskal Wallis")
    axis(1, at = c(0,max(x)), labels = c(0,max(x)), tick = TRUE,lwd=1)
    axis(2, at = c(0,max(y)), labels =c(0,max(y)), tick = TRUE,lwd=2)
  }

}

#' @import stats
#' @import graphics
#' @title Density Plot from the Concordance Coefficient and Kruskal Wallis Normalized Statistics
#' @description This function performs the graphical visualization of the density distribution of Concordance Coefficient and Kruskal Wallis Statistic.
#' @param C_freq Frecuency and probability distribution of Concordance Coefficient obtained from the functions: \code{\link{CT_Distribution}}, \code{\link{CT_Distribution_Sim}}.
#' @param H_freq Frecuency and probability distribution of Kruskal Wallis statistic obtained from the functions: \code{\link{CT_Distribution}}, \code{\link{CT_Distribution_Sim}}.
#' @examples
#' Sample_Sizes <- c(5,4)
#' Distributions <-  CT_Distribution(Sample_Sizes, H = 1)
#' C_freq <- Distributions$C_freq
#' H_freq <- Distributions$H_freq
#' CT_Density_Plot(C_freq, H_freq)
#'
#'
#' Sample_Sizes <- c(5,5,5)
#' Distributions <-  CT_Distribution_Sim(Sample_Sizes, Num_Sim = 1000, H = 1)
#' C_freq <- Distributions$C_freq
#' H_freq <- Distributions$H_freq
#' CT_Density_Plot(C_freq, H_freq)
#' @export
CT_Density_Plot<-function(C_freq = NULL, H_freq = NULL){

  if ((!is.null(C_freq)) && (!is.matrix(C_freq) || (ncol(C_freq) != 4)))
    stop("'C_freq' must be a matrix obtained from the function CT_Distribution or CT_Distribution_Sim")

  if ((!is.null(H_freq)) && (!is.matrix(H_freq) || (ncol(H_freq) != 3)))
    stop("'H_freq' must be a matrix obtained from the function CT_Distribution or CT_Distribution_Sim")

  KT <- NULL
  KW <- NULL

  if(length(C_freq)>1){
    frecuency_Total <- 10*length(C_freq[,3])/sum(C_freq[,3])
    for(i in 1:length(C_freq[,3])) {
      KT <- c(KT,rep(C_freq[i,2],frecuency_Total*C_freq[i,3]))
    }
  }

  if(length(H_freq)>1){
    frecuency_Total <- 10*length(H_freq[,2])/sum(H_freq[,2])
    for(i in 1:length(H_freq[,2])) {
      KW <- c(KW,rep(H_freq[i,1],frecuency_Total*H_freq[i,2]))
    }
  }

  if (is.numeric(KT) && is.numeric(KW)) {

    y_lim <- c(0, max(max(density(KT)$y),max(density(KW/max(KW))$y)))
    xy_legend <- y_lim + c(0.8,.0)

    plot(density(KT),type="l",pch=1.2, cex=20., lwd=1.95,lty=1, frame="false",axes="false",main="",ylab=list(""),xlab=list(""),
         ylim=y_lim,col="black",xlim=c(-0.25,1.55))
    lines(density(KW/max(KW)),pch=1.2, cex=1., xlab="", ylab="", col="black", lwd=1.95,lty=2)
    axis(1, at = c(0,.2,.4,.6,.8,1), labels = c(0,.2,.4,.6,.8,1), tick = TRUE,lwd=1)
    legend(xy_legend[1],xy_legend[2] , legend=c("Concordance", "Kruskal Wallis"),
           col=c("black", "black"), cex=.81, lty=c(1,3),lwd=2)

  }
  else {
    if (is.numeric(KT)) {

      y_lim <- c(0, max(density(KT)$y))
      xy_legend <- y_lim + c(0.8,.0)

      plot(density(KT),type="l",pch=1.2, cex=20., lwd=1.95,lty=1, frame="false",axes="false",main="",ylab="",xlab="",
           ylim=y_lim,col="black",xlim=c(-0.25,1.55))
      axis(1, at = c(0,.2,.4,.6,.8,1), labels = c(0,.2,.4,.6,.8,1), tick = TRUE,lwd=1)
      legend(xy_legend[1],xy_legend[2] , legend=c("Concordance"),
             col="black", cex=.81, lty=1,lwd=2)

    }
    else {
      if (is.numeric(KW)) {

        y_lim <- c(0, max(density(KW/max(KW))$y))
        xy_legend <- y_lim + c(0.8,.0)

        plot(density(KW/max(KW)),type="l",pch=1.2, cex=20., lwd=1.95,lty=2,frame="false",axes="false",main="",ylab="",xlab="",
             ylim=y_lim,col="black",xlim=c(-0.25,1.55))
        axis(1, at = c(0,.2,.4,.6,.8,1), labels = c(0,.2,.4,.6,.8,1), tick = TRUE,lwd=1)
        legend(xy_legend[1],xy_legend[2] , legend=c("Kruskal Wallis"),
               col="black", cex=.81, lty=3,lwd=2)
      }
    }
  }
}


#SIMULATED

#' @import stats
#' @title Simulated Probability Distribution of Concordance Coefficient and Kruskal Wallis Statistic
#' @description This function computes by simulation the probability distribution tables of Concordance Coefficient and Kruskal Wallis Statistic.
#' @param Sample_Sizes Numeric vector (\code{n1},...,\code{nk}) containing the number of repetitions of each element, i.e., the size of each sample in the experiment.
#' @param Num_Sim Number of simulations in order to obtain the probability distribution of the statistics. The default is 10000.
#' @param H 0 by default. If set to 1, the probability distribution table of Kruskal-Wallis Statistic is also calculated and returned.
#' @return The function returns a list with the following elements:
#' \enumerate{
#'  \item{ \code{C_freq}: Matrix with the probability distribution of Concordance Coefficient. Each row in the matrix contains the disorder, the value of the statistic, the frequency and its probability.}
#'  \item{ \code{H_freq}: Matrix with the probability distribution of Kruskal Wallis Statistic. Each row in the matrix contains the value of the statistic, the frequency and its probability (only if H = 1).}
#' }
#'@examples
#' Sample_Sizes <- c(5,4)
#' CT_Distribution_Sim(Sample_Sizes, Num_Sim = 1000)
#' CT_Distribution_Sim(Sample_Sizes, Num_Sim = 1000, H = 1)
#' @export
CT_Distribution_Sim<-function(Sample_Sizes, Num_Sim = 10000, H = 0){

  if(!is.numeric(Sample_Sizes))
    stop("Some elements of 'Sample_Sizes' are not numeric")

  N_g <- length(Sample_Sizes)
  if (N_g < 2)
    stop("'Sample_Sizes' must be a numeric vector with at least 2 elements")

  Num_Sim <- as.integer(Num_Sim)
  if (Num_Sim < 2){
    warning("'Num_Sim' must be an integer greater than one. Default value is used (10000)")
    Num_Sim <- 10000
  }

  if ((H != 0) && (H != 1))
    warning("The probability distribution table of Kruskal Wallis Statistic is not calculated. Set H = 1 for it")

  N_Samples <- length(Sample_Sizes)
  N_permutations <- factorial(sum(Sample_Sizes))
  for(i in 1:N_Samples) N_permutations <- N_permutations/factorial(Sample_Sizes[i])

  N <- sum(Sample_Sizes)
  b <- 0
  GP <- 0

  N_T <- 0
  for(i in 1:(N_g-1)){
    for (j in (i+1):N_g){
      N_T<- N_T + Sample_Sizes[i]*Sample_Sizes[j]
    }
  }

  for(i in 1:N_g)if((Sample_Sizes[i])%%2!=0)b<-b+1

  N_min <- 0
  for(i in 1:(N_g-1)){
    for (j in (i+1):N_g){
      N_min <- N_min + trunc(Sample_Sizes[i]*Sample_Sizes[j]/2)
    }
  }

  GB <- 0
  if(b>1){
    if(b%%2==0){
      l <- b/2
      GB <- l*(3*l-1)/2
    }else{
      l <- (b-1)/2
      GB <- l*(3*l+1)/2
    }
  }
  if(b==1){GB<-0}
  N_min <- N_min + GB

  opt <- c(rep(0,Num_Sim))
  rho <- c(rep(0,Num_Sim))
  KW <- c(rep(0,Num_Sim))

  for(n in 1:Num_Sim)
  {
    x <- rank(runif(N, 1, 1000000))

    for(j in 1:N)
    {
      temp <- 0
      i_temp <- 0
      for(i in 1:N_g){
        temp <- temp+Sample_Sizes[i]
        if((x[j]<=temp)&&(i_temp==0)){
          x[j] <- i
          i_temp <- 1;
        }
      }
    }

    pref <- matrix(rep(0,N_g*N_g), nrow=N_g, ncol=N_g)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        pref[x[i],x[j]] <- pref[x[i],x[j]] + 1
      }
    }

    opt[n] <- LOP(pref)$obj_val
    opt[n] <- N_T-opt[n]
    rho[n] <- 1 - opt[n]/(N_T-N_min)

    if(n==1) message("Simulations (x1000):")
    if(n%%1000==0) message(c(n/1000,"  "), appendLF = FALSE)

    if (H==1) {
      KW[n] <- sum(tapply(1:N, x, "sum")^2 / tapply(1:N, x, "length"))
      KW[n] <- 12*KW[n]/(N*(N+1)) -3*(N+1)
    }
  }
  message("")

  if (any(sapply(rho, is.na))) {
    C_freq <-  cbind(NA,NA,NA,NA)
    colnames(C_freq)<-c("disorder", "Concordance Coefficient", "Frequency", "Probability")

    H_freq <- cbind(NA,NA,NA)
    colnames(H_freq)<-c("H Statistic", "Frequency", "Probability")
  } else {

    rho <- round(sort(rho),2)
    count <- 1
    count2 <- 1
    rho2 <- 0
    rho2[1] <- rho[1]
    C_freq <- 0
    C_freq[1] <- count2

    for(i in 2:length(rho))
    {
      if(rho[i]==rho[i-1]){ count <- count +1}
      if(rho[i]>rho[i-1]){  C_freq[count2] <- count; count2 <- count2+1; rho2[count2] <- rho[i]; count <- 1 }
    }
    C_freq[count2] <- count

    if (H==1) {
      KW <- round(sort(KW),2)
      count <- 1
      count2 <- 1
      KW2 <- 0
      KW2[1] <- KW[1]
      H_freq <- 0
      H_freq[1] <- count2

      for(i in 2:length(KW))
      {
        if(KW[i]==KW[i-1]){ count <- count +1}
        if(KW[i]>KW[i-1]){  H_freq[count2]<- count; count2 <- count2+1; KW2[count2] <- KW[i]; count <- 1 }
      }
      H_freq[count2] <- count
      sum(H_freq)

      H_freq  <- round(H_freq*N_permutations/Num_Sim,0)
      H_freq <- cbind(KW2,H_freq,round((H_freq/sum(H_freq)),4))
      colnames(H_freq)<-c("H Statistic","Frequency","Probability")
    }

    C_freq  <- round(C_freq*N_permutations/Num_Sim,0)
    C_freq <- cbind(round(((N_T-N_min)*(1-rho2)),0),rho2,C_freq,round((C_freq/sum(C_freq)),4))
    colnames(C_freq)<-c("disorder"," Concordance Coefficient","Frequency","Probability")
  }

  if (H==1) return(list(C_freq=C_freq,H_freq=H_freq))
  else return(list(C_freq=C_freq))
}

#' @title Critical Values of Concordance and Kruskal Wallis Tests. Simulation Method
#' @description This function computes the critical values of the Concordance and Kruskal Wallis Tests. Simulated p-value of desired significance levels of .10, .05 and .01.
#' @param Sample_Sizes Numeric vector (\code{n1},...,\code{nk}) containing the number of repetitions of each element, i.e., the size of each sample in the experiment.
#' @param Num_Sim Number of simulations in order to obtain the probability distribution of the statistic. The default is 10000.
#' @param H 0 by default. If set to 1, the critical values of the Kruskal-Wallis Test are also calculated and returned.
#' @return The function returns a list with the following elements:
#' \enumerate{
#'  \item{ \code{C_results}: Concordance Coefficient results. Critical values and approximate p-values for a desired significance levels of 0.1, .05 and .01. }
#'  \item{ \code{H_results}: Kruskal Wallis results. Critical values and approximate p-values for a desired significance levels of 0.1, .05 and .01 (only if H = 1). }
#' }
#'@examples
#' Sample_Sizes <- c(4,4,2)
#' CT_Critical_Values_Sim(Sample_Sizes, Num_Sim = 1000)
#' CT_Critical_Values_Sim(Sample_Sizes, Num_Sim = 1000, H = 1)
#' @export
CT_Critical_Values_Sim<-function(Sample_Sizes, Num_Sim = 10000, H = 0){

  if(!is.numeric(Sample_Sizes))
    stop("Some elements of 'Sample_Sizes' are not numeric")

  if (length(Sample_Sizes) < 2)
    stop("'Sample_Sizes' must be a numeric vector with at least 2 elements")

  Num_Sim <- as.integer(Num_Sim)
  if (Num_Sim < 2){
    warning("'Num_Sim' must be an integer greater than one. Default value is used (10000)")
    Num_Sim <- 10000
  }

  if ((H != 0) && (H != 1)){
    warning("The critical values of the Kruskal Wallis Test are not calculated. Set H = 1 for it")
    H <- 0
  }

  distributions <- CT_Distribution_Sim(Sample_Sizes, Num_Sim, H)

  C_freq <- distributions$C_freq

  if (is.na(C_freq[1,1])) {
    C_results <- matrix(c(rep(NA,9)), nrow=3, byrow=TRUE)
    colnames(C_results)<-c("|  disorder", "|  Concordance Coefficient ","|  p-value")
    rownames(C_results)<-c("Sig level .10", "Sig level .05","Sig level .01")
    H_results <- matrix(c(rep(NA,6)), nrow=3, byrow=TRUE)
    colnames(H_results)<-c("|  H Statistic", "| p-value")
    rownames(H_results)<-c("Sig level .10", "Sig level .05", "Sig level .01")
  }
  else {

    C_result <- c(rep(0,9))
    for(k in 1:3){
      if(k==1)  alpha <- 0.1
      if(k==2)  alpha <- 0.05
      if(k==3)  alpha <- 0.01

      sum <- 0
      for(i in 1:(length(C_freq[,4])-1)){
        sum <- sum+C_freq[i,4]
        if( sum > (1-alpha) ){
          C_result[1+(k-1)*3] <- C_freq[i+1,1]
          C_result[2+(k-1)*3] <- C_freq[i+1,2]
          C_result[3+(k-1)*3] <- 1-sum
          alpha<- -1000
        }
      }
      i <- length(C_freq[,3])

      if(C_result[1+(k-1)*3]==0){
        if(k==1)  alpha <- 0.1
        if(k==2)  alpha <- 0.05
        if(k==3)  alpha <- 0.01
        if(C_freq[i,4]< alpha){
          C_result[1+(k-1)*3]<- C_freq[i,1]
          C_result[2+(k-1)*3]<- C_freq[i,2]
          C_result[3+(k-1)*3]<- C_freq[i,4]
        }
        else {
          C_result[1+(k-1)*3]<- NA
          C_result[2+(k-1)*3]<- NA
          C_result[3+(k-1)*3]<- NA
        }
      }
    }

    if (H == 1) {

      H_freq  <- distributions$H_freq

      H_result <- c(rep(0,6))
      for(k in 1:3){
        if(k==1)  alpha <- 0.1
        if(k==2)  alpha <- 0.05
        if(k==3)  alpha <- 0.01

        sum <- 0
        for(i in 1:(length(H_freq[,3])-1)){
          sum<-sum+H_freq[i,3]
          if( sum > (1-alpha) ){
            H_result[1+(k-1)*2] <- H_freq[i+1,1]
            H_result[2+(k-1)*2] <- 1-sum
            alpha<- -1000
          }
        }

        i <- length(H_freq[,3])

        if(H_result[1+(k-1)*2] == 0){
          if(k==1)  alpha <- 0.1
          if(k==2)  alpha <- 0.05
          if(k==3)  alpha <- 0.01
          if(H_freq[i,3]< alpha){
            H_result[1+(k-1)*2] <- H_freq[i,1]
            H_result[2+(k-1)*2] <- H_freq[i,3]
          }
          else {
            H_result[1+(k-1)*2] <- NA
            H_result[2+(k-1)*2] <- NA
          }
        }
      }

      H_results <- matrix(H_result, nrow=3, byrow = TRUE)
      colnames(H_results)<-c("|  H Statistic", "|  p-value")
      rownames(H_results)<-c("Sig level .10", "Sig level .05","Sig level .01")
    }

    C_results <- matrix(C_result, nrow=3, byrow = TRUE)
    colnames(C_results)<-c("|  disorder", "|  Concordance Coefficient", "|  p-value")
    rownames(C_results)<-c("Sig level .10", "Sig level .05", "Sig level .01")
  }

  if (H == 1) return(list(C_results=C_results, H_results=H_results))
  else return(list(C_results=C_results))
}

#' @import stats
#' @title Hypothesis test for testing whether samples originate from the same distribution
#' @description This function performs the hypothesis test for testing whether samples originate from the same distribution.
#' @param Sample_List List of numeric data vectors with the elements of each sample.
#' @param Num_Sim The number of used simulations. The default is 10000.
#' @param H 0 by default. If set to 1, the Kruskal-Wallis Test is also performed.
#' @return The function returns a list with the following elements:
#' \enumerate{
#' \item{ \code{results}: Table with the statistics and the signification levels.}
#'  \item{ \code{C_p-value}: Concordance test signification level. }
#'  \item{ \code{H_p-value}: Kruskal Wallis test signification level (only if H = 1). }
#' }
#'@examples
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#' Sample_List <- list(x, y, z)
#' CT_Hyphotesis_Test(Sample_List, Num_Sim = 1000, H = 1)
#'
#'
#' ## Example
#' A <- c(12,13,15,20,23,28,30,32,40,48)
#' B <- c(29,31,49,52,54)
#' C <- c(24,26,44)
#' Sample_List <- list(A, B, C)
#' CT_Hyphotesis_Test(Sample_List, Num_Sim = 1000, H = 1)
#'
#' ## Example with ties
#' A <- c(12,13,15,20,24,29,30,32,40,49)
#' B <- c(29,31,49,52,54)
#' C <- c(24,26,44)
#' Sample_List <- list(A, B, C)
#' CT_Hyphotesis_Test(Sample_List, Num_Sim = 1000, H = 1)
#' @export
CT_Hyphotesis_Test<-function(Sample_List, Num_Sim = 10000, H = 0){

  if(is.list(Sample_List)) {

    if (!all(sapply(Sample_List, is.numeric)))
      stop("Some elements of 'Sample_List' are not numeric")

    N_g <- length(Sample_List)
    Sample_Sizes <- lengths(Sample_List)

    if (N_g < 2)
      stop("'Sample_List' must be a list with at least 2 elements")

    if (any(Sample_Sizes == 0))
      stop("Some of the samples is empty")

    Num_Sim <- as.integer(Num_Sim)
    if (Num_Sim < 2){
      warning("'Num_Sim' must be an integer greater than one. Default value is used (10000)")
      Num_Sim <- 10000
    }

    if ((H != 0) && (H != 1))
      warning("The Kruskal Wallis Statistic is not calculated. Set H = 1 for it")

    groups <- rep(seq_len(N_g), Sample_Sizes)
    samples <- unlist(Sample_List)

    order_elements <- groups[order(samples)]
    sort_samples <- sort(samples)

    if (H == 1) {
      results <- matrix(rep(0,4),nrow = 2)
      results <- data.frame(results)
      rownames(results) <- c("Concordance Coefficient","Kruskal Wallis")
    } else {
      results <- matrix(rep(0,2),nrow = 1)
      results <- data.frame(results)
      rownames(results) <- c("Concordance Coefficient")
    }
    colnames(results) <- c("Statistic","p-value")

    N <- sum(Sample_Sizes)
    x <- c(1:N)
    b <- 0
    GP <- 0

    N_T <- 0
    for(i in 1:(N_g-1)){
      for (j in (i+1):N_g){
        N_T <- N_T + Sample_Sizes[i]*Sample_Sizes[j]
      }
    }

    for(i in 1:N_g) if((Sample_Sizes[i])%%2!=0) b <- b+1

    N_min <- 0
    for(i in 1:(N_g-1)){
      for (j in (i+1):N_g){
        N_min <- N_min + trunc(Sample_Sizes[i]*Sample_Sizes[j]/2)
      }
    }

    GB <- 0
    if(b>1){
      if(b%%2==0){
        l <- b/2
        GB <- l*(3*l-1)/2
      }else{
        l <- (b-1)/2
        GB <- l*(3*l+1)/2
      }
    }
    if(b==1){GB <- 0}
    N_min <- N_min + GB

    opt <- c(rep(0,Num_Sim))
    rho <- c(rep(0,Num_Sim))
    KW <- c(rep(0,Num_Sim))

    for(n in 1:Num_Sim){
      x <- rank(runif(N, 1, 1000000))

      for(j in 1:N){
        temp <- 0
        i_temp <- 0
        for(i in 1:N_g){
          temp <- temp+Sample_Sizes[i]
          if( (x[j]<=temp) && (i_temp==0) ){
            x[j] <- i
            i_temp <- 1;
          }
        }
      }
      pref <- matrix(rep(0,N_g*N_g), nrow=N_g, ncol=N_g)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          pref[x[i],x[j]] <-  pref[x[i],x[j]] + 1
        }
      }

      opt[n] <- LOP(pref)$obj_val

      if(n==1) message("Simulations (x1000):")
      if(n%%1000==0) message(c(n/1000,"  "), appendLF = FALSE)

      if (H == 1) {
         KW[n] <- sum(tapply(1:N, x, "sum")^2 / tapply(1:N, x, "length"))
         KW[n] <- 12*KW[n]/(N*(N+1)) -3*(N+1)
      }
    }
    message("")

    x <- order_elements
    opt_SUCESO <- 0
    pref <- matrix(rep(0,N_g*N_g), nrow=N_g, ncol=N_g)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        if (sort_samples[i] == sort_samples[j]) {
          pref[x[i],x[j]] <-  pref[x[i],x[j]] + .5
          pref[x[j],x[i]] <-  pref[x[j],x[i]] + .5
        } else pref[x[i],x[j]] <-  pref[x[i],x[j]] + 1
      }
    }
    opt_SUCESO <- LOP(pref)$obj_val

    CC <- 1 - (N_T-opt_SUCESO)/(N_T-N_min)
    results[1,1] <- round(CC,3)

    if (H == 1) {
      ties <- table(unlist(Sample_List))
      KW_SUCESO <- sum(tapply(rank(samples), groups, "sum")^2 / tapply(rank(samples), groups, "length"))
      KW_SUCESO <- (12*KW_SUCESO/(N*(N+1)) -3*(N+1))/(1 - (sum(ties^3 - ties)/(N^3 - N)))
      results[2,1] <- round(KW_SUCESO,3)

      COUNT_KW <- 0
      for(i in 1:Num_Sim)
        if(KW_SUCESO > KW[i]) COUNT_KW <- COUNT_KW +1

      p_value_KW <- 1 - COUNT_KW/Num_Sim

      results[2,2] <- p_value_KW
    }

    COUNT_KT <- 0
    for(i in 1:Num_Sim)
      if(opt_SUCESO > opt[i]) COUNT_KT <- COUNT_KT + 1

    p_value_KT <- 1 - COUNT_KT/Num_Sim

    results[1,2] <- p_value_KT

    if (H == 1) return(list(results=results, C_p_value=p_value_KT, H_p_value=p_value_KW))
    else return(list(results=results, C_p_value=p_value_KT))
  }
  else stop("'Sample_List' must be a list with at least 2 elements")
}


#' @title Concordance Coefficient and Kruskal Wallis Statistic
#' @description This function computes the Concordance Coefficient and Kruskal Wallis Statistic.
#' @param Sample_List List of numeric data vectors with the elements of each sample.
#' @param H 0 by default. If set to 1, the Kruskal Wallis Statistic is also calculated and returned.
#' @return The function returns a list with the following elements:
#' \enumerate{
#' \item{ \code{Sample_Sizes}: Numeric vector of sample sizes.}
#'  \item{ \code{order_elements}: Numeric vector containing the elements order. }
#' \item{ \code{disorder}: Disorder of the permutation given by \code{order_elements}.}
#'  \item{ \code{Concordance_Coefficient}: 1-relative disorder of permutation given by \code{order_elements}. }
#'  \item{ \code{H_Statistic}: Kruskal Wallis Statistic (only if H = 1). }
#' }
#'@examples
#' ## Example
#' A <- c(12,13,15,20,23,28,30,32,40,48)
#' B <- c(29,31,49,52,54)
#' C <- c(24,26,44)
#' Sample_List <- list(A, B, C)
#' CT_Coefficient(Sample_List)
#' CT_Coefficient(Sample_List, H = 1)
#'
#' ## Example with ties
#' A <- c(12,13,15,20,24,29,30,32,40,49)
#' B <- c(29,31,49,52,54)
#' C <- c(24,26,44)
#' Sample_List <- list(A, B, C)
#' CT_Coefficient(Sample_List, H = 1)
#' @export
CT_Coefficient<-function(Sample_List, H = 0){

  if(is.list(Sample_List)) {

    if (!all(sapply(Sample_List, is.numeric)))
      stop("Some elements of 'Sample_List' are not numeric")

    N_g <- length(Sample_List)
    Sample_Sizes <- lengths(Sample_List)

    if (N_g < 2)
      stop("'Sample_List' must be a list with at least 2 elements")

    if (any(Sample_Sizes == 0))
      stop("Some of the samples is empty")

    if ((H != 0) && (H != 1))
      warning("The Kruskal Wallis Statistic is not calculated. Set H = 1 for it")

    groups <- rep(seq_len(N_g), Sample_Sizes)
    samples <- unlist(Sample_List)

    order_elements <- groups[order(samples)]
    sort_samples <- sort(samples)

    N <- sum(Sample_Sizes)
    x <- c(1:N)
    b <- 0
    GP <- 0

    N_T <- 0
    for(i in 1:(N_g-1)){
      for (j in (i+1):N_g){
        N_T <- N_T + Sample_Sizes[i]*Sample_Sizes[j]
      }
    }

    for(i in 1:N_g) if((Sample_Sizes[i])%%2!=0) b <- b+1

    N_min <- 0
    for(i in 1:(N_g-1)){
      for (j in (i+1):N_g){
        N_min <- N_min + trunc(Sample_Sizes[i]*Sample_Sizes[j]/2)
      }
    }

    GB <- 0
    if(b>1){
      if(b%%2==0){
        l <- b/2
        GB <- l*(3*l-1)/2
      }else{
        l <- (b-1)/2
        GB <- l*(3*l+1)/2
      }
    }
    if(b==1){GB <-0}
    N_min <- N_min + GB

    x <- order_elements
    opt_SUCESO <- 0
    pref <- matrix(rep(0,N_g*N_g), nrow=N_g, ncol=N_g)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        if (sort_samples[i] == sort_samples[j]) {
          pref[x[i],x[j]] <-  pref[x[i],x[j]] + .5
          pref[x[j],x[i]] <-  pref[x[j],x[i]] + .5
        } else pref[x[i],x[j]] <-  pref[x[i],x[j]] + 1
      }
    }
    opt_SUCESO <- LOP(pref)$obj_val

    disorder <- (N_T-opt_SUCESO)

    rel_disorder <- 1 - (N_T-opt_SUCESO)/(N_T-N_min)

    if (H == 1){
      ties <- table(unlist(Sample_List))
      KW <- sum(tapply(rank(samples), groups, "sum")^2 / tapply(rank(samples), groups, "length"))
      KW <- (12*KW/(N*(N+1))-3*(N+1))/(1 - (sum(ties^3 - ties)/(N^3 - N)))

      return(list(Sample_Sizes=Sample_Sizes, order_elements=order_elements,
                            disorder=disorder, Concordance_Coefficient=rel_disorder,
                            H_Statistic=KW))
    }
    else return(list(Sample_Sizes=Sample_Sizes, order_elements=order_elements,
                     disorder=disorder, Concordance_Coefficient=rel_disorder))
  }
  else  stop("'Sample_List' must be a list with at least 2 elements")
}



