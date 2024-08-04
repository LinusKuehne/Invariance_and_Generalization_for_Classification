# this script contains utility functions for the experiments on synthetic data


# generate uniformly dist. samples from [-max, -min] union [min, max]
# n: number of samples
# min: non-negative
# max: non-negative, satisfies max > min
runifstrong <- function(n, min, max){
  
  # sample the absolute value
  abs.samples <- runif(n, min, max)
  
  # randomly sample the sign
  sign <- sample(c(-1,1), size = n, replace = T)
  
  return(sign*abs.samples)
}



# computes binary cross-entropy
# y: label in {0,1} 
# y.hat: probability in (0,1) (probability that Y=1 given some predictors)
BCE <- function(y, y.hat){
  
  y.hat.norm <- y.hat
  
  # cap extreme values such that BCE is still defined
  y.hat.norm[y.hat == 1] <- 0.99999999999
  y.hat.norm[y.hat == 0] <- 0.00000000001
  
  nbce <- y*log(y.hat.norm) + (1-y)*log(1-y.hat.norm)
  
  return(-mean(nbce))
}



# computes weighted binary cross-entropy
# y: label in {0,1} 
# y.hat: probability in (0,1) (probability that Y=1 given some predictors)
BCE.weighted <- function(y, y.hat){

  # proportion of 1's
  ones <- mean(y) 
  
  # proportion of 0's
  zeros <- 1-ones 
  
  # compute weights
  w1 <- 1/(2*ones)
  w0 <- 1/(2*zeros)
  
  # cap extreme values such that the log is defined
  y.hat.norm <- y.hat
  y.hat.norm[y.hat > 0.99999999999] <- 0.99999999999
  y.hat.norm[y.hat < 0.00000000001] <- 0.00000000001
  
  nbce <- w1*y*log(y.hat.norm) + w0*(1-y)*log(1-y.hat.norm)
  
  return(-mean(nbce))
}





# compute linear combination of columns in dat using the vector coef as weights
# dat: matrix whose columns should be linearly combined
# coef: coefficients used for the linear combination
lincomb <- function(dat, coef){
  res <- numeric(nrow(dat))
  
  for(j in 1:length(coef)){
    res <- res + coef[j]*dat[,j]
  }
  return(res)
}






#-------------------------------------------------------------------------------
# methods for ICP
#-------------------------------------------------------------------------------


# compute output of ICP (intersection of all sets for which invariance can't be rejected)
# pvals: vector of p-values corresponding to the sets in the list "sets" (global variable)
# alpha: level of test (here: tuning parameter)
calc.ICP.output <- function(pvals, alpha = 0.05){
  
  inv.sets <- sets[pvals > alpha]

  # Compute the intersection of all vectors in the list
  intersection.result <- Reduce(intersect, inv.sets)
  
  return(intersection.result)

}




#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# surrogate measures for level+power 
#-------------------------------------------------------------------------------


# returns 1 if s.hat \subseteq parents, and 0 otherwise
# s.hat: vector of variables (ICP output)
# parents: vector of variables (parents of response)
FWE <- function(s.hat, parents){
  r <- all(s.hat %in% parents)
  return(as.numeric(!r))
}




# returns Jaccard(s.hat, parents)
# s.hat: vector of variables (ICP output)
# parents: vector of variables (parents of response)
jaccard.ICP <- function(s.hat, parents){
  
  # deal with the cases of empty sets
  if(length(s.hat) == 0){
    if(length(parents) == 0){
      return(1)
    } else{
      return(0)
    }
  }
  
  # compute size of intersection 
  intersection <- length(intersect(s.hat, parents))
  
  # compute size of union
  union <- length(parents) + length(s.hat) - intersection
  
  # return Jaccard index
  return(intersection/union)
}



# jaccard index specifically for the standard DAG
# pvalues: vector of length 7 with entries corresp. to the sets
#         list(c(1,2,3), c(1,3), c(1,2), c(2,3), c(1), c(2), c(3), c(0))
# thresh: threshold such that for a pvalue above thresh, we say set is "invariant"
jaccard.standard <- function(pvalues, thresh = 0.05){

  # ground truth invariant sets
  gt <- c(2,5)
  
  # predicted invariant subsets
  pred <- which(pvalues > thresh)
  
  # compute size of the intersection
  intersection <- length(intersect(gt, pred))
  
  # compute size of the union
  union <- length(gt) + length(pred) - intersection
  
  # return Jaccard index
  return(intersection/union)
}






# compute jaccard score for general DAG
# pvalues: vector with entries corresp. to the sets powerSet(1:num.covariates) (using the rje library)
# thresh: threshold such that for a pvalue above thresh, we say the set is "invariant"
# stable.sets: sets which satisfy the d-separation corresponding to invariance (subset of powerSet("X1", "X2", ...))
# num.covariates are the number of covariates (how many X-variables we have)
jaccard <- function(pvalues, thresh = 0.05, stable.sets, num.covariates){

  # generate vector of variable names
  Xnames <- rep("A", num.covariates)
  for(w in 1:num.covariates){
    number <- as.character(w)
    name <- paste("X", number, sep="") 
    Xnames[w] <- name  
  }
  
  # list of all subsets
  ps <- powerSet(Xnames)
  
  # get the indices for the stable sets
  gt <- (1:length(ps))[ps %in% stable.sets]
  
  # the predicted invariant subsets
  pred <- which(pvalues > thresh)
  
  # compute the size of the intersection
  intersection <- length(intersect(gt, pred))
  
  # if A = B = empty, then J(A,B) := 1
  if((length(gt)==0) && (length(pred)==0)){
    return(1)
  }
  
  if(intersection==0){
    return(0)
  }
  
  # compute the size of the union 
  union <- length(gt) + length(pred) - intersection
  
  # return the Jaccard index
  return(intersection/union)
  
}




#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# stuff relating to DAGs
#-------------------------------------------------------------------------------


# generate a vector of ordered variable names with the covariates (and Y) first and then the interventions
# d: number of covariates (including Y)
# Y.ind: index for the covariate which should be Y
# num.int: number of interventions
generate.varnames <- function(d, Y.ind, num.int){
  
  # total number of variables
  tot.length <- d + num.int
  
  var.names <- rep("A", (d+num.int))
  
  cov.number <- 0
  for(i in 1:d){
    if(i==Y.ind){
      var.names[Y.ind] <- "Y"
    } else{
      cov.number <- cov.number + 1
      number <- as.character(cov.number)
      name <- paste("X", number, sep="")
      var.names[i] <- name
    }
  }
  
  
  for(i in (d+1):tot.length){
    number <- as.character(i-d)
    name <- paste("I",number, sep="")
    var.names[i] <- name
  }
  return(var.names)
}






# computes the stable blanket
# dag.full: dagitty dag including I's
# dag.cov: dagitty dag without I's
# num.int: number of interventions
# number of covariates (including Y)
stableBlanket <- function(dag.full, dag.cov, num.int, d){
  
  
  # covariate names (just "X1", "X2", ...)
  cov.names <- rep("A", d-1)
  for(w in 1:(d-1)){
    number <- as.character(w)
    name <- paste("X", number, sep="") 
    cov.names[w] <- name  
  }
  
  # find variables which have been intervened on 
  int.vars <- NULL
  for(w in 1:num.int){
    number <- as.character(w)
    name <- paste("I", number, sep="")
    int.vars <- append(int.vars, children(x=dag.full, v=name))
  }
  int.vars <- unique(int.vars)
  
  # children of Y
  CH.Y <- unique(children(x=dag.cov, v="Y"))
  
  # children of Y which have been intervened on
  CH.Y.int <- intersect(CH.Y, int.vars)
  
  # complement of Nint (i.e. all children of Y which have been intervened on and all descendants of such children)
  complement.N.int <- descendants(x=dag.cov, v=CH.Y.int)
  
  # Nint
  N.int <- setdiff(x=cov.names, y = complement.N.int)
  
  # we now construct the three sets of the stable blanket
  set1 <- parents(x=dag.cov, v="Y")
  set2 <- intersect(N.int, CH.Y)
  set3 <- setdiff(parents(x=dag.cov, v=set2), c("Y"))
  
  # return stable blanket (ordered)
  SB.Y <- mixedsort(union(set1, union(set2, set3)))
  
  return(SB.Y)
}





# compute list of all intervention stable sets (sets which satisfy the d-separation associated with invariance)
# dag.full: dagitty dag object including all covariates, the response, and the interventions
# d: number of covariates including Y
# num.int: number of interventions
int.stable.sets <- function(dag.full, d, num.int){
  
  # list of possible variables
  cov.names <- rep("A", d-1)
  for(k in 1:(d-1)){
    number <- as.character(k)
    name <- paste("X", number, sep="") 
    cov.names[k] <- name  
  }
  
  # vector of interventions
  I.vec <- rep("A", num.int)
  for(k in 1:num.int){
    number <- as.character(k)
    name <- paste("I",number, sep="")
    I.vec[k] <- name
  }
  
  # compute all possible subsets of covariates
  sets <- powerSet(x=cov.names)
  
  # initialize list of intervention stable sets
  int.stable <- list()
  
  n.sets <- length(sets)
  
  stab.counter <- 1
  for(k in 1:n.sets){
    
    # check whether I dsep Y | S
    stable <- dseparated(x = dag.full, X = I.vec, Y = "Y", Z = sets[[k]])
    
    if(stable){
      # add subset to list
      int.stable[[stab.counter]] <- sets[[k]]
      
      # increase counter
      stab.counter <- stab.counter + 1
    }
  }
  
  return(int.stable)
}


# compute Markov blanket
# dag.cov: dagitty dag object including all covariates and the response
markovBlanket <- function(dag.cov){
  
  # compute the Markov blanket with dagitty, and order it according to our needs
  re <- mixedsort(dagitty::markovBlanket(x=dag.cov, v="Y"))
  
  return(re)
}



# returns oracle ICP output using formula by Mogensen, Thams, Peters (2022) from Invariant Ancestry Search paper
# dag.cov: dagitty dag object including all covariates and the response
# dag.full: dagitty dag object including all covariates, the response, and the interventions
# num.int: number of interventions
ICP.oracle <- function(dag.cov, dag.full, num.int){
  
  
  # int names (just "I1", "I2", ...)
  int.names <- rep("A", num.int)
  for(w in 1:num.int){
    number <- as.character(w)
    name <- paste("I", number, sep="") 
    int.names[w] <- name  
  }
  
  
  # parents of Y
  pa.Y <- parents(x = dag.cov, v = "Y")
  
  # children of E
  ch.E <- children(x = dag.full, v = int.names)
  
  # ancestors of Y
  an.Y <- ancestors(x = dag.cov, v = "Y", proper = T)
  
  # intersection of the ancestors of Y and the children of E
  an.ch <- intersect(an.Y, ch.E)
  
  # get the parents of the nodes an.ch
  pa.an.ch <- parents(x = dag.cov, v = an.ch)
  
  # compute the union of ch.E and pa.an.ch
  un <- union(ch.E, pa.an.ch)
  
  # oracle ICP output
  s.ICP <- intersect(pa.Y, un)
  
  return(s.ICP)
}





#-------------------------------------------------------------------------------
