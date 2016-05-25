# Originally based on Binomial Trees from http://www.theresearchkitchen.com/archives/738
# Greatly expanded to actually price options.
# Oliver Frolovs, 2016

# Generate a binomial lattice for option pricing problem
# as specified by the function's arguments.
# The Payoff argument expects a function of the form f(Asset, Strike).
# The default argument values match the example 11.1 in (Hull, 2007) page 253.
genlattice <- function(Asset=810, Volatility=0.2, IntRate=0.05, DividendRate=0.02, Strike=800, Expiry=0.5, NoSteps=2, Payoff=euro_call_payoff) {
  
  # The number of tree nodes to process.
  count <- sum(1 : (NoSteps+1))
  
  # This data frame will store asset and option prices.
  # The mapping from tree node (i,j) to linear index
  # inside the data frame will have to be computed.
  X <- data.frame(matrix(NA, nrow=count, ncol=2))
  names(X) <- c("asset", "option")
  
  # Time between price movements.
  dt = Expiry / NoSteps
  
  # Option price discount factor.
  DiscountFactor <- exp(-IntRate * dt)
  
  # The up and down jump factors and
  # corresponding (synthetic) probabilities using 
  # Cox, Ross, and Rubinstein (1979) method.
  u = exp(Volatility * sqrt(dt))
  d = 1/u
  a = exp((IntRate - DividendRate) * dt)
  p = (a - d) / (u - d)
  
  # Compute the asset and option prices, starting 
  # from the last node of the tree, which is
  # its bottom right corner when viewed as a graph.
  # Work up and backwards. Backwards, comrades!
  for (i in NoSteps:0) {
    for (j in i:0) {
      X$asset[count] <- Asset * u^(i-j) * d^j
      
      # Compute the payoff directly for the last step's nodes,
      # otherwise use a formula.
      if (i == NoSteps) {
        X$option[count] <- Payoff(X$asset[count], Strike)
      } else {
        up   <- X$option[sum(1:(i+1), j, 1)]
        down <- X$option[sum(1:(i+1), j+1, 1)]
        X$option[count] <- DiscountFactor * (p * up + (1-p) * down)
      }
      
      count <- count - 1
    }
  }
  
  return(X)
}

# generate a graph specification that can be fed into graphviz:
dotlattice <- function(S, digits=2, labels=TRUE) {
  shape <- ifelse(labels == TRUE, "plaintext", "point")
  
  cat("digraph G {", "\n", sep="")
  cat("node[shape=",shape,"];","\n", sep="")
  cat("rankdir=LR;","\n")
  
  cat("edge[arrowhead=none];","\n")
  
  # Create a dot node for each element in the lattice
  for (i in 1:nrow(S)) {
    x <- round(S$asset[i], digits=digits)
    y <- round(S$option[i], digits=digits)
    cat("node", i, "[label=\"", x, ", ", y, "\"];", "\n", sep="")
  }
  
  # The number of levels in a binomial lattice of length N
  # is `$\frac{\sqrt{8N+1}-1}{2}$`
  L <- ((sqrt(8*nrow(S)+1)-1)/2 - 1)
  
  k<-1
  for (i in 1:L) {
    tabs <- rep("\t",i-1)
    j <- i
    while(j>0) {
      cat("node",k,"->","node",(k+i),";\n",sep="")
      cat("node",k,"->","node",(k+i+1),";\n",sep="")
      k <- k + 1
      j <- j - 1
    }
  }
  
  cat("}", sep="")
}

# x <- capture.output(dotlattice(genlattice(N=8, u=1.1, d=0.9)))
# cat(x, file="lattice1.dot")

#
# A fistful of option payoff functions of different types.
# 

# Vanilla European call payoff
euro_call_payoff <- function(Asset, Strike) {
  return( max(0, Asset - Strike) )
}

# Vanilla European put payoff
euro_put_payoff <- function(Asset, Strike) {
  return( max(0, Strike - Asset) )
}
