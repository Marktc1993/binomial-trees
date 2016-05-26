# Originally based on Binomial Trees from http://www.theresearchkitchen.com/archives/738
# Greatly expanded to actually price options.
# Oliver Frolovs, 2016

# The four valuators, the public interface to this library.
genlattice.vanilla.european.call <- function(Asset, Volatility, IntRate, DividendRate, Strike, Expiry, NoSteps) {
  return (genlattice.european(Asset=Asset, Volatility=Volatility, IntRate=IntRate, DividendRate=DividendRate, Strike=Strike, Expiry=Expiry, NoSteps=NoSteps, Payoff=payoff.vanilla.call))
}

genlattice.vanilla.european.put <- function(Asset, Volatility, IntRate, DividendRate, Strike, Expiry, NoSteps) {
  return (genlattice.european(Asset=Asset, Volatility=Volatility, IntRate=IntRate, DividendRate=DividendRate, Strike=Strike, Expiry=Expiry, NoSteps=NoSteps, Payoff=payoff.vanilla.put))
}

genlattice.vanilla.american.call <- function(Asset, Volatility, IntRate, DividendRate, Strike, Expiry, NoSteps) {
  return (genlattice.american(Asset=Asset, Volatility=Volatility, IntRate=IntRate, DividendRate=DividendRate, Strike=Strike, Expiry=Expiry, NoSteps=NoSteps, Payoff=payoff.vanilla.call))
}

genlattice.vanilla.american.put <- function(Asset, Volatility, IntRate, DividendRate, Strike, Expiry, NoSteps) {
  return (genlattice.american(Asset=Asset, Volatility=Volatility, IntRate=IntRate, DividendRate=DividendRate, Strike=Strike, Expiry=Expiry, NoSteps=NoSteps, Payoff=payoff.vanilla.put))
}

#
# Some examples, uncommented to be easy to copy.
#

if (FALSE) {
  # (Hull, 7th ed, p. 253), Example 11.1
  x <- genlattice.vanilla.european.call(Asset=810, Volatility=0.2, IntRate=0.05, DividendRate=0.02, Strike=800, Expiry=0.5, NoSteps=2)
  
  # (Hull, 7th ed, p. 250), Figure 11.10
  x <- genlattice.vanilla.american.put(Asset=50, Volatility = 0.3, IntRate = 0.05, DividendRate = 0, Strike = 52, Expiry = 2, NoSteps = 2)
  
  # Rendering and saving
  y <- capture.output(dotlattice(x, digits=4))
  cat(y, file="lattice.dot")
}

#
# The following functions are private.
#

# Generate a binomial lattice for European option.
genlattice.european <- function(Asset, Volatility, IntRate, DividendRate, Strike, Expiry, NoSteps, Payoff, Type) {
  
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

# Generate a binomial lattice for American option.
genlattice.american <- function(Asset, Volatility, IntRate, DividendRate, Strike, Expiry, NoSteps, Payoff) {
  
  # The number of tree nodes to process.
  count <- sum(1 : (NoSteps+1))
  
  # This data frame will store asset and option prices.
  # The mapping from tree node (i,j) to linear index
  # inside the data frame will have to be computed.
  # The early exercise flag is also stored.
  X <- data.frame(matrix(NA, nrow=count, ncol=3))
  names(X) <- c("asset", "option", "exercise")
  
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
      AssetCurrent <- Asset * u^(i-j) * d^j
      X$asset[count] <- AssetCurrent
      
      # Compute the payoff directly for the last step's nodes,
      # otherwise use a formula. 
      # FIXME do we care about "early exercise" flag at final nodes?
      if (i == NoSteps) {
        X$option[count]   <- Payoff(AssetCurrent, Strike)
        X$exercise[count] <- (X$option[count] > 0)
      } else {
        up   <- X$option[sum(1:(i+1), j, 1)]
        down <- X$option[sum(1:(i+1), j+1, 1)]
        
        # Possible option values when discounted or when early exercise is applied
        V <- DiscountFactor * (p * up + (1-p) * down)
        V.early <- Payoff(AssetCurrent,Strike)
        
        # The greatest of two possible values is stored
        X$option[count] <- max(V, V.early)
        
        # Should the option be exercised early?
        X$exercise[count]  <- (V.early > V)
      }
      
      count <- count - 1
    }
  }
  
  return(X)
}

# Generates a graph specification that can be fed into graphviz.
# Input: the binomial lattice produced by one of genlattice family functions.
dotlattice <- function(S, digits=2) {
  
  shape <- "plaintext"
  
  cat("digraph G {", "\n", sep="")
  cat("node[shape=",shape,"];","\n", sep="")
  cat("rankdir=LR;","\n")
  
  cat("edge[arrowhead=none];","\n")
  
  # Create a dot node for each element in the lattice
  for (i in 1:nrow(S)) {
    x <- round(S$asset[i], digits=digits)
    y <- round(S$option[i], digits=digits)
    
    # Detect the American tree and draw accordingly
    early.exercise <- ""
    if (("exercise" %in% colnames(S)) && S$exercise[i]) {
      early.exercise <- "shape=oval,"
    }
    
    cat("node", i, "[", early.exercise, "label=\"", x, ", ", y, "\"];", "\n", sep="")
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



#
# A fistful of option payoff functions of different types.
# 

payoff.vanilla.call <- function(Asset, Strike) {
  return( max(0, Asset - Strike) )
}

payoff.vanilla.put <- function(Asset, Strike) {
  return( max(0, Strike - Asset) )
}

