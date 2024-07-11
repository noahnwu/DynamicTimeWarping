DTW.Forrest = function(x,y, DistFun = proxy::dist, 
               dist.only = T, ...) { 
  
  #compute distance matrix
  #note that the DistFun specified must return a numberic matrix with dimensions x,y
  #where x and y are the length of the input vectors.
  #overall performance will be better if DistFun is vectorized
  dist.mat = DistFun(x,y, ...)
  nm = sum(dim(dist.mat))
  #add a row and a column with infinite values. 
  #this is not strictly necessary. Might function out of the box without 
  #this step
  dist.mat=cbind(Inf, dist.mat)
  dist.mat = rbind(Inf, dist.mat)
  #initialize empty cost matrix
  cost.mat = matrix(nrow = nrow(dist.mat), 
                    ncol = ncol(dist.mat))
  rownames(dist.mat) = NULL
  colnames(dist.mat) = NULL
  cost.mat[1, ] = Inf
  cost.mat[, 1] = Inf
  cost.mat[1, 1] = 0
  dir.track = c()
  counter = 1
  for(i in 2:nrow(dist.mat)) {
    for(j in 2:ncol(dist.mat)) {
      if((i == 2 & j == 2)) {
        #do not need to weight the diagonal
        #on the first step (first points in sequence will always align)
        neighbors = c(
          "oblique" = cost.mat[i-1, j-1] + dist.mat[i, j],
          "left" = cost.mat[i-1, j] + dist.mat[i, j],
          "up" = cost.mat[i, j-1] + dist.mat[i, j]
        ) 
      } else { 
        neighbors = c(
          #oblique movements through the cost matrix are penalized 
          #('symmetric2' step pattern in dtw package)
          #that allows for normalization of output distance
          "oblique" = cost.mat[i-1, j-1] + 2*dist.mat[i, j],
          "left" = cost.mat[i-1, j] + dist.mat[i, j],
          "up" = cost.mat[i, j-1] + dist.mat[i, j]
        ) 
      }
      
      cost.mat[i, j] = min(neighbors)
      if(dist.only == F) { 
        dir.track[counter] = names(neighbors)[which.min(neighbors)]
        counter = counter + 1
      }
    }
    
    
  }
  if(dist.only == F) { 
    warp.mat = matrix(dir.track, 
                      nrow = nrow(cost.mat)-1, 
                      byrow = T)
    index_1 = c()
    index_2 = c()
    counter = 1
    cost.mat = cost.mat[-1, -1]
    i = nrow(cost.mat)
    j = ncol(cost.mat)
    while((i+j) > 1) { 
      index_1[counter] = i
      index_2[counter] = j
      counter = counter + 1
      if(warp.mat[i,j] == "up") { 
        j = j-1
      } else if(warp.mat[i,j] == "left") { 
        i = i-1   
      } else {
        i = i-1
        j= j-1
      }
    }
    
    return(list("distance" = cost.mat[nrow(cost.mat), ncol(cost.mat)], 
                "normalized.distance" = cost.mat[nrow(cost.mat), ncol(cost.mat)]/nm, 
                index_1 = rev(index_1), 
                index_2 = rev(index_2), 
                "cost.mat" = cost.mat, 
                "dist.mat" = dist.mat[-1, -1]))
    
  } else { 
    return(list(
      "distance" = cost.mat[nrow(cost.mat), ncol(cost.mat)],
      "normalized.distance" = cost.mat[nrow(cost.mat), ncol(cost.mat)]/nm, 
      "cost.mat" = cost.mat[-1, -1]
    ))
  }
  
  
} 

