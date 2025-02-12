# Various -----------------------------------------------------------------
get_Yij0 <- function(i, j, Y){
  Yij0 <- sum(Y[((i - 1)*mi + j), ])
  return(Yij0)
}

get_Bij <- function(i, j, B){
  Bij <- B[((i - 1)*mi + j), ]
  return(Bij)
}

get_Z_ijl <- function(i, j, l, Z){
  Z_ijl <- Z[i + (j - 1), (l + 1)]
  return(Z_ijl)
}
