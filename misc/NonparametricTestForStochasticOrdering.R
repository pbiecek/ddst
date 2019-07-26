# Calculation of Md, QT and other quantities

l.j = function(k, i, z) {
  a.j = (2 * i - 1) / 2 ^ (k + 1)
  left = -sqrt((1 - a.j) / a.j)
  right = sqrt(a.j / (1 - a.j))
  score = ifelse(z < a.j, left, right)
  return(score)
}

L.j = function(k, i, x, y, m, n) {
  N = m + n
  H = ecdf(c(x, y))
  rx = H(x) - 1 / (2 * N)
  ry = H(y) - 1 / (2 * N)
  cx = sum(l.j(k, i, rx)) / m
  cy = sum(l.j(k, i, ry)) / n
  score = (cy - cx) * sqrt(m * n / N)
  return(score)
}


test.M.d = function(x, y, m, n, k.N) {
  N = m + n
  d.N = 2 ^ (k.N + 1) - 1
  L.vec = matrix(0, 1, d.N)
  for (k in 0:k.N) {
    for (i in 1:(2 ^ k)) {
      j = 2 ^ k - 1 + i
      L.vec[1, j] = L.j(k, i, x, y, m, n)
    }
  }
  M.d = min(L.vec[1, ])
  result = c(M.d, L.vec)
  return(result)
}

test.Q = function(x, y, m, n, k.N, t) {
  N = m + n
  d.N = 2 ^ (k.N + 1) - 1
  L.vec = matrix(0, 1, d.N)
  L.vec.trun = matrix(0, 1, d.N)
  for (k in 0:k.N) {
    for (i in 1:(2 ^ k)) {
      j = 2 ^ k - 1 + i
      L.vec[1, j] = L.j(k, i, x, y, m, n)
      L.vec.trun[1, j] = max(-L.vec[1, j], 0)
    }
  }
  Q.d.vec = matrix(0, 1, k.N + 1)
  D.vec = 2 ^ (0:k.N + 1) - 1
  for (k in 1:(k.N + 1)) {
    Q.d.vec[1, k] = L.vec.trun[1, 1:D.vec[k]] %*% L.vec.trun[1, 1:D.vec[k]]
  }
  S = which.max(Q.d.vec[1, ] - D.vec * log(N))
  M = which.max(Q.d.vec[1, ])
  if (max(-L.vec[1, ]) <= sqrt(t * log(N))) {
    T = S
  } else{
    T = M
  }
  Q.S = Q.d.vec[1, S]
  Q.T = Q.d.vec[1, T]
  result = c(Q.T, Q.S, T, S)
  return(result)
}



