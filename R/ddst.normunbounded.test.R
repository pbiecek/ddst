#
#
#
#
#
#
# d.n = 20
# pols = hermite.h.polynomials(d.n)
# H.j.stan <- function(x, j) {
#   fpol = as.function(pols[[j + 1]])
#   2 ^ (-j / 2) * fpol(x / sqrt(2)) / sqrt(factorial(j))
# }
#
# delta.k = function(x, n, k) {
#   x.s = sort(x)
#   ps = 0:n / n
#   phi = dnorm(qnorm(ps))
#   h.j = H.j.stan(qnorm(ps), k)
#   h.j[1] = 0
#   h.j[n + 1] = 0
#   w.k = matrix(0, 1, n)
#   for (j in 1:n) {
#     w.k[1, j] = phi[j] * h.j[j] - phi[j + 1] * h.j[j + 1]
#   }
#   score = (w.k[1,] %*% x.s) / sqrt(k + 1)
#   return(score)
# }
# C.vec = function(x, n, d.n) {
#   score = matrix(0, 1, d.n)
#   for (j in 1:d.n)
#     score[1, j] = delta.k(x, n, j) * sqrt((j + 2) * n) / delta.k(x, n, 0)
#   return(score[1, ])
# }
# n = 100
# nr = 10000
# C.0 = matrix(0, d.n, nr)
# e.0 = matrix(0, 1, d.n)
# v.0 = matrix(0, 1, d.n)
#
# for (i in 1:nr) {
#   x = rnorm(n)
#   C.0[, i] = C.vec(x, n, d.n)
# }
#
# for (j in 1:d.n) {
#   e.0[1, j] = mean(C.0[j, ])
#   v.0[1, j] = sd(C.0[j, ]) * sqrt((n - 1) / n)
# }
#
# cal.R = function(x, n) {
#   S.sq = (sum((x - mean(x)) ^ 2)) / n
#   score = n * (1 - delta.k(x, n, 0) ^ 2 / S.sq)
#   return(score)
# }
#
# cal.N.test = function(x, n, d.n, r, e.0, v.0, s, so) {
#   cal.C.vec = (C.vec(x, n, d.n) - e.0) / v.0
#   cal.C.vec.sq = cal.C.vec ^ 2
#   cal.N.d.n = cumsum(cal.C.vec.sq)
#   d.n.o = ceiling(d.n / 2)
#   cal.Co.vec.sq = cal.C.vec.sq[2 * 1:d.n.o]
#   cal.No.d.n = cumsum(cal.Co.vec.sq)
#   dim = 1:d.n
#   dim.o = 1:d.n.o
#   cal.A.s = which.max(cal.N.d.n - dim * s)
#   cal.Ao.so = which.max(cal.No.d.n - dim.o * so)
#   cal.A.2 = which.max(cal.N.d.n - dim * 2)
#   cal.Ao.2 = which.max(cal.No.d.n - dim.o * 2)
#   r.n = cal.R(x, n)
#   if (r.n <= r) {
#     cal.A = cal.A.s
#     cal.Ao = cal.Ao.so
#   } else {
#     cal.A = cal.A.2
#     cal.Ao = cal.Ao.2
#   }
#   cal.N.cal.A.s = cal.N.d.n[cal.A.s]
#   cal.N.cal.A = cal.N.d.n[cal.A]
#   cal.N.cal.A.2 = cal.N.d.n[cal.A.2]
#   cal.No.d.n.cal.Ao.so = cal.No.d.n[cal.Ao.so]
#   cal.No.d.n.cal.Ao = cal.No.d.n[cal.Ao]
#   cal.No.d.n.cal.Ao.2 = cal.No.d.n[cal.Ao.2]
#   result = c(
#     cal.N.cal.A.s,
#     cal.N.cal.A,
#     cal.N.cal.A.2,
#     cal.A.s,
#     cal.A,
#     cal.A.2,
#     cal.No.d.n.cal.Ao.so,
#     cal.No.d.n.cal.Ao,
#     cal.No.d.n.cal.Ao.2,
#     cal.Ao.so,
#     cal.Ao,
#     cal.Ao.2,
#     cal.C.vec
#   )
#   return(result)
# }
