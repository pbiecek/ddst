#' Data Driven Smooth Test for Upward Trend Alternatives
#'
#' Performs data driven smooth test for upward trend in k-sample problem.
#' Suppose that we have random samples from k distributions F_i where i = 1, ..., k.
#' The null hypothesis is that there is lack of trend,
#' i.e. F_1 >= ... >= F_k and F_i != F_j for some i and j.
#' The alternative is that there is a trend
#' i.e. F_1 >= ... >= F_k and F_i != F_j for some i and j.
#' This test is implemented as a special case of an umbrella test.
#'
#' @param x a list of k (non-empty) numeric vectors of data
#' @param r.N a (k-1)-dimensional vector specifying the levels of complexity of the grids considered, only for advanced users
#' @param alpha a significance level
#' @param t.p an alpha-dependent (k-1)-dimensional vector of the tunning parameters in the penalties in the model selection rules T.o
#' @param t.n an alpha-dependent (k-1)-dimensional vector of the tunning parameters in the penalties in the model selection rules T.tilde
#' @param t an alpha-dependent tunning parameter in the penalty in the model selection rule
#' @param B an integer specifying the number of runs for a p-value and a critical value computation if any
#' @param compute.cv a logical value indicating whether to compute a critical value corresponding to the significance level alpha or not
#'
#' @references An automatic test for the umbrella alternatives. Wylupek (2016) \url{https://onlinelibrary.wiley.com/doi/abs/10.1111/sjos.12231}
#' @export
#' @examples
#' set.seed(7)
#' # H0 is true
#' x = runif(80)
#' y = runif(80) + 0.2
#' z = runif(80) + 0.4
#' t <- ddst.upwardtrend.test(list(x, y, z), t.p = 2.2, t.n = 2.2)
#' t
#' plot(t)
#'
#' # H0 is false
#' # known fixed alternative
#' x1 = rnorm(80)
#' x2 = rnorm(80) + 2
#' x3 = rnorm(80) + 4
#' x4 = rnorm(80) + 3
#' t <- ddst.upwardtrend.test(list(x1, x2, x3, x4), t.p = 2.2, t.n = 2.2)
#' t
#' plot(t)
#'
#' @keywords htest
`ddst.upwardtrend.test` <-
  function(x,
           r.N  = rep(4, length(x)-1),
           alpha = 0.05,
           t.p, t.n,
           B = 10000, compute.cv = FALSE) {   # tlh.p = 2.2, tl.n = 2.2
##    if (is.list(x)) {
      # x is list with coordinates
      x.vector <- unlist(x)
      n <- sapply(x, length)
##    } else {
##      # x is a vector, other vectors are in ...
##      x.vector <- unlist(c(x, list(...)))
##      n <- c(length(x),
##             sapply(list(...), length))
##    }
    p = length(n)

    # test for trend is implemented as a special case of the umbrella test
    coord = ddst.umbrella.Nk(x.vector, n,
                             tlh.p = t.p, tl.n = t.n,
                             p = p)

    scoresflat <- c(coord$score.mat)
    scoresnames <- paste0(rep(1:length(n), times = length(n)),
                          "-",
                          rep(1:length(n), each = length(n)))
    scoresnames <- scoresnames[!is.na(scoresflat)]
    scoresflat  <- scoresflat[!is.na(scoresflat)]
    names(scoresflat) <- scoresnames


    l = coord$U.T
    attr(l, "names") = "U.T"
    t = p
    attr(t, "names") = "p"
    result = list(statistic = l,
                  parameter = t,
                  coordinates = scoresflat,
                  method = "Data Driven k-Sample Upward Trend Test")
    result$data.name = paste(paste(as.character(substitute(x)), collapse = ""),
                             ", t.p: ",
                             t.p,
                             ", t.n: ",
                             t.n,
                             sep = "")
    class(result) = c("htest", "ddst.test", "ddst.upward.test")

    result
  }



# l.j = function(r,i,v) {
#   a.j = (2*i-1)/2^(r+1)
#   left = -sqrt((1-a.j)/a.j)
#   right = sqrt(a.j/(1-a.j))
#   score = ifelse(v < a.j, left, right)
#   return(score)
# }
#
# L.j.l = function(r,i,xl,xl1,nl,nl1){
#   Nl = nl + nl1
#   Hl = ecdf(c(xl,xl1))
#   rxl = Hl(xl) - 1/(2*Nl)
#   rxl1 = Hl(xl1) - 1/(2*Nl)
#   cxl = mean( l.j(r,i,rxl) )
#   cxl1 = mean( l.j(r,i,rxl1) )
#   score = (cxl1 - cxl)*sqrt(nl*nl1/Nl)
#   return(score)
# }
# r.Nl = 4
# d.Nl = 2^(r.Nl+1)-1
# test.Ul = function(xl,xl1,nl,nl1,tl.p,tl.n,alpha.l) {
#   Nl = nl+nl1
#   L.vec.l = matrix(0,1,d.Nl)
#   L.vec.l.p = matrix(0,1,d.Nl)
#   L.vec.l.n = matrix(0,1,d.Nl)
#   for(r in 0:r.Nl){
#     for(i in 1:(2^r)){
#       j = 2^r - 1 + i
#       L.vec.l[1,j] = L.j.l(r,i,xl,xl1,nl,nl1)
#       L.vec.l.p[1,j] = max( L.vec.l[1,j], 0 )
#       L.vec.l.n[1,j] = max( -L.vec.l[1,j], 0 )
#     }
#   }
#   Q.d.vec.l.p = matrix(0,1,r.Nl+1)
#   Q.d.vec.l.n = matrix(0,1,r.Nl+1)
#   D.vec.l = 2^(0:r.Nl+1) - 1
#   for(r in 1:(r.Nl+1)){
#     Q.d.vec.l.p[1,r] = L.vec.l.p[1,1:D.vec.l[r]] %*% L.vec.l.p[1,1:D.vec.l[r]]
#     Q.d.vec.l.n[1,r] = L.vec.l.n[1,1:D.vec.l[r]] %*% L.vec.l.n[1,1:D.vec.l[r]]
#   }
#   Sl.p = which.max( Q.d.vec.l.p[1,] - D.vec.l*log(Nl) )
#   Ml.p = which.max( Q.d.vec.l.p[1,] )
#   if( max( L.vec.l.p[1,] ) <= sqrt(tl.p*log(Nl)) ){
#     Tl.p = Sl.p
#   }else{
#     Tl.p = Ml.p
#   }
#   Q.Sl.p = Q.d.vec.l.p[1,Sl.p]
#   Q.Tl.p = Q.d.vec.l.p[1,Tl.p]
#   Sl.n = which.max( Q.d.vec.l.n[1,] - D.vec.l*log(Nl) )
#   Ml.n = which.max( Q.d.vec.l.n[1,] )
#   if( max( L.vec.l.n[1,] ) <= sqrt(tl.n*log(Nl)) ){
#     Tl.n = Sl.n
#   }else{
#     Tl.n = Ml.n
#   }
#   z.alpha.l = qnorm(alpha.l)
#   Z.Sl.n = sign( min(L.vec.l[1,1:D.vec.l[Sl.n]] - z.alpha.l) )
#   Z.Tl.n = sign( min(L.vec.l[1,1:D.vec.l[Tl.n]] - z.alpha.l) )
#   U.Sl = Z.Sl.n * Q.Sl.p
#   U.Tl = Z.Tl.n * Q.Tl.p
#   result = c(U.Tl, Z.Tl.n, Q.Tl.p, U.Sl, Z.Sl.n, Q.Sl.p)
#   return(result)
# }
# test.U = function(x,n,tl.p,t.n,alpha){
#   k = length(n)
#   N = sum(n)
#   n.cum = c(0,cumsum(n))
#   p = matrix(0,1,k-1)
#   Z.T.n.vec = matrix(0,1,k-1)
#   Z.S.n.vec = matrix(0,1,k-1)
#   Q.T.p.vec = matrix(0,1,k-1)
#   Q.S.p.vec = matrix(0,1,k-1)
#   n1 = n[1]
#   nk = n[k]
#   coeff = (k-1)/(2*N - (n1+nk))
#   for(l in 1:(k-1)) {
#     nl = n[l]
#     nl1 = n[l+1]
#     Nl = nl + nl1
#     p[1,l] = Nl/N
#     xl = x[(n.cum[l]+1):n.cum[l+1]]
#     xl1 = x[(n.cum[l+1]+1):n.cum[l+2]]
#     tl.n = t.n[l]
#     wl = coeff*Nl
#     alpha.l = wl*alpha
#     statistics = test.Ul(xl,xl1,nl,nl1,tl.p,tl.n,alpha.l)
#     Z.T.n.vec[1,l] = statistics[2]
#     Q.T.p.vec[1,l] = statistics[3]
#     Z.S.n.vec[1,l] = statistics[5]
#     Q.S.p.vec[1,l] = statistics[6]
#   }
#   U.T = min(Z.T.n.vec)*sum(p*Q.T.p.vec)
#   U.S = min(Z.S.n.vec)*sum(p*Q.S.p.vec)
#   score = c(U.T,U.S)
#   return(score)
# }
