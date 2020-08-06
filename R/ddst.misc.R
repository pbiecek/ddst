`ddst.IIC` <-
function(coord, n, c=2.4) {
 IIC = as.numeric(max(coord[1],diff(coord)) < c*log(n))
 which.max(coord - (IIC*log(n) + (1-IIC)*2)*(1:length(coord)))
}

`ddst.base.cos` <-
  function(x, j) {
    if (j==0)
      return(1)
    sqrt(2) * cos(pi*j*x)
  }

`ddst.phi` <-
  function(x, j, base = ddst.base.legendre) {
    mean(base(x,j))
  }

`ddst.base.legendre` <-
  function(x, j) {
    ddst.polynomial.fun[[j+1]](x)
  }

