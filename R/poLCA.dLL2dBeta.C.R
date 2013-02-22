poLCA.dLL2dBeta.C <-
function(rgivy,prior,x) {
    classes <- dim(prior)[2]
    numx <- dim(x)[2]
    ret <-  .C("d2lldbeta2",
                as.double(t(rgivy)),
                as.double(t(prior)),
                as.double(t(x)),
                as.integer(dim(x)[1]),
                as.integer(classes),
                as.integer(numx),
                grad = double((classes-1)*numx),
                hess = double(((classes-1)*numx)^2)                
            )
    return(list(grad=ret$grad,hess=-matrix(ret$hess,ncol=((classes-1)*numx),byrow=TRUE)))
}
