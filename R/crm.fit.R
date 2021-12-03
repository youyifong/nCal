# fit.4pl=FALSE; var.model="power"; robust="mean"; method="drm"; max.iter=20; reltol=1e-3; verbose=TRUE; log.both.sides=F; gof.threshold=0.2
crm.fit=function (formula, data, fit.4pl=FALSE, var.model=c("constant","power"),
    robust="mean", method=c("gls-pl","gnls","mle"), max.iter=50, reltol=1e-3, 
    gof.threshold=0.2, log.both.sides=FALSE, verbose=FALSE) 
{
    
    var.model=match.arg(var.model)
    method=match.arg(method)
    
    # if there is well_role, make sure only standard is used
    if ("well_role" %in% names(data)) {
        data=data[data$well_role %in% c("Standard","standard"),]
        if (nrow(data)==0) stop("there is not standard data")
    }
    
    # reorder, b/c drm 5pl is sensitive to the ordering of rows
    predictor.coln=all.vars(formula)[2]    
    data=data[order(data[[predictor.coln]]),] 
    
    if (log.both.sides) {
        if (var.model=="power") stop("log both sides only work with constant variance")    
        # log-both-sides is implemented by bcVal=0
        fit = drm.fit (formula, data, fit.4pl=fit.4pl, robust=robust, verbose=verbose, bcVal = 0, gof.threshold=gof.threshold)
        fit$sigmasq.hat = getVarComponent.drc(fit) 
        fit$log.both.sides=TRUE       
    
    } else if (var.model=="constant") {    
        fit = drm.fit (formula, data, fit.4pl=fit.4pl, robust=robust, verbose=verbose, gof.threshold=gof.threshold)
        fit$sigmasq.hat = getVarComponent.drc(fit)        
    
    } else if (var.model=="power" & method=="gnls") {    
        fit = gnls.fit (formula, data, fit.4pl=fit.4pl, verbose=verbose)
    
    } else if (var.model=="power" & method=="gls-pl") {
        fit = glspl (formula, data, max.iter=max.iter, reltol=reltol, fit.4pl=fit.4pl, gof.threshold=gof.threshold, verbose=verbose)
        
    } else if (var.model=="power" & method=="mle") {
        fit.tmp = glspl (formula, data, max.iter=1, reltol=reltol, fit.4pl=fit.4pl, gof.threshold=gof.threshold, verbose=verbose)
        start=cla2gh(coef(fit.tmp))
        start=c(start, g.hat=fit.tmp$var.power, sigmasq.hat=fit.tmp$sigmasq.hat)
        fit = mle.dr(formula, data, start=start, max.iter=max.iter, reltol=reltol, fit.4pl=fit.4pl, verbose=T)
    }
    
    class(fit)=c("crm",class(fit))
    return (fit)
    
}


deviance.crm=function(object, ...){
    if ("drc" %in% class(object)) {
        fitte=fitted(object)
        gamma=object$var.power
        if (is.null(gamma)) gamma=0
        y=object$dataList$resp # diff between origResp and resp is caused by bcVal not log() in the formula
        sigmasq.hat=object$sigmasq.hat
        resid.=resid(object)
        sum( gamma*log(fitte) + log(sigmasq.hat) + resid.^2/(fitte^gamma * sigmasq.hat) )
    
    } else {
        return (NA)
    } 
}

print.crm=function(x,..., digits=3) {

    ###################################################################################
    # has to copy print.drc here b/c I can't figure out how to use callNextMethod()
    object <- x
    classList <- class(object)
    cat(paste("\n", "A 'drc' model.", "\n", sep = ""))
    cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
    if (length(coef(object)) > 0) {
        cat("Coefficients:\n")
        print.default(format(coef(object), digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else {
        cat("No coefficients\n")
    }
    cat("\n")
    ###################################################################################
    
    cat("crm parameters.\n")
    cat("Sigma squared: ", signif(x$sigmasq.hat,2), "\n", sep="")
    if (!is.null(x$var.power)) cat("Gamma: ", signif(x$var.power,2), "\n", sep="")
}

lines.crm=function(x, ...) {
    lines5PL (coef(x), xlim=range(x$dataList$dose), ...)
}

coef.crm=function(object, parameterization=c("cla","gh","ed50b","edb50"), ...) {
    p=match.arg(parameterization)
    out=object$coefficients
    if (p=="cla") {
        out
    } else if (p=="gh") {
        cla2gh(out)        
    } else if (p=="ed50b") {
        cla2ed50b(out)        
    } else if (p=="ed50") {
        cla2ed50(out)        
    }
}

## plot.crm is difficult to write because plot.drc is not exported and it is legal to use ::: in a package. callNextMethod does not seem to pass the argument. 
#plot.crm=function(x, ...) {
#    
#    callNextMethod()
#}
