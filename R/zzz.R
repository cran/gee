if(R.version$major == 0 && R.version$minor < 99) {
### need to know the variance function
quasi <-
function (link = "identity", variance = "constant")
{
    linktemp <- substitute(link)
    ##this is a function used in  glm()
    ##it holds everything personal to the family
    ##converts link into character string
    if (is.expression(linktemp))
        linktemp <- eval(linktemp)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link")
            linktemp <- eval(link)
    }
    stats <- make.link(linktemp)
    ##converts variance into character string
    variancetemp <- substitute(variance)
    if (!is.character(variancetemp)) {
        variancetemp <- deparse(variancetemp)
        if (linktemp == "variance")
            variancetemp <- eval(variance)
    }
    # end switch(.)
    switch(variancetemp, "constant" = {
        variance <- function(mu) rep(1, length(mu))
        dev.resids <- function(y, mu, wt) wt * ((y - mu)^2)
        validmu <- function(mu) TRUE
    }, "mu(1-mu)" = {
        variance <- function(mu) mu * (1 - mu)
        validmu <- function(mu) all(mu > 0) && all(mu < 1)
        dev.resids <- function(y, mu, wt) 2 * wt * (y * log(ifelse(y ==
            0, 1, y/mu)) + (1 - y) * log(ifelse(y == 1, 1, (1 -
            y)/(1 - mu))))
    }, "mu" = {
        variance <- function(mu) mu
        validmu <- function(mu) all(mu > 0)
        dev.resids <- function(y, mu, wt) 2 * wt * (y * log(ifelse(y ==
            0, 1, y/mu)) - (y - mu))
    }, "mu^2" = {
        variance <- function(mu) mu^2
        validmu <- function(mu) all(mu != 0)
        dev.resids <- function(y, mu, wt) pmax(-2 * wt * (log(y/mu) -
            (y - mu)/mu), 0)
    }, "mu^3" = {
        variance <- function(mu) mu^3
        validmu <- function(mu) all(mu > 0)
        dev.resids <- function(y, mu, wt) wt * ((y - mu)^2)/(y *
            mu^2)
    }, stop(paste(variancetemp, "not recognised, possible variances",
        "are \"mu(1-mu)\", \"mu\", \"mu^2\", \"mu^3\" and \"constant\"")))
    initialize <- expression({
        n <- rep(1, nobs)
        mustart <- y
    })
    aic <- function(y, n, mu, wt, dev) NA
    structure(list(family = "quasi", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
        validmu = validmu, valideta = stats$valideta, varfun = variancetemp),
              class = "family")
}
}
