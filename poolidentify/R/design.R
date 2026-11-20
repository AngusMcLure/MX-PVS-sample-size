# Sharp limits of the independence fit given the mean and individual ICC.
naive_bias_bounds <- function(prevalence, icc, pool_size) {
    args <- list(prevalence, icc, pool_size)
    if (!all(vapply(args, function(x)
        is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))) ||
        prevalence < 0 || prevalence > 1 || icc < 0 || icc > 1 ||
        pool_size < 1 || pool_size != floor(pool_size))
        stop("Specify scalar prevalence and ICC in [0,1] and a positive integer pool size.")
    if (prevalence == 0 || prevalence == 1 || pool_size == 1 || icc == 0)
        return(c(lower = 0, upper = 0))
    a <- prevalence * (1 - icc)
    b <- prevalence + icc * (1 - prevalence)
    positive <- (prevalence / b) * (-expm1(pool_size * log1p(-b)))
    # At high positivity, subtracting the positive probability from one
    # erases rare negative pools that still carry prevalence information.
    log_negative <- if (positive < 0.5) log1p(-positive) else
        log(icc * (1 - prevalence) / b +
            (prevalence / b) * exp(pool_size * log1p(-b)))
    worst <- -expm1(log_negative / pool_size)
    best <- -expm1((log1p(-prevalence) + (pool_size - 1) * log1p(-a)) / pool_size)
    c(lower = worst - prevalence, upper = best - prevalence)
}
