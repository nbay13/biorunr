# All code from Dr. Harold Pimentel's blog post
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

#' @export count.to.tpm
count.to.tpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

#' @export count.to.fpkm
count.to.fpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

#' @export fpkm.to.tpm
fpkm.to.tpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#' @export count.to.effcounts
count.to.effcounts <- function(counts, len, effLen)
{
    counts * (len / effLen)
}