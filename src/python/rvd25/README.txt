RVD2.5 is an extension of RVD that includes a prior on the position-specific parameters.

The model is
The position is j = 1,...,J. The replicate is i = 1,...,N. The model only addresses error/reference positions and does not model individual nucleotide frequencies.


\mu_j | \mu_0, \sigma_0 ~ Gaussian(\mu_0, \sigma_0)
\theta_{ij} | \mu_j, M ~ Beta(\mu_j, M)
r_{ij} | \theta_{ij}, n_{ij} ~ Binomial(\theta_{ij}, n_{ij})

where n_{ij} is the total counts at position j in replicate i. 

The error read count at position j in replicate i is modeled by the binomial random variable r_{ij}. The probability of an error at position j in replicate i is \theta_{ij}. The error probability has a prior beta distribution with position-specific rate parameter \mu_j and precision M. The precision is constant across all position, but the prior error rate is position-specific. The position error rate, \mu_j, has as a prior distribution a Log-Normal with parameters \mu_0 and \sigma_0. This is to ensure that that error rate is non-negative. The precision parameter for the Beta distribution has an improper prior.

Initial experiments suggest that this model works well is the read depth,n, is large enough. But when n is small, the precision parameter, M, blows up. I think this is because it doesn't have a prior constraining it. When n is small, there isn't enough data to keep that parameter in check and the model wants to shrink all of the theta's to the location mean. In RVD2.6 I'm going to make M location-specific and put a Gamma prior on it across locations.

