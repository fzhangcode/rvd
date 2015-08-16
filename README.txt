RVD2.7 is an extension of RVD that includes a prior on the position-specific parameters.

The model is
The position is j = 1,...,J. The replicate is i = 1,...,N. The model only addresses error/reference positions and does not model individual nucleotide frequencies.

\mu_j | \mu_0, M_0 ~ Beta(\mu_0, M_0)
\theta_{ij} | \mu_j, M_j ~ Beta(\mu_j, M_j)
r_{ij} | \theta_{ij}, n_{ij} ~ Binomial(\theta_{ij}, n_{ij})

where n_{ij} is the total counts at position j in replicate i. 

The error read count at position j in replicate i is modeled by the binomial random variable r_{ij}. The probability of an error at position j in replicate i is \theta_{ij}. The error probability has a prior beta distribution with position-specific rate parameter \mu_j and precision M_j. The position error rate, \mu_j, has as a Beta prior distribution  with parameters \mu_0 and M_0. This is to ensure that that error rate is between 0 and 1. The precision parameter M_j has an improper prior. This is useful for situations when there is a significant minor allele.


