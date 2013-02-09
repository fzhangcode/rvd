RVD2.3 is an extension of RVD that includes a prior on the position-specific parameters.

The model is
The position is j = 1,...,M. The replicate is i = 1,...,N. The nucleotide is k={A,C,T,G}.

\alpha_j | a,b ~ Gamma(a,b)
\theta_{ij} | \alpha_j ~ Dirichlet(\alpha_j)
r_{ij} | \theta_{ij}, n_{ij} ~ Multinomial(\theta_{ij}, n_{ij})

where n_{ij} is the total counts at position j in replicate i. 