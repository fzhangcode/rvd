RVD2.4 is an extension of RVD that includes a prior on the position-specific parameters.

The model is
The position is j = 1,...,M. The replicate is i = 1,...,N. The nucleotide is k={A,C,T,G}.

This model is similar to RVD2.3 except it has nucleotide-specfic Gamma parameters. This allows the model to handle the baseline distribution for reference and error bases independently.

\alpha_j^k | a_k,b_k ~ Gamma(a_k,b_k)
\theta_{ij} | \alpha_j ~ Dirichlet(\alpha_j)
r_{ij} | \theta_{ij}, n_{ij} ~ Multinomial(\theta_{ij}, n_{ij})

where n_{ij} is the total counts at position j in replicate i. 