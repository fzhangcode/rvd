% SRC
%
% Files
%   bb_loglik_comp        - Beta-Binomial complete-data log-likelihood
%   beta_bino_em          - Expectation-Maximization Algorithm for Beta-Binomial
%   combineRawReads       - Combines forward/reverse and pairs while filtering
%   compileRawDepth       - Converts raw depths to reference and error depth counts.
%   createDepthChart      - Creates depth chart from BAM file and filters low quality reads
%   exportCallTable       - Filters calls by resolution threshold and exports call tables for BB model
%   filter_strand_by_diff - Filter strand reads based on fwd/rev differential
%   filterStrandError     - Filter strands based on error rate differential
%   findsecondbase        - Identify the second most frequent base by read depth
%   fitBB                 - Estimate the Beta-Binomial parameters using reference 
%   getRVDConfig          - Load configuration file for RVD
%   importSampleInfo      - Import and save sample information for dataset
%   importSequenceData    - Reads sequence data from BAM files to data matrix for analysis 
%   RVD                   - Process depth charts to variant calls using the RVD algorithm
%   summary_var           - Outputs a summary table of distribution statistics of samples
%   testbb                - Normal test for Beta-Binomial Model
%   testBBAll             - Test for high error rates compared to null model for all samples

