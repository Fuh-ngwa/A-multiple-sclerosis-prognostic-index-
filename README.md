# A-multiple-sclerosis-prognostic-index-
Our inability to reliably predict disease outcomes in multiple sclerosis remains an issue for clinicians and clinical trialists. This study aims to create, from available clinical, genetic and environmental factors; a clinical–environmental–genotypic prognostic index to predict the probability of new relapses and disability worsening. The analyses cohort included prospectively assessed multiple sclerosis cases (N = 253) with 2858 repeated observations measured over 10 years. N = 219 had been diagnosed as relapsing-onset, while N = 34 remained as clinically isolated syndrome by the 10th-year review.

These are the R-codes for the results presented in the paper: "Developing a clinical–environmental–genotypic prognostic index for relapsing-onset multiple sclerosis and clinically isolated syndrome"
This paper is available at https://doi.org/10.1093/braincomms/fcab288

## Part I: Analysis of worsening of disability (WoD)

1-Compute the CEPI, GPI, and CEGPI using the R-script "1-Analysis of WoD_Develope the prognostic index.R"

2-Perform dynamic landmark prediction on the prognostic indices using the R-script "2-Analysis of WoD_Dynamic landmark predicitons.R"

## Part II: Analysis of recurrent relapsing events (RRE)

3-Compute the CEPI, GPI, and CEGPI using the R-script "1-Analysis of RRE_develope the prognostic index.R"

4-Perform dynamic landmark prediction on the prognostic indices using the R-script "2-Analysis of RRE_Dynamic landmark predicitons.R"

## Part III: Analysis of relapses and/or worsening of disability (RWoD)

-To perform this analysis, combined the EDSS status (0,1) with relapse status (0,1) and repeat the analysis in Part I and II.
