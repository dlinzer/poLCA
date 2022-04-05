# poLCA

### Polytomous Variable Latent Class Analysis

[poLCA][] is a software package for the estimation of latent class models and latent class regression models for polytomous outcome variables, implemented in the [R][] statistical computing environment.

Latent class analysis (also known as latent structure analysis) can be used to identify clusters of similar "types" of individuals or observations from multivariate categorical data, estimating the characteristics of these latent groups, and returning the probability that each observation belongs to each group. These models are also helpful in investigating sources of confounding and nonindependence among a set of categorical variables, as well as for density estimation in cross-classification tables. Typical applications include the analysis of opinion surveys; rater agreement; lifestyle and consumer choice; and other social and behavioral phenomena.

The basic latent class model is a finite mixture model in which the component distributions are assumed to be multi-way cross-classification tables with all variables mutually independent. The model stratifies the observed data by a theoretical latent categorical variable, attempting to eliminate any spurious relationships between the observed variables. The latent class regression model makes it possible for the researcher to further estimate the effects of covariates (or "concomitant" variables) on predicting latent class membership.

poLCA uses expectation-maximization and Newton-Raphson algorithms to find maximum likelihood estimates of the parameters of the latent class and latent class regression models.


## Package authors

[Drew A. Linzer](https://votamatic.org/about-me/)  

[Jeffrey Lewis](https://www.sscnet.ucla.edu/polisci/faculty/lewis/)  


## Installation

To install the package directly through [R][], type

```R
install.packages("poLCA", dependencies = TRUE)
```

and select a CRAN mirror.  Once the installation is complete, enter

```R
library(poLCA)
```

to load the package into memory for use.

poLCA is distributed through the Comprehensive R Archive Network, [CRAN](https://cran.r-project.org).  The compiled package source and MacOS and Windows binary files can be downloaded from https://cran.r-project.org/web/packages/poLCA.

The poLCA package appears in CRAN Task Views for [Cluster Analysis & Finite Mixture Models](https://CRAN.r-project.org?view=Cluster), and [Psychometric Models and Methods](https://CRAN.r-project.org?views=Psychometrics). poLCA is provided free of charge, subject to version 2 of the GPL or any later version. 


## Documentation

[Download user's manual (PDF)](inst/doc/poLCA-manual-1-4.pdf?raw=true). The package is also documented internally upon installation.  For help in [R][], type

```R
?poLCA
```


## Citation

Users of poLCA are requested to cite the software package as:

Linzer, Drew A. and Jeffrey Lewis. 2022. "poLCA: Polytomous Variable Latent Class Analysis." R package version 1.6. https://dlinzer.github.com/poLCA.

and

Linzer, Drew A. and Jeffrey Lewis. 2011. "poLCA: an R Package for Polytomous Variable Latent Class Analysis." Journal of Statistical Software. 42(10): 1-29. https://www.jstatsoft.org/v42/i10


## Contact 

Please direct all inquiries, comments, and reports of bugs to drew@votamatic.org.


[poLCA]: https://dlinzer.github.io/poLCA/
[R]: https://cran.r-project.org
