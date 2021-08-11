# IITA 2021 Genomic Prediction and Mate Selection

[**Start here**](https://wolfemd.github.io/IITA_2021GS/)

[**Jump to the Results here!**](https://wolfemd.github.io/IITA_2021GS/07-Results.html)

-   I recently (Aug 3, 2021) completed thorough testing to develop a protocol for genomic *mate* selection in NextGen Cassava and will implement it here for the first time in practice! Check [here](https://wolfemd.github.io/implementGMSinCassava/) to see the implementation tests / code-base development documented.
-   Subsequently, I have developed my entire code base into an R package [genomicMateSelectR](https://wolfemd.github.io/genomicMateSelectR/index.html). The package is fully documented (in rough draft version), but doesn't yet include a tutorial / vignette.
-   The IITA DArTseqLD report (DCas21_6038), which contains the GS C5 (i.e. progeny of crosses made in 2020, "TMS20F"), was recently (July 19, 2021) released.
-   The imputations of DCas21_6038, genomic predictions and mate selection that follows will leverage [genomicMateSelectR](https://wolfemd.github.io/genomicMateSelectR/index.html) functions and will form the principal example of their use. Install it e.g. `devtools::install_github("wolfemd/genomicMateSelectR", ref = 'master')` .

