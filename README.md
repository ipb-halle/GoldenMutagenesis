
<!-- README.md is generated from README.Rmd. Please edit that file -->
GoldenMutagenesis
=================

The Golden Gate cloning technique has been proven to be a highly efficient toolbox for a variety of cloning setups. Based on its modular concept it is particularly suitable for the use in multiple-site mutagenesis approaches. In this technical note we developed a protocol termed Golden Mutagenesis for the rapid, easy, reliable and cheap formation of mutagenesis libraries. One to five positions could be altered in parallel or simultaneously within two days. To facilitate the implementation of this technique, this R-library has been developed for the automated primer design and the graphical evaluation of sequencing results to determine the quality of the library.
This library is currently still under development and will be enhanced by an R shiny web-application.

Installation
------------

You can install this development prerelease version of GoldenMutagenesis from github with:

``` r
# install.packages("devtools")
devtools::install_github("ipb-halle/GoldenMutagenesis", ref="prerelease-dev")
```

Example
-------

To start with this library we recommend to have a look on our vignettes!  
[Point Mutagenesis](https://github.com/ipb-halle/GoldenMutagenesis/blob/master/vignettes/Point_Mutagenesis.md)  
[Multiple Site Directed Mutagenesis - Example1](https://github.com/ipb-halle/GoldenMutagenesis/blob/master/vignettes/MSD.md)   
[Multiple Site Directed Mutagenesis - Example2](https://github.com/ipb-halle/GoldenMutagenesis/blob/master/vignettes/MSD2.md)

Interactive Example
-------

You can try out our interactive user notebook for Multiple Site Directed Mutagenesis: [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesis/binder?filepath=notebooks%2FMSD_USER.ipynb)
