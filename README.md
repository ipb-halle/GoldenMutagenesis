README
================
Chris Ulpinnis & Pascal Püllmann
2018-10-11

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build Status](https://travis-ci.org/ipb-halle/GoldenMutagenesis.svg?branch=master)](https://travis-ci.org/ipb-halle/GoldenMutagenesis)

<img src=https://raw.githubusercontent.com/ipb-halle/GoldenMutagenesis/master/inst/gm_biocsticker.png height="150px"></img>

# GoldenMutagenesis

The Golden Gate cloning technique has been proven to be a highly
efficient toolbox for a variety of cloning setups. Based on its modular
concept it is particularly suitable for the use in multiple-site
mutagenesis approaches. In this technical note we developed a protocol
termed Golden Mutagenesis for the rapid, easy, reliable and cheap
formation of mutagenesis libraries. One to five positions could be
altered in parallel or simultaneously within two days. To facilitate the
implementation of this technique, this R-library has been developed for
the automated primer design and the graphical evaluation of sequencing
results to determine the quality of the library.  
This library is currently still under development and will be enhanced
by an R shiny web-application.

## Installation

You can install GoldenMutagenesis from github with:

``` r
# install.packages("devtools")
devtools::install_github("ipb-halle/GoldenMutagenesis")
```

## Example

To start with this library we recommend to have a look on our
vignettes\!  
[Point
Mutagenesis](https://ipb-halle.github.io/GoldenMutagenesis/articles/Point_Mutagenesis.html)  
[Multiple Site Directed Mutagenesis -
Example1](https://ipb-halle.github.io/GoldenMutagenesis/articles/MSD.html)  
[Multiple Site Directed Mutagenesis -
Example2](https://ipb-halle.github.io/GoldenMutagenesis/articles/MSD2.html)  
[Multiple Site Directed Mutagenesis -
Example3](https://ipb-halle.github.io/GoldenMutagenesis/articles/MSD3.html)  
[Quick Quality
Control](https://ipb-halle.github.io/GoldenMutagenesis/articles/QQC.html)  

## Interactive Example

## Webtool Beta Version
**You can try our webtool at http://msbi.ipb-halle.de/GoldenMutagenesisWeb or on binderhub at [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesisWeb/master?urlpath=shiny).** 


The Webtools git repository is available at: https://github.com/ipb-halle/GoldenMutagenesisWeb 


*Please note that the webtool is still in beta phase. There are still errors and also the user interface is not the final user experience. If you encounter any problems please send an email to me.* 


## You can try out our interactive user notebooks for Mutagenesis:

### Multiple Site Directed

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesis/binder?filepath=notebooks%2FMSD_USER.ipynb)

### Single Point
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesis/binder?filepath=notebooks%2FSPM_USER.ipynb)

### Domestication only
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesis/binder?filepath=notebooks%2FDomesticate_only_USER.ipynb)
