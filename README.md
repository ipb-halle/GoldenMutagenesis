README
================
Chris Ulpinnis & Pascal Püllmann
2020-03-05

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build Status](https://travis-ci.org/ipb-halle/GoldenMutagenesis.svg?branch=master)](https://travis-ci.org/ipb-halle/GoldenMutagenesis)
[![GitHub license](https://img.shields.io/github/license/ipb-halle/GoldenMutagenesis)](https://github.com/ipb-halle/GoldenMutagenesis/blob/master/LICENSE)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ipb-halle/GoldenMutagenesis)](https://github.com/ipb-halle/GoldenMutagenesis/releases)
[![GitHub issues](https://img.shields.io/github/issues/ipb-halle/GoldenMutagenesis)](https://github.com/ipb-halle/GoldenMutagenesis/issues)

# GoldenMutagenesis
<img src=https://raw.githubusercontent.com/ipb-halle/GoldenMutagenesis/master/inst/gm_biocsticker.png height="150px"></img>

## Abstract
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

The publication is available at https://rdcu.be/bMfta. 

## Web-interface
A web-based user-interface is available at https://msbi.ipb-halle.de/GoldenMutagenesisWeb.

The GitHub page of the webinterface can be found at https://github.com/ipb-halle/GoldenMutagenesisWeb.

A docker container containing the web interface can be found at https://hub.docker.com/r/sneumann/goldenmutagenesisweb.

## Installation

You can install GoldenMutagenesis from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ipb-halle/GoldenMutagenesis")
```

## Documentation

The documentation can be accessed at https://ipb-halle.github.io/GoldenMutagenesis/.

## Examples

### Vignettes
To start with this library we recommend to have a look at our
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

### Jupiter Notebooks at binder

You can also use predefined Jupiter notebooks on binder to use the library.

#### Multiple Site Directed Mutagenesis
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesis/binder?filepath=notebooks%2FMSD_USER.ipynb)

#### Single Point Mutagenesis
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesis/binder?filepath=notebooks%2FSPM_USER.ipynb)

#### Domestication only
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ipb-halle/GoldenMutagenesis/binder?filepath=notebooks%2FDomesticate_only_USER.ipynb)
