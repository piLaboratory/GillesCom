# GillesCom
[![Travis-CI Build Status](https://travis-ci.org/piLaboratory/GillesCom.svg?branch=master)](https://travis-ci.org/piLaboratory/GillesCom)

## An R package for Gillespie Simulation of Ecological Communities Dynamics
### Paulo I. Prado & Andre Chalom

This package is being developed to provide simulation tools for generating ecological community data (such as 
  species abundance distributions) using the Gillespie continuous time method.

### Installation
GillesCom is currently only supported in Linux environments.
```R
library(devtools)
install_github("piLaboratory/GillesCom", build_vignettes = TRUE)
```

See the examples on the package vignette with

```R
library(GillesCom)
vignette("GillesCom")
```
