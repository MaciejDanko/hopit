---
output: github_document
---

[![Travis-CI Build Status](https://travis-ci.org/MaciejDanko/hopit.svg?branch=master)](https://travis-ci.org/MaciejDanko/hopit)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/MaciejDanko/hopit?branch=master&svg=true)](https://ci.appveyor.com/project/MaciejDanko/hopit)

# Analysis of the reporting styles using generalized ordered probit models

## Installation

1. Make sure you have the most recent version of R
2. Run the following code in your R console 

   ```R
   install.packages("hopit")
   ```

## Updating to the latest version of `hopit` package

You can track (and contribute to) the development of `hopit` at https://github.com/MaciejDanko/hopit. To install it:

1. Install the release version of `devtools` from CRAN with `install.packages("devtools")`.

2. Make sure you have a working development environment.
    * **Windows**: Install [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/).
    * **Mac**: Install `Xcode` from the Mac App Store.
    * **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).

3. Install the development version of `hopit`:

   ```R
   devtools::install_github("MaciejDanko/hopit")
   ```
   or optionally

   ```R
   devtools::install_github("MaciejDanko/hopit", build_vignettes=TRUE)
   ```

## Intro
Get started with `hopit` by checking the [vignette](https://github.com/MaciejDanko/hopit/blob/master/vignettes/introduction_to_hopit.pdf) 
 ```R
 browseVignettes(package = "hopit") 
 ```

## Contributing
This software is an academic project. We welcome any issues and pull requests.
* If `hopit` is malfunctioning, please report the case by submitting an issue on GitHub.
