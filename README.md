[![Travis-CI Build Status](https://travis-ci.org/MaciejDanko/hopit.svg?branch=master)](https://travis-ci.org/MaciejDanko/hopit)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/MaciejDanko/hopit?branch=master&svg=true)](https://ci.appveyor.com/project/MaciejDanko/hopit)

# R-package *hopit*: Hierarchical ordered probit models with application to reporting heterogeneity.

The *hopit* package provides R functions to fit and analyze ordered response data in the context of reporting heterogeneity.

### Installation
1. Make sure you have the most recent version of R

2. Run the following code in your R console 

   ```R
   install.packages("hopit") 
   ```

### Updating to the latest version of `hopit` package
You can track (and contribute to) the development of `hopit` at https://github.com/MaciejDanko/hopit. To install it:

1. Install the release version of `devtools` from CRAN with `install.packages("devtools")`.

2. Make sure you have a working development environment.
    * **Windows**: Install [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/).
    * **Mac**: Install `Xcode` from the Mac App Store.
    * **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).

3. To install the development version of `hopit` run:
   ```R
   devtools::install_github("MaciejDanko/hopit")
   ```
   
### Introduction
Get started with `hopit` by checking the [vignette](https://github.com/MaciejDanko/hopit/blob/master/vignettes/introduction_to_hopit.pdf) or run:

   ```R
   browseVignettes(package = "hopit") 
   ```

### Contributing
This software is an academic project. Any issues and pull requests are welcome.
* If `hopit` is malfunctioning, please report the case by submitting an issue on GitHub.

### References
Jurges H (2007). “True health vs response styles: exploring cross-country differences in self-reported health.” Health Economics, 16(2), pp. 163-178. doi: 10.1002/hec.1134.

Oksuzyan A, Danko MJ, Caputo J, Jasilionis D and Shkolnikov VM (2019). “Is the story about sensitive women and stoical men true? Gender differences in health after adjustment for reporting behavior.” Social Science & Medicine, 228, pp. 41-50. doi: 10.1016/j.socscimed.2019.03.002. 
