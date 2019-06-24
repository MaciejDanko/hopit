[![Travis-CI Build Status](https://travis-ci.org/MaciejDanko/hopit.svg?branch=master)](https://travis-ci.org/MaciejDanko/hopit)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/MaciejDanko/hopit?branch=master&svg=true)](https://ci.appveyor.com/project/MaciejDanko/hopit)
[![status](http://joss.theoj.org/papers/73b926670df79a6cfa48cffa7d4775a6/status.svg)](http://joss.theoj.org/papers/73b926670df79a6cfa48cffa7d4775a6)
[![CRAN_Download_Badge1](https://cranlogs.r-pkg.org/badges/grand-total/hopit)](https://CRAN.R-project.org/package=hopit)
[![CRAN_Download_Badge2](https://cranlogs.r-pkg.org/badges/hopit)](https://CRAN.R-project.org/package=hopit)
[![license](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://github.com/MaciejDanko/hopit/blob/master/LICENSE)

# R-package *hopit*: Hierarchical ordered probit models with application to reporting heterogeneity.

Self-reported health, happiness, attitudes, and other statuses or perceptions are often the subject of biases that may come from different sources. For example, the evaluation of an individual’s own health may depend on previous medical diagnoses, functional status, and symptoms and signs of illness; as on well as life-style behaviors, including contextual social, gender, age-specific, linguistic and other cultural factors (Oksuzyan et al. 2019). The **hopit** package offers versatile functions for analyzing different self-reported ordinal variables, and for helping to estimate their biases. Specifically, the package provides the function to fit a generalized ordered probit model that regresses original self-reported status measures on two sets of independent variables (King et al. 2004; Jurges 2007; Oksuzyan et al. 2019). included in the regression are individual statuses and characteristics that are directly related to the self-reported variable. In the case of self-reported health, these could be chronic conditions, mobility level, difficulties with daily activities, performance on grip strength tests, anthropometric measures, and lifestyle behaviors. The second set of independent variables (threshold variables) is used to model cut-points between adjacent self-reported response categories as functions of individual characteristics, such as gender, age group, education, and country (Oksuzyan et al. 2019). The model helps to adjust for specific socio-demographic and cultural differences in how the continuous latent health is projected onto the ordinal self-rated measure. The fitted model can be used to calculate an individual predicted latent status variable, a latent index, and standardized latent coefficients; and makes it possible to reclassify a categorical status measure that has been adjusted for inter-individual differences in reporting behavior.

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
Get started with `hopit` by checking the [vignette](https://github.com/MaciejDanko/hopit/blob/master/vignettes/vig_hopit.pdf) or run:

   ```R
   browseVignettes(package = "hopit") 
   ```

### Contributing
This software is an academic project. Any issues and pull requests are welcome.
* If `hopit` is malfunctioning, please report the case by submitting an issue on GitHub.

### References
King, GCh, Murray JL, Salomon JA, and Tandon A. (2004). “Enhancing the Validity and Cross-Cultural Comparability of Measurement in Survey Research.” American Political Science Review 98 (1). Cambridge University Press: 191–207. doi: 10.1017/S000305540400108X.

Jurges H (2007). “True health vs response styles: exploring cross-country differences in self-reported health.” Health Economics, 16(2), pp. 163-178. doi: 10.1002/hec.1134.

Oksuzyan A, Danko MJ, Caputo J, Jasilionis D, and Shkolnikov VM (2019). “Is the story about sensitive women and stoical men true? Gender differences in health after adjustment for reporting behavior.” Social Science & Medicine, 228, pp. 41-50. doi: 10.1016/j.socscimed.2019.03.002. 
