---
title: '`hopit`: An R Package for Analysis of Reporting Behavior Using Generalized Ordered Probit Models'
tags:
- self-reported status 
- reporting heterogeneity
- reporting behavior
- latent index
- latent health
- disability weights
authors:
- name: Maciej J. Dańko
  affiliation: 1
  orcid: 0000-0002-7924-9022
affiliations:
- name: Max Planck Institute for Demographic Research, Rostock, Germany
  index: 1
date: "3 May 2019"
bibliography: paper.bib
link-citations: true
# header-includes: \usepackage{float}
---
# Summary
``hopit`` is an open source software library written in the R [@Rteam2019] and C++ [@RcppEigen2013; @Rcpp2011] programming languages. The ``hopit`` package provides versatile methods for fitting and analyzing ordered response data in the context of heterogeneity in self reporting behavior.

The ordered response data classify a measure of interest into ordered categories collected during a survey. For example, if the dependent variable is a happiness rating, then a respondent typically answers a question such as: “Taking all things together, would you say you are ... ?$"$; and then selects from response options such as: "very happy", "pretty happy", "not too happy", and "very unhappy" [@Liao2005]. Similarly, if interviewees are asked to evaluate their health in general (e.g., “Would you say your health is ... ?”), they may choose among several categories, such as "very good", "good", "fair", "bad", and "very bad" [@King2004; @Jurges2007; @Rebelo2014; @OKSUZYAN2019]. In political science, a respondent may be asked for an opinion about recent legislation (e.g. “Rate your feelings about the proposed legislation.") and asked to choose among categories like "strongly oppose", "mildly oppose", "indifferent", "mildly support", and "strongly support" [@GreeneHensher2010]. It is easy to imagine other multi-level ordinal variables that might be used during a survey and to which the methodology described below could be applied.

In practice, it is assumed that when responding to a survey question about their general happiness, health, feelings, attitudes or other status, participants are assessing their true value of this unobserved continuous variable, and project it onto the discrete scale provided. The thresholds that individuals use to categorize their true status by selecting a specific response option may be affected by the reference group chosen, their earlier life experiences, and cross-cultural differences in using scales. Thus, the responses of individuals may differ depending on their gender, age, cultural background, education, and personality traits; among other factors.

From the perspective of reporting behavior modeling, one of the main tasks researchers face is to compute this continuous estimate of the underlying, latent measures of individuals based on several specific characteristics of the responses considered (e.g., health variables or happiness variables), and to account for variations in reporting across socio-demographic and cultural groups. More specifically, to build a latent, underlying measure, a generalized hierarchical ordered threshold model is fitted that regresses the reported status/attitude/feeling on two sets of independent variables [@Boes2006; @Green2014]. When the dependent reported ordered variable is self-rated health status, then the first set of variables – i.e., health variables – assess specific aspects of individuals’ health, such as measures of chronic conditions, mobility, difficulties with a range of daily activities, grip strength, anthropometric characteristics, and lifestyle behaviors. Using the second set of independent variables (threshold variables), the model also adjusts for differences across socio-demographic and cultural groups, such as differences in cultural background, gender, age, and education [@King2004; @Jurges2007; @OKSUZYAN2019; but see @Rebelo2014].

The ``hopit`` package delivers functions and methods to fit (```hopit```), summarize (e.g., ```summary```), check (e.g., ```profile```), and compare (e.g., ```AIC``` and ```anova```) fitted models. The latent and threshold formulas are defined separately. The interactions can be specified both within and between these formulas. Depending on how an interactions between latent and threshold variables is interpreted, it can be added to either the latent or the threshold formula. The package has also an option to include a survey design using the ``survey`` package [@Lumley2004; @Lumley2019]. 

Once the model is fitted, the model estimates are used to determine reporting behavior; i.e., how the continuous latent measure is projected onto the categorical response. In practice, this is done by comparing actual categorical ordered responses with theoretical responses that are adjusted for heterogeneity in reporting behaviors, and are, therefore, more comparable across individuals.

<!--
\begin{figure}[H]
{\centering\includegraphics{HI-2.pdf}

}
\caption[Health index vs. self-reported health for a model fitted to an exemplary data.]{Health index vs. self-reported health for a model fitted to an exemplary data.}\label{fig:HI-2}
\end{figure}-->
![Health index vs. self-reported health for a model fitted to an exemplary data.](HI-2.pdf)

One of the first steps of the analysis is the standardization of the latent measure to obtain the latent index. In the self-rated health example the latent health index is a proxy for the true underlying health of an individual, and varies from 0, representing the (model-based) worst health state in the sample, to 1, representing the (model-based) best health state in the sample (Fig. 1).

The predicted latent measure obtained from the model is also used to standardize the latent variable coefficients. In the self-rated health example the standardized coefficients are called disability weights [@Jurges2007; @OKSUZYAN2019] and are calculated for each health variable to provide information about the impact of a specific health measure on the latent index (Fig. 2). The disability weight for a health variable is equal to the ratio of the corresponding health coefficient and the difference between the lowest and the highest values of predicted latent health. In other words, the disability weight reduces health index by some given amount or percentage; i.e. the health index of every individual is reduced by the same amount if the person had a heart attack or other heart problems)[@Jurges2007; @OKSUZYAN2019]. While the latent index is intended to reflect the underlying health, happiness or other status across individuals, the standardized coefficients (disability weights for health case), are computed for an average individual in the study population. 

<!-- 
\begin{figure}[H]
{\centering\includegraphics{D-2.pdf}

}
\caption[Disability weights for the fitted model.]{Disability weights for the fitted model.}\label{fig:D-2}
\end{figure} -->
![Disability weights for the fitted model.](D-2.pdf) 

Reporting behavior analysis is based on the reclassification of individuals into new response categories. There are two methods of reclassification: (1) @Jurges2007 percentile method [see also @Rebelo2014; @OKSUZYAN2019] and (2) reclassification directly based on model-estimated thresholds. In the first method, the classification is based on the calculated latent index which is adjusted for inter-individual differences in reporting behavior. This method is based on the original distribution of the categorical response variable (see also @OKSUZYAN2019). 

The package offers functions ```boot_hopit``` and ```percentile_CI``` for calculating the confidence intervals for any measure derived from the model using parametric bootstrap methods. In each of the bootstrap repetitions, a set of new model coefficients is drawn from the multivariate normal distribution, assuming the originally estimated model coefficients as a mean and using the model estimated variance-covariance matrix. The drawn coefficients are then used to calculate the measure of interest via a user defined function. In the example presented in Fig. 3, the confidence intervals of the difference between the original and the adjusted frequencies of bad health are calculated. The bad health is determined by the presence of "```Poor```" or "```Fair```" self-rated health categories.

The results (Fig.3) show that men tend to over-report bad health at ages (50,60] and (50,70], whereas women tend to over-report bad health at ages [70,80); and that both sexes at ages above 80 tend to under-report bad health. See also @OKSUZYAN2019 for similar analyses done on true SHARE data.

<!-- 
\begin{figure}[H]
{\centering\includegraphics{B-2.pdf}

}
\caption[Differences between original and adjusted prevalances of bad health for the fitted model. The confidence intervals were calculated using percentile bootstrap method.]{Differences between original and adjusted prevalances of bad health for the fitted model. The confidence intervals were calculated using percentile bootstrap method.}\label{fig:B-2}
\end{figure} -->
![Differences between original and adjusted prevalances of bad health for the fitted model. The confidence intervals were calculated using percentile bootstrap method.](B-2.pdf) 

# Acknowledgements

I thank Anna Oksuzyan, Christian Dudel, Marius Pascariu, Laszlo Nemeth, and Oskar Burger for their comments and suggestions. I also thank the Max-Planck Institute for Demographic Research for all their support.

In all examples presented above I use ```healthsurvey```, which is a completely artificial data set that is simulated using the distributions of some major health and socio-demographic characteristics. The distributions and the data structure is roughly based on the WAVE1 SHARE database (DOIs: 10.6103/SHARE.w1.600) see @Borsch2013 for technical details about SHARE database. None of the records represent any part of the true data. 

The SHARE data collection has been primarily funded by the European Commission through FP5 (QLK6-CT-2001-00360), FP6 (SHARE-I3: RII-CT-2006-062193, COMPARE: CIT5-CT-2005-028857, SHARELIFE: CIT4-CT-2006-028812) and FP7 (SHARE-PREP: N°211909, SHARE-LEAP: N°227822, SHARE M4: N°261982). Additional funding from the German Ministry of Education and Research, the Max Planck Society for the Advancement of Science, the U.S. National Institute on Aging (U01_AG09740-13S2, P01_AG005842, P01_AG08291, P30_AG12815, R21_AG025169, Y1-AG-4553-01, IAG_BSR06-11, OGHA_04-064, HHSN271-201300071C) and from various national funding sources is gratefully acknowledged (see www. share-project.org).

# References

