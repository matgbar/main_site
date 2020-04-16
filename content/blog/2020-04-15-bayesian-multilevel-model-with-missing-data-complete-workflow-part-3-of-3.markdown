---
title: Bayesian Multilevel Model with Missing Data Complete Workflow (Part 3 of 3)
author: ''
date: '2020-04-15'
slug: bayesian-multilevel-model-with-missing-data-complete-workflow-part-3-of-3
categories:
  - Modeling
tags:
  - Bayesian
  - Missing Data
  - Stan
---



# Overview:
This is the third post in a three-part blog series I put together. You can find part 1 [here](../../../../../blog/2018/07/05/2018-07-05-bayesian-multilevel-model-with-missing-data-complete-work-flow/) and part 2 is [here](../../../../../blog/2018/07/05/2018-07-05-bayesian-multilevel-model-with-missing-data-complete-work-flow/). Part 1 provies a quick overview of exploring missing data properties and sets the table for creating a simple missingness model. Part 2, covers the creation of an imputation model and explores the output of that model. This third post will walkthrough running the analysis using `brms` (CITE). 

Even for a toy problem, because I am using Bayesian estimation for the model parameters, the size of the output can get pretty large. Multiplied across 10 or 20 imputed datasets and you can see how the file sizes start to grow. In practice, I tend to save these outputs in a separate external hard drive, and load them when I need to do more extensive analyses.

## Setting up the model: 
Here will run a simple model in which we predict the association between dispostional negativity and momentary positive affect, assuming a random intercept for individual's momentary positive mood ratings. For more context about the data and study design review the previous posts (linked above).

First the basic model assuming a Gaussian prior for the outcome distribution. 

```r
PosAff_DN_form <- 
  brmsformula(
    PosAff ~ 1 + c.DN + c.Best + (1 + c.Best|ID)
  ) + gaussian()
```

If you are familiar with `lme4` syntax (CITE + and link), this should model formulation should be fairly recognizable. We can translate this back to the underlying mathematical formulation. But first a quick conceptual overview of the model. At the between-subjects level, we have an effect that captures the association between individual dispostional negativity scores and average levels of momentary positive effect. At the momentary level, there is a component of the model that characterizes the association between ratings of recent positive events and momentary positive mood. 

Now for the model. 

## Within-Subjects Portion of the Model

$$
PosAff_{ti} \sim Normal(\pi_{0i} + \pi_{1i}(PosEvnt_{ti}-\bar{PosEvnt_{i}}), \sigma_{ti})
$$
## Between-Subjects Portion of the Model

$$
`\begin{equation}
  \begin{bmatrix}
    \pi_{0i} \\
    \pi_{1i}
  \end{bmatrix}
\sim Normal \Big(
  \begin{bmatrix}
    \beta_{00} + \beta_{01}(DN_i) \\
    \beta_{10}
  \end{bmatrix}
, \Sigma \Big)
\end{equation}`
$$
where the `\(\beta\)`s represent the "fixed" or population effects in the model. 

## References
