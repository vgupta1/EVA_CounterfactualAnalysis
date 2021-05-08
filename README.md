# Off-Policy and Counterfactual Analysis

In the interest of reproducibility of research, this repository provides all code necessary to reproduce the off-policy evaluation and counterfactual analysis in the paper [Deploying an Artificial Intelligence System for COVID-19 Testing at the Greek Border](https://dx.doi.org/10.2139/ssrn.3789038)  

If you find it hepful please consider citing:
```
@article{bastani2021deploying,
  title={Deploying an Artificial Intelligence System for Covid-19 Testing at the Greek Border},
  author={Bastani, Hamsa and Drakopoulos, Kimon and Gupta, Vishal and Vlachogiannis, Jon and Hadjicristodoulou, Christos and Lagiou, Pagona and Magiorkinis, Gkikas and Paraskevis, Dimitrios and Tsiodras, Sotirios}, 
  year={2021},
  month={Feb.},
  url={http://dx.doi.org/10.2139/ssrn.3789038}, 
  journal={SSRN}
}
```

For code related to the actual EVA system and how it allocated tests, see our [partner repository](https://github.com/vgupta1/EvaTargetedCovid19Testing).



## Data Availability

:warning: **All data files in this repository are realistic but SYNTHETIC data.** :warning: 

They do not represent actual prevalence rates in Summer 2020.  

The actual data that support the findings of the above paper are available from the Ministry of
Digital Governance but restrictions apply to the availability of these data, which were used under license are not publicly available. Data are however available from the authors upon reasonable request and with permission from the Ministry of Digital Governance.

## Repository Structure
### Synthetic Data Files
* sample_data.csv
* sample_test_results.csv
>The first gives the data format and structure of the PLF database and which (anonymized) passengers arrived at each point of >time to each port of entry.  The second file contains allocation decisions by our bandit algorithm.

### OtherData
This folder contains non-confidential static data used by the bandit allocation algorithm and hence the off-policy evaluation.

### Source Code
* helpers.r
* mse_helpers.r
>These two files contain small helper functions used across all scripts for core computations like fitting and upating our ?>empirical bayes model.

* genHistoricalEB_TS.R
>This file uses all of the available data to create a time series of estimated prevalence for each type.  In contrast to the >bandit allocation, this analysis is in hindsight, and hence free to use data from time periods after time t in order to >estimate prevalence at time t.  In particular i) to estimate prevalence at time t, we use data from [t-3, t+3] days, and ii) we always define types at the level of (country, color_designation) (i.e. without cities/states), for consistency.  

* smoothingEBPriors.R
>The outputs of the moment-matching with only 7 days of data [t-3, t+3] is a bit choppy.  This file performs a cubic spline to >smooth the prior estimates (over time) and then updates the appropriate priors.  

* buildOPEdatabase.R
>Does the lion's share of work for the off-policy-evaluation to assess the value of EVA's targeting relatie to random >surveillance testing.  Outputs a complete analysis of estimates by day and type for random surveillance and EVA.  These can >be fed into other scripts to produce plots/summary analysis etc.  

### Outputs

* grey_arrivs.csv
* grey_eb_preds.csv
>These two file contain our estimates of i) the number of arrivals from and ii) the Covid-19 prevalence for every country as a >daily time series.  These estimates are computed under the counterfactual assumption that we had *not* grey-listed the >country. 

##### OPE Outputs 
This folder contains outputs generated in the course of the off-policy analysis. Specifically,  
* hist_eb_timeseries_TRUE_Window_3_MinTest_30_UnWeighted.csv
> Contains hindisight, empirial bayes estimates of prevalence for each country

* hist_eb_timeseries_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv
> Contains updated estimates obtained by smoothing the priors in the previous estimates using cubic splines.  We find these estimates to have slightly better out-of-sample mse and hence prefer them in our analysis.

* ope_dat_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv
> Contains the daily estimated number of infections caught by our policy and random-surveillance, and random-surveillance without the benefit of our grey-listing procedure. 


## Executing the Code
To execute the code, one should follow the following steps:
1. Run the "genHistoricalEB_TS.R" to generate the time-series of estimated prevalence for each country. 
2. Run the "smoothingEBPriors.R" to smooth these estimates for improved performance.
3. @Hamsa.  pelase fill in what we do next to get the grey counterfactuals.
4. Run buildOPEdatabase.R to generate all off-policy evaluationa nd counterfactual analysis.  Summary statistics and plots can easily be created from the resulting dumped files.  

