# Off-Policy and Counterfactual Analysis

In the interest of reproducibility of research, this repository provides all code necessary to reproduce the off-policy evaluation and counterfactual analysis in the paper [Efficient and Targeted COVID-19 Border Testing via Reinforcement Learning](https://www.nature.com/articles/s41586-021-04014-z).  A previous version of this paper was entitled "Deploying an Artificial Intelligence System for COVID-19 Testing at the Greek Border."  



If you find it hepful please consider citing:
```
@article{bastani2021Efficient,
  title={Efficient and Targeted COVID-19 Border Testing via Reinforcement Learning},
  author={Bastani, Hamsa and Drakopoulos, Kimon and Gupta, Vishal and Vlachogiannis, Jon and Hadjicristodoulou, Christos and Lagiou, Pagona and Magiorkinis, Gkikas and Paraskevis, Dimitrios and Tsiodras, Sotirios}, 
  year={2021},
  month={Sept.},
  url={https://doi.org/10.1038/s41586-021-04014-z}, 
  journal={Nature}
}
```

For code related to the actual EVA system and how it allocated tests, see our [partner repository](https://github.com/vgupta1/EvaTargetedCovid19Testing).

## Data Availability

:warning: **All data files in this repository are realistic but SYNTHETIC data.** :warning: 

They do not represent actual prevalence rates in Summer 2020 but do represent the structure of data used for our analysis.  

Actual, aggregated, anonymized data are available at [this repository](https://github.com/kimondr/EVA_Public_Data). These data aggregate
passenger arrival and testing information over pairs of consecutive days, country of origin, and point of entry. 

The finer granularity data that support the (exact) findings of this study are protected by GDPR. These finer granularity data are available from the Greek Ministry of Civil Protection but restrictions apply to the availability of these data, which were used under license for the current study, and so are not publicly available. Access to these data can only be granted by the Greek Ministry of Civil Protection (info@gscp.gr) for research that is conducted in the public interest for public health (GDPR Recital 159) and scientific purposes (GDPR Article 89). 

Finally, the population-level epidemiological metrics used in our analysis can be obtained freely from the [Our World In Data](https://github.com/owid/covid-19-data/tree/master/public/data).


## Repository Structure
### Synthetic Data Files
* sample_data.csv
* sample_test_results.csv
>The first gives the data format and structure of the PLF database and which (anonymized) passengers arrived at each point of time to each port of entry.  The second file contains allocation decisions by our bandit algorithm.

### OtherData
This folder contains non-confidential static data used by the bandit allocation algorithm and hence the off-policy evaluation.  Many of these data files are generated by files in "src" described below.  The exceptions are:

* first_viable_date_port.csv
>Documents the first date when scanners and border agents were deployed for each port of entry.  Data was not available before this date.

* grey_list_start_end.csv
>Documents when various countries were added/removed from grey-list.  See our paper for details on how these decisions were made.

* short_to_iso.csv
>Dataset to convert 2 letter ISO Codes to 3 Letter ISO Codes.

* testScalings.csv
>As discussed in our paper, some of our benchmark policies based on public data do not fully allocate all tests.  This file shows the fraction of tests allocated.  It is used for scaling when presenting results to allow for "apples-to-apples" comparisons to Eva. 

### Source Code
* helpers.r
* mse_helpers.r
>These two files contain small helper functions used across all scripts for core computations like fitting and updating our empirical bayes model.

* genHistoricalEB_TS.R
>This file uses all of the available data to create a time series of estimated prevalence for each type.  In contrast to the bandit allocation, this analysis is in hindsight, and hence free to use data from time periods after time t in order to estimate prevalence at time t.  In particular i) to estimate prevalence at time t, we use data from [t-3, t+3] days, and ii) we always define types at the level of (country, color_designation) (i.e. without cities/states), for consistency.  

* smoothingEBPriors.R
>The outputs of the moment-matching with only 7 days of data [t-3, t+3] is a bit choppy.  This file performs a cubic spline to smooth the prior estimates (over time) and then updates the appropriate priors.  

* buildOPEdatabase.R
>Does the lion's share of work for the off-policy-evaluation to assess the value of EVA's targeting relative to random surveillance testing.  Outputs a complete analysis of estimates by day and type for random surveillance and EVA.  These can be fed into other scripts to produce plots/summary analysis etc.  

* compare_to_epi.R
>Generates a time series of risk scores for each country based on its epidemiological metrics (cases per capita, deaths per capita, positivity rates). It then generates propensity scores for testing each type of passenger -- had we tested roughly proportionally to these epidemiological metrics -- for each day based on arrivals, testing budgets. It also computes scalings that are recorded manually in testScalings.csv (see paper for details).

* clusteringLags.ipynb
>This is Julia notebook that takes as input the file "country_lag_auc_profile.csv.csv" produced by produce_lag_profiles.py and solves a linear binary optimization problem to cluster countries based on the number of days their public information lags real-time data. 

* grey_ebpred.r
* grey_arrivals.r
> These two files produce counterfactual predictions of the prevalence and arrival rates from greylisted countries had they not been greylisted. They take as input OtherData/Xall.csv (which summarizes publicly reported metrics) and create the files grey_eb_preds.csv and grey_arrivs.csv respectively.

* produce_public_vs_private_raw_data.py
> This file prepares the timeseries of publicly posted epidemiological metrics Xall.csv and Xall_new.csv. It produces for every country and date pair the timeseries of these metrics for the interval [t-20,t+19]. Set the variable "raw" to FALSE to obtain "smoothed" estimates (corresponding to Xall_new.csv) and TRUE to obtain "raw" values (corresponding to Xall.csv).


* publicdata_efficacy.r
> This file tries to predict privately observed prevalence rates using GBM on publicly reported metrics (summarized in OtherData/Xall.csv). It tests 5 models with varying features and reports the resulting ROC curves.




### Outputs

* grey_arrivs.csv
* grey_eb_preds.csv
>These two file contain our estimates of i) the number of arrivals from and ii) the Covid-19 prevalence for every country as a daily time series.  These estimates are computed under the counterfactual assumption that we had *not* grey-listed the country. 

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
3. Run the "produce_public_vs_private_raw_data.py" file to generate public data timeseries for each country and date.  Run this file twice, once with "raw" set to TRUE and once with "raw" set to FALSE as described above.
4. Run "grey_ebpred.r" and "grey_arrivals.r" to obtain counterfactual estimates of prevalence and arrivals had a country not been greylisted.
5. Run buildOPEdatabase.R to generate all off-policy evaluation and counterfactual analysis.  Summary statistics and plots can easily be created from the resulting dumped files.  

