import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.ensemble import GradientBoostingClassifier


#--This code evaluates how well publicly available data (14-day case data)
#--predicts the risk of a country, i.e.  whether its true prevalence
#-- exceeds its median true prevalence over the summer. Intuitively, if a 
#--country's public data lags the true prevalence of travelers 
#--by l days, then we should be able to predict it's risk status using case 
#--data l days into the future. This code produces for every country and every
#--lag the achieved AUC. 



#--The train_evaluate function expects as input the 14 day case data as well
#--as the target variables, trains a model on the training data and returns 
#--the off-sample AUC. 
 
def train_evaluate(X, y):
  #--defining the chosen model--
  model=GradientBoostingClassifier()
  #--splitting data into training and testing. Split is done by time since 
  #--we are working with timeseries data. 
  X_train, X_test, y_train, y_test = train_test_split( X, y, test_size=0.5, shuffle=False)
  #--training model on training data
  clf=model.fit(X_train, y_train)
  #--evaluating model on test data
  y_pred=clf.predict_proba(X_test)
  #--obtaining false positive and false negative rate
  fpr, tpr, threshold = metrics.roc_curve(y_test, y_pred[:,0], pos_label=False)
  #--returning achieved AUC for the given X,y pair. 
  return metrics.auc(fpr, tpr)




#Download publicly reported case  data from OWID
url="https://covid.ourworldindata.org/data/owid-covid-data.csv"
df = pd.read_csv(url)
df=df.drop_duplicates()



#--Transforming 3 digit ISO codes to 2 digit equivalents
url="../OtherData/short_to_iso.csv"
short_to_iso=pd.read_csv(url,error_bad_lines=False)
df=pd.merge(df,short_to_iso,how='left', left_on='iso_code', right_on='alpha-3')


#--Loading Estimates of true Prevalence----
estimates=pd.read_csv('../OPE_outputs/ope_dat_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv');
estimates=estimates.sort_values('eb_type').groupby(['date_entry','country']).tail(1)


#--Adding prevalence estimates on public data-----------
df=pd.merge(df,estimates,how='left',left_on=['alpha-2','date'] ,\
            right_on=['country','date_entry'])
df=df.sort_values(by=['alpha-2','date'])
df['date']=pd.to_datetime(df['date'])
df=df.drop_duplicates()


#--Producing List of all countries for which data is available
countries=df['alpha-2']
countries=countries.drop_duplicates()

#--Reading list of Non-Blacklisted countries
url="../OtherData/countries_allowed.csv"
allow=pd.read_csv(url, header=None)[0].values.tolist()
grey=['ES','BE','MT','BG','RO','GR','CZ','SE','NL']

#--Focusing on allowed countries that are not greylisted, since for the latter
#--measured prevalence is inaccurate during the period when they were greylisted
considered_countries=set(countries)-set(grey)

#--Initializing lags vector and result dataframe
lags=[];
info=pd.DataFrame(columns=['Country','Lag','AUC'])



maximizers=[];
country_list=[]
lags=[]
crap=[]
lowess = sm.nonparametric.lowess
maxxx=0          
decent=[]


#--For each country 
for country in considered_countries:
    res=[] #--KIMON check
    this_country=df[df['alpha-2']==country]
    #--for each country we require that we have enough data to produce 
    #--prevalence estimates for at least 70 days. This is in order to exclude
    #--countries with very rare arrivals. 
    valid_eb=this_country[this_country['eb_prev'].notnull()]
    enough_tests=(len(valid_eb)>70);
    if enough_tests:
        #--for each lag up to 20 days
        for lag in range(0,20):
            #--we produce the X: cases y:risk status 
            X=[];
            y=[];
            for index,row in valid_eb.iterrows():
                date=row['date']
                moved_dates=pd.date_range(date-pd.Timedelta(14-lag,unit='D'),\
                                          date+pd.Timedelta(lag,unit='D'))
                cases=this_country[np.isin(this_country['date'], moved_dates)]\
                    ['new_cases_smoothed_per_million'].values.tolist()
                X.append(cases)
                y.append(valid_eb[valid_eb['date']==date].iloc[0]['eb_prev'])
            X=np.array(X)
            y=np.array(y)
            #--A country is defined to be risky on a given day if it's 
            #--prevalence is higher it's median during the summer
            y=(y>np.median(y))
            
            #--For the given country and lag we train and evaluate a model that
            #--classifies risky or not and adding resulting AUC to result. 
            country_lag_auc=train_evaluate(X, y)
            info=info.append({'Country':country, 'Lag':lag, \
               'AUC': max(country_lag_auc,1-country_lag_auc)},ignore_index=True)
            

#--Producing Country Lag AUC profile to be used in clustering
info.to_csv("../OtherData/country_lag_auc_profile.csv")






















    