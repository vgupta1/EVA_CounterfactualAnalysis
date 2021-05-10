import numpy as np

import pandas as pd



history_start=20
future_end=20
window_f=7;
window_e=3;




#Downloading actual case and death data from owid and drop duplicates
url="https://covid.ourworldindata.org/data/owid-covid-data.csv"
df = pd.read_csv(url)
df=df.drop_duplicates()

#---Specifying correct date format
df['date']=pd.to_datetime(df['date'], format='%Y-%m-%d')

#Filtering for after July 15th. Estimates are valid after 08/02, so window of +- 20 days.
df=df[(df['date']>=pd.to_datetime('07/01/2020')) \
      & (df['date']<pd.to_datetime('11/22/2020'))]
    
#----Keeping only relevant fields: cases, deaths, tests and dealing with
#----Countries that make corrections in reported data (reporting many times
#----for same date)                                             
df=df[['iso_code','date','location','new_cases_smoothed_per_million',\
       'new_deaths_smoothed_per_million','new_tests_per_thousand']]
df=df.sort_values(by=['iso_code','date'])
df=df.loc[df.groupby(['iso_code','date']).date.idxmax()]


#--Transforming three digit ISO code to 2 digit ISO code using conversion table
#--found in other data 
url="../OtherData/short_to_iso.csv"
short_to_iso=pd.read_csv(url,error_bad_lines=False)
df=pd.merge(df,short_to_iso,how='left', left_on='iso_code', right_on='alpha-3')

#----Remove duplicates-----
df=df.drop_duplicates()

#-------Reading greylisting and prevalence data from OPE estimation outputs 
estimates=\
    pd.read_csv('../OPE_Outputs/ope_dat_TRUE_Window_3_MinTest_30_SmoothPrior_TRUE_2001_0.9.csv');
estimates['date']=pd.to_datetime(estimates['date_entry'])
estimates=estimates.sort_values('eb_type').groupby(['date','country']).tail(1)


#Adding prevalence estimates and greylisting status on public  data-----------
df=pd.merge(df,estimates,how='left',left_on=['alpha-2','date'] \
            ,right_on=['country','date'])
df=df.sort_values(by=['alpha-2','date'])
df=df.drop_duplicates()



#----Preparing column names for cases deaths and tests
#----Index specifies offset with respect to the row's reference date
cases=['cases'+str(i) for i in range(-history_start,future_end)]
deaths=['deaths'+str(i) for i in range(-history_start,future_end)]
tests=['tests'+str(i) for i in range(-history_start,future_end)]

#----Initializing new dataframe with formed column names----
X_dataframe=pd.DataFrame(columns=['country','date','grey','prev']\
                         +cases+deaths+tests)


#----Producing all X-y combos------------------------
for country in df['alpha-2'].drop_duplicates(): #---for every country
    #---and every date ------------
    for date in df[df['country']==country].date.drop_duplicates(): 
        #---generate window values from current_date-looking back 
        #by history_start to current date looking forward by future_end
        data_window=df[(df['alpha-2']==country) \
             & (df['date']>=date-pd.Timedelta(history_start,unit='day'))\
             & (df['date']<date+pd.Timedelta(future_end,unit='day'))] 
        #---produce entries for dataframe checking for appropriate length
        if len(data_window)==history_start+future_end:
                    data_window=data_window.sort_values(by=['date'])
                    #----adding greylisting status of current date
                    grey_status=data_window[data_window['date']==date]\
                        ['isCtryGrey'].values.tolist()
                     #----adding prevalence estimate of current date
                    prevalence=data_window[data_window['date']==date]\
                        ['eb_prev'].values.tolist()
                    #----adding public data of window around reference date
                    X_dataframe.loc[len(X_dataframe),\
                        ['country','date','grey','prev']+cases+deaths+tests]\
                        =np.array([country]+[date.strftime('%Y-%m-%d')]+ \
                        grey_status+prevalence+ \
                        data_window['new_cases_smoothed_per_million'].values.tolist()\
                        +data_window['new_deaths_smoothed_per_million'].values.tolist() \
                        +data_window['new_tests_per_thousand'].values.tolist())


#-----Store for Greylisting Counterfactuals and Public Data Evaluation----
X_dataframe.to_csv('../OtherData/Xall.csv')



