import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import boxcox

df = pd.read_csv('TimeSeries/data/raw/AirPassengers.csv' , index_col='Month' , parse_dates=True)


df['#Passengers'].plot(figsize=(20,16))


np.sqrt(df['#Passengers']).plot(figsize=(20,10))
np.log(df['#Passengers']).plot(figsize=(20,10))

data , lam = boxcox(df['#Passengers'])
df['boxcox'] = data
df['boxcox'].plot(figsize=(20,16))
