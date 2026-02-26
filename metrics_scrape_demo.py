import pandas as pd 
import datetime
from db_lat_long_scraper import *

loc_df = pd.read_excel('sample_loc_df.xlsx')
month_offset = 1

loc_df['future_date'] = pd.to_datetime(loc_df['date']) + pd.DateOffset(months=month_offset)

loc_df['date'] = loc_df['date'].astype(str)
loc_df['future_date'] = loc_df['future_date'].astype(str)

loc_df['date_final'] = loc_df['date'] + "/" + loc_df['future_date']

loc_df['ndvi'] = loc_df.apply(lambda loc_df: retrieve_ndvi(loc_df.latitude, loc_df.longitude, loc_df.date_final), axis=1)

loc_df[['min_tmp','max_tmp','prcp']] = loc_df.apply(lambda loc_df: retrieve_tmp_prcp(loc_df.latitude, loc_df.longitude, loc_df.date, loc_df.future_date), axis=1)
