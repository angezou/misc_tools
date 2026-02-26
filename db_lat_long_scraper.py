from pystac_client import Client
from odc.stac import load
from datetime import datetime


from openmeteopy import OpenMeteo
from openmeteopy.hourly import HourlyHistorical
from openmeteopy.daily import DailyHistorical
from openmeteopy.options import HistoricalOptions
from openmeteopy.utils.constants import *


def retrieve_ndvi(point_lat, point_lon, time_of_interest):
    """
    Given latitude, longtitude, and time period, returns ndvi values 

    Args:
        point_lat (float): latitude
        point_lon (float): longitude
        time_of_interest (time period): YYYY-MM-DD/YYYY-MM-DD
    """
    sentinel_search_url = "https://earth-search.aws.element84.com/v1"
    sentinel_stac_client = Client.open(sentinel_search_url)

    search = sentinel_stac_client.search(

           collections=["sentinel-2-l2a"], 
           bbox=[point_lon - 0.0005, point_lat - 0.0005, point_lon + 0.0005, point_lat + 0.0005],
           datetime = time_of_interest)
    
    if len(search.item_collection()) == 0:
        return("NaN")

    selected_item=min(search.items(),key=lambda item: item.properties["eo:cloud_cover"])
    selected_bands = ["nir", "red"]
    data = load([selected_item], bands=selected_bands, bbox = [point_lon - 0.0005, point_lat - 0.0005, point_lon + 0.0005, point_lat + 0.0005]).isel(time=0)
    sentinel_df = data.to_dataframe()
    sentinel_df["ndvi"] = (sentinel_df["nir"].astype(float) - sentinel_df["red"].astype(float))/(sentinel_df["nir"].astype(float) + sentinel_df["red"].astype(float))
    avg_ndvi = sum(sentinel_df['ndvi'])/len(sentinel_df['ndvi'])
    return(avg_ndvi)


def retrieve_tmp_prcp(point_lat, point_lon, start_date, end_date):
    """
    Given latitude, longtitude, and time period, returns minimum temperature, maximum temperature, and total precipitation 

    Args:
        point_lat (float): latitude
        point_lon (float): longitude
        time period start date: YYYY-MM-DD
        time period end date: YYYY-MM-DD
    """
    
    hourly = HourlyHistorical()
    daily = DailyHistorical()
    options = HistoricalOptions(point_lat, point_lon, start_date=start_date, end_date=end_date)

    mgr = OpenMeteo(options, hourly.all(), daily.all())
    # Download data
    meteo = mgr.get_pandas()
    min_tmp = meteo[1][['temperature_2m_min']].min()
    max_tmp = meteo[1][['temperature_2m_max']].max()
    sum_prcp = meteo[1][['precipitation_sum']].sum()
    climate_metrics = pd.Series([min_tmp.iloc[0], max_tmp.iloc[0], sum_prcp.iloc[0]])
    return(climate_metrics)


