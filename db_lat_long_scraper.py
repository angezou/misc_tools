from pystac_client import Client
import stackstac
from odc.stac import load


import matplotlib.pyplot as plt


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
    sentinel_df["ndvi"] = (sentinel_df["nir"] - sentinel_df["red"])/(sentinel_df["nir"] + sentinel_df["red"])
    avg_ndvi = sum(sentinel_df['ndvi'])/len(sentinel_df['ndvi'])
    return(avg_ndvi)


