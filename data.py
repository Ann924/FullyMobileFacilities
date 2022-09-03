import random
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def read_file(county_name: str):
    if county_name!= "albemarle" and county_name!="charlottesville_city":
        print("Invalid County Name")
    else:
        filename = f"usa_va_{county_name}_adult_activity_location_assignment_week.csv"
        df = pd.read_csv(filename)
        return df

def clean_data(df_day: pd.DataFrame):
    # TODO: drop meaningless columns
    
    HOME_SHIFT=1000000000

    home_activities = df_day[df_day["activity_type"]==1]["lid"]
    df_day.loc[home_activities.index, "lid"] += HOME_SHIFT
    
    return df_day

def get_location_directions(df: pd.DataFrame):
    return df[["lid","longitude", "latitude"]].groupby("lid").mean().to_dict('index')

def find_potential_facilities(df: pd.DataFrame):
    return set(df[df["activity_type"]!=1].lid)

def pid_hour_breakdown(df: pd.DataFrame, day: int, start_hour: int, end_hour: int):
    
    HR = 3600
    
    hour_start = day*24*HR + start_hour*HR
    hour_end = day*24*HR + end_hour*HR
    
    df["end_time"] = df["start_time"]+df["duration"]
    
    df = df[(df["start_time"]<hour_end) & (df["end_time"]>hour_start)].copy()
    
    df["combined_loc"] = df[["lid", "longitude", "latitude", "start_time", "end_time"]].apply(
        lambda x: (int(x["lid"]), (x["longitude"], x["latitude"]), (x["start_time"], x["end_time"])), axis=1)
    
    df_pid = (df.groupby("pid")["combined_loc"].apply(list)).apply(
        lambda x: {(i-day*24*HR)//HR: [loc for loc in x if (loc[2][0]<=i and loc[2][1]>=i+HR)] for i in range(hour_start, hour_end, HR)})
    
    return df_pid.to_frame()

def random_filter(df:pd.DataFrame, random_state=42):
    df["condensed"] = df["combined_loc"].apply(lambda x: [(h, l[0]) for h, loc in x.items() for l in loc])
    df["condensed_len"] = df["condensed"].apply(lambda x: len(x))

    # remove clients without any location-hour assignments
    df = df.drop(list(df[df["condensed_len"]==0].index))

    df_selected = df["condensed"].apply(lambda x: x[random.randint(0, len(x)-1)])

    return df_selected.to_dict({})

def get_data(county_name: str, day:int, start_hour: int, end_hour: int):

    df_full = read_file(county_name)
    df_clean = clean_data(df_full)

    # Must be separate since some of the activity locations are recorded as home visitations (but are not potential facility locations)
    potential_facilities = find_potential_facilities(df_clean)
    location_directory = get_location_directions(df_clean)

    df_pid = pid_hour_breakdown(df_clean, day, start_hour, end_hour)
    pid_assignment = random_filter(df_pid)

    return potential_facilities, location_directory, pid_assignment

#potential_facilities, location_directory, pid_assignment = get_data("charlottesville_city", 5, 6, 20)
potential_facilities, location_directory, pid_assignment = get_data("albemarle", 5, 6, 20)
print(len(potential_facilities), len(location_directory), len(pid_assignment))

"""

TODO:

Store processed data in csv for easier access and quicker runtime
Speed up dataframe computations (and storage/space requirements)
Comment methods
Potentially change datastructure

"""