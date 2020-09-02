#! /usr/bin/python3.8

"""
Author: Pritha Ghosh
Date: 1 September 2020
Place: Warsaw, Poland

Script for plotting rainfall distribution
Usage:
./rainfall_distribution.py

Input data was retrieved from https://tinyurl.com/y2sdu6fz
"""

import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("district.csv")

cities = ["KOLKATA", "NEW DELHI", "MUMBAI CITY", "CHENNAI", "BANGALORE URB", "HYDERABAD"]
colours = ["r", "b", "m", "g", "y", "c"]

df = data[data["DISTRICT"].isin(cities)].reset_index(drop = True)

labels = df.iloc[:, 2:14].columns

for i, row in df.iterrows():
    legend_added = False
    for j, col in enumerate(df.iloc[:, 2:14]):
        name = row["DISTRICT"]
        color = colours[cities.index(name)]
        if legend_added == True:
            name = None
        else:
            legend_added = True
        plt.bar(j, row[col], alpha = 0.2, width = 1, color = color, label = name)

title = "District Rainfall Normal (in mm): Data Period 1951-2000"
ylab = "Rainfall (in mm)"

plt.xticks(range(len(labels)), labels)
plt.legend()
plt.title(title)
plt.ylabel(ylab)
plt.show()