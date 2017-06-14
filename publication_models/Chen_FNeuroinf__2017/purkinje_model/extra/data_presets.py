import sys
import os
import random
import numpy as np
from numpy import arange, sin, pi

from constants import AVOGADRO, E_CHARGE

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def getMeanData(data):
    sums = {}
    n_twin = 0
    first_data = data.values()[0]
    for s in first_data.keys():
        sums[s] = 0.0
    for t in data.keys():
        n_twin += 1
        current_data = data[t]
        for s in current_data.keys():
            sums[s] += current_data[s]

    for s in sums.keys():
        sums[s] = sums[s] / n_twin
    return sums

def getInfluxRate(neuron_curr, surf_area, vol):
    # Data file now store values in SI units
    sum_curr = neuron_curr * surf_area
    # For Ca2+, I = Q/t = N2e/t, so N/t = I/2e, where e is elementary charge,
    # devided by volume (in liter = vol * 1e3) and AVOGADRO
    influx_mol_per_sec_per_liter = sum_curr / 2.0 / E_CHARGE / vol / AVOGADRO / 1e3
    return influx_mol_per_sec_per_liter

def readData(data_file):
    dataset = {}
    file = open(data_file, 'r')
    for line in file:
        line_secs = line.split()
        # entry info
        if line_secs[0] == '#Entries:':
            dataset["Entries"] = line_secs[1:]
            dataset["Data"] = []
        else:
            data_line = [float(value) for value in line_secs]
            dataset["Data"].append(data_line)
    return dataset

def SI2NEURON(dataset):
    new_dataset = {}
    new_dataset["Entries"] = dataset["Entries"]
    new_dataset["Data"] = []
    for data_line in dataset["Data"]:
        new_data_line = []
        time = data_line[0] * 1000
        new_data_line.append(time)
        for v in data_line[1:]:
            new_data_line.append(v*1e6)
        new_dataset["Data"].append(new_data_line)
    return new_dataset

def genCaInfluxProfile(dataset, sur_areas, vols, start_time, end_time, time_win):
    roi_entries = dataset["Entries"][1:]
    n_entries = len(roi_entries)
    curr_start_time = start_time
    curr_end_time = start_time + time_win
    n_win = 0
    sum_curr = [0.0] * n_entries
    influx_profile = {}
    influx_profile["Entries"] = dataset["Entries"]
    influx_profile["Data"] = []
    for data in dataset["Data"]:
        time = data[0]
        
        curr_data = data[1:]
        if start_time > time: continue
        
        if isclose(time, curr_end_time):
            influx_data = []
            influx_data.append(curr_start_time)
            for i in range(n_entries):
                mean_curr = sum_curr[i] / n_win
                influx_rate = getInfluxRate(mean_curr, sur_areas[roi_entries[i]], vols[roi_entries[i]])
                influx_data.append(influx_rate)
            influx_profile["Data"].append(influx_data)
            sum_curr = [0.0] * n_entries
            n_win = 0
            curr_start_time = curr_end_time
            curr_end_time = curr_start_time + time_win
            
        if isclose(time, end_time): break
        for i in range(n_entries):
            sum_curr[i] += abs(curr_data[i])
        n_win += 1
    return influx_profile