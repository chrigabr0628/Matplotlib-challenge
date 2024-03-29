#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Dependencies and Setup
import os
import csv
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from pathlib import Path


# In[2]:


# Specify the new directory path
new_directory = 'C:\\Users\\chris\\Documents\\Module 5 Challenge'

# Change the current working directory
os.chdir(new_directory)

# Print the new current working directory
print("New Current Working Directory:", os.getcwd())


# In[3]:


# Study data files
mouse_metadata_csv = Path("data/Mouse_metadata.csv")
study_results_csv = Path("data/Study_results.csv")


# In[4]:


# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_csv)
study_results = pd.read_csv(study_results_csv)


# In[ ]:


# Combine the data into a single DataFrame
merged_df = pd.merge(mouse_metadata_csv, study_results_csv, how="left", on=["Mouse ID"])

# Display the data table for preview
merged_df.head()


# In[ ]:




