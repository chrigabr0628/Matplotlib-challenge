#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.

# Analysis:
# 
# 
# Extraction, Loading, and Transformation:
# 
# Set up dependencies. Created new directory path. Loaded in csv files. Merged datasets into one single data frame. Calculated mouse count and number of duplicate mouse ID. Cleaned data by removing duplicated mouse ID. Created a summary table pf the cleaned/merged data. Created graphs and charts to visualize the data.
# 
# 
# 
# Outcome:
# 
# There were more male mice than female mice 
# Larger mice appear to have had the largest tumor volum .
# According to the box plot, the drugs Capomulin and Ramicane appear to have had the greatest positive effect on tumor volume as they both had the lowest, final tumor volume. Ramicane appears to have had a slightly greater effec. 
# Initially, the drug Capomulin appears to have increased the tumor volume, significantly decreased the tumor volume around 20 days, raised the tumor volume again around 25 days, so on and so forth. It appears the Capomulin drug was volet le.
# The drug Propriva appears to have had the least positive effect on tumor v.e.ct on tumor volume.
#  
# 

# In[1]:


#PLEASE NOTE: I recieved code source from AskBCS, Xpert Learning Assitant, and tutors.


# In[2]:


# Dependencies and Setup
import os
import csv
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from pathlib import Path
from scipy.stats import linregress


# In[3]:


# Specify the new directory path
new_directory = 'C:\\Users\\chris\\Documents\\Module 5 Challenge\\Pymaceuticals-challenge'

# Change the current working directory
os.chdir(new_directory)

# Print the new current working directory
print("New Current Working Directory:", os.getcwd())


# In[4]:


# Study data files
mouse_metadata_csv = Path("data/Mouse_metadata.csv")
study_results_csv = Path("data/Study_results.csv")


# In[5]:


# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_csv)
study_results = pd.read_csv(study_results_csv)


# In[6]:


# Combine the data into a single DataFrame
merged_df = pd.merge(mouse_metadata, study_results, on=["Mouse ID"])

# Display the data table for preview
merged_df.head()


# In[7]:


mice_count = mouse_metadata["Mouse ID"].count()
mice_count


# In[8]:


# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint

duplicate_mice = merged_df[merged_df.duplicated(["Mouse ID", "Timepoint"])]
duplicate_mice_data = duplicate_mice["Mouse ID"].unique()
duplicate_mice_data


# In[9]:


# Optional: Get all the data for the duplicate mouse ID. 

mice_data = merged_df[merged_df["Mouse ID"].isin(duplicate_mice_data)]
mice_data


# In[10]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.

clean_df = merged_df[~merged_df["Mouse ID"].isin(duplicate_mice_data)]
clean_df


# In[11]:


# Checking the number of mice in the clean DataFrame.

len(clean_df["Mouse ID"].unique())


# Summary Statistics

# In[12]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.


summary_statistics = clean_df.groupby("Drug Regimen")
mean = summary_statistics["Tumor Volume (mm3)"].mean()
median = summary_statistics["Tumor Volume (mm3)"].median()
variance = summary_statistics["Tumor Volume (mm3)"].var()
standard_deviation = summary_statistics["Tumor Volume (mm3)"].std()
SEM = summary_statistics["Tumor Volume (mm3)"].sem()

summary_statistics = pd.DataFrame({"Mean": mean, 
                                   "Median": median,
                                   "Variance": variance,
                                   "Standard Deviation": standard_deviation,
                                   "SEM" : SEM})
summary_statistics


# In[13]:


# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)
# Using the aggregation method, produce the same summary statistics in a single line

clean_df.groupby("Drug Regimen").agg({"Tumor Volume (mm3)":["mean", "median", "var", "std", "sem"]})


# Bar and Pie Charts

# In[14]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.

clean_df["Drug Regimen"].value_counts().plot.bar()

plt.ylabel("# of Observed Mouse Timepoints")
plt.show()


# In[15]:


# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.

mouse_pyplot = clean_df["Drug Regimen"].value_counts()

x = mouse_pyplot.index

y = mouse_pyplot.values


plt.bar(x, y)
plt.xticks(rotation = 90)

plt.ylabel("# of Observed Mouse Timepoints")
plt.show()


# In[16]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas

plot = clean_df["Sex"].value_counts().plot.pie(subplots=True, autopct='%1.1f%%', figsize=(11, 6))
plt.ylabel("Sex")
plt.show()


# In[17]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot

mouse_piechart = clean_df["Sex"].value_counts()

x = mouse_piechart.index

y = mouse_piechart.values


plt.pie(y, labels = x, autopct='%1.1f%%')
plt.ylabel("Sex")
plt.show()


# Quartiles, Outliers and Boxplots

# In[18]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
greatest_timepoint = clean_df.groupby(["Mouse ID"])["Timepoint"].max()

# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_df = pd.merge(greatest_timepoint, clean_df, how="inner", on=["Mouse ID", "Timepoint"])
merged_df.head()


# In[19]:


# Put treatments into a list for the for loop
treatments = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']

# Create empty list to fill with tumor volume data (for plotting)
tumor_vol_data = []

# Iterate over the treatments list
for treatment in treatments:
    # Filter the DataFrame by the current treatment
    tumor_data = merged_df.loc[merged_df["Drug Regimen"] == treatment, "Tumor Volume (mm3)"]
    
    # Append tumor volume data for the current treatment
    tumor_vol_data.append(tumor_data)
    
    # Calculate quartiles and IQR for the current treatment
    quartiles = tumor_data.quantile([0.25, 0.5, 0.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq - lowerq
    lower_bound = lowerq - (1.5 * iqr)
    upper_bound = upperq + (1.5 * iqr)
    outliers = tumor_data.loc[(tumor_data < lower_bound) | (tumor_data > upper_bound)]
    
    print(f"{treatment}'s potential outliers: {outliers}")


# In[20]:


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.

fig, ax = plt.subplots()
ax.boxplot(tumor_vol_data, labels=treatments, flierprops=dict(markerfacecolor='r', markersize=12))
ax.set_ylabel('Final Tumor Volume (mm3)')
plt.show()


# Line and Scatter Plots

# In[21]:


# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin

capomulin_data = clean_df.loc[(clean_df["Drug Regimen"] == "Capomulin")  & (clean_df['Mouse ID'] == 'l509')]

plt.plot(capomulin_data['Timepoint'], capomulin_data['Tumor Volume (mm3)'])
plt.title("Capomulin treatment of mouse l509")
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.show()


# In[22]:


# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen

capomulin_data2 = clean_df[clean_df["Drug Regimen"] == "Capomulin"]
capomulin_data2 = capomulin_data2[["Mouse ID", "Weight (g)", "Tumor Volume (mm3)"]]
avg_tumor_vol = capomulin_data2.groupby("Mouse ID").mean()
plt.scatter(avg_tumor_vol["Weight (g)"], avg_tumor_vol["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()


# Correlation and Regression

# In[26]:


# Calculate the correlation coefficient and a linear regression model
# for mouse weight and average observed tumor volume for the entire Capomulin regimen

correlation = capomulin_data2.groupby("Mouse ID").agg({"Tumor Volume (mm3)": 'mean', 'Weight (g)': 'mean'})
correlation_final = st.pearsonr(correlation["Weight (g)"], correlation["Tumor Volume (mm3)"])[0]
(slope, intercept, rvalue, pvalue, stderr) = linregress(correlation["Weight (g)"], correlation["Tumor Volume (mm3)"])
regression_values = correlation["Weight (g)"] * slope + intercept

print(f"The correlation between mouse weight and the average tumor volume is {correlation_final:.2f}")

plt.scatter(correlation["Weight (g)"], correlation["Tumor Volume (mm3)"])
plt.plot(correlation["Weight (g)"], regression_values, "r-")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()

