# The following code is designed to check for overlap between GS array and sequencing panel

# import libraries
import pandas as pd
import numpy as np
from tqdm import tqdm

# import data
probeSeqData = pd.read_csv('NGS580_targets.csv')
gsaData = pd.read_csv('GSA-24v1-0_A1.annotated.txt',delimiter='\t')

# All chromosomes need to be character values
gsaData['Chr'] = gsaData['Chr'].astype(str)
probeSeqData['Chromosome'] = probeSeqData['Chromosome'].astype(str)
probeSeqData['Chromosome'] = probeSeqData['Chromosome'].str.replace("chr", "")
## Now that we have the data imported, we need to setup the location range columns

# Steps:
#
# 1. Ensure that the gsaData contains unique rows.
# 2. Construct the interval [start,stop] for each chromosome
# 3. Check if a given value maps to the interval [start,stop] in probeSeqData
# 4. Record only if True or False

# First, what are all the chromosomes we need to check
allChrs = gsaData.Chr.unique().tolist()

## Checking scheme:

# 1. Iterate over each chromosome
# 2. Get GSA data and probe data subsets for this chromosome
# 3. Build the temporary range value variables
# 4. Check if the chromosome + mapInfo combination exists in the probe data

# this is where we will record the result
gsaMappedData = pd.DataFrame(columns=['Chromosome','Location','isInProbe'])

# iterate over the list allChrs
# NOTE: tqdm is just for the nice progress bar
for chromo in tqdm(allChrs):
    # get the gsa data for this chromo
    tempGSA = gsaData[gsaData.Chr == chromo]
    # keep only columns we need
    tempGSA = tempGSA[['Chr','MapInfo']]
    # rename the columns
    tempGSA.columns = ['Chromosome','Location']
    # create a flag to say if the location is in probe data
    tempGSA['isInProbe'] = False
    # drop duplicate rows
    tempGSA.drop_duplicates(inplace=True)
    
    # get the probe data to create ranges
    tempProbe = probeSeqData[probeSeqData.Chromosome == chromo]
    
    # Extract the range values from probeSeqData
    ## Look up on how to extract a single column from a pandas dataframe as a list
    startVals = tempProbe.Start.tolist()
    stopVals = tempProbe.Stop.tolist()
    
    # get all the locations from tempGSA
    locations = tempGSA.Location.unique().tolist()
    
    # Get number of ranges
    N = len(startVals)
    
    # iterate over each location
    for location in locations:
        # iterate over each range
        for i in range(0,N):
            tempRange = range(startVals[i],stopVals[i]+1)
            if location in tempRange:
                tempGSA = tempGSA.set_value(tempGSA.Location == location,'isInProbe',True)
                
    # take the data for this chromosome and add it to the output dataset
    gsaMappedData = pd.concat([gsaMappedData,tempGSA])

# print a report
locationsMapped = gsaMappedData.isInProbe.value_counts()[True]
totalLocations = gsaMappedData.shape[0]
percentMatch = locationsMapped/totalLocations * 100
print('Total',str(locationsMapped),'out of',str(totalLocations),'locations were mapped')
print('We have a',str(percentMatch),'% match')

# save the mapped data
gsaMappedData.to_csv('gsaMappedData.csv',index=False,index_label=False)