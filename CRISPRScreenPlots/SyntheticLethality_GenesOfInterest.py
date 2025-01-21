#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Use the CRISPRScreenAnalysis.R output files of a synthetic lethality screen, adjust the settings, and run the script in Spyder.
# This script can normalize the median log2 fold change to the essential and non-essential genes of a synthetic lethality screen, and plots T1control/T0 against T1treated/T0. 
# This normalization can improve the comparison treated - control if the treated arm did not have equal cell divisions, however the separation between the essentials and 
# non-essentials will not be improved. Synthetic lethal genes will be located around the lower half of the vertical 0 axis.
# Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
"""
##############################################################################################################

# Folder screen data:
folder = "C:/BioWin/CRISPRscreen/Z138"

# Cell line:
cellID = "Z138"

# Is there a T1drug/T1control comparison, 0: no, 1: yes
t2t1com = 1

# Size graph
size = 5

# Genes to emphasize, 0: all significant (from T1drug/T1control comparison), 1: specific genes 
allsignif = 0

GenesOfInterest = []
if allsignif==1:
  # If specific genes, Which ones?
  GenesOfInterest = ["BTK", "SYK", "PIK3R1"]

# Show all gene symbols, 0: no, 1: yes
GeneSymbol = 0

# Axes labels:
xlab = "Control (Relative log2 median fold change)"
ylab = "Venetoclax (Relative log2 median fold change)"

# # BCR-controlled adhesion screens:  
# xlab = 'PMA (Log2 median fold change)'
# ylab = r'$\alpha$IgM (Log2 median fold change)'

# Normalize to essentials and non-essentials (only for lethality), 0: no, 1: yes
NormalizeX = 1
NormalizeY = 1

# Normalize to (log)mean or median of the (non-)essentials, 0: mean, 1: median
meanOrmedian = 1

# Axes limit, 0 = automatic, 1: custom
Axlim = 0
# If automatic, Equal X and Y axes, 0 = no, 1: yes
XYequal = 1

if Axlim==1:
# Custom axes limits:    
    xmin, xmax, xticks=-3, 1, 1
    ymin, ymax, yticks=-3, 1, 1

# Colors:
call = 'silver'
cpos = 'red'
cneg = 'blue'
chit = 'black'

##############################################################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

df1 = pd.read_csv(folder+'/DESeq2 T0vsT1 Genes.csv', sep=',')
if NormalizeX == 0:
    df1['Control'] = np.log2(df1['MedianFoldChange'])
if NormalizeX == 1:
    if meanOrmedian==0:
        mCP = np.mean(np.log2(df1['MedianFoldChange'][df1['Type']=='p']))
        mCN = np.mean(np.log2(df1['MedianFoldChange'][df1['Type']=='n']))
    if meanOrmedian==1:
        mCP = np.median(np.log2(df1['MedianFoldChange'][df1['Type']=='p']))
        mCN = np.median(np.log2(df1['MedianFoldChange'][df1['Type']=='n']))
    df1['Control'] = (np.log2(df1['MedianFoldChange'])-mCN)/np.abs(mCP-mCN)
    
df2 = pd.read_csv(folder+'/DESeq2 T0vsT2 Genes.csv', sep=',')
if NormalizeY == 0:
    df2['Treated'] = np.log2(df2['MedianFoldChange'])
if NormalizeY == 1:
    if meanOrmedian==0:
         mCP = np.mean(np.log2(df2['MedianFoldChange'][df2['Type']=='p']))
         mCN = np.mean(np.log2(df2['MedianFoldChange'][df2['Type']=='n']))
    if meanOrmedian==1:
         mCP = np.median(np.log2(df2['MedianFoldChange'][df2['Type']=='p']))
         mCN = np.median(np.log2(df2['MedianFoldChange'][df2['Type']=='n']))
    df2['Treated'] = (np.log2(df2['MedianFoldChange'])-mCN)/np.abs(mCP-mCN)

if t2t1com == 1:
    df12 = pd.read_csv(folder+'/DESeq2 T1vsT2 Genes.csv', sep=',')
    df12['Stat'] = df12.loc[:,['rhoDepleted','rhoEnriched']].min(axis=1)
    df12['fdr'] = df12.loc[:,["fdrDepleted","fdrEnriched"]].min(axis=1)
    if allsignif==0:
      GenesOfInterest = df12['GeneSymbol'][df12['fdr']<0.1].tolist()

df = df1[['Type','GeneSymbol','Control']].merge(df2[['GeneSymbol','Treated']])
if t2t1com == 1:
    df = df.merge(df12[['GeneSymbol', 'Stat']])

dfpos = df[df['Type']=='p']
dfneg = df[df['Type']=='n']
dfhits = df[df['GeneSymbol'].isin(GenesOfInterest)].reset_index()

if Axlim==0:
    # Calculate axis limits:    
    xmin, xmax=round(df['Control'].min(),2)-0.3, round(df['Control'].max(),2)+0.3
    ymin, ymax=round(df['Treated'].min(),2)-0.3, round(df['Treated'].max(),2)+0.3
    if XYequal==1:
        xmin,ymin=min(xmin,ymin),min(xmin,ymin)
        xmax,ymax=max(xmax,ymax),max(xmax,ymax) 
    xticks=round((xmax-xmin)/5.1,2)
    yticks=round((ymax-ymin)/5.1,2)

# Calculate prediction interval from essentials/non-essentials
dfcon = pd.concat((dfpos,dfneg))
model = smf.ols('Treated ~ Control', dfcon)
results = model.fit()
regpoints = pd.Series(np.linspace(xmin,xmax,10), name='Control')
predictions = results.get_prediction(regpoints).summary_frame(.05)

# Make figure
plt.figure(figsize=(size,size))
ax = plt.subplot2grid((5,5), (1,0), colspan=4, rowspan=4)  
plt.suptitle(cellID, fontweight="bold") 

plt.fill_between(regpoints, predictions['obs_ci_lower'], predictions['obs_ci_upper'], alpha=.25, label='_nolegend_', color="silver") # pred int
plt.fill_between(regpoints, predictions['mean_ci_lower'], predictions['mean_ci_upper'], alpha=.15, label='_nolegend_',  color="silver") # conf int
plt.plot(regpoints, predictions['mean'], label='_nolegend_', color="silver", linestyle='--', linewidth=0.5)

if t2t1com==0:
    plt.scatter(df['Control'], df['Treated'], marker='o', color=call, s=1) 
    plt.scatter(dfpos['Control'], dfpos['Treated'], marker='o', color=cpos, s=1)
    plt.scatter(dfneg['Control'], dfneg['Treated'], marker='o', color=cneg, s=1)
if t2t1com==1:
    plt.scatter(df['Control'], df['Treated'], marker='o', color=call, s=-2*np.log10(df['Stat']))
    plt.scatter(dfpos['Control'], dfpos['Treated'], marker='o', color=cpos, s=-2*np.log10(dfpos['Stat']))
    plt.scatter(dfneg['Control'], dfneg['Treated'], marker='o', color=cneg, s=-2*np.log10(dfneg['Stat']))
if GeneSymbol==1:
    for i, txt in enumerate(df['GeneSymbol']):
        plt.annotate(txt, (df['Control'][i]+(xmax-xmin)/100, df['Treated'][i]), rotation=22.5, fontsize=7, color="gray", ha='left')
if len(GenesOfInterest)>0:
    if t2t1com==0:
        plt.scatter(dfhits['Control'], dfhits['Treated'], marker='o', color=chit,  s=1)
    if t2t1com==1:
        plt.scatter(dfhits['Control'], dfhits['Treated'], marker='o', color=chit,  s=-2*np.log10(dfhits['Stat']))
    for i, txt in enumerate(dfhits['GeneSymbol']):
        plt.annotate(txt, (dfhits['Control'][i]+(xmax-xmin)/100, dfhits['Treated'][i]), rotation=22.5, fontsize=7, color="black", ha='left')  

lgnd = plt.legend(('All genes', 'Essentials', 'Non-essentials'), fontsize=8, loc=2)
for handle in lgnd.legendHandles:
    handle.set_sizes([10])

plt.axvline(x=0, color=cneg, linestyle='--', linewidth=0.5)
plt.axhline(y=0, color=cneg, linestyle='--', linewidth=0.5)
if NormalizeX == 1:
    plt.axvline(x=-1, color=cpos, linestyle='--', linewidth=0.5)
if NormalizeY == 1:
    plt.axhline(y=-1, color=cpos, linestyle='--', linewidth=0.5)
plt.plot([-10, 10], [-10, 10], color='gray', linestyle='--', linewidth=0.5)
plt.xlabel(xlab, fontsize=10)
plt.ylabel(ylab, fontsize=10)
plt.xticks(np.arange(xmin, xmax, xticks))
plt.xlim(xmin,xmax)
plt.yticks(np.arange(ymin, ymax, yticks))
plt.ylim(ymin,ymax)

plt.subplot2grid((5,5), (1,4), rowspan=4)
sns.kdeplot(y=dfpos['Treated'], fill= True, linewidth= 1, color=cpos)
sns.kdeplot(y=dfneg['Treated'], fill= True, linewidth= 1, color=cneg)
plt.ylim(ymin,ymax)
plt.axis('off')  

plt.subplot2grid((5,5), (0,0), colspan=4)
sns.kdeplot(x=dfpos['Control'], fill= True, linewidth= 1, color=cpos)
sns.kdeplot(x=dfneg['Control'], fill= True, linewidth= 1, color=cneg)
plt.xlim(ymin,ymax)
plt.axis('off')  

plt.tight_layout(pad=0, w_pad=-24, h_pad=0)
plt.subplots_adjust(wspace=0)
plt.subplots_adjust(hspace=0)
plt.savefig(folder+'/CRISPR_SL_'+cellID+'_PY.pdf', bbox_inches='tight', transparent=True)