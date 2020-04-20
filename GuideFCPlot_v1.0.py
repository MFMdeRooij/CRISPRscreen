# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 13:41:02 2020

@author: Martin F.M. de Rooij, PhD
"""
# After performing the DESeq2 script sucessfully, you can produce publishing-grade GuideFC plots
###############################################################################################################
#                                                      SETTINGS

# Copy-paste the required folder (use / instead of \ )
folder = "H:/BioWin/Screens/Namalwa"

# Gene symbols (if there are subset, cluster the subsets)
genes = ['BTK', 'SYK', 'PIK3R1', 'CSK', 'PRKCE', 'PRKCB', 'ACTR2', 'GUK1', 'MAP2K1', 'MAP2K2', 'AKT1', 'AKT2', 'AKT3', 'LYN']
# #Alternatively, take all significant genes:
# dfhits = pd.read_csv('H:/BioWin/Screens/Namalwa"DESeq2 T1vsT2 Genes.csv', sep=',')
# genes = dfhits[dfhits['fdrDepleted']<0.1]['GeneSymbol']

# Are there subsets of genes: 0 = no, 1 = yes
geneSubsets = 1

if geneSubsets==1:
    # number of genes per subset:
    subset1 = 5
    subset2 = 3
    subset3 = 6 # For 2 subsets this should be 0

# Titles:
title1 = 'PMA/input'
title2 = r'$\alpha$IgM/input'
title3 = r'$\alpha$IgM/PMA'  

# X-axis limits:
xmin, xmax, xticks = -1.25, 1, 0.25

# RRA score: 0=depletion, 1=enrichment
rrascore = 0
###############################################################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

genenr = np.arange(len(genes))+1
genedf = pd.DataFrame({'nr':genenr, 'GeneSymbol': genes})

# Tables
dfname = glob.glob(folder+'/*T0vsT1 Guides.csv')[0]  
dfP = pd.read_csv(dfname, sep=',')
dfP['l2fc']=np.log(dfP['FoldChange'])/np.log(2)
dfname = glob.glob(folder+'/*T0vsT1 Genes.csv')[0]  
dfPrra = pd.read_csv(dfname, sep=',')

dfname = glob.glob(folder+'/*T0vsT2 Guides.csv')[0]  
dfM = pd.read_csv(dfname, sep=',')
dfM['l2fc']=np.log(dfM['FoldChange'])/np.log(2)
dfname = glob.glob(folder+'/*T0vsT2 Genes.csv')[0] 
dfMrra = pd.read_csv(dfname, sep=',')

dfname = glob.glob(folder+'/*T1vsT2 Guides.csv')[0]  
dfC = pd.read_csv(dfname, sep=',')
dfC['l2fc']=np.log(dfC['FoldChange'])/np.log(2)
dfname = glob.glob(folder+'/*T1vsT2 Genes.csv')[0] 
dfCrra = pd.read_csv(dfname, sep=',')

i=0
for df,dfrra in [(dfP,dfPrra),(dfM,dfMrra),(dfC,dfCrra)]:  
    dfg=df[df['GeneSymbol'].isin(genes)].iloc[:,[0,8]]  
    dfg = pd.merge(dfg,genedf, on='GeneSymbol')  
    if geneSubsets==1:
        dfg1 = dfg[dfg['GeneSymbol'].isin(genes[0:subset1])]
        dfg2 = dfg[dfg['GeneSymbol'].isin(genes[subset1:subset1+subset2])]
        dfg3 = dfg[dfg['GeneSymbol'].isin(genes[subset1+subset2:subset1+subset2+subset3])]
    
    dfrra=pd.merge(dfrra[dfrra['GeneSymbol'].isin(genes)],genedf, on='GeneSymbol')
    dfrra['color']='black'
    if rrascore==0:
        dfrra['10lrra']=-1*np.log(dfrra['rhoDepleted'])/np.log(10)
        dfrra.loc[dfrra['fdrDepleted']<0.1,'color']='red'
    elif rrascore==1:
        dfrra['10lrra']=-1*np.log(dfrra['rhoEnriched'])/np.log(10)
        dfrra.loc[dfrra['fdrEnriched']<0.1,'color']='red'
       
    plt.figure(1, figsize=(30,7))
    ax = plt.subplot2grid((5,24), (0,i*8), colspan=6, rowspan=1) 
    if i==2:
        sns.distplot(df['l2fc'], hist=False, kde_kws = {'shade': False, 'linewidth': 2}, color='magenta', label='All guides')
        sns.distplot(df[df['Type']=='n']['l2fc'], hist=False, kde_kws = {'shade': False, 'linewidth': 2}, color='orange', label='Non-essential genes') 
        sns.distplot(df[df['padj']<0.1]['l2fc'], hist=False, kde_kws = {'shade': False, 'linewidth': 2}, color='darkorchid', label='Significant (FDR < 0.1)')
    else:
        sns.distplot(df['l2fc'], hist=False, kde_kws = {'shade': False, 'linewidth': 2}, color='magenta', label=None)
        sns.distplot(df[df['Type']=='n']['l2fc'], hist=False, kde_kws = {'shade': False, 'linewidth': 2}, color='orange', label=None) 
        sns.distplot(df[df['padj']<0.1]['l2fc'], hist=False, kde_kws = {'shade': False, 'linewidth': 2}, color='darkorchid', label=None)      
    plt.xlim(xmin,xmax)
    if i==0:
        plt.title(title1)  
    elif i==1:
        plt.title(title2) 
    elif i==2:
        plt.title(title3)    
    plt.axis('off')   
    
    plt.subplot2grid((5,24), (1,i*8), colspan=6, rowspan=4)  
    plt.scatter(df['l2fc'], np.repeat(7,len(df)), marker="|", s=1000000, c='silver', alpha=0.05)
    if geneSubsets==1:
        plt.scatter(dfg1['l2fc'], dfg1['nr'], marker="|", s=300,  c='navy')
        plt.scatter(dfg2['l2fc'], dfg2['nr'], marker="|", s=300, c='firebrick')
        plt.scatter(dfg3['l2fc'], dfg3['nr'], marker="|", s=300, c='green')
    else:
       plt.scatter(dfg['l2fc'], dfg['nr'], marker="|", s=300,  c='navy') 
    plt.xlim(xmin,xmax)   
    plt.yticks(np.arange(1,1+len(genes)), genes)
    plt.ylim(len(genes)+0.5,0.5)    
    plt.xlabel('$^{2}$log fold change')
    plt.axvline(x=0, color='black', linewidth=1.5)
    if geneSubsets==1:
        plt.axhline(y=len(pd.unique(dfg1['GeneSymbol']))+0.5, color='black', linewidth=1.5)
        plt.axhline(y=len(pd.unique(dfg1['GeneSymbol']))+len(pd.unique(dfg2['GeneSymbol']))+0.5, color='black', linewidth=1.5)   
    for a in np.arange(len(genes)-1)+1:
        plt.axhline(y=a+0.5, color='black', linewidth=0.5)

    plt.subplot2grid((5,24), (1,i*8+6), colspan=1, rowspan=4)
    h=plt.barh(dfrra['nr'], dfrra['10lrra'], align='center', color=dfrra['color'])
    plt.ylim(len(genes)+0.5,0.5)
    if rrascore==0:
        plt.text(0,0,'-$^{10}$log RRA$_{dep}$ score', fontsize=8)
    elif rrascore==1:
        plt.text(0,0,'-$^{10}$log RRA$_{enr}$ score', fontsize=8)
    plt.axis('off') 
    plt.text(dfrra['10lrra'].max(),len(genes)-0.5, 'red: FDR < 0.1', rotation=-90) 
    i+=1

plt.subplots_adjust(hspace=0)
plt.subplots_adjust(wspace=0)
plt.savefig(folder+'/GuideFCPlotTest.pdf', bbox_inches='tight', transparent=True)
