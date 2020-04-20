# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 09:26:48 2018

@author: Martin F.M. de Rooij, PhD
"""

# After performing the DESeq2 script sucessfully, you can produce publishing-grade MA plots 
##########################################################################################################
#                                                 SETTINGS

# Copy-paste the required folder (use / instead of \ )
folder = 'H:/BioWin/Namalwa'

# Which gene:
gene = 'BTK'

# Save plots as PDF (0=no, 1=yes)
save = 1

# Show significant guides as triangles (0=no, 1=yes)
guidesignificant = 1

# Titles:
maintitle = 'NAMALWA'
xtitle = '$^{10}$log average read counts (input)'
ytitle = '$^{2}$log fold change (CXCL12/input)'
#0 and 1 subscript: '$^{2}$log fold change (t$_1$/t$_0$)'
#greek alpha letter: r'$^{2}$log fold change ($\alpha$IgM/PMA)'   ->include the upstream r for greek letters 

# X- and Y-axis limits:
xmin=0
xmax=4.2
xticks=1

ymin=-1.25
ymax=1
yticks=0.25

# Font size of X- and Y-axis titles (the main title is 10 bigger)
fs = 15

# Colors (see https://matplotlib.org/1.5.1/examples/color/named_colors.html)
call = 'green'
cpos = 'firebrick'
cneg = 'blue'
chit = 'black'
##########################################################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

files = glob.glob(folder+'/*Guides.csv')
for file in files:
    df = pd.read_csv(file, sep=',')

    df['GeneSymbol'] = df['GeneSymbol'].str.upper()
    df['l10rc'] = np.log10(df['BaseMeanA']+1)
    df['l2fc'] = np.log2(df['FoldChange'])
    df['l2fc'].fillna(0, inplace=True)
    df['l2fc'].replace([np.inf, -np.inf], 0, inplace=True)
    df.fillna(1, inplace=True) 
       
    df = df[df['Type']=='x'].append(df[df['Type']!='x'].sample(frac=1, random_state=10))
    
    #MA plots
    plt.figure(1, figsize=(10,10))
    ax = plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=2)   
           
    for i, row in df.iterrows():    
        label_pos_x = row['l10rc']
        label_pos_y = row['l2fc']
        # create and add the dots and text annotation to the scatterplot
        plt.scatter(
            label_pos_x,
            label_pos_y,
            color=cpos if row['Type']=='p' else  cneg if row['Type']=='n' else call,
            marker = "v" if row['padj'] < 0.1 and row['l2fc'] < 0 and guidesignificant==1 else "^" if row['padj'] < 0.1 and row['l2fc'] > 0 and guidesignificant else "o",
            s=25
        )   
    for i, row in df[df['GeneSymbol']==gene].iterrows():    
        label_pos_x = row['l10rc']
        label_pos_y = row['l2fc']
        plt.scatter(
            label_pos_x,
            label_pos_y,
            color=chit, 
            marker = "v" if row['padj'] < 0.1 and row['l2fc'] < 0 and guidesignificant==1 else "^" if row['padj'] < 0.1 and row['l2fc'] > 0 and guidesignificant else "o",
            s=50
        ) 

    plt.xticks(np.arange(xmin, xmax, xticks))
    plt.yticks(np.arange(ymin, ymax, yticks))
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    
    plt.xlabel(xtitle, fontsize=fs)
    plt.ylabel(ytitle, fontsize=fs)
    plt.title(maintitle, fontsize=fs+10, y=1.05)
    
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    
    from matplotlib.lines import Line2D
    if guidesignificant==1:
        legend_elements = [Line2D([0], [0], marker='o', markerfacecolor=call, color=(0,0,0,0), label='All guides', markersize=fs-5),
                       Line2D([0], [0], marker='o', markerfacecolor=cpos, color=(0,0,0,0), label='Essentials', markersize=fs-5),
                       Line2D([0], [0], marker='o', markerfacecolor=cneg, color=(0,0,0,0), label='Non-essentials', markersize=fs-5),
                       Line2D([0], [0], marker='o', markerfacecolor=chit, color=(0,0,0,0), label=gene, markersize=fs-5),
                       Line2D([0], [0], marker='^', markerfacecolor='white', color=(0,0,0,0), label='Significantly enriched', markersize=fs-5),
                       Line2D([0], [0], marker='v', markerfacecolor='white', color=(0,0,0,0), label='Significantly depleted', markersize=fs-5)]
    elif guidesignificant==0: 
        legend_elements = [Line2D([0], [0], marker='o', markerfacecolor=call, color=(0,0,0,0), label='All guides', markersize=fs-5),
                       Line2D([0], [0], marker='o', markerfacecolor=cpos, color=(0,0,0,0), label='Essentials', markersize=fs-5),
                       Line2D([0], [0], marker='o', markerfacecolor=cneg, color=(0,0,0,0), label='Non-essentials', markersize=fs-5),
                       Line2D([0], [0], marker='o', markerfacecolor=chit, color=(0,0,0,0), label=gene, markersize=fs-5)]
        
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.6, 1), fontsize=fs, frameon=False)
          
    plt.axhline(0, color='black', linewidth=1)
    plt.subplot2grid((3,3), (0,2), rowspan=2)
    sns.distplot(df['l2fc'], hist=False, kde_kws = {'shade': True, 'linewidth': 1}, color=call, vertical=True, axlabel=False)
    sns.distplot(df[df['Type']=='n']['l2fc'], hist=False, kde_kws = {'shade': True, 'linewidth': 1}, color=cneg, vertical=True, axlabel=False)
    sns.distplot(df[df['Type']=='p']['l2fc'], hist=False, kde_kws = {'shade': True, 'linewidth': 1}, color=cpos, vertical=True, axlabel=False)
    sns.distplot(df[df['GeneSymbol']==gene]['l2fc'], hist=False, kde_kws = {'shade': True, 'linewidth': 1}, color=chit, vertical=True, axlabel=False)
    plt.ylim(ymin,ymax)
    plt.axis('off')  
    plt.tight_layout(pad=0, w_pad=-24, h_pad=0)
    plt.subplots_adjust(wspace=0)
    if save==0:
        plt.show()
    elif save == 1:    
        plt.savefig(file+'MAplot'+gene+'.pdf', bbox_inches='tight', transparent=True)
    plt.close()
