# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:46:40 2019

@author: Martin F.M. de Rooij, PhD
"""

# After performing the DESeq2 script sucessfully, you can produce publishing-grade Volcano plots 
##########################################################################################################
#                                                 SETTINGS

# Copy-paste the required folder (use / instead of \ )
folder = "H:/BioWin/Screens/Namalwa"

# Which genes to highlight:
gene = ['BTK', 'SYK', 'PIK3R1']

# Save plots as PDF (0=no, 1=yes)
save = 1

# Show significant genes as triangles (0=no, 1=yes)
genesignificant = 1

# Titles:
maintitle = 'NAMALWA'

x1title = '$^{2}$log fold change (PMA/input)'
x2title = r'$^{2}$log fold change ($\alpha$IgM/input)'
x3title = r'$^{2}$log fold change ($\alpha$IgM/PMA)'

#0 and 1 subscript: '$^{2}$log fold change (t$_1$/t$_0$)'
#greek alpha letter: r'$^{2}$log fold change ($\alpha$IgM/PMA)'   -> include the upstream r for greek letters 

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
import glob
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

n = 0
files = glob.glob(folder+'/*Genes.csv')
for file in files:
    n+=1
    dfg = pd.read_csv(file, sep=',')
    
    dfg['GeneSymbol'] = dfg['GeneSymbol'].str.upper()
    dfg['l2mfc'] = np.log2(dfg['MedianFoldChange'])
    dfg['l10rdep'] = np.log10(dfg['rhoDepleted']+10**-15)*-1
    dfg['l10renr'] = np.log10(dfg['rhoEnriched']+10**-15)*-1
    dfg = dfg[dfg['Type']=='x'].append(dfg[dfg['Type']!='x'].sample(frac=1, random_state=10))
 
    dfgene = dfg[dfg['GeneSymbol'].isin(gene)]
    
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10,5))
    # remove the underlying axes
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05)
    gs = axs[0].get_gridspec()
    axbig = fig.add_subplot(gs[0, 0:2])

    ymin=-1
    ymax=16
    xmin=min(dfg['l2mfc']-0.5)
    xmax=max(dfg['l2mfc']+0.5)
    
    if n==1:
        xtitle = x1title
    elif n==2:
        xtitle = x2title
    elif n==3:
        xtitle= x3title
                    
    axbig.set_xticks([])
    axbig.set_yticks([])
    axbig.spines['bottom'].set_linewidth(0)
    axbig.spines['top'].set_linewidth(0)
    axbig.spines['left'].set_linewidth(0)
    axbig.spines['right'].set_linewidth(0)
    axbig.text(0.5, 1.1, maintitle, fontsize=fs+10, ha='center')
    axbig.text(0.5, -0.25, xtitle, fontsize=fs, ha='center')
    axbig.set_facecolor('None')

    axs[0].set_xlim(xmin,xmax)
    axs[0].set_ylim(ymin,ymax)
    axs[0].set_ylabel(r'- $^{10}$log $\rho$ ($\alpha$RRA$_{depleted}$)', fontsize=fs)
    axs[0].axvline(0, color='black', linewidth=1)  
    axs[0].spines['bottom'].set_linewidth(1.5)
    axs[0].spines['top'].set_linewidth(1.5)
    axs[0].spines['left'].set_linewidth(1.5)
    axs[0].spines['right'].set_linewidth(1.5)                            
    for i, row in dfg.iterrows():    
        label_pos_x = row['l2mfc']
        label_pos_y = row['l10rdep']
        # create and add the dots and text annotation to the scatterplot
        axs[0].scatter(
           label_pos_x,
           label_pos_y,
           color=cpos if row['Type']=='p' else  cneg if row['Type']=='n' else call,
           marker = "v" if row['fdrDepleted'] < 0.1 and genesignificant==1 else "^" if row['fdrEnriched'] < 0.1 and genesignificant==1 else "o",
           s=25
        )    

    for i, row in dfgene.iterrows(): 
        label = row['GeneSymbol']
        label_pos_x = row['l2mfc']
        label_pos_y = row['l10rdep']
        # create and add the dots and text annotation to the scatterplot
        axs[0].scatter(
            label_pos_x,
            label_pos_y,
            color=chit,
            marker = "v" if row['fdrDepleted'] < 0.1 and genesignificant==1 else "^" if row['fdrEnriched'] < 0.1 and genesignificant==1 else "o" ,
            s=25
        ) 
        axs[0].text(
            label_pos_x+0.01,
            label_pos_y+0.3,
            label,
            color=chit,
            rotation=10,
            fontsize=fs
        )

    axs[1].yaxis.tick_right()
    axs[1].yaxis.set_label_position('right')
    axs[1].set_xlim(xmin,xmax)
    axs[1].set_ylim(ymin,ymax)
    axs[1].set_ylabel(r'- $^{10}$log $\rho$ ($\alpha$RRA$_{enriched}$)', rotation=-90, fontsize=fs, labelpad=20)
    axs[1].axvline(0, color='black', linewidth=1)
    axs[1].spines['bottom'].set_linewidth(1.5)
    axs[1].spines['top'].set_linewidth(1.5)
    axs[1].spines['left'].set_linewidth(1.5)
    axs[1].spines['right'].set_linewidth(1.5) 
    for i, row in dfg.iterrows():    
            label_pos_x = row['l2mfc']
            label_pos_y = row['l10renr']
            # create and add the dots and text annotation to the scatterplot
            axs[1].scatter(
                label_pos_x,
                label_pos_y,
               color=cpos if row['Type']=='p' else  cneg if row['Type']=='n' else call,
               marker = "v" if row['fdrDepleted'] < 0.1 and genesignificant==1 else "^" if row['fdrEnriched'] < 0.1 and genesignificant==1 else "o",
               s=25
            )    
    for i, row in dfgene.iterrows(): 
        label = row['GeneSymbol']
        label_pos_x = row['l2mfc']
        label_pos_y = row['l10renr']
        # create and add the dots and text annotation to the scatterplot
        axs[1].scatter(
            label_pos_x,
            label_pos_y,
            color=chit,
            marker = "v" if row['fdrDepleted'] < 0.1 and genesignificant==1 else "^" if row['fdrEnriched'] < 0.1 and genesignificant==1 else "o",
            s=25
        ) 
        axs[1].text(
            label_pos_x+0.01,
            label_pos_y+0.3,
            label,
            color=chit,
            rotation=10,
            fontsize=fs
        )

    if save==0:
        plt.show()
    elif save == 1:    
        plt.savefig(file+'Volcano'+maintitle+'.pdf', bbox_inches='tight', transparent=True)
    plt.close()
