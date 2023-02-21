# -*- coding: utf-8 -*-
"""
If one or both drug are not toxic you cannot use Chou-Talalay, but you can use Bliss. However, due to that most drug-concentration curves are S-curves, 
you can also find Biss-synergy with combinations of the same drug.
Combine 2 drugs in different concentrations in a matrix format, normalize the control (no drugs) to 100%, and save it to a CSV file. 
Omit the drug names and concentration from the CSV file, but fill them in the variable here below. The number of concentrations of drug 1 and 2 should 
match with the number of columns and row in the CSV file. 
You can test this script with BlissSynergyCell1.csv and BlissSynergyCell2.csv.
Author: M.F.M. de Rooij PhD, Amsterdam UMC, Spaargaren Lab, 2023, info: m.f.derooij@amsterdamumc.nl
"""
import os
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
#####################################################################################################################################
# VARIABLES
    
# Folder
os.chdir("H:/BioWin/Bliss/")

# Find csv files automatically (Name the files after the cells)
Cells = []
for file in glob.glob("*.csv"):
    Cells.append(file) 

# Drug names and concentrations:
# Drug1 = Columns
Drug1 = "Drug1 ($\mu$M)"
Conc1 = [0,0.5,1,2,4]

# Drug2 = Rows
Drug2 = "Drug2 ($\mu$M)"
Conc2 = [0,1,2,4,8]

#####################################################################################################################################
for nCell,Cell in enumerate(Cells):
    CellName = Cell.split(".csv")[0]
    mObs = np.genfromtxt(Cell, delimiter=",")/100
    
    mExp = mObs[0,]*mObs[:,0][0:, np.newaxis]
    mdBliss = mExp-mObs
    mrBliss=mObs/mExp
      
    plt.figure(1, figsize=(7*len(Cells),25))

    colors = LinearSegmentedColormap.from_list("Custom", ["red", "darkblue"], N=100) 
    
    plt.subplot2grid((5,len(Cells)),(0,nCell))
    ax = sns.heatmap((mObs*100).round(0), annot=True, fmt='g', annot_kws={"size": 20}, square=True, vmin=0, vmax=100, cmap=colors,
                     cbar_kws={"shrink": 0.5, "ticks": [0,25,50,75,100], "extend": "both"}, linewidths=2, linecolor="white", 
                     xticklabels=Conc1, yticklabels=Conc2)
    plt.title(CellName,fontsize=25, fontweight="bold")
    plt.xlabel(Drug1, fontsize=18)
    plt.xticks(fontsize=16)
    plt.ylabel(Drug2, fontsize=18)
    plt.yticks(fontsize=16)
    cbar = ax.collections[0].colorbar
    cbar.set_label("Viability (%)", size=16, rotation=-90, labelpad=25)
    
     
    colors = LinearSegmentedColormap.from_list("Custom", ["red", "lightgray", "lime"], N=100)
    
    plt.subplot2grid((5,len(Cells)),(1,nCell))
    ax = sns.heatmap((mdBliss*100).round(0), annot=True, fmt='g', annot_kws={"size": 20}, square=True, vmin=-25, center=0, vmax=50, 
                     cmap=colors, cbar_kws={"shrink": 0.5, "ticks": [-25,0,25,50], "extend": "both"}, linewidths=2, 
                     linecolor="white", xticklabels=Conc1, yticklabels=Conc2)
    plt.title(CellName,fontsize=25, fontweight="bold")
    plt.xlabel(Drug1, fontsize=18)
    plt.xticks(fontsize=16)
    plt.ylabel(Drug2, fontsize=18)
    plt.yticks(fontsize=16)
    cbar = ax.collections[0].colorbar
    cbar.set_label("$\Delta$Bliss", size=16, rotation=-90, labelpad=25)
    
    plt.subplot2grid((5,len(Cells)),(2,nCell))
    ax = sns.heatmap((mrBliss*100).round(0), annot=True, fmt='g', annot_kws={"size": 20}, square=True, vmin=0, center=100, vmax=200, 
                     cmap=colors, cbar_kws={"shrink": 0.5, "ticks": [0,100,200], "extend": "both"}, linewidths=2, 
                     linecolor="white", xticklabels=Conc1, yticklabels=Conc2)
    plt.title(CellName,fontsize=25, fontweight="bold")
    plt.xlabel(Drug1, fontsize=18)
    plt.xticks(fontsize=16)
    plt.ylabel(Drug2, fontsize=18)
    plt.yticks(fontsize=16)
    cbar = ax.collections[0].colorbar
    cbar.set_label("Relative Bliss", size=16, rotation=-90, labelpad=25)

    
    def colorFader(c1,c2,mix=0): # make a colorrange
        c1=np.array(mpl.colors.to_rgb(c1))
        c2=np.array(mpl.colors.to_rgb(c2))
        return mpl.colors.to_hex((1-mix)*c1 + mix*c2)
    
    colors= list(["red"])
    for i in range(mrBliss.shape[0]-1):
        colors.append(colorFader('lightblue','darkblue', i/(mrBliss.shape[0]-2)))
        
    plt.subplot2grid((5,len(Cells)),(3,nCell))
    for i in range(mrBliss.shape[0]):
        plt.plot(range(mrBliss.shape[1]), mrBliss[i,:]*100, color=colors[i])
    plt.title(CellName,fontsize=25, fontweight="bold")
    plt.xlabel(Drug1, fontsize=18)
    plt.xticks(ticks=range(mrBliss.shape[1]), labels=Conc1, fontsize=16)
    plt.ylabel("Relative Viability (%)", fontsize=18)
    plt.yticks(fontsize=16)
    plt.ylim(ymin=0)
    plt.legend(Conc2, title=Drug2, bbox_to_anchor=(1.04, 0.5), loc="center left")
    
    colors= list(["red"])
    for i in range(mrBliss.shape[1]-1):
        colors.append(colorFader('lightblue','darkblue', i/(mrBliss.shape[1]-2)))
    
    plt.subplot2grid((5,len(Cells)),(4,nCell))
    for i in range(mrBliss.shape[1]):
        plt.plot(range(mrBliss.shape[0]), mrBliss[:,i]*100, color=colors[i])
    plt.title(CellName,fontsize=25, fontweight="bold")
    plt.xlabel(Drug2, fontsize=18)
    plt.xticks(ticks=range(mrBliss.shape[0]), labels=Conc2, fontsize=16)
    plt.ylabel("Relative Viability (%)", fontsize=18)
    plt.yticks(fontsize=16)
    plt.ylim(ymin=0)
    plt.legend(Conc1, title=Drug1, bbox_to_anchor=(1.04, 0.5), loc="center left")
    
plt.tight_layout()
plt.savefig("BlissSynergyPY.pdf")