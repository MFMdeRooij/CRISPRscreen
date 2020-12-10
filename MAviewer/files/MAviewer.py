#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 09:26:48 2018

@author: Martin F.M. de Rooij, PhD
"""
###################################################################################################
# Use MAviewer:
# Open with the sh (Linux) or bat (Windows) file.
# To add new screens: copy-paste the CRISPRScreenAnalysis.R output to the files/DataCRISPR03New folder. For one screen make one folder. If your screen is a synthetic lethal screen with 3 
# comparisons, make 3 folders. Extent the index at line 459 , you can rename the screen experiment at line 508/595, and 600-611.
# When you have more than 12 screens, you can add an extra page, like page 2 and 3, add also the new page to line 482.
###################################################################################################
#General
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
import sys
import glob
import webbrowser
from tkinter import messagebox
import warnings
warnings.filterwarnings('ignore')
###############################################################################
#Option1
def viewGenes(headDir, dirs, gene, tumortype=0, numberCor=10, MAorVul=1, scale=1, legend=1, guideID=1, save=1, correlategene=1, log=1):
    gene = gene.replace(" ", "")
    #View genes in all MA plots
    if re.match('CRISPR', headDir):
        #View genes in all MA plots of CRISPR screens
        nameGene = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[0]+'/*Genes.csv')[0]           
        dfgenes = pd.read_csv(nameGene, sep=',')
        dfgenes['GeneSymbol'] = dfgenes['GeneSymbol'].str.upper()
        present = any(dfgenes.GeneSymbol==gene)
        if present:
            rep=0
            if MAorVul==1:
                #MA plots
                for cell in dirs:
                    rep+=1
                    nameGuide = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[rep-1]+'/*Guides.csv')[0]
                    df = pd.read_csv(nameGuide, sep=',')
                    nameGene = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[rep-1]+'/*Genes.csv')[0]
                    dfg = pd.read_csv(nameGene, sep=',')
                    dfg.fillna(1, inplace=True) 
                    
                    df['GeneSymbol'] = df['GeneSymbol'].str.upper()
                    dfg['GeneSymbol'] = dfg['GeneSymbol'].str.upper()
                    df['l10rc'] = np.log10(df['BaseMeanA']+1)
                    df['l2fc'] = np.log2(df['FoldChange'])
                    df['l2fc'].fillna(0, inplace=True)
                    df['l2fc'].replace([np.inf, -np.inf], 0, inplace=True)
                    df.fillna(1, inplace=True) 
    
                    dfpos = df[df['Type']=='p']
                    dfneg = df[df['Type']=='n']
                    dfhit = df[df['GeneSymbol'] == gene]
            
                    padjD = dfg[dfg['GeneSymbol'] == gene]['fdrDepleted'].min()
                    padjD1  = round(float(padjD),3)
                    padjD2 = str(padjD1)
                    if padjD1 < 0.001:
                        padjD2 = '<0.001'               
                    padjE = dfg[dfg['GeneSymbol'] == gene]['fdrEnriched'].min()
                    padjE1  = round(float(padjE),3)
                    padjE2 = str(padjE1)
                    if padjE1 < 0.001:
                        padjE2 = '<0.001'
             
                    #MA plots
                    plt.figure(1, figsize=(30,12))
                    plt.subplot2grid((3,16), (0,(rep-1)*4) if rep < 5 else (1,(rep-5)*4) if rep < 9 else (2,(rep-9)*4), colspan=2)   
                    plt.scatter(df['l10rc'],df['l2fc'],c='green', s=3)
                    plt.scatter(dfpos['l10rc'],dfpos['l2fc'],c='orangered', s=2)
                    plt.scatter(dfneg['l10rc'],dfneg['l2fc'],c='cornflowerblue', s=2)
                    plt.scatter(dfhit['l10rc'],dfhit['l2fc'],c='black', s=10)         
                    if scale==0:
                        xmin=min(df['l10rc'])
                        xmax=max(df['l10rc'])
                        ymin=min(df['l2fc'])
                        ymax=max(df['l2fc'])
                    if scale==1:
                        xmin=-1
                        xmax=5
                        ymin=-7
                        ymax=3
                    xmin=xmin-2
                    xmax=xmax+1
                    ymin=ymin-1
                    ymax=ymax+1
                    plt.xlim(xmin,xmax)
                    plt.ylim(ymin,ymax)
                    plt.xlabel('$^{10}$log read counts (t$_0$)', fontsize=10)
                    plt.ylabel('$^{2}$log fold change (t$_1$/t$_0$)', fontsize=10)
                    plt.title(cell.split( "csv")[0], fontsize=12)
                    if guideID==1:                       
                        for i, row in dfhit.iterrows():   
                            label = row['Guide']
                            label_pos_x = row['l10rc']
                            label_pos_y = row['l2fc']
                            plt.text(
                                label_pos_x+0.01,
                                label_pos_y+0.01,
                                label,
                                color='black',
                                fontsize=6
                                )
                    if legend==1:
                        plt.legend(('All Guides', 'Essentials', 'Non-Essentials', gene), loc=3)
                    plt.text(xmin+(xmax-xmin)/2, ymax-0.05*(ymax-ymin), '$\itP$$_{adj}$ = '+padjD2+ r' ($\alpha$RRA$_{dep}$) / $\itP$$_{adj}$ = '+padjE2+ r' ($\alpha$RRA$_{enr}$)', fontsize=6, horizontalalignment='center')
                    plt.axhline(0, color='black', linewidth=1)
                    plt.subplot2grid((3,16), (0,(rep-1)*4+2) if rep < 5 else (1,(rep-5)*4+2) if rep < 9 else (2,(rep-9)*4+2))
                    sns.kdeplot(y=df['l2fc'], linewidth=1, fill=True, color='green')
                    sns.kdeplot(y=dfpos['l2fc'], linewidth=1, fill=True, color='red')
                    sns.kdeplot(y=dfneg['l2fc'], linewidth=1, fill=True, color='blue')
                    sns.kdeplot(y=dfhit['l2fc'], linewidth=1, fill=True, color='black')
                    plt.ylim(ymin,ymax)
                    plt.axis('off')
                plt.tight_layout(pad=0, w_pad=-4, h_pad=3, rect=[0.05, 0.05, 0.95, 0.95])
                plt.subplots_adjust(wspace=0)
            
            if MAorVul==2:
                # Vulcano plots

                fig, axs = plt.subplots(ncols=12, nrows=3, figsize=(30,12))
                # remove the underlying axes
                plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.7)
                for ax in axs[0:, 2]:
                    ax.remove()
                for ax in axs[0:, 5]:
                    ax.remove()
                for ax in axs[0:, 8]:
                    ax.remove()
                for ax in axs[0:, 11]:
                    ax.remove()
                   
                gs = axs[0,0].get_gridspec()
                for cell in dirs:
                    rep+=1
                    nameGene = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[rep-1]+'/*Genes.csv')[0]
                    dfg = pd.read_csv(nameGene, sep=',')
                    dfg.fillna(1, inplace=True) 
                    dfg['l2mfc'] = np.log2(dfg['MedianFoldChange'])
                    dfg['l10rdep'] = np.log10(dfg['rhoDepleted']+10**-15)*-1
                    dfg['l10renr'] = np.log10(dfg['rhoEnriched']+10**-15)*-1
                    
                    dfhitdep = dfg[dfg['fdrDepleted'] < 0.1]
                    dfhitenr = dfg[dfg['fdrEnriched'] < 0.1]
                    
                    dfpos = dfg[dfg['Type']=='p']
                    dfneg = dfg[dfg['Type']=='n']
                    
                    ymin=-1
                    ymax=16
                    xmin=min(dfg['l2mfc']-0.5)
                    xmax=max(dfg['l2mfc']+0.5)
           
                    # Vulcano plots
                    
                    # Empty plot for title
                    axbig = fig.add_subplot(gs[0, (rep-1)*3:(rep-1)*3+2]) if rep < 5 else fig.add_subplot(gs[1, (rep-5)*3:(rep-5)*3+2]) if rep < 9 else fig.add_subplot(gs[2, (rep-9)*3:(rep-9)*3+2])
                        
                    axbig.set_xticks([])
                    axbig.set_yticks([])
                    axbig.spines['bottom'].set_linewidth(0)
                    axbig.spines['top'].set_linewidth(0)
                    axbig.spines['left'].set_linewidth(0)
                    axbig.spines['right'].set_linewidth(0)
                    axbig.text(0.5, 1.1, cell.split( "csv")[0], fontsize=12, ha='center')
                    axbig.text(0.5, -0.25, '$^{2}$log median fold change (t$_1$/t$_0$)', fontsize=10, ha='center')
                    axbig.set_facecolor('None')
                                        
                    # Left vulcano plot
                    xcoor = (rep-1)*3 if rep < 5 else (rep-5)*3 if rep < 9 else (rep-9)*3
                    ycoor = 0 if rep < 5 else 1 if rep < 9 else 2
                                        
                    axs[ycoor,xcoor].scatter(dfg['l2mfc'],dfg['l10rdep'],c='green', s=3)
                    axs[ycoor,xcoor].scatter(dfhitdep['l2mfc'],dfhitdep['l10rdep'], facecolors='none', edgecolors='black', s=10, alpha=0.5) 
                    axs[ycoor,xcoor].scatter(dfpos['l2mfc'],dfpos['l10rdep'],c='orangered', s=1)
                    axs[ycoor,xcoor].scatter(dfneg['l2mfc'],dfneg['l10rdep'],c='cornflowerblue', s=1)

                    axs[ycoor,xcoor].set_xlim(xmin,xmax)
                    axs[ycoor,xcoor].set_ylim(ymin,ymax)
                    axs[ycoor,xcoor].set_ylabel(r'- $^{10}$log $\rho$ ($\alpha$RRA$_{dep}$)', fontsize=10)
                                           
                    df_geneX = dfg.loc[dfg['GeneSymbol']==gene]
                    for i, row in df_geneX.iterrows():    
                        label = row['GeneSymbol']
                        label_pos_x = row['l2mfc']
                        label_pos_y = row['l10rdep']
                        # create and add the dots and text annotation to the scatterplot
                        axs[ycoor,xcoor].scatter(
                            label_pos_x,
                            label_pos_y,
                            c='yellow',
                            s=12
                        ) 
                        axs[ycoor,xcoor].text(
                            label_pos_x+0.05,
                            label_pos_y+0.05,
                            label,
                            color='black',
                            weight ='bold',
                            rotation=10,
                            fontsize=10
                        )       
                    
                    # Right vulcano plot
                    xcoor += 1        

                    axs[ycoor,xcoor].scatter(dfg['l2mfc'],dfg['l10renr'],c='green', s=3)
                    axs[ycoor,xcoor].scatter(dfhitenr['l2mfc'],dfhitenr['l10renr'], facecolors='none', edgecolors='black', s=10, alpha=0.5) 
                    axs[ycoor,xcoor].scatter(dfpos['l2mfc'],dfpos['l10renr'],c='orangered', s=1)
                    axs[ycoor,xcoor].scatter(dfneg['l2mfc'],dfneg['l10renr'],c='cornflowerblue', s=1)

                    axs[ycoor,xcoor].yaxis.tick_right()
                    axs[ycoor,xcoor].yaxis.set_label_position('right')
                    axs[ycoor,xcoor].set_xlim(xmin,xmax)
                    axs[ycoor,xcoor].set_ylim(ymin,ymax)
                    axs[ycoor,xcoor].set_ylabel(r'- $^{10}$log $\rho$ ($\alpha$RRA$_{enr}$)', rotation=-90, fontsize=10, labelpad=20)
                    
                    df_geneX = dfg.loc[dfg['GeneSymbol']==gene]
                    for i, row in df_geneX.iterrows():    
                        label = row['GeneSymbol']
                        label_pos_x = row['l2mfc']
                        label_pos_y = row['l10renr']
                        # create and add the dots and text annotation to the scatterplot
                        axs[ycoor,xcoor].scatter(
                            label_pos_x,
                            label_pos_y,
                            c='yellow',
                            s=12
                        )         
                        axs[ycoor,xcoor].text(
                            label_pos_x+0.05,
                            label_pos_y+0.05,
                            label,
                            color='black',
                            weight ='bold',
                            rotation=10,
                            fontsize=10
                        ) 
                        
                # remove empty plots 
                maxplot=12-len(dirs)
                if maxplot>0:
                    for repm in range(len(dirs),12):
                        xcoor = (repm)*3 if repm < 4 else (repm-4)*3 if repm < 8 else (repm-8)*3
                        ycoor = 0 if repm < 4 else 1 if repm < 8 else 2   
                        axs[ycoor,xcoor].axis('off')
                        axs[ycoor,xcoor+1].axis('off')

            
            if save==0:
                plt.show()
            elif save==1:
                if MAorVul==1:
                    plt.savefig(cwd+'/SavedMAplots/'+headDir+gene+'_MA.pdf', bbox_inches='tight')
                if MAorVul==2:
                    plt.savefig(cwd+'/SavedMAplots/'+ 'MA'+headDir+gene+'_VP.pdf', bbox_inches='tight')
                plt.close()
        else:
            messagebox.showwarning("Warning", "This gene annotation is not present in the used CRISPR library.\nCheck your spelling, or check for official HGNC gene symbols.")    


    elif re.match('RNAseq', headDir):
        nameGene = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[0]+'/*NorTPM.csv')[0]           
        dfgenes = pd.read_csv(nameGene, sep=',')
        dfgenes = dfgenes.drop_duplicates()
        dfgenes.loc[dfgenes['hgnc_symbol'].isnull(), 'hgnc_symbol'] = dfgenes.loc[dfgenes['hgnc_symbol'].isnull(), 'ensembl_gene_id']  
        dfgenes['hgnc_symbol'] = dfgenes['hgnc_symbol'].str.upper()
        dfgenes['hgnc_symbol'][dfgenes['hgnc_symbol'].duplicated()] = dfgenes['hgnc_symbol'][dfgenes['hgnc_symbol'].duplicated()] + "_2"

        # Tumor colors
        dftype = pd.DataFrame({'cell':dfgenes.columns[2:len(dfgenes.columns)], 'type': tumortype}) 
        my_color=["black", "blue"]*100
        my_palette = dict(zip(dftype['type'].unique(), my_color[0:len(dftype['type'].unique())]))    
        tumor_colors = dftype['type'][0:len(dftype['type'])].map(my_palette)  

        genes = gene.split(",")
        numberofgenes = len(genes)
        
        if numberofgenes==1:
            present = any(dfgenes.hgnc_symbol==gene)          
            if present:
                if correlategene==0:
                    ax=plt.gca()
                    plt.barh(np.arange(len(dfgenes.columns)-2), dfgenes[dfgenes['hgnc_symbol']==gene].iloc[0,2:len(dfgenes.columns)], align='center')
                    plt.yticks(np.arange(len(dfgenes.columns)-2), dfgenes.columns[2:len(dfgenes.columns)])
                    plt.gca().invert_yaxis()
                    plt.ylabel('Cell line', fontsize=10)    
                    plt.xlabel('Gene expression (normalized TPM)', fontsize=10) 
                    plt.title(gene)   
                   
                    # Tumor colors
                    for ytick, color in zip(ax.get_yticklabels(), tumor_colors):
                        ytick.set_color(color) 
                    plt.show()

                elif correlategene==1:
                    dfexpr = dfgenes
                    dfexpr.index = dfexpr['hgnc_symbol']
                    dfexpr = dfexpr.drop(["ensembl_gene_id",'hgnc_symbol'], axis=1) 

                    if log==0:
                        dfexpr = np.log(dfexpr+0.1)
                        
                    corr = dfexpr.loc[gene,]
                    rowmeans=dfexpr.mean(axis=1)
                    if log==0:
                        rowmeansSel = rowmeans.loc[rowmeans>-2]
                    else: rowmeansSel = rowmeans.loc[rowmeans>1]
                    dfexpr = dfexpr.reindex(rowmeansSel.index)
                    
                    corrSer = dfexpr.corrwith(corr, axis=1)
                    corrSer = corrSer.sort_values(ascending = False) 
                    dfexprCorr = dfexpr.reindex(corrSer.index)
                    genes1 = list(dfexprCorr.index[0:numberCor+1])
                    genes2 = list(dfexprCorr.index[len(dfexprCorr)-numberCor:len(dfexprCorr)])
                    genes = (genes1+genes2)
                    numberofgenes = len(genes)
                      
            else:
                messagebox.showwarning("Warning", "This gene annotation is not present in the human genome (hg38).\nCheck your spelling, or check for official HGNC gene symbols.")   

        if numberofgenes>1:
            dfexpr = dfgenes[dfgenes['hgnc_symbol'].isin(genes)]
            dfexpr.index = dfexpr['hgnc_symbol']
            dfexpr = dfexpr.reindex(genes)
            dfexpr = dfexpr.drop(["ensembl_gene_id",'hgnc_symbol'], axis=1)           
    
            fig, ax = plt.subplots()
            ax = plt.gca()       
            # Plot the heatmap
            rowmeans=dfexpr.mean(axis=1)
            rowsd=dfexpr.std(axis=1)
            data=dfexpr.subtract(rowmeans, axis=0)
            data=data.divide(rowsd, axis=0)
            im = ax.imshow(data, cmap="coolwarm")     
            # Create colorbar
            cbar = ax.figure.colorbar(im, ax=ax)
            cbar.ax.set_ylabel("Gene expression (Color scale:Z score per gene, values: normalized TPM)", rotation=-90, va="bottom")       
            # We want to show all ticks...
            ax.set_xticks(np.arange(data.shape[1]))
            ax.set_yticks(np.arange(data.shape[0]))
            # ... and label them with the respective list entries.
            fs = 4 if numberofgenes>25 else 6 if numberofgenes>15 else 10
            ax.set_xticklabels(list(dfexpr.columns.values), fontsize=fs)
            ax.set_yticklabels(list(dfexpr.index.values), fontsize=fs)       
            # Let the horizontal axes labeling appear on top.
            ax.tick_params(top=True, bottom=False,
                           labeltop=True, labelbottom=False)       
            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
                     rotation_mode="anchor")       
            # Turn spines off and create white grid.
            for edge, spine in ax.spines.items():
                spine.set_visible(False)       
            ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
            ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
            lw= 1 if numberofgenes>25 else 2 if numberofgenes>15 else 3
            ax.grid(which="minor", color="w", linestyle='-', linewidth=lw)
            ax.tick_params(which="minor", bottom=False, left=False)               
    
            data = dfexpr.values
            # Set default alignment to center, but allow it to be
            # overwritten by textkw.
            kw = dict(horizontalalignment="center",verticalalignment="center")
            #kw.update(textkw)
            # Get the formatter in case a string is supplied
            import matplotlib
            valfmt="{x:.1f}"
            valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)
            # Loop over the data and create a `Text` for each "pixel".
            # Change the text's color depending on the data.
            texts = []
            fs = 0 if numberofgenes>25 else 5 if numberofgenes>15 else 7
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    kw.update(color="black", fontsize = fs)
                    text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                    texts.append(text)       

            plt.ylabel('Gene', fontsize=10)
            ax.xaxis.set_label_position('top')
            plt.xlabel('Cell line', fontsize=10) 

            # Tumor colors
            for xtick, color in zip(ax.get_xticklabels(), tumor_colors):
                xtick.set_color(color) 

            try:
                # Prepare a vector of color mapped to the correlations 
                listCor1 = list(corrSer.values[0:numberCor+1])
                listCor2 = list(corrSer.values[len(dfexprCorr)-numberCor:len(dfexprCorr)])
                listCor = (listCor1+listCor2)  
                listCor = [round(float(elem), 2) for elem in listCor]
                dfexpr['r'] = listCor
            
                def rgb(minimum, maximum, value):
                    minimum, maximum = float(minimum), float(maximum)
                    ratio = 2 * (value-minimum) / (maximum - minimum)
                    r = (int(max(0, 255*(1 - ratio))))/255
                    g = (int(max(0, 255*(ratio - 1))))/255
                    b = 0 #(255 - r - g)/255
                    a = 1
                    return r, g, b, a
            
                calcol = lambda x: rgb(-1,1,x)
                corr_colors=dfexpr['r'].map(calcol)            

                for ytick, color in zip(ax.get_yticklabels(), corr_colors):
                    ytick.set_color(color)

            except NameError:
                pass

            fig.tight_layout()
            plt.show()


    return
###############################################################################
#Option2
def viewGuides(headDir, dirs, rep):
    #View guides in an MA plot
    from matplotlib.widgets import Button
    from matplotlib.text import Annotation
    global colID
    colID=0
    if re.match('CRISPR', headDir):
        #View guides in an MA plots of CRISPR screens
        # Click on MA plot, and view guides
        nameGuide = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[rep-1]+'/*Guides.csv')[0]
        df = pd.read_csv(nameGuide, sep=',')
        df['GeneSymbol'] = df['GeneSymbol'].str.upper()
        df['l10rc'] = np.log10(df['BaseMeanA']+1)
        df['l2fc'] = np.log2(df['FoldChange'])
        df['l2fc'].fillna(0, inplace=True)
        df['l2fc'].replace([np.inf, -np.inf], 0, inplace=True)
        df.fillna(1, inplace=True) 
        df['col'] =  'green'
        df['size'] =  2
        df.loc[df['Type']=='p', 'col'] = 'orangered'
        df.loc[df['Type']=='p', 'size'] = 10
        df.loc[df['Type']=='n', 'col'] = 'cornflowerblue'
        df.loc[df['Type']=='n', 'size'] = 10
        
        axis_values_x = df['l10rc']
        axis_values_y = df['l2fc']
        instances_colors = df['col']
        size = df['size']
        # draw a scatter-plot of the generated values
        fig = plt.figure(figsize=(20, 20))
        ax = plt.subplot() 
        
        # extract the scatterplot drawing in a separate function so we ca re-use the code
        def draw_scatterplot():
            ax.scatter(
                axis_values_x,
                axis_values_y,
                c=instances_colors,
                s=size,
                picker=True
            )
            ax.set_title(dirs[rep-1].split( "csv")[0], fontsize=20)
            ax.set_xlabel('$^{10}$log read counts (t$_0$)', fontsize=15)
            ax.set_ylabel('$^{2}$log fold change (t$_1$/t$_0$)', fontsize=15)
            
        # draw the initial scatterplot
        draw_scatterplot()
    
        # create and add an annotation object (a text label)
        def annotate(axis, text, x, y):
            text_annotation = Annotation(text, xy=(x, y), fontsize=10, xycoords='data')
            axis.add_artist(text_annotation)
        
        # define the behaviour -> what happens when you pick a dot on the scatterplot by clicking close to it
        def onpick(event):
            ind = event.ind
            geneX = df['GeneSymbol'].loc[ind]
            geneX = geneX.tolist() 
            df_geneX = df.loc[df['GeneSymbol'].isin(geneX)]
            global colID
            for i, row in df_geneX.iterrows():    
                color = (('black', 'gray', 'orange', 'purple', 'brown')*(colID//5+1))[colID]
                label = row['Guide']
                label_pos_x = row['l10rc']
                label_pos_y = row['l2fc']
                # create and add the dots and text annotation to the scatterplot
                ax.scatter(
                    label_pos_x,
                    label_pos_y,
                    c=color,
                    s=15
                ) 
                ax.text(
                    label_pos_x+0.01,
                    label_pos_y+0.01,
                    label,
                    color=color
                )
                # force re-draw
                ax.figure.canvas.draw_idle()
            colID+=1
        # connect the click handler function to the scatterplot
        fig.canvas.mpl_connect('pick_event', onpick)
        # create the "clear all" button, and place it somewhere on the screen
        ax_clear_all = plt.axes([0.15, 0.15, 0.05, 0.05])
        button_clear_all = Button(ax_clear_all, 'Clear')
        # define the "clear all" behaviour
        def onclick(event):
            global colID
            colID=0
            # clear all artist object of the scatter plot
            ax.cla()
            # re-populate the scatterplot only with the dots not the labels
            draw_scatterplot()
            #  force re-draw
            ax.figure.canvas.draw_idle()
        # link the event handler function to the click event on the button
        button_clear_all.on_clicked(onclick)
        # initial drawing of the scatterplot
        plt.plot()
        # present the scatterplot
        plt.show()
    return
###############################################################################
#Option2A
def viewHits(headDir, dirs, rep, tumortype=0, number=25, direction=1, rhoFC=1, screenOrRnaseq=1):   
    #View gene top hits
    nameGene = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[rep-1]+'/*Genes.csv')[0]  
    dfg = pd.read_csv(nameGene, sep=',')    
    if re.match('CRISPR', headDir):
        dfg['ml2fc']=np.log2(dfg['MedianFoldChange'])
        dfg['ml2fc'].fillna(0, inplace=True)
        dfg['ml2fc'].replace([np.inf, -np.inf], 0, inplace=True)  
        dfg['col']='black' 
        dfg.loc[dfg['Type']=='p', 'col'] = 'red'
        dfg.loc[dfg['Type']=='n', 'col'] = 'blue' 
    if re.match('Microarray', headDir):
        dfg['ml2fc']=dfg['log2MedianFC']
        dfg['ml2fc'].fillna(0, inplace=True)
        dfg['ml2fc'].replace([np.inf, -np.inf], 0, inplace=True)  
        dfg['col']='black'    
    dfg.fillna(1, inplace=True)  
    if rhoFC==1:
        if direction==1:
            dfg['logrho']=np.log10(dfg['rhoDepleted'])
        elif direction==0:
            dfg['logrho']=np.log10(dfg['rhoEnriched'])
        dfg=dfg.sort_values(by=['logrho'])  
        if screenOrRnaseq==1:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.barh(np.arange(number), dfg['logrho'][0:number], align='center')
            plt.yticks(np.arange(number), dfg['GeneSymbol'][0:number])
            if re.match('CRISPR', headDir):
                for ytick, color in zip(ax.get_yticklabels(), dfg['col'][0:number]):
                    ytick.set_color(color)     
            plt.gca().invert_yaxis()
            plt.ylabel('Gene Symbol', fontsize=10)    
            if direction==1:
                plt.xlabel(r'$^{10}$log $\rho$ ($\alpha$RRA$_{dep}$)', fontsize=10)
                rhocutoff = dfg[dfg['fdrDepleted']>=0.1]['rhoDepleted']
                plt.text(np.log10(min(rhocutoff)), number/2, r'$\itP$$_{adj}$ < 0.1 ($\alpha$RRA$_{dep}$)', rotation=-90, va= 'center', ha='right', color='red' , fontsize=8)
            elif direction==0:
                plt.xlabel(r'$^{10}$log $\rho$ ($\alpha$RRA$_{enr}$)', fontsize=10)
                rhocutoff = dfg[dfg['fdrEnriched']>=0.1]['rhoEnriched']
                plt.text(np.log10(min(rhocutoff)), number/2, r'$\itP$$_{adj}$ < 0.1 ($\alpha$RRA$_{enr}$)',rotation=-90, va= 'center', ha='right', color='red' , fontsize=8)
            plt.axvline(np.log10(min(rhocutoff)), color='red', linewidth=1) 
            plt.title('Gene Top '+str(number)+' ('+dirs[rep-1].split( "csv")[0]+')')   
            plt.show()
    if rhoFC==0:    
        if direction==1:
            dfg=dfg.sort_values(by=['ml2fc'])  
        elif direction==0:
            dfg=dfg.sort_values(by=['ml2fc'], ascending=False)   
        if screenOrRnaseq==1:    
            fig = plt.figure()
            ax = fig.add_subplot(111) 
            plt.barh(np.arange(number), dfg['ml2fc'][0:number], align='center')
            plt.yticks(np.arange(number), dfg['GeneSymbol'][0:number])
            if re.match('CRISPR', headDir):
                for ytick, color in zip(ax.get_yticklabels(), dfg['col'][0:number]):
                    ytick.set_color(color)
            plt.gca().invert_yaxis()
            plt.ylabel('Gene Symbol', fontsize=10)
            plt.xlabel(r'$^{2}$log Median Fold Change', fontsize=10)
            plt.title('Gene Top '+str(number)+' ('+dirs[rep-1].split( "csv")[0]+')') 
            plt.show()
    if screenOrRnaseq==0:
        gene=str(dfg['GeneSymbol'][0:number].to_csv(header=None, index=False).strip('\n').split('\n')).replace("'","").replace("[","").replace("]","").replace("\\r","")
        viewGenes(headDir='RNAseq', dirs=['01_CellLines'], tumortype=tumortype, gene=gene, correlategene=0, log=0) 
      
    return
###############################################################################
#Option2B
def viewScreendetails(headDir, dirs, rep):
    if re.match('CRISPR', headDir):
        import webbrowser 
        screenFile = glob.glob(cwd+'/files/Data'+headDir+'/'+dirs[rep-1]+'/*Plots.pdf')[0]
        webbrowser.open_new(screenFile)
    return   
###############################################################################
#Data sets
cwd=os.getcwd()
def dataSet(cwd, dataset):
    #Use choosen dataset
    def listdirs(folder):
        # Find directories
        if os.path.exists(folder):
             return [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]
        else:
             return []
    global dirs
    global headDir

    # RNAseq tumortypes
    global tumortype
    tumortype=[1,1,1,1,2,3,3,3,3,3,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,7,7,7,7]
      
    # CRISPR Screens
    if dataset==2: 
        order = [3,4,8,9,2,6,7,10,0,1,5]
        headDir = 'CRISPR02Phelan' 
    elif dataset==3:
        order = [0]
        headDir = 'CRISPR03New'

    # RNAseq
    elif dataset==100:
        order = [0]
        headDir = 'RNAseq'  
        
    dirS = listdirs(cwd+"/files/Data"+headDir)
    dirS.sort()
    dirs = [dirS[i] for i in order]
    return(headDir,dirs)

###############################################################################
#Menu
import tkinter as tk                
from tkinter import font  as tkfont 
from tkinter.ttk import Entry
from tkinter import *
colbg='lightskyblue'
class SampleApp(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        container = tk.Frame(self)
        container.pack(side="left", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.frames = {}
        for F in (StartPage, PageOne, PageTwo, PageThree, PageTwenty, PageOneHundred):
            page_name = F.__name__
            frame = F(parent=container, controller=self)
            self.frames[page_name] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame("StartPage")
    def show_frame(self, page_name):
        '''Show a frame for the given page name'''
        frame = self.frames[page_name]
        frame.tkraise()
        frame.configure(background=colbg)

class StartPage(tk.Frame): 
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        tk.Label(self, text='Welcome at MAviewer 2019\n MArtin de Rooij, Amsterdam UMC, Spaargaren Lab', bg='white' , font = tkfont.Font(family='Comic Sans MS', size=15)).pack(side="top", fill="x", pady=20)
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text='MAIN MENU', bg=colbg, font = tkfont.Font(family='Times New Roman', size=20, slant="italic", weight="bold")).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack() 
        tk.Label(self, text='Which dataset do you want to analyze?', bg=colbg, font = tkfont.Font(family='Times New Roman', size=15)).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        global startselected
        startselected = tk.IntVar(value=2)       
        tk.Label(self, text='CRISPR screens:', bg=colbg, font = tkfont.Font(family='Times New Roman', size=15)).pack(anchor = 'w')
        tk.Radiobutton(self,text='Phelan et al Nature 2018 - Staudt (Brunello - DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=2, variable=startselected).pack(anchor = 'w')
        tk.Radiobutton(self,text='New Screens', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=3, variable=startselected).pack(anchor = 'w')      
        
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()      
        tk.Label(self, text='RNAseq:', bg=colbg, font = tkfont.Font(family='Times New Roman', size=15)).pack(anchor = 'w')
        tk.Radiobutton(self,text='Panel of B cell lines (https://www.ncbi.nlm.nih.gov/sra)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=100, variable=startselected).pack(anchor = 'w')
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()

        
        
        
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text="What do you want to visualize?", bg=colbg, font = tkfont.Font(family='Times New Roman', size=15)).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Button(self, text="Genes", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: [controller.show_frame("PageOne" if startselected.get() != 100 else "PageOneHundred"), dataSet(cwd=cwd, dataset=startselected.get())]).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Button(self, text="Guides", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: [clickedGuides(startselected.get()), dataSet(cwd=cwd, dataset=startselected.get())]).pack()       
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=350)).pack()
        tk.Button(self, text="Exit", font=('Times New Roman', '15'), fg='white', bg='red', command=lambda:app.destroy()).pack(side='right')
        
        def clickedGuides(startselected):
           if startselected==2:
               controller.show_frame("PageTwo")
           elif startselected==3:
               controller.show_frame("PageThree")
           elif startselected==100:
               messagebox.showwarning("Warning", "You have to choose Genes.\n(Guides are normally not expressed in cell lines, crazy bitch!)")


#MA or Volcano plots whole datasets               
class PageOne(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text='Which gene do you want to visualize?', bg=colbg, font=('Times New Roman', '15')).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text='Gene:', bg=colbg, font=('Times New Roman', '15')).pack()
        e1 = Entry(self, width=30)
        e1.pack()
        scale=tk.IntVar(value=1)
        legend=tk.IntVar(value=1)
        guideID=tk.IntVar(value=1)       
        MAorVul = tk.IntVar(value=1)
        tk.Radiobutton(self,text='MA plots (guide level)   ', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=1, variable=MAorVul).pack()
        tk.Checkbutton(self, text = "Equal scales    ", font=('Times New Roman', '10'), bg=colbg, onvalue=1, offvalue=0, variable=scale, height=1, width = 20).pack()
        tk.Checkbutton(self, text = "Show legends ", font=('Times New Roman', '10'), bg=colbg, onvalue=1, offvalue=0, variable=legend, height=1, width = 20).pack()
        tk.Checkbutton(self, text = "Show guide ID", font=('Times New Roman', '10'), bg=colbg, onvalue=1, offvalue=0, variable=guideID, height=1, width = 20).pack()
        tk.Radiobutton(self,text='Volcano plot (gene level)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=2, variable=MAorVul).pack()        
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Button(self, text="View Plots", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: viewGenes(headDir=headDir, dirs=dirs, gene=e1.get().upper(), scale=scale.get(), MAorVul=MAorVul.get(), legend=legend.get(), guideID=guideID.get(), save=0)).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        #tk.Button(self, text="Save MA/VOL plots as PDF", font=('Times New Roman', '12'), fg='black', bg='yellow', command=lambda: viewGenes(headDir=headDir, dirs=dirs, gene=e1.get().upper(), scale=scale.get(), MAorVul=MAorVul.get(), legend=legend.get(), guideID=guideID.get(), save=1)).pack() 
        tk.Label(self, text="You can change titles and axis scales using the edit button\n\nYou can save the figure in any format using the save button", bg=colbg, font = tkfont.Font(family='Times New Roman', size=12)).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=50)).pack()
        tk.Button(self, text="Main Menu", font=('Times New Roman', '15'), fg='white', bg="red", command=lambda: controller.show_frame("StartPage")).pack(side='bottom')          

# CRISPR Screens        
class PageTwo(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        tk.Label(self, text='CRISPR screens Phelan et al Nature 2018 - Staudt (Brunello - DLBCL)', bg=colbg, font=('Comic Sans MS', '15')).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text='Which cell do you want to analyze?', bg=colbg, font=('Times New Roman', '15')).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        selected = tk.IntVar(value=1)
        tk.Radiobutton(self,text='HBL1 (ABC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=1, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='HLY1 (ABC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=2, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='TMD8 (ABC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=3, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='U2932 (ABC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=4, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='DOHH2 (GC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=5, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='SUDHL4 (GC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=6, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='SUDHL5 (GC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=7, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='WSUDLCL2 (GC-DLBCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=8, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='ARP1C (MM)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=9, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='DEL (ALCL)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=10, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='KMS11 (MM)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=11, variable=selected).pack(anchor = 'w')
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack() 
        tk.Button(self, text="Explore MA plot", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: viewGuides(headDir=headDir, dirs=dirs, rep=selected.get())).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=30)).pack() 
        tk.Button(self, text="View Gene Top:", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: viewHits(headDir=headDir, dirs=dirs, rep=selected.get(), tumortype=tumortype, number=int(number.get()), direction=direction.get(), rhoFC=rhoFC.get(), screenOrRnaseq=screenOrRnaseq.get())).pack() 
        number = Scale(self, from_=1, to=100, bg='white')
        number.set(25)
        number.pack(anchor=CENTER)
        direction=tk.IntVar(value=1)
        tk.Checkbutton(self, text = "Depletion (on) / Enrichment (off)                        ", font=('Times New Roman', '12'), bg=colbg, onvalue=1, offvalue=0, variable=direction, height=1, width = 35).pack()
        rhoFC=tk.IntVar(value=1)
        tk.Checkbutton(self, text = u"\u03B1RRA (on) / Median Fold Change (off)             ", font=('Times New Roman', '12'), bg=colbg, onvalue=1, offvalue=0, variable=rhoFC, height=1, width = 35).pack()
        screenOrRnaseq=tk.IntVar(value=1)
        tk.Checkbutton(self, text = "Screen Results (on) / RNAseq B Cell Panel (off)", font=('Times New Roman', '12'), bg=colbg, onvalue=1, offvalue=0, variable=screenOrRnaseq, height=1, width = 35).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=30)).pack()       
        tk.Button(self, text="View Screen Details", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: viewScreendetails(headDir=headDir, dirs=dirs, rep=selected.get())).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()       
        tk.Button(self, text="Search Guide Sequence", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: controller.show_frame("PageTwenty")).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Button(self, text="Main Menu", font=('Times New Roman', '15'), fg='white', bg="red", command=lambda: controller.show_frame("StartPage")).pack(side='bottom')    
    
class PageThree(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        tk.Label(self, text='New Screens', bg=colbg, font=('Comic Sans MS', '15')).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text='Which cell do you want to analyze?', bg=colbg, font=('Times New Roman', '15')).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        selected = tk.IntVar(value=1)
        tk.Radiobutton(self,text='Cell 1 (Brunello Kinome Screen)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=1, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 2', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=2, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 3', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=3, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 4', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=4, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 5', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=5, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 6', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=6, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 7', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=7, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 8', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=8, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 9', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=9, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 10', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=10, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 11', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=11, variable=selected).pack(anchor = 'w')
        tk.Radiobutton(self,text='Cell 12', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=12, variable=selected).pack(anchor = 'w')
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack() 
        tk.Button(self, text="Explore MA plot", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: viewGuides(headDir=headDir, dirs=dirs, rep=selected.get())).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=30)).pack() 
        tk.Button(self, text="View Gene Top:", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: viewHits(headDir=headDir, dirs=dirs, rep=selected.get(), tumortype=tumortype, number=int(number.get()), direction=direction.get(), rhoFC=rhoFC.get(), screenOrRnaseq=screenOrRnaseq.get())).pack() 
        number = Scale(self, from_=1, to=100, bg='white')
        number.set(25)
        number.pack(anchor=CENTER)
        direction=tk.IntVar(value=1)
        tk.Checkbutton(self, text = "Depletion (on) / Enrichment (off)                        ", font=('Times New Roman', '12'), bg=colbg, onvalue=1, offvalue=0, variable=direction, height=1, width = 35).pack()
        rhoFC=tk.IntVar(value=1)
        tk.Checkbutton(self, text = u"\u03B1RRA (on) / Median Fold Change (off)             ", font=('Times New Roman', '12'), bg=colbg, onvalue=1, offvalue=0, variable=rhoFC, height=1, width = 35).pack()
        screenOrRnaseq=tk.IntVar(value=1)
        tk.Checkbutton(self, text = "Screen Results (on) / RNAseq B Cell Panel (off)", font=('Times New Roman', '12'), bg=colbg, onvalue=1, offvalue=0, variable=screenOrRnaseq, height=1, width = 35).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=30)).pack()       
        tk.Button(self, text="View Screen Details", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: viewScreendetails(headDir=headDir, dirs=dirs, rep=selected.get())).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()       
        tk.Button(self, text="Search Guide Sequence", font=('Times New Roman', '15'), fg='black', bg='yellow', command=lambda: controller.show_frame("PageTwenty")).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Button(self, text="Main Menu", font=('Times New Roman', '15'), fg='white', bg="red", command=lambda: controller.show_frame("StartPage")).pack(side='bottom')    
    
#Find guide sequence
class PageTwenty(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        tk.Label(self, text='From which gene do you want to search guides sequences?', bg=colbg, font=('Times New Roman', '15')).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack() 
        tk.Label(self, text="Gene:", bg=colbg, font = tkfont.Font(family='Times New Roman', size=15)).pack() 
        e1 = Entry(self, width=10)
        e1.pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack() 
        tk.Button(self, text="Search", font=('Times New Roman', '12'), fg='black', bg='yellow', command=lambda: searchGuideSequence(Tex=Tex, startselected=startselected.get(), gene=e1.get().upper())).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=30)).pack()
        tk.Label(self, text="To order Gibson cloning oligos, copy-paste the right sequences in oligo-ordering", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack() 
        Tex = tk.Text(self, height=9, width=125, bg='cornsilk') 
        Tex.pack()         
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=30)).pack()
        tk.Label(self, text="To find genomic locations, copy-paste this in BLAT", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack() 
        Texblat = tk.Text(self, height=25, width=25, bg='cornsilk') 
        Texblat.pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Button(self, text="Open BLAT", font=('Times New Roman', '12'), fg='black', bg='yellow', command=lambda: webbrowser.open('https://genome.ucsc.edu/cgi-bin/hgBlat?command=start')).pack()   
        tk.Button(self, text="Back", font=('Times New Roman', '15'), fg='white', bg="red", command=lambda:  [clickedGuides(startselected.get()), dataSet(cwd=cwd, dataset=startselected.get())]).pack(side='bottom')       
        
        def clickedGuides(startselected):
           if startselected==2:
               controller.show_frame("PageTwo")
           elif startselected==3:
               controller.show_frame("PageThree")

        def searchGuideSequence(Tex, startselected, gene):
           if startselected==2:
               library = pd.read_csv(cwd+'/files/Libraries/LibraryWholeGenomeBrunelloLentiGuideHGNC.csv', sep=',')
           elif startselected==3:
               library = pd.read_csv(cwd+'/files/Libraries/LibraryKinomeBrunelloLentiGuideHGNC.csv', sep=',')
           if startselected != 11:
               present = any(library['Gene']==gene)
               if present:               
                   seq = library[library['Gene']==gene]
                   pd.options.display.max_colwidth = 100
                   seq['Gibson cloning oligo for LentiGuide'] = 'TATCTTGTGGAAAGGACGAAACACCG'+seq['Sequence']+'GTTTTAGAGCTAGAAATAGCAAGTTAAAA'
                   Tex.delete(1.0, tk.END)
                   Tex.insert(tk.END, seq.to_string(index=False))
                   Tex.pack()
                   seq.index = range(len(seq.index))
                   seqblat = ""
                   for i in range(len(seq.index)):
                       seqblat=seqblat+">"+seq['sgRNA'][i]+'\n'
                       seqblat=seqblat+seq['Sequence'][i]+ 'NGG\n\n'
                   Texblat.delete(1.0, tk.END)
                   Texblat.insert(tk.END, seqblat)
                   Texblat.pack()
                                  
               else:
                   messagebox.showwarning("Warning", "This gene annotation is not present in the used CRISPR library.\nCheck your spelling, or check for official HGNC gene symbols.")  
           else:
                   seq = "No sequences available"
                   Tex.delete(1.0, tk.END)
                   Tex.insert(tk.END, seq)
                   Tex.pack()  
                   Texblat.delete(1.0, tk.END)
                   Texblat.insert(tk.END, seq)
                   Texblat.pack()

# RNAseq       
class PageOneHundred(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text='Which gene(s) do you want to visualize?', bg=colbg, font=('Times New Roman', '15')).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Label(self, text='Gene(s):', bg=colbg, font=('Times New Roman', '15')).pack()
        e1 = Entry(self, width=30)
        e1.pack()
        tk.Label(self, text='(Separate mutiple genes by a comma)', bg=colbg, font=('Times New Roman', '10')).pack()       
        correlategene = tk.IntVar(value=0)
        tk.Radiobutton(self,text='View expression of selected gene(s)                      ', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=0, variable=correlategene).pack()
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Radiobutton(self,text='View expression of correlated genes (enter 1 gene)', bg=colbg, font = tkfont.Font(family='Times New Roman', size=12), value=1, variable=correlategene).pack()
        log=tk.IntVar(value=1)
        tk.Checkbutton(self, text = "Correlation: linear (on) / logaritmic (off)", font=('Times New Roman', '12'), bg=colbg, onvalue=1, offvalue=0, variable=log, height=1, width = 30).pack()
        tk.Label(self, text='How many genes to correlate?', bg=colbg, font=('Times New Roman', '12')).pack()
        numberCor = Scale(self, from_=1, to=50, bg='white')
        numberCor.set(10)
        numberCor.pack(anchor=CENTER)
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=10)).pack()
        tk.Button(self, text="View Gene Expression", font=('Times New Roman', '12'), fg='black', bg='yellow', command=lambda: viewGenes(headDir=headDir, dirs=dirs, tumortype=tumortype, gene=e1.get().upper(), numberCor=numberCor.get(), correlategene=correlategene.get(), log=log.get())).pack() 
        tk.Label(self, text="", bg=colbg, font = tkfont.Font(family='Times New Roman', size=50)).pack()
        tk.Button(self, text="Main Menu", font=('Times New Roman', '15'), fg='white', bg="red", command=lambda: controller.show_frame("StartPage")).pack(side='bottom')
                   
if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
