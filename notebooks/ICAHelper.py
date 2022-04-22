#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: fujita
"""

import matplotlib.cm as cm
from sklearn.decomposition import FastICA
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
import seaborn as sns
import scipy.stats as sts
import numpy as np
import pandas as pd

class ICA_class:


    def __init__(self):

        self.file_dir = ''
        self.save_dir=''
        self.timepointlist=[0,2,4,6,8,12,16,24]
        
        




        
    def CEV(self,X, Xmodel):
        return 1 - self.TSS(X - Xmodel)/self.TSS(X)


    
    def mkBar(self,List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize):
        x = np.linspace( 1, len(List1), len(List1) )
        
        fig = plt.figure(figsize=(5,3))
        ax1 = fig.add_axes((0, 1, 1, 1))
    
        ax1.bar( [0+0.1*i for i in range(len(xticks))], List1,width=0.05, color=Color,tick_label=xticks,linewidth=0.5,ec=Color)
        
        xmin, xmax, ymin, ymax = ax1.axis() 
        
        p = ax1.hlines([0], xmin, xmax, "black", linestyles='solid',linewidth=1)     # hlines
        
        axis=['top','bottom','left','right']
        line_width=[1,1,1,1]
        
        for a,w in zip(axis, line_width):  # change axis width
            ax1.spines[a].set_linewidth(w)
    
        ax1.set_xlim(xmin, xmax)
        ax1.set_title(Title,fontsize=Titlesize)
        ax1.set_xlabel(xlabel,fontsize=xsize)
        ax1.set_ylabel(ylabel,fontsize=size)
        #plt.xticks(x_tick)
        ax1.set_xticklabels(labels=xticks,rotation=270,fontsize=xsize) 
        #ax1.set_xticklabels(fontsize=size)    
        xmin, xmax, ymin, ymax = ax1.axis() 
        #plt.plot([xmin, xmax],[1, 1], "black", linestyle='dashed') # normal way
        #plt.legend(handles=ax1,loc='best')

    def mkhistList(List1,save_dir,xlabel,filename,bins,title):
            fig = plt.figure(figsize=(6.4,4.8))
            plt.hist(List1,bins=bins,ec='black',align='left' )#,cumulative=1,normed=True)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)  

            plt.xlabel(xlabel,fontsize=20); plt.ylabel('Frequency',fontsize=20);plt.tick_params(labelsize=20);plt.title(title,size=20)
            ##plt.savefig(save_dir +  filename +'Hist.pdf',bbox_inches="tight")  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])

    
    def plotPCACovRatio(self,cov_ratio):

      num_var = len(cov_ratio)
      x_tick = np.arange(1,num_var+1)
      plt.bar(x_tick, cov_ratio)
      plt.plot(x_tick, np.cumsum(cov_ratio),'-o', mfc='none',mec='b',mew=2,linewidth=3)
      plt.xticks(x_tick, fontsize=40)#20
      plt.yticks(fontsize=40)#20
    
      plt.axis([1-0.4, num_var+0.5, 0,1])
      plt.xlabel("Number of PC", fontsize=40)#10
      plt.ylabel("Explained variance ratio", fontsize=40)#10
      plt.rcParams['axes.linewidth'] = 1.5
    
    def TSS(self,Y):
        return np.sum(np.array(Y)**2)        


    def ICAsort(self,X,metagenes,metasamples):
        k = metagenes.shape[1]
        ICAexplain = pd.Series([self.CEV(X, np.outer(metagenes.iloc[:,i], metasamples.iloc[:,i])) for i in range(k)])
        sorted_meta_index = ICAexplain.sort_values(ascending=False).index
        ICAexplain = ICAexplain[sorted_meta_index]
        ICAexplain.index = metagenes.columns
        metagenes = metagenes.iloc[:,sorted_meta_index].set_axis(metagenes.columns,axis=1,inplace=False)
        metasamples = metasamples.iloc[:,sorted_meta_index].set_axis(metasamples.columns,axis=1,inplace=False)
        return metagenes, metasamples, ICAexplain  

    def ICA4ZscoredGeneSampleMatrix(self,E, k, seed=None):
        decomposer = FastICA(n_components=k, random_state=seed,tol=0.00001,max_iter=10000)#max_iter=200, tol=0.0001,デフォルト
        decomposer.fit(E)
        metagenes = pd.DataFrame(decomposer.transform(E), index=E.index, columns=["metagene"+str(i) for i in range(1,k+1)])
        metasamples = pd.DataFrame(decomposer.mixing_, index=E.columns, columns=["metasample"+str(i) for i in range(1,k+1)])
        return metagenes, metasamples
    
    
    def mkTablegenePathway(self,pvals):
        import re
        Signcomp = list(pvals.columns[(pvals=='*').any()])
        IndList=[];NumList=[]
        for i in Signcomp:
            jj = re.sub(r"\D", "",i)
            try:
                enrich=pd.read_excel(save_dir+'/Enrich_'+jj+'_KEGG_2019_Mouse.xlsx',header=0,index_col=0)
                IndList += list(enrich[enrich['Adjusted P-value']<0.1]['Term'])
                NumList+=[jj]
            except:
                pass
        tempDF = pd.DataFrame(data=None,index=list(set(IndList)),columns=NumList)
 
            
        for i in Signcomp:
            
            jj = re.sub(r"\D", "",i)
            try:
                enrich=pd.read_excel(save_dir+'/Enrich_'+jj+'_KEGG_2019_Mouse.xlsx',header=0,index_col=0)
                Idx= list(enrich[enrich['Adjusted P-value']<0.1]['Term']) 
                tempDF[jj][Idx] = '*'
            except:                pass

   
    
    def robust_ica(self,X, n_components=100, n_repeat=10, **kwargs):
        from joblib import Parallel,delayed
        import seaborn as sns
        print(n_components,n_repeat)
        print("Getting", n_components, "ICs *", n_repeat, "times in parallel.")
        results = Parallel(n_jobs=-1)(delayed(self.ICA4ZscoredGeneSampleMatrix)(X, k=n_components, seed=i, **kwargs) for i in range(n_repeat))
        
        AllMetagenes = np.concatenate([results[i][0].values for i in range(n_repeat)],axis=1)
        AllMetasamples = np.concatenate([results[i][1].values for i in range(n_repeat)],axis=1)
    
        print("Calculating distances b/w ICs.")
    
        CorMat = pd.DataFrame(AllMetagenes).corr().values
        DisMat = 1-abs(CorMat)
    
        # Run DBSCAN
        from sklearn.cluster import DBSCAN
        
        min_samples = int(n_repeat/2)
        dbscan = DBSCAN(eps=0.1, min_samples=min_samples, metric='precomputed')
        label_assign = dbscan.fit_predict(DisMat)
        label_set = sorted(set(label_assign)-{-1})
        #    n_clusters = len(set(label_assign))
        n_clusters = len(label_set)
    
        # Coalescing metagene indices for each cluster label
        lab_ind = {label:pd.Series(label_assign)[label_assign==label].index for label in label_set}
    
        print('DBSCAN identified',n_clusters,'clusters of ICs.')
    
        # Rotate all clusters with respect to the 0th and then collect them.
        ClusteredMetagenes = {lab:[] for lab in label_set}
        for lab in label_set:
            ClusteredMetagenes[lab].append(AllMetagenes[:,lab_ind[lab][0]])
            for j in range(1, len(lab_ind[lab])):
                if CorMat[lab_ind[lab][0], lab_ind[lab][j]]>0: #CorMat[lab_ind[lab]].T[lab_ind[lab]][0][j] > 0:
                    ClusteredMetagenes[lab].append(AllMetagenes[:,lab_ind[lab][j]])
                else:
                    ClusteredMetagenes[lab].append(-AllMetagenes[:,lab_ind[lab][j]])
            ClusteredMetagenes[lab] = pd.DataFrame(ClusteredMetagenes[lab]).T
    
        # Calculate cluster averages
        ClusterMeanMetagenes = pd.DataFrame([ClusteredMetagenes[lab].mean(axis=1) for lab in label_set]).T
        ClusterMeanMetagenes.columns = ["metagene"+str(i+1) for i in range(ClusterMeanMetagenes.shape[1])] 
        ClusterMeanMetagenes.index = X.index
    
        # Perform ortho-normalization
        W = ClusterMeanMetagenes.values
    
        W = W / np.sqrt(np.linalg.norm(W.T@W,np.inf))
        for i in range(100):
            W = (3/2) * W - (1/2) * W@W.T@W
        W = pd.DataFrame(W, index=ClusterMeanMetagenes.index,columns=ClusterMeanMetagenes.columns)
    
        print("Orthogonality of ICs is", (np.linalg.det(W.T@W)*100).round(0),"%")
        combined = pd.concat([ClusterMeanMetagenes, W],axis=1)

        sns.heatmap(combined.corr(),vmin=-1,vmax=1,cmap="bwr") 
        ##plt.savefig(save_dir+'Combinedetagen.pdf',bbox_inches="tight")   
        # Calculate and sort metasample
        S = W.T.dot(X).T
        S.columns = [s.replace("metagene","metasample") for s in S.columns]
    
        return W, S    