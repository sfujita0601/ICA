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
        
        
    def Analeachcomponents(self,SampleDF,compmetaboDF,save_dir):
        import matplotlib.cm as cm

        timepoint=self.timepointlist
        subjectnum = int((len(SampleDF.index)/len(timepoint)))

        OptionDict=dict()
        #BHQ=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20211117/QvalueBH.xlsx',header=0,index_col=0)    
        #BHQ=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20211201/QvalueBH.xlsx',header=0,index_col=0)    
        #ICA_newest
        #BHQ=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20211208/QvalueBH.xlsx',header=0,index_col=0)    
        #BHQ=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20211223/QvalueBH.xlsx',header=0,index_col=0)    
        #SignIndex=BHQ[(BHQ<0.1).sum(axis=1)>1].index
        
        BHQ =pd.read_excel( '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20220113/result_adjR-squared_Individual_Time.xlsx',header=0,index_col=0) 
        #miceOGTTorgansICA_T
        #BHQ =pd.read_excel( '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20220227/row_MultiorganT_ICA/result_adjR-squared_Individual_Time.xlsx',header=0,index_col=0) 
        SignIndex=BHQ.index#[BHQ>0.7].dropna().index

        #OptionDict=dict()
        #全成分をplotする
        ColorList=['darkblue','darkgreen','darkred','darkmagenta']
        #self.plotMultComp(comptimeDF[[0,1]],save_dir,'comptime',ColorList)
        #TagHel.plotMultComp(comptimeDF[['X3','X4']],save_dir,'comptimet',colorList=ColorList)
        #self.plotComp(compsubjectDF[[0,1,2,3]],save_dir,'cBHQompsubject')#各成分をplotする
        
        OptionDict['xlabel']='Component  1';OptionDict['ylabel']='Component 2';OptionDict['Annotate']=0;OptionDict['title']='Compindividuals';OptionDict['calcR'] =''
        #ColorList = ['black' for i in range(len(compsubjectDF.index))]
        #OptionDict['Label'] = list(compmetaboDF.index);##OptionDict['markersize']=5 ;OptionDict['LabelSize']=0.1
        #GraphHel.mkScatterWHist(compsubjectDF[0],compsubjectDF[1],save_dir,ColorList,OptionDict)
        
        DF =SampleDF
        for ij in range(len(SignIndex)):
        ### 個人ごとの順番から時点ごとに帰る
         #   for i in range(len(timepoint)):                    
                #時点ごと
          #      if i==0:
           #         DF=SampleDF.iloc[[i + len(timepoint)*(ii-1) for ii in range(1, subjectnum+1)],:]
                    #print(str([i + len(timepoint)*(ii-1) for ii in range(1, subjectnum+1)])+'_'+str(i))
            #    else:
                    
             #       tempDF=SampleDF.iloc[[i + len(timepoint)*(ii-1) for ii in range(1, subjectnum+1)],:]
              #      DF=pd.concat([DF,tempDF],axis=0)
            #ColorList = [cm.jet(i/20)*14 for i in range(20)]
            #ColorList = list(itertools.chain.from_iterable([[cm.jet(i/20)]*14 for i in range(20)]))
            self.mkBar(list(DF['metasample'+str(SignIndex[ij])]), 'Compsample'+str(SignIndex[ij]), 'Component 1', '', 'black', list(DF.index), 10, 2,10,save_dir)
        #GraphHel.mkBar(list(DF['metasample'+str(SignIndex[ij])]), 'Compsample'+str(SignIndex[ij]), 'Component 1', '', 'black', list(DF.index), 10, 2,10,save_dir)

        #GraphHel.mkBar(list(compsubjectDF[1]), 'Compsubject2', 'Component 2', '', 'black', [str(i+1) for i in range(3)], 10, 10,10,save_dir)
        try:
            self.mkBar(list(OptionDict['compconditionDF1']), 'Compcondition1', 'Component 1', '', 'black', [str(i+1) for i in range(6)], 10, 10,10,save_dir)
            self/mkBar(list(OptionDict['compconditionDF2']), 'Compcondition2', 'Component 2', '', 'black', [str(i+1) for i in range(6)], 10, 10,10,save_dir)
        except:
            pass
        
    def ANOVA_Time_Indiv_Fasting(self,sampleDF,ICAexplain,save_dir):
        import seaborn as sns
        from statsmodels.stats.outliers_influence import variance_inflation_factor
        from statsmodels.api import OLS
        from sklearn import preprocessing
        import statsmodels.api as sm
        import category_encoders as ce   
        import statsmodels.formula.api as smf
        import statsmodels.stats as sts
        import scikit_posthocs as sp

        #各時点に対して、ダミー変数を与える
        timepoint=self.timepoint
        metasamples=sampleDF
        data4anova=sampleDF.copy()
        ResultDF = pd.DataFrame(data=None, index=range(len(sampleDF.columns)),columns= ['intercept']+timepoint[1:] + ['WT'])
        ResultparamDF= pd.DataFrame(data=None, index=range(len(sampleDF.columns)),columns=   ['intercept']+ timepoint[1:]+ ['WT'])
        ResultRDF= pd.DataFrame(data=None, index=range(len(sampleDF.columns)),columns=['rsquared_adj'])
        ResulVIFDF= pd.DataFrame(data=None, index=range(len(sampleDF.columns)),columns=   ['intercept']+ timepoint[1:]+  ['WT'])

        Organs = ['Liver']
        TimeList_New =   [[Organs[i]]*int((len(sampleDF.index)/len(Organs))) for i in range(len(Organs))]
        data4anova['Organs'] = list(itertools.chain.from_iterable(TimeList_New ))
        #sampleDF_Encoder = ce.OneHotEncoder(cols=['Individual'])

        #sampleDF_Encoder_onehot = sampleDF_Encoder.fit_transform(sampleDF['Individual']) #実行
        #print(pd.concat([sampleDF['Individual'], sampleDF_Encoder_onehot], axis=1)) #確認
        
        #X = sampleDF_Encoder_onehot 
        #sampleDF=sampleDF.drop('Individual',axis=1)

        #ResultDF = pd.DataFrame(data=None, index=range(len(sampleDF.columns)),columns=timepoint)
        #ResultparamDF= pd.DataFrame(data=None, index=range(len(sampleDF.columns)),columns=timepoint)
        #ResultRDF= pd.DataFrame(data=None, index=range(len(sampleDF.columns)),columns=[0])
        
        TimeList_New =   [[timepoint[i]]* 5 for i in range(len(timepoint))]*2
        data4anova['Time'] =  list(itertools.chain.from_iterable(TimeList_New ))
        #data4anova=sampleDF.drop('Time',axis=1) 
        #sampleDF_Encoder = ce.OneHotEncoder(cols=['Time'])

        #sampleDF_Encoder_onehot = sampleDF_Encoder.fit_transform(sampleDF['Time']) #実行
        #X = pd.concat([X,sampleDF_Encoder_onehot ],axis=1)
        #data4anova=sampleDF.drop('Time',axis=1)     
        TimeList_New = [0]* len(timepoint) * 5 + [1]* len(timepoint) * 5
        data4anova['WTOB'] = TimeList_New 
        
        factors =['Time','WTOB']
        # カテゴリカル因子の定義
        categorical_factors = ' + '.join([ 'C(Q("'+ factor +'"))' for factor in factors])
    
        # モデルの構築(やや時間かかる)
        #categorical_factors+str(-1), 

        cat_model = {metasample:smf.ols(formula = metasample + ' ~ ' + categorical_factors,  data = data4anova).fit() for metasample in metasamples.columns}
        
        # 分散分析
        anova_res = {metasample:sm.stats.anova_lm(cat_model[metasample],typ=2) for metasample in metasamples.columns}
        # 一般線形モデル＆分散分析
    
        for i in range(len(metasamples.columns)):
            ResultRDF.iloc[i,0]= cat_model[metasamples.columns[i]].rsquared_adj       

            ResultDF.iloc[i,:]=list(cat_model[metasamples.columns[i]].pvalues)
            ResultparamDF.iloc[i,:]=list(cat_model[metasamples.columns[i]].params)
 
            variables = cat_model[metasamples.columns[i]].model.exog
            ResulVIFDF.iloc[i,:]=[variance_inflation_factor(variables, i) for i in range(variables.shape[1])]
            #vif 
        
        # 寄与率を計算する
        def factor_contribution_ratio(metasample):
            return (anova_res[metasample]["sum_sq"]#[:-1]
                     ) / anova_res[metasample]["sum_sq"].sum()
        #残差引かなくていいのでは？
        #- anova_res[metasample].loc["Residual","sum_sq"]            
        #'mean_sq'なのか？
                    #* anova_res[metasample].loc["Residual","mean_sq"] ) / anova_res[metasample]["sum_sq"].sum()
        
        # 主効果（）の平均によって符号を決める
        
        def mean_effect_sign(metasample, factor):
            return np.sign(cat_model[metasample].params[cat_model[metasample].params.index.str.contains(factor)].values.mean())
        
        # 寄与率をヒートマップにして出力する
        factors = factors
        conts = pd.DataFrame([factor_contribution_ratio(metasample) for metasample in metasamples.columns]).T
        conts.index = factors +['Residual']
        conts.columns = metasamples.columns
        
        signed_conts = pd.DataFrame([[mean_effect_sign(metasample, factor)
                               for metasample in metasamples.columns]
                               for factor in factors +['Residual']])
        ###
        signed_conts.iloc[2,:]=[1 for i in range(len(signed_conts.columns))]
        ###        
        signed_conts.loc[:,:] = conts.values*signed_conts.values
        signed_conts.index = factors+['Residual']
        signed_conts.columns = [metasample + " (" + '{:.2%}'.format(ICAexplain["metagene"+str(i+1)]) + ")" for i,metasample in enumerate(metasamples.columns)]

        signed_conts.to_excel(save_dir+'signed_conts.xlsx')


        GraphHel.mkhistList(list(ResultDF.values.flatten()),save_dir,'pvalues','pvalues',10,'')#1つのリストのヒストグラム
        GraphHel.mkhistList(list(ResultparamDF.values.flatten()),save_dir,'coefficient','coefficient',10,'')#1つのリストのヒストグラム
        GraphHel.mkhistList(list(ResultRDF.values.flatten()),save_dir,'adjR-squared','adjR-squared',10,'')#1つのリストのヒストグラム
                    
        ResultDF.index=list(ResultDF.index+1);ResultparamDF.index=list(ResultDF.index);ResultRDF.index=list(ResultDF.index);
        ResultDF.to_excel(save_dir+'result_pvalues_Individual_Time.xlsx')
        ResultparamDF.to_excel(save_dir+'result_coefficient_Individual_Time.xlsx')
        ResultRDF.to_excel(save_dir+'result_adjR-squared_Individual_Time.xlsx')
        ResulVIFDF.to_excel(save_dir+'result_VIF.xlsx')

        # 分散分析の q valueに応じて*を出す
        alpha = 10**-1
        pvals = pd.DataFrame([anova_res[metasample]["PR(>F)"][:-1] for metasample in metasamples.columns]).T
        _,qvals,_,_ = sts.multitest.multipletests(pvals.values.reshape(-1), alpha=0.1,method='fdr_bh')
        pvals.to_excel(save_dir+'pvals.xlsx')

        pvals.loc[:,:] = np.where(qvals.reshape(pvals.shape) < alpha, "*", "")
              
        s = pd.Series(['' for i in range(len(pvals.columns))], index=pvals.columns, name='Residual')
        pvals = pvals.append(s)
        pvals.index = factors +['Residual']
        pvals.columns = metasamples.columns
        ResultRDF.to_excel(save_dir+'ResultRDF.xlsx')
        conts.to_excel(save_dir+'conts.xlsx')
        
        _,qvals,_,_ = sts.multitest.multipletests(ResultDF.values.reshape(-1), alpha=0.1,method='fdr_bh')
        ResultDF.loc[:,:] = np.where(qvals.reshape(ResultDF.shape) < alpha, "*", "")
        ResultDF.to_excel(save_dir+'result_qvalues_Individual_Time.xlsx')
 
        for j in ['Time','WTOB']:#,'TwoHGlucose','TwoHInsulin']:
            GraphHel.mkhistList(list(np.sort(conts.loc[j].values.flatten())),save_dir,j,j,9,'')#1つのリストのヒストグラム
            print(j +'_'+ str(SC.threshold_otsu(list(conts.loc[j].values.flatten()))))

        plt.figure(figsize=(1+conts.shape[1]/3, conts.shape[0]/3))
        sns.heatmap(signed_conts,cmap="PuOr_r", vmin=-1, vmax=1, annot=pvals, fmt = '',annot_kws={'size':14, 'weight':'bold'});
        plt.xticks(size='13',rotation=270);plt.yticks(size='10',rotation=0)

        plt.savefig(save_dir+'signed_conts.pdf',bbox_inches="tight")  
        Col=list(metasamples.columns)

        Signcomp = list(pvals.columns[pvals.loc['Time']=='*'])
        #Signcomp = [0,1,2,3,5,6,7,8,10,12,13,14,15,18,19,20,22,23,24,24,25,28,30,32,34,35,39,40]
        #Signcomp = [0,2,3,5,8,9,11]
        #for i in [1,2,5,6,8,9,12,17,26,28,44]:
        ScheffePvalDict = {metasample:sp.posthoc_scheffe(data4anova,val_col=metasample,group_col='Time') for metasample in Signcomp}
        for  count,i in enumerate(Signcomp):
            if count==0:
                pvalsTime=ScheffePvalDict[i]
            else:
                tempDF = ScheffePvalDict[i]
                pvalsTime=pd.concat([pvalsTime,tempDF])
        print(pvalsTime)
        _,qvals,_,_ = sts.multitest.multipletests(pvalsTime.values.reshape(-1), alpha=0.1,method='fdr_bh')
        pd.DataFrame(pvalsTime.to_excel(save_dir+'Scheffe_pvalues_Time.xlsx'))

        pd.DataFrame(qvals.reshape(pvalsTime.shape)).to_excel(save_dir+'Scheffe_qvalues_Time.xlsx')

        Signcomp_WTOB = list(pvals.columns[pvals.loc['WTOB']=='*'])
        
       # Signcomp = [0,1,3,4,6,7,8,10,12,13,14,15,19,20,22,23,24,24,28,30,32,34,41]
        #Signcomp = [0,3,4,8,9,11]

        ScheffePvalDict = {metasample:sp.posthoc_scheffe(data4anova,val_col=metasample,group_col='WTOB') for metasample in Signcomp_WTOB}
        for  count,i in enumerate(Signcomp_WTOB):
            if count==0:
                pvalsWTOB=ScheffePvalDict[i]
            else:
                tempDF = ScheffePvalDict[i]
                pvalsWTOB=pd.concat([pvalsWTOB,tempDF])

        _,qvals,_,_ = sts.multitest.multipletests(pvalsWTOB.values.reshape(-1), alpha=0.1,method='fdr_bh')
        pd.DataFrame(pvalsWTOB.to_excel(save_dir+'Scheffe_pvalues_WTOB.xlsx'))

        pd.DataFrame(qvals.reshape(pvalsWTOB.shape)).to_excel(save_dir+'Scheffe_qvalues_WTOB.xlsx')

        #Signcomp = [0,1,2,3,5,6,7,8,10,12,13,14,15,18,19,20,22,23,24,24,25,28,30,32,34,35,39,40]
        #Signcomp = [0,2,3,5,8,9,11]
        ##複数因子でtestしたFDRを補正
        pvalsComb = np.concatenate([pvalsTime.values.reshape(-1),pvalsWTOB.values.reshape(-1)])
        _,qvals,_,_ = sts.multitest.multipletests(pvalsComb, alpha=0.1,method='fdr_bh')
        qvalsTime=qvals[:len(pvalsTime.values.reshape(-1))]
        qvalsTimeDF=pd.DataFrame(qvalsTime.reshape(pvalsTime.shape),index=[metasample+'_'+ str(timepoint[i]) for metasample in Signcomp for i in range(len(timepoint)) ], columns=[str(timepoint[i]) for i in range(len(timepoint)) ])
        qvalsTimeDF.to_excel(save_dir+'Scheffe_qvalues_Time_comb.xlsx')
        #Signcomp = [0,1,3,4,6,7,8,10,12,13,14,15,19,20,22,23,24,24,28,30,32,34,41]
        #Signcomp = [0,3,4,8,9,11]
   
        qvalsWTOB=qvals[len(pvalsTime.values.reshape(-1)):len(pvalsTime.values.reshape(-1))+len(pvalsWTOB.values.reshape(-1))]
        WTOB = ['WT', 'OB']
        qvalsWTOBDF=pd.DataFrame(qvalsWTOB.reshape(pvalsWTOB.shape),index=[metasample+'_'+ WTOB[i]  for metasample in Signcomp_WTOB for i in range(2)], columns=[WTOB[i] for i in range(2) ])
        qvalsWTOBDF.to_excel(save_dir+'Scheffe_qvalues_WTOB_comb.xlsx')
        
        ### Steel Dwass
        Signcomp = list(pvals.columns[pvals.loc['Time']=='*'])
        SDPvalDict = {metasample:sp.posthoc_dscf(data4anova,val_col=metasample,group_col='Time') for metasample in Signcomp}
        for  count,i in enumerate(Signcomp):
            if count==0:
                pvalsTime=SDPvalDict[i]
            else:
                tempDF = SDPvalDict[i]
                pvalsTime=pd.concat([pvalsTime,tempDF])
        print(pvalsTime)
        _,qvals,_,_ = sts.multitest.multipletests(pvalsTime.values.reshape(-1), alpha=0.1,method='fdr_bh')
        pd.DataFrame(pvalsTime.to_excel(save_dir+'SD_pvalues_Time.xlsx'))

        pd.DataFrame(qvals.reshape(pvalsTime.shape)).to_excel(save_dir+'Scheffe_qvalues_Time.xlsx')

        Signcomp_WTOB = list(pvals.columns[pvals.loc['WTOB']=='*'])
        
       # Signcomp = [0,1,3,4,6,7,8,10,12,13,14,15,19,20,22,23,24,24,28,30,32,34,41]
        #Signcomp = [0,3,4,8,9,11]

        SDPvalDict = {metasample:sp.posthoc_dscf(data4anova,val_col=metasample,group_col='WTOB') for metasample in Signcomp_WTOB}
        for  count,i in enumerate(Signcomp_WTOB):
            if count==0:
                pvalsWTOB=SDPvalDict[i]
            else:
                tempDF = SDPvalDict[i]
                pvalsWTOB=pd.concat([pvalsWTOB,tempDF])

        _,qvals,_,_ = sts.multitest.multipletests(pvalsWTOB.values.reshape(-1), alpha=0.1,method='fdr_bh')
        pd.DataFrame(pvalsWTOB.to_excel(save_dir+'SD_pvalues_WTOB.xlsx'))

        pd.DataFrame(qvals.reshape(pvalsWTOB.shape)).to_excel(save_dir+'Scheffe_qvalues_WTOB.xlsx')

        #Signcomp = [0,1,2,3,5,6,7,8,10,12,13,14,15,18,19,20,22,23,24,24,25,28,30,32,34,35,39,40]
        #Signcomp = [0,2,3,5,8,9,11]
        ##複数因子でtestしたFDRを補正
        pvalsComb = np.concatenate([pvalsTime.values.reshape(-1),pvalsWTOB.values.reshape(-1)])
        _,qvals,_,_ = sts.multitest.multipletests(pvalsComb, alpha=0.1,method='fdr_bh')
        qvalsTime=qvals[:len(pvalsTime.values.reshape(-1))]
        qvalsTimeDF=pd.DataFrame(qvalsTime.reshape(pvalsTime.shape),index=[metasample+'_'+ str(timepoint[i]) for metasample in Signcomp for i in range(len(timepoint)) ], columns=[str(timepoint[i]) for i in range(len(timepoint)) ])
        qvalsTimeDF.to_excel(save_dir+'SD_qvalues_Time_comb.xlsx')
        #Signcomp = [0,1,3,4,6,7,8,10,12,13,14,15,19,20,22,23,24,24,28,30,32,34,41]
        #Signcomp = [0,3,4,8,9,11]
   
        qvalsWTOB=qvals[len(pvalsTime.values.reshape(-1)):len(pvalsTime.values.reshape(-1))+len(pvalsWTOB.values.reshape(-1))]
        WTOB = ['WT', 'OB']
        qvalsWTOBDF=pd.DataFrame(qvalsWTOB.reshape(pvalsWTOB.shape),index=[metasample+'_'+ WTOB[i]  for metasample in Signcomp_WTOB for i in range(2)], columns=[WTOB[i] for i in range(2) ])
        qvalsWTOBDF.to_excel(save_dir+'SD_qvalues_WTOB_comb.xlsx')
        return(pvals)



        
    def CEV(self,X, Xmodel):
        return 1 - self.TSS(X - Xmodel)/self.TSS(X)


    
    def mkBar(self,List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize,save_dir):#任意のリストを代入して棒グラフを描画する
        x = np.linspace( 1, len(List1), len(List1) )
        
        fig = plt.figure(figsize=(5,3))
        ax1 = fig.add_axes((0, 1, 1, 1))
    
        ax1.bar( [0+0.1*i for i in range(len(xticks))], List1,width=0.05, color=Color,tick_label=xticks,linewidth=0.5,ec=Color)
        
        #ax1.set_ylim([-1.0,1.0])
        
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
        plt.savefig(save_dir+Title+'_.pdf',bbox_inches="tight")
        plt.close()    
    
    def plotPCACovRatio(self,cov_ratio):
      # 分散の説明率プロット
      num_var = len(cov_ratio)
      x_tick = np.arange(1,num_var+1)
      plt.bar(x_tick, cov_ratio)
      plt.plot(x_tick, np.cumsum(cov_ratio),'-o', mfc='none',mec='b',mew=2,linewidth=3)
      plt.xticks(x_tick, fontsize=40)#20
      plt.yticks(fontsize=40)#20
    
      plt.axis([1-0.4, num_var+0.5, 0,1])#左に寄せる
      plt.xlabel("Number of PC", fontsize=40)#10
      plt.ylabel("Explained variance ratio", fontsize=40)#10
      plt.rcParams['axes.linewidth'] = 1.5# 軸の線幅edge linewidth。囲みの太さ
    
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

    def ICA4ZscoredGeneSampleMatrix(self,E, k, seed=None):#Noneだと乱数の種もランダム
        decomposer = FastICA(n_components=k, random_state=seed,tol=0.00001,max_iter=10000)#max_iter=200, tol=0.0001,デフォルト
        decomposer.fit(E)
        metagenes = pd.DataFrame(decomposer.transform(E), index=E.index, columns=["metagene"+str(i) for i in range(1,k+1)])
        metasamples = pd.DataFrame(decomposer.mixing_, index=E.columns, columns=["metasample"+str(i) for i in range(1,k+1)])
        return metagenes, metasamples
    
    
    def mkTablegenePathway(self,pvals,save_dir):
        import re
        Signcomp = list(pvals.columns[(pvals=='*').any()])
        IndList=[];NumList=[]
        for i in Signcomp:
            jj = re.sub(r"\D", "",i)
            #enrich=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20220305/Enrich_'+jj+'.xlsx',header=0,index_col=0)
            try:
                enrich=pd.read_excel(save_dir+'/Enrich_'+jj+'_KEGG_2019_Mouse.xlsx',header=0,index_col=0)
                IndList += list(enrich[enrich['Adjusted P-value']<0.1]['Term'])
                NumList+=[jj]
            except:
                pass
        tempDF = pd.DataFrame(data=None,index=list(set(IndList)),columns=NumList)
 
            
        for i in Signcomp:
            
            jj = re.sub(r"\D", "",i)
            #enrich=pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20220305/Enrich_'+jj+'.xlsx',header=0,index_col=0)
            try:
                enrich=pd.read_excel(save_dir+'/Enrich_'+jj+'_KEGG_2019_Mouse.xlsx',header=0,index_col=0)
                Idx= list(enrich[enrich['Adjusted P-value']<0.1]['Term']) 
                tempDF[jj][Idx] = '*'
            except:                pass
        tempDF.to_excel(save_dir+'genePathwayenrich.xlsx')    

    def PypeR(Pvalue,SwitchDict):
        import pyper
        import pandas as pd
        
        # Python で CSV のデータを読み出す
        #wine = pd.read_csv("wine.csv")
        
        # R のインスタンスを作る
        r = pyper.R(use_numpy = 'True', use_pandas='True',RCMD ='/usr/local/bin/R')
        
        # Python のオブジェクトを R に渡す
        r.assign("data", Pvalue)
        #r.assign("filename", SwitchDict['Time'])
        r.assign("filename", SwitchDict['EngLabel'])
        data=r.get("data")
        #print(data)
        r("'+' <- function(e1, e2) { \
      if (is.character(c(e1, e2))) { \
        paste(e1, e2, sep = '') \
      } else { \
        base::'+'(e1, e2) \
      } \
    }")
        r("file_name1 <- filename")
        r('file_name2 <- "c"')
        r("today <- Sys.Date()")
        r("Today = strftime(today, format='%Y%m%d')")
        
        r("save_dir <- '/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/' + Today + '/' + file_name1 + '/'")
        r('if (!file.exists(save_dir)) \
    { \
      dir.create(save_dir) \
    }')
        r("Pvalue <- data")
        #print(Pvalue)
        #rownames(Pvalue)<-  data$X__1
        #r("Pvalue <- Pvalue[,-1]")
        r('library(qvalue)')
        r('library(xlsx)')
        r('Pvalue <- as.matrix((data))')
        r('Pvalue<-as.numeric(Pvalue)')
            
        r('  q.values <- p.adjust(Pvalue, method = "BH")')
        result = r.get("q.values")
        print(result)
        #1行しかないとき
        if len(Pvalue.index) < 2:
            r('  QValue <- matrix(q.values,ncol=length(names(data)), nrow=1)')
        #そうでない時
        else:
            r('  QValue <- matrix(q.values,ncol=length(colnames(data)), nrow=length(rownames(data)))')
            r('  colnames(QValue) <- names(data)')
            r('rownames(QValue) <- rownames(data)')
        
        r('write.csv(QValue,save_dir + "/QvalueBH" + file_name1 +"_" + file_name2 + ".csv",append=T, quote=F, col.names=F,fileEncoding = "CP932")')
        r('write.xlsx(QValue, file=save_dir + "/QvalueBH"+ file_name1 +"_" + file_name2 + ".xlsx", sheetName="sheet1", row.names=T)')
        QvalueBH = r.get("QValue")
    
        
    
    
        r('qobj <- qvalue(p = Pvalue, fdr.level=0.1)')
        if len(Pvalue.index) < 2:
            r('  colnames(qobj$qvalues) <- names(data)')
            r('rownames(qobj$qvalues) <- rownames(data)')
        else:
            r('  QValue <- matrix(qobj$qvalues,ncol=length(colnames(data)), nrow=length(rownames(data)))')
    
            r('  colnames(QValue) <- colnames(data)')
            r('rownames(QValue) <- rownames(data)')        
        r('  write.csv(QValue,save_dir + "/QvalueStorey" + file_name1 +"_" + file_name2 + ".csv",append=T, quote=F, col.names=F,fileEncoding = "CP932")')
        r('  write.xlsx(QValue, file=save_dir + "/QvalueStorey"+ file_name1 +"_" + file_name2 + ".xlsx", sheetName="sheet1", row.names=T)')
        QvalueStorey = r.get("QValue")
    
        r('  xlab="λ"')
        r('ylab="π0(λ)"')
        r('pdf(save_dir + "Qvalueplot.pdf")')
        r('plot(qobj$lambda, qobj$pi0.lambda,xlab=xlab,ylab=ylab)')
        r('title(main = bquote(π0==.(qobj$pi0)))')
        r('dev.off()')
      
        #r("  #pvalue vs qvalue \
        r('pdf(save_dir + "PvaluevsQvalueplot.pdf")')
        r("plot(qobj$pvalues,qobj$qvalues,xlab='Pvalue',ylab='Qvalue')")
        r('dev.off()')  
        
          #Qvalueヒストグラム
        r('pdf(save_dir + "QvalueHits.pdf")')
        r('WOnanQvalue <- qobj$qvalues[!is.na(qobj$qvalues)]')
        r('histtitle <- length(WOnanQvalue[WOnanQvalue<0.1])')
          #  \
        r('hist(qobj$qvalues,xlim=c(0,1),xlab="Qvalue",main="")')
        r('title(main = bquote(num(q<0.1)==.(histtitle)))')
        r('dev.off()')
          # \')
        r('pdf(save_dir + "QvalueAll.pdf")')
        r('plot(qobj)')
        r('dev.off()') 
        
        
        result = r.get("Pvalue")
        #print(result)
        #r("pdf(save_dir + 'PvalueHist.pdf')") #pythonで吐いてるからいらない
        #r("hist(Pvalue,breaks=20)")
        #r("dev.off()")
            # R のソースコードを実行する
        r("source(file='../../R/AnalParamAll.R')")
        
        return(QvalueStorey,QvalueBH)    
    
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
    
        # クラスタのラベルごとにmetageneのindexをまとめる
        lab_ind = {label:pd.Series(label_assign)[label_assign==label].index for label in label_set}
    
        print('DBSCAN identified',n_clusters,'clusters of ICs.')
    
        # 全部のクラスタを0番目に対して回転させてから集める
        ClusteredMetagenes = {lab:[] for lab in label_set}
        for lab in label_set:
            ClusteredMetagenes[lab].append(AllMetagenes[:,lab_ind[lab][0]])
            for j in range(1, len(lab_ind[lab])):
                if CorMat[lab_ind[lab][0], lab_ind[lab][j]]>0: #CorMat[lab_ind[lab]].T[lab_ind[lab]][0][j] > 0:
                    ClusteredMetagenes[lab].append(AllMetagenes[:,lab_ind[lab][j]])
                else:
                    ClusteredMetagenes[lab].append(-AllMetagenes[:,lab_ind[lab][j]])
            ClusteredMetagenes[lab] = pd.DataFrame(ClusteredMetagenes[lab]).T
    
        # クラスター平均を出す
        ClusterMeanMetagenes = pd.DataFrame([ClusteredMetagenes[lab].mean(axis=1) for lab in label_set]).T
        ClusterMeanMetagenes.columns = ["metagene"+str(i+1) for i in range(ClusterMeanMetagenes.shape[1])] 
        ClusterMeanMetagenes.index = X.index
    
        # 正規直交化を行う
        W = ClusterMeanMetagenes.values
    
        W = W / np.sqrt(np.linalg.norm(W.T@W,np.inf))
        for i in range(100):
            W = (3/2) * W - (1/2) * W@W.T@W
        W = pd.DataFrame(W, index=ClusterMeanMetagenes.index,columns=ClusterMeanMetagenes.columns)
    
        print("Orthogonality of ICs is", (np.linalg.det(W.T@W)*100).round(0),"%")
        combined = pd.concat([ClusterMeanMetagenes, W],axis=1)

        sns.heatmap(combined.corr(),vmin=-1,vmax=1,cmap="bwr") 
        #plt.savefig(save_dir+'Combinedetagen.pdf',bbox_inches="tight")   
        # metasampleを計算してソートする
        S = W.T.dot(X).T
        S.columns = [s.replace("metagene","metasample") for s in S.columns]
    
        return W, S    