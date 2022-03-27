#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: fujita
"""


import matplotlib.cm as cm
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import statsmodels.stats.multitest as multi
import seaborn as sns
import scipy.stats as sts
import statsmodels.stats as stsmodel
import matplotlib.pyplot as plt
import itertools

import numpy as np
import pandas as pd

class Anal_class:
    def __init__(self):
        self.timepointlist=[0,2,4,6,8,12,16,24]

    def Analeachcomponents(self,SampleDF,compmetaboDF):
        import matplotlib.cm as cm

        timepoint=self.timepointlist
        subjectnum = int((len(SampleDF.index)/len(timepoint)))

        OptionDict=dict()

        SignIndex=SampleDF.columns

        for ij in range(len(SignIndex)):
            self.plotBar(list(SampleDF[str(SignIndex[ij])]), 'Compsample'+str(ij+1), 'Compsample'+str(ij+1), '', 'black', list(SampleDF.index), 10, 2,10)

    def ANOVA_Time_Indiv_Fasting(self, metasamples,ICAexplain):
        import seaborn as sns
        from statsmodels.stats.outliers_influence import variance_inflation_factor

        from statsmodels.api import OLS
        from sklearn import preprocessing
        import statsmodels.api as sm
        import category_encoders as ce   
        import statsmodels.formula.api as smf
        import statsmodels.stats as sts
        import scikit_posthocs as sp

        timepoint=self.timepointlist
        data4anova=metasamples.copy()
        ResultDF = pd.DataFrame(data=None, index=range(len(metasamples.columns)),columns= ['intercept']+timepoint[1:] + ['WT'])
        ResultparamDF= pd.DataFrame(data=None, index=range(len(metasamples.columns)),columns=   ['intercept']+ timepoint[1:]+ ['WT'])
        ResultRDF= pd.DataFrame(data=None, index=range(len(metasamples.columns)),columns=['rsquared_adj'])
        
        
        TimeList=   [[timepoint[i]]* 5 for i in range(len(timepoint))]*2
        data4anova['Time'] =  list(itertools.chain.from_iterable(TimeList ))
    
        WTOBList = [0]* len(timepoint) * 5 + [1]* len(timepoint) * 5
        data4anova['WTOB'] = WTOBList
        

        factors =['Time','WTOB']
        # Definition of categorical factors
        categorical_factors = ' + '.join([ 'C(Q("'+ factor +'"))' for factor in factors])
    
        # Building the model (somewhat time consuming)
        cat_model = {metasample:smf.ols(formula = metasample + ' ~ ' + categorical_factors,  data = data4anova).fit() for metasample in metasamples.columns}
        
        # ANOVA
        anova_res = {metasample:sm.stats.anova_lm(cat_model[metasample],typ=2) for metasample in metasamples.columns}
    
        for i in range(len(metasamples.columns)):
            ResultRDF.iloc[i,0]= cat_model[metasamples.columns[i]].rsquared_adj       
            ResultDF.iloc[i,:]=list(cat_model[metasamples.columns[i]].pvalues)
            ResultparamDF.iloc[i,:]=list(cat_model[metasamples.columns[i]].params)
            variables = cat_model[metasamples.columns[i]].model.exog

        
        # Calculate the contribution ratio
        def factor_contribution_ratio(metasample):
            return (anova_res[metasample]["sum_sq"]#[:-1]
                     ) / anova_res[metasample]["sum_sq"].sum()

        # The sign is determined by the average of the main effects
        def mean_effect_sign(metasample, factor):
            return np.sign(cat_model[metasample].params[cat_model[metasample].params.index.str.contains(factor)].values.mean())
        
        # Output a heatmap of contribution rates
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

        #signed_conts.to_excel(save_dir+'signed_conts.xlsx')


        self.mkhistList(list(ResultDF.values.flatten()),'pvalues','pvalues',10,'')
        self.mkhistList(list(ResultparamDF.values.flatten()),'coefficient','coefficient',10,'')
        self.mkhistList(list(ResultRDF.values.flatten()),'adjR-squared','adjR-squared',10,'')
                    
        #ResultDF.index=list(ResultDF.index+1);ResultparamDF.index=list(ResultDF.index);ResultRDF.index=list(ResultDF.index);
        #ResultDF.to_excel(save_dir+'result_pvalues_Individual_Time.xlsx')
        #ResultparamDF.to_excel(save_dir+'result_coefficient_Individual_Time.xlsx')
        #ResultRDF.to_excel(save_dir+'result_adjR-squared_Individual_Time.xlsx')

        # Produces * according to the q value of ANOVA
        alpha = 10**-1
        pvals = pd.DataFrame([anova_res[metasample]["PR(>F)"][:-1] for metasample in metasamples.columns]).T
        _,qvals,_,_ = sts.multitest.multipletests(pvals.values.reshape(-1), alpha=0.1,method='fdr_bh')
        #pvals.to_excel(save_dir+'pvals.xlsx')

        pvals.loc[:,:] = np.where(qvals.reshape(pvals.shape) < alpha, "*", "")
              
        s = pd.Series(['' for i in range(len(pvals.columns))], index=pvals.columns, name='Residual')
        pvals = pvals.append(s)
        pvals.index = factors +['Residual']
        pvals.columns = metasamples.columns
        #ResultRDF.to_excel(save_dir+'ResultRDF.xlsx')
        #conts.to_excel(save_dir+'conts.xlsx')
        
        _,qvals,_,_ = sts.multitest.multipletests(ResultDF.values.reshape(-1), alpha=0.1,method='fdr_bh')
        ResultDF.loc[:,:] = np.where(qvals.reshape(ResultDF.shape) < alpha, "*", "")
        #ResultDF.to_excel(save_dir+'result_qvalues_Individual_Time.xlsx')
 
        for j in ['Time','WTOB']:
            self.mkhistList(list(np.sort(conts.loc[j].values.flatten())),j,j,9,'')

        plt.figure(figsize=(1+conts.shape[1]/3, conts.shape[0]/3))
        sns.heatmap(signed_conts,cmap="PuOr_r", vmin=-1, vmax=1, annot=pvals, fmt = '',annot_kws={'size':14, 'weight':'bold'});
        plt.xticks(size='13',rotation=270);plt.yticks(size='10',rotation=0)
        #plt.savefig(save_dir+'signed_conts.pdf',bbox_inches="tight")  
        
        ###       Scheffe 
        Col=list(metasamples.columns)
        Signcomp = list(pvals.columns[pvals.loc['Time']=='*'])

        ScheffePvalDict = {metasample:sp.posthoc_scheffe(data4anova,val_col=metasample,group_col='Time') for metasample in Signcomp}
        for  count,i in enumerate(Signcomp):
            if count==0:
                pvalsTime=ScheffePvalDict[i]
            else:
                tempDF = ScheffePvalDict[i]
                pvalsTime=pd.concat([pvalsTime,tempDF])

        _,qvals,_,_ = sts.multitest.multipletests(pvalsTime.values.reshape(-1), alpha=0.1,method='fdr_bh')
        #pd.DataFrame(pvalsTime.to_excel(save_dir+'Scheffe_pvalues_Time.xlsx'))
        #pd.DataFrame(qvals.reshape(pvalsTime.shape)).to_excel(save_dir+'Scheffe_qvalues_Time.xlsx')

        Signcomp_WTOB = list(pvals.columns[pvals.loc['WTOB']=='*'])
        ScheffePvalDict = {metasample:sp.posthoc_scheffe(data4anova,val_col=metasample,group_col='WTOB') for metasample in Signcomp_WTOB}
        for  count,i in enumerate(Signcomp_WTOB):
            if count==0:
                pvalsWTOB=ScheffePvalDict[i]
            else:
                tempDF = ScheffePvalDict[i]
                pvalsWTOB=pd.concat([pvalsWTOB,tempDF])

        _,qvals,_,_ = sts.multitest.multipletests(pvalsWTOB.values.reshape(-1), alpha=0.1,method='fdr_bh')

        ##Corrects FDRs tested by multiple factors
        pvalsComb = np.concatenate([pvalsTime.values.reshape(-1),pvalsWTOB.values.reshape(-1)])
        _,qvals,_,_ = sts.multitest.multipletests(pvalsComb, alpha=0.1,method='fdr_bh')
        qvalsTime=qvals[:len(pvalsTime.values.reshape(-1))]
        qvalsTimeDF=pd.DataFrame(qvalsTime.reshape(pvalsTime.shape),index=[metasample+'_'+ str(timepoint[i]) for metasample in Signcomp for i in range(len(timepoint)) ], columns=[str(timepoint[i]) for i in range(len(timepoint)) ])
        #qvalsTimeDF.to_excel(save_dir+'Scheffe_qvalues_Time_comb.xlsx')
        print('Scheffe_qvalues_Time')
        display(qvalsTimeDF)

   
        qvalsWTOB=qvals[len(pvalsTime.values.reshape(-1)):len(pvalsTime.values.reshape(-1))+len(pvalsWTOB.values.reshape(-1))]
        WTOB = ['WT', 'OB']
        qvalsWTOBDF=pd.DataFrame(qvalsWTOB.reshape(pvalsWTOB.shape),index=[metasample+'_'+ WTOB[i]  for metasample in Signcomp_WTOB for i in range(2)], columns=[WTOB[i] for i in range(2) ])
        #qvalsWTOBDF.to_excel(save_dir+'Scheffe_qvalues_WTOB_comb.xlsx')
        print('Scheffe_qvalues_WTOB')
        display(qvalsWTOBDF)

        
        ### Steel Dwass
        Signcomp = list(pvals.columns[pvals.loc['Time']=='*'])
        SDPvalDict = {metasample:sp.posthoc_dscf(data4anova,val_col=metasample,group_col='Time') for metasample in Signcomp}
        for  count,i in enumerate(Signcomp):
            if count==0:
                pvalsTime=SDPvalDict[i]
            else:
                tempDF = SDPvalDict[i]
                pvalsTime=pd.concat([pvalsTime,tempDF])

        _,qvals,_,_ = sts.multitest.multipletests(pvalsTime.values.reshape(-1), alpha=0.1,method='fdr_bh')
        #pd.DataFrame(pvalsTime.to_excel(save_dir+'SD_pvalues_Time.xlsx'))
        #pd.DataFrame(qvals.reshape(pvalsTime.shape)).to_excel(save_dir+'Scheffe_qvalues_Time.xlsx')

        Signcomp_WTOB = list(pvals.columns[pvals.loc['WTOB']=='*'])
        SDPvalDict = {metasample:sp.posthoc_dscf(data4anova,val_col=metasample,group_col='WTOB') for metasample in Signcomp_WTOB}
        for  count,i in enumerate(Signcomp_WTOB):
            if count==0:
                pvalsWTOB=SDPvalDict[i]
            else:
                tempDF = SDPvalDict[i]
                pvalsWTOB=pd.concat([pvalsWTOB,tempDF])

        _,qvals,_,_ = sts.multitest.multipletests(pvalsWTOB.values.reshape(-1), alpha=0.1,method='fdr_bh')
        #pd.DataFrame(pvalsWTOB.to_excel(save_dir+'SD_pvalues_WTOB.xlsx'))
        #pd.DataFrame(qvals.reshape(pvalsWTOB.shape)).to_excel(save_dir+'Scheffe_qvalues_WTOB.xlsx')

        ##Corrects FDRs tested by multiple factors
        pvalsComb = np.concatenate([pvalsTime.values.reshape(-1),pvalsWTOB.values.reshape(-1)])
        _,qvals,_,_ = sts.multitest.multipletests(pvalsComb, alpha=0.1,method='fdr_bh')
        qvalsTime=qvals[:len(pvalsTime.values.reshape(-1))]
        qvalsTimeDF=pd.DataFrame(qvalsTime.reshape(pvalsTime.shape),index=[metasample+'_'+ str(timepoint[i]) for metasample in Signcomp for i in range(len(timepoint)) ], columns=[str(timepoint[i]) for i in range(len(timepoint)) ])
        #qvalsTimeDF.to_excel(save_dir+'SD_qvalues_Time_comb.xlsx')

   
        qvalsWTOB=qvals[len(pvalsTime.values.reshape(-1)):len(pvalsTime.values.reshape(-1))+len(pvalsWTOB.values.reshape(-1))]
        WTOB = ['WT', 'OB']
        qvalsWTOBDF=pd.DataFrame(qvalsWTOB.reshape(pvalsWTOB.shape),index=[metasample+'_'+ WTOB[i]  for metasample in Signcomp_WTOB for i in range(2)], columns=[WTOB[i] for i in range(2) ])
        #qvalsWTOBDF.to_excel(save_dir+'SD_qvalues_WTOB_comb.xlsx')        

        return(pvals)
    
    def enrichr(self,gene_list, gene_set_library):#TRANSFAC and JASPAR PWMs
        import json
        import requests
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
        genes_str = '\n'.join(gene_list)
        description = 'Example gene list'
        payload = {
            'list': (None, genes_str),
            'description': (None, description)
        }
    
        response = requests.post(ENRICHR_URL, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')
    
        data = json.loads(response.text)
    
        # enrich
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
        query_string = '?userListId=%s&backgroundType=%s'
        user_list_id = data["userListId"]
        response = requests.get(
            ENRICHR_URL + query_string % (user_list_id, gene_set_library)
         )
        if not response.ok:
            raise Exception('Error fetching enrichment results')
        data = pd.DataFrame(json.loads(response.text)[gene_set_library], columns=["Index","Term","P-value","Odds Ratio","Combined Score", "Genes", "Adjusted P-value", "Old P-value", "Old Adjusted P-value"]).set_index("Index")
        return data   
    
    def metagenes_topBH(self,metagenes, qvalcutoff, namedict,LiverNewDF, OptionDict):
        from statsmodels.stats.multitest import local_fdr
        if 'ENS' in list(metagenes.index)[0]:
            import mygene
            mg = mygene.MyGeneInfo()
            ens = list(metagenes.index)#['ENSG00000148795', 'ENSG00000165359', 'ENSG00000150676']
            #ginfo = mg.querymany(ens, scopes='ensembl.gene')
            geneSyms = mg.querymany(ens , scopes='ensembl.gene', fields='symbol', species='mouse')
            Name=[]
            for j in range(len(geneSyms)):
                try:
                        Name += [ geneSyms[j]['symbol'] ]
                except:
                    Name += [ ens[j]]
            metagenes.index = Name
            LiverNewDF.index = Name
            OptionDict['std'].index = Name
            
        else:          
            Name = list(metagenes.index)
            #[namedict[i] for i in list(metagenes.index)] ## ID　convert
            metagenes.index = Name
            OptionDict['std'].index =list(metagenes.index)
            LiverNewDF.index = list(metagenes.index)
        
        modules = dict()
        metagenes = self.mkZscore(metagenes,list(metagenes.index),list(metagenes.columns), 'col')
        for ii, metagene in enumerate(np.array(metagenes.T)):
        
            Pvalue = sts.chi2.sf(metagenes.iloc[:,ii].astype(np.float32)**2,1)
            _,QvalueBH,_,_ = stsmodel.multitest.multipletests(Pvalue, alpha=0.1,method='fdr_bh')

            QvalueBH=pd.DataFrame(QvalueBH)
            QvalueBH.index= list(metagenes.index)
            genelist = list(QvalueBH[QvalueBH[0]<qvalcutoff].index)#[namedict[i] for i in list(QvalueBH.index)]
            modules[metagenes.columns[ii]] = genelist

            for j in ['KEGG_2019_Mouse','TRANSFAC_and_JASPAR_PWMs','MSigDB_Computational','Disease_Perturbations_from_GEO_down','Disease_Perturbations_from_GEO_up','GO_Biological_Process_2021']:
                try:
                    print(self.enrichr(genelist,j))
                except:
                    print('error_' +str(ii+1)+j)
        return( modules)        
            #TensorClss.Reconstruction_Timecourse_miceFasting(LiverNewDF.loc[genelist],OptionDict,'_'+str(ii+1))

    
    def mkhistList(self,List1,xlabel,filename,bins,title):#1つのリストのヒストグラム
            fig = plt.figure(figsize=(6.4,4.8))
            plt.hist(List1,bins=bins,ec='black',align='left' )#,cumulative=1,normed=True)#plt.hist(list1[0:50],bins = 5,normed=True)#, bins = 20)  
            #alin 各棒の中心を X 軸目盛上のどの横位置で出力するか。 ‘left‘,‘mid‘,‘right‘ から選択。デフォルト値: ‘mid’
            #plt.xlim([0,1]) #モデルパラメタのとき
            plt.xlabel(xlabel,fontsize=20); plt.ylabel('Frequency',fontsize=20);plt.tick_params(labelsize=20);plt.title(title,size=20)
    ### temp 縦軸をlog10に
            #plt.yscale('log');#plt.xticks(np.arange(0, 110 + 1, 10)) #np.arange(-110, 20 + 1, 20) np.arange(0, 110 + 1, 10)
            #plt.savefig(save_dir +  filename +'Hist.pdf',bbox_inches="tight")  ;plt.close()#print(ConList2[4][~np.isnan(ConList2[4])])
    

    def mkZscore(self,DF,IndLabel,ColLabel,axis):
        if axis=='col':##zscore in column direction
            ax=0
        else:##zscore in row direction
            ax=1
        meanNP = np.array(DF) -np.nanmean(np.array(DF),axis=ax,keepdims=True)
        stdNP = np.nanstd(np.array(DF),axis=ax)
        ZscoredDF = pd.DataFrame(index=IndLabel,columns=ColLabel)
        if axis=='col':
            for i in range(len(ColLabel)): 
                ZscoredNP = meanNP[:,i]/stdNP[i]
                ZscoredDF[ColLabel[i]] = ZscoredNP  
        else:    
            for i in range(len(IndLabel)): 
                ZscoredNP = meanNP[i,:]/stdNP[i]
                ZscoredDF.loc[IndLabel[i]] = ZscoredNP
        return(ZscoredDF) 
    def plotBar(self,List1, Title, xlabel, ylabel, Color, xticks, size, xsize,Titlesize) : 
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

        ax1.set_xticklabels(labels=xticks,rotation=270,fontsize=xsize)    
        xmin, xmax, ymin, ymax = ax1.axis() 


            
    def PypeR_MetabanalystEnrich(self,pcorr_mtr, SwitchDict):
        import pyper        
        # R のインスタンスを作る
        r = pyper.R(use_numpy = 'True', use_pandas='True')
        
        r.assign("data", pcorr_mtr)
        r.assign("filename", SwitchDict['file_name'])
        data=r.get("data")
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
        result = r.get("save_dir")
        print(result)
        r('cov_select<-function(model=x, cov=y, n.obs=0, AIC=0, cov_orig=NULL){ \
          require(ggm) \
          if(is.null(cov_orig)){cov_orig<-cov} \
          #  print(sprintf("AIC= %.4f", AIC ),quote=F) \
          # 偏相関行列を作成 \
          pmat<- (-1)*solve(cov) / sqrt(diag(solve(cov)) %*% t(diag(solve(cov)))) \
          diag(pmat)<- 1 \
          #偏相関係数を絶対値に変換\
          amat<-abs(pmat) \
          #モデルの係数が0の箇所を無限大に設定 \
          amat[which(model==0)]<-Inf \
          #偏相関の絶対値が最小の要素を0にした修正モデルを作成\
          model_post<-rep(1,nrow(cov))%*%t(rep(1,nrow(cov)))\
          model_post[which.min(amat)]<-0 \
          model_post<-model_post * t(model_post) * model \
          # モデルのフィットとAICの算出 \
          f<-fitConGraph(model_post,cov_orig,n.obs) \
          AIC_post<-f$dev-2*f$df \
          # モデルの適合度が最大になるまで反復 \
          if (AIC_post<AIC){ \
            Recall(model_post,f$Shat,n.obs,AIC=AIC_post,cov_orig=cov_orig) \
          } \
          #最終的に得られたモデルを描画 & 偏相関行列を表示 \
          else{ \
            diag(pmat)<-1 \
            pmat[which(model==0)]<-0 \
            model.0<-model*0 \
            f<-fitConGraph(model,cov_orig,n.obs) \
            f.0<-fitConGraph(model.0,cov_orig,n.obs) \
            # GFIの算出 \
            S<-cov2cor(cov_orig) \
            H<-cov2cor(f$Shat) \
            p<-nrow(cov_orig) \
            num<-sum(diag((solve(H)%*%(S-H))%*%(solve(H)%*%(S-H)))) \
            den<-sum(diag((solve(H)%*%S)%*%(solve(H)%*%S))) \
            GFI<-1-(num/den) \
            AGFI<-1-(p*(p+1)*(1-GFI))/(2*f$df) \
            RMSEA<-sqrt(max(((f$dev-f$df)/(f$df*(nrow(cov_orig)-1))),0)) \
            CFI<-(f.0$dev*f.0$df-f$dev*f$df)/(f.0$dev*f.0$df)\
            res<-f \
            res<-c(f,GFI=GFI, AGFI=AGFI,RMSEA=RMSEA,CFI=CFI) \
            return(list(fit=res,model=model, covmat=pmat)) \
          } \
        }')   
        
        r('rating_cov<-data')
        r('rating<-data')
        r('ncov <- ncol(rating)')
        r('model_post <- matrix(rep(1, ncov^2), nrow=ncov, ncol=ncov)')
        r('diag(model_post) <- 0')
        r('colnames(model_post) <- names(rating)')
        r('rownames(model_post) <- names(rating)')
        r('model <- cov_select(model_post, rating_cov, ncov)')
        
        r('  Selectedmtr <- matrix(model$covmat,ncol=length(colnames(rating)), nrow=length(rownames(rating)))')
        r('  colnames(Selectedmtr) <- names(rating)')
        r('rownames(Selectedmtr) <- names(rating)')
        result = r.get("Selectedmtr")
        print(result)
    
        r('write.csv(Selectedmtr,save_dir + "/Selectedmtr" + file_name1 + ".csv",append=T, quote=F, col.names=F,fileEncoding = "CP932")')
        r('write.xlsx(Selectedmtr, file=save_dir + "/Selectedmtr"+ file_name1 + ".xlsx", sheetName="sheet1", row.names=T)')
    
        r('qgraph(model$covmat, edge.labels=T, minimum=.2,theme = "Borkulo")')
        
        r('pdf(save_dir + "CovarianceSelection_Network.pdf")')
        #r('plot(qobj$lambda, qobj$pi0.lambda,xlab=xlab,ylab=ylab)')
        #r('title(main = bquote(π0==.(qobj$pi0)))')
        r('dev.off()')
        r("source(file='../../R/CovarianceSelection.R')")
        
    def Reconstruction_Timecourse_miceFasting(self,aaDF,OptionDict,save_dir):
        #NewhumanOGTTDF = pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20200727/NewhumanOGTTDF.xlsx',header=0,index_col=0)   
    
        SubjectName = ['1']*2;#sorted(set(SubjectNameSet), key=SubjectNameSet.index)
        Timepoint = [0,2,4,6,8,12,16,24]
        
        ####
        #aaDF=DF.loc[:,(DF!=0).any()]
        #aaDF = DF.loc[list(tensorDF.sort_index().index)]
        #OptionDict['std']=OptionDict['std'].loc[list(tensorDF.sort_index().index)]
        ###
        #aaDF = aaDF.T
        OptionDict['std'].index=[str(OptionDict['std'].index[i]) for i in range(len(OptionDict['std'].index))]
        ###
        TCTClass = TensorClassTimeCourse()
        
        #aaDF,aaMolList = TCTClass.metabolome_time(tempMolDF)#metabolomeデータ、時間成分で分ける
        #aaDF,aaMolList = TCTClass.metabolome_subject(tempMolDF)#metabolomeデータ、個人成分で分ける
        #aaDF,aaMolList = TCTClass.metabolome_all(tempMolDF)
       ### aaDF,aaMolList = TCTClass.metabolome_all([])#metabolomeデータ、どの都度設定
                    #d = r.search(j)     
        #r.search(jj).group(3)
                        #r = re.compile("(.*)(_)(.*)") 
                    #d = r.search(j) 
    
        #aaDF = aaDF.drop('name',axis=1);aaDF = aaDF.drop('Label',axis=1);
        #OptionDict={}
        OptionDict['switch']='miceFasting_mean'#'BC':#BCの順に描画なら 'BC_all'#まとめ　のやつ、'metabolome'Bolus20 'metabolome_ColorComp':Bolu20, 成分で色分け _Absで絶対値 'miceOGTT_mean', 'miceOGTT'
        #OptionDict['Compvalue'] = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20201213/metabolome/compsubject.xlsx',header=0,index_col=0)['X6'])
        ###最新版_20210314
        #OptionDict['Compvalue'] = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210314/metabolome/compsubject.xlsx',header=0,index_col=0)['X2'])
        #OptionDict['Compvalue'] = (OptionDict['Compvalue'] - np.min(OptionDict['Compvalue'])) / (np.max(OptionDict['Compvalue'])-np.min(OptionDict['Compvalue']))
        #OptionDict['CompvalueTrue'] = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20201213/metabolome/compsubject.xlsx',header=0,index_col=0)['X6'])
         ###最新版_20210314
        #OptionDict['CompvalueTrue'] = list(pd.read_excel('/Users/fujita/Google ドライブ/Kuroda lab/Research/Metabolome/result/Property/20210314/metabolome/compsubject.xlsx',header=0,index_col=0)['X2'])
        
        aaDF['sort'] = list(aaDF.index)
        aaDF=aaDF.sort_values(by='sort',ascending=True) 
        #OptionDict['CompvalueTrue'] = list(np.sort(OptionDict['CompvalueTrue']))
        #OptionDict['Compvalue']= list(np.sort(OptionDict['Compvalue']))
        aaDF=aaDF.drop('sort',axis=1)
        
        OptionDict['Label'] = list(aaDF.columns)
        ylim_l, ylim_h, TitleDF = TCTClass.calcylim(aaDF,OptionDict)   
        aaDF.to_excel(save_dir+'allconcatDF.xlsx')
        try:
            OptionDict['std'].to_excel(save_dir+'allconcatDF_std.xlsx')
        except:
            pass
        TCTClass.pltoTimeCourse(aaDF,Timepoint,SubjectName,ylim_l,ylim_h,[],OptionDict,save_dir)
        
        