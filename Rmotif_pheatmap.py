import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import re,os,sys
#from sklearn import preprocessing
from matplotlib.backends.backend_pdf import PdfPages

#def Standarscaler(data):
#    data=(data-data.mean())/data.std()
#    return data
#def minmaxscale(data):
#    data=(data-data.min())/(data.max()-data.min())
#    return data
#min_max_scaler = preprocessing.MinMaxScaler()

#np.random.seed(0)

def Rmotifplot(matrixfile,Rmotifmatrix,pdfout):
    pp = PdfPages(pdfout)
    data=pd.read_table(matrixfile,sep='\t',header=0)
    singal=pd.read_table(Rmotifmatrix,sep='\t',header=0)

    #uniform_data=min_max_scaler.fit_transform(data)
    #data1=minmaxscale(data['p50Ob'])
    #data2=minmaxscale(data['p50U1'])
    uniform_data=data
    #cmap = sns.cubehelix_palette(start = 1.5, rot = 3, gamma=0.8, as_cmap = True)

    #f, (ax4,ax1,ax2,ax3) = plt.subplots(figsize = (3.5,5),ncols=4)
    plt.figure(figsize=(3.5,5))
    grid = plt.GridSpec(1, 6, wspace=0.5, hspace=0.5)
    ax1=plt.subplot(grid[0,1:4])
    ax2=plt.subplot(grid[0,4])
    #sns.heatmap(uniform_data, ax = ax1,cmap='YlOrRd',cbar_kws={"orientation": "vertical", "shrink":0.40, "aspect":40})
    sns.heatmap(uniform_data, ax = ax1,cmap='YlOrRd',cbar_kws = {"use_gridspec":False, "location":"left","shrink":0.40})

    ax1.set_title('expression')
    #ax1.set_xlabel('')
    #ax1.set_xlabel('picture1')
    #ax1.set_xticklabels([]) #设置x轴图例为空值
    #ax1.set_ylabel('kind')
    sns.heatmap(singal, ax = ax2,cmap='PiYG',vmin=-2,vmax=2,cbar=False)
    #ax2.set_title('R_motif')
    ax2.set_xlabel('region')
    #ax2.set_ylabel('kind')
    ax2.set_ylabel('')
    ax2.set_yticklabels([])
    ax2.set_yticks([])
    ax1.set_yticks([])
    label_ax1 = ax1.get_xticklabels() 
    plt.setp(label_ax1, rotation=45, horizontalalignment='right')
    label_ax2 = ax2.get_xticklabels() 
    plt.setp(label_ax2, rotation=45, horizontalalignment='right')

    pp.savefig()
    pp.close()


directory1=sys.argv[1]
for eachfile in os.listdir(directory1):
    if re.search(r'\.order\.matrix$',eachfile):
        matrixfile=eachfile
        Rmotifmatrix=eachfile+'.Rmotif'
        pdfout=re.sub(r'.normalize.order.matrix','',eachfile)+'_Rmotif.pdf'
        Rmotifplot(matrixfile,Rmotifmatrix,pdfout)


