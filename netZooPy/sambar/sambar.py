import numpy as np
import pandas as pd
import time
import os
import networkx as nx
from scipy.spatial.distance import pdist,cosine,squareform
from scipy.cluster.hierarchy import linkage,cut_tree
import pkg_resources
"""
Description:
    Python implementation of the Subtyping Agglomerated Mutations By Annotation Relations (SAMBAR) method as implemented in R https://github.com/mararie/SAMBAR.
    SAMBAR, or Subtyping Agglomerated Mutations By Annotation Relations, is a method to identify subtypes based on somatic mutation data. SAMBAR was used to identify mutational subtypes in 23 cancer types from The Cancer Genome Atlas (Kuijjer ML, Paulson JN, Salzman P, Ding W, Quackenbush J, *British Journal of Cancer* (May 16, 2018), doi: 10.1038/s41416-018-0109-7, https://www.nature.com/articles/s41416-018-0109-7, BioRxiv, doi: https://doi.org/10.1101/228031).

Usage:
    To use the package you can import it using: ```import pysambar```. And then access the different functions implemented with ```pysambar.function()```.
    As an example you can find mutation data of Uterine Corpus Endometrial Carcinoma (UCEC) primary tumpor samples from The Cancer Genome Atlas. This data is in the ToyData folder as well as the MSigDb "Hallmark" gene sets. 
    The program will compute the pathway mutation scores and clustering for *k*=2-4 (by default) and output the corrected mutation scores, the pathway mutation scores, and the clustering table. 

Example:
    Run the SAMBAR method with the ToyData from the UCEC mutation data:
    ``` 
    import pysambar as sm
    pathways, groups = sm.sambar("/ToyData/mut.ucec.csv","ToyData/esizef.csv",'ToyData/genes.txt','ToyData/h.all.v6.1.symbols.gmt')
    ``` 
    The output of this command will be four files:
    ``` pt_out.csv```  -> Pathway mutation score matrix.
    ``` mt_out.csv```  -> Processed gene mutation score matrix.
    ``` dist_matrix.csv ``` -> Distance matrix with binomial distance in numpy condensed format.
    ``` clustergroups.csv```  -> Matrix of pertinence to a cluster.
    The function also returns the pathway matrix dataframe and the cluster group dataframe as python variables.

Notes:
    Flags by default in the sambar function:
    ``` normPatient=True```  -> Normalizes the mutation data by number of mutations in a sample.
    ``` kmin=2,kmax=4```  -> Cut-offs of the cluster tree.
    ``` gmtMSigDB=True```  -> If the signature file comes from MSigDB the second element of each line is removed to process the format available at the database. This can be toggled off if using a custom signature file.
    ``` subcangenes=True```  -> Subset to cancer associated genes. By default uses the file provided in ToyData.

Functions:
            This package includes the functions ```sambar```,```desparsify```,```corgenelength```,```convertgmt```,```clustering``` as well as an implementation of the binomial distance (Millar dissimilarity from the package vegdist from R. To see the full description of each of this functions use ```help(pysambar.function)```.
"""
## Default toydata files
esize = pkg_resources.resource_filename('netZooPy', 'ToyData/esizef.csv')
genes = pkg_resources.resource_filename('netZooPy', 'ToyData/genes.txt')
sign  = pkg_resources.resource_filename('netZooPy', 'ToyData/h.all.v6.1.symbols.gmt')
mut   = pkg_resources.resource_filename('netZooPy', 'ToyData/mut.ucec.csv')

def corgenelength(mut,cangenes,esize,normbysample=True,subcangenes=True):
    """
    Description:
        Function to normalize gene mutation scores by gene length.
        mut should be a dataframe of mutation scores with genes as columns and samples as rows. 
        (VERY IMPORTANT, IF OTHERWISE MATRIX SHOULD BE TRANSPOSED OR IT WON'T WORK!!!)
    
    Inputs:
        mut         : Mutation scores.
        cangenes    : A set of cancer associated genes.
        esize       : A dataframe of gene lengths.
        normbysample: True : Normalizes the gene mutation scores in a sample by the total mutations within the sample. 
                      False: Deactivate normalization.
        subcangenes : True: Subsets mutation data to cancer-associated genes.
                      False: Takes all genes.

    Outputs:
        mut         : Mutation scores normalized by gene length.
    """
    #Subsets mutation data to cancer-associated genes.
    if(subcangenes):
    	mut_cangenes = cangenes.intersection(set(mut.columns))
    	mut = mut[list(mut_cangenes)]
    else: mut_cangenes = set(mut.columns)
    #Intersection of genes in the mutation data and the length data.
    mut_esize = list(mut_cangenes.intersection(set(esize.columns)))
    mut = mut[mut_esize]
    esize = esize[mut_esize]
    
    #Divide each number of mutations by gene size
    mut = mut.div(esize.iloc[0])
    
    #Drop samples (rows) without mutations (each column equals 0)
    mut = mut.drop(list(mut[(mut == 0).all(1)].index))
    
    #Normalize mutation score by patient number of mutations
    if(normbysample):
        mut["sum"]=mut.sum(axis=1) #Adds auxiliar column with the sum of columns (sum of mutation scores within a sample)
        mut = mut.div(mut["sum"],axis=0) #Corrects gene mutation scores in a sample by the total mutation score in the sample.
        mut = mut.drop(columns=["sum"])
    
    #Returns mut matrix of gene mutation scores by sample adjusted by length and total number of mutations within sample.
    #Sorted columns for niceness
    assert mut.shape!=(0,0),"Are you sure your mutation matrix has the proper orientation? Genes as columns, samples as rows."
    return mut[sorted(mut.columns)]

def convertgmt(gmtfile, cangenes,gmtMSigDB=True,subcangenes=True):
    """
    Description:
        This function takes as input the name of a gmt file containing lists of genes associated to pathways. 
        It outputs an adjacency matrix of genes and pathways. 
        It also subsets the genes to a list of cancer-associated genes. 

    Inputs:
        gmtfile    : Path the gmt file.
        cangenes   : A set of cancer associated genes.
        gmtMSigDB  : True : gmt file from MSigDB 
                     False: file not form MSigDB
        subcangenes: True : Subsets mutation data to cancer-associated genes.
                     False: Takes all genes.

    Outputs:
        sign_matrix: Adjacency matrix of genes and pathways.
    """
    file = open(gmtfile) # Loads the gmt file
    pw = file.readlines() #Reads file into list of strings

    a = dict() #Sets up dictionary
    for entry in pw:
        u = entry[:-1].split("\t",-1) #Process each string separating genes and removing final \n.
        if(gmtMSigDB):
        	a[u[0]] = u[2:] #Sets dictionary key the pathway id, value the list of genes removing the link entry.
        else: a[u[0]] = u[1:]
    #Set of all genes in the dictionary
    b=set()
    for i in a.keys():
        b = b.union(set(a[i]))

    C = sorted(list(b)) #List of genes in the file (sorted)
    P = sorted(a.keys()) #List of pathways in the file (sorted)
    
    G = nx.Graph(a) #Creates a bipartite graph from the dictionary
    A = nx.adjacency_matrix(G,nodelist=(C+P)).todense() #Creates the adjacency matrix of the graph G in the proper order. 
    #Note that this matrix is square with size len(C)+len(P), the biadjacency matrix caused issues with big files.
    
    #Builds a dataframe from the adjacency matrix and tags the rows and columns.
    sign_matrix = pd.DataFrame(A)
    sign_matrix.columns=C+P
    sign_matrix.index=C+P
    
    #Subsets the dataframe so it has as columns the genes and as rows the pathways. Orders for niceness.
    sign_matrix = sign_matrix[sorted(C)]
    sign_matrix = sign_matrix.loc[sorted(P)]
 
    #Subset pathways to cancer-associated genes
    if(subcangenes):
    	sign_cangenes = cangenes.intersection(set(sign_matrix.columns))
    	sign_matrix = sign_matrix[list(sign_cangenes)]
    
    return sign_matrix[sorted(sign_matrix.columns)]

def desparsify(mutdata,exonsize,gmtfile,cangenes,normMut=True,gmtMSigDB=True,subcangenes=True):
    """
    Description:
        Applies the sambar method to de-sparcify the mutation data using the pathway signatures in the gmtfile.
        
    Inputs:
        mutdata     : Mutation scores.
        exonsize    : A dataframe of gene lengths.
        gmtfile     : Path the gmt file.
        cangenes    : A set of cancer associated genes.
        normbysample: True : Normalizes the gene mutation scores in a sample by the total mutations within the sample. 
                      False: Deactivate normalization.
        gmtMSigDB   : True : gmt file from MSigDB 
                      False: file not form MSigDB
        subcangenes : True : Subsets mutation data to cancer-associated genes.
                      False: Takes all genes.

    Outputs:
        mt            : Genes in both pathways and mutation data.
        pathway_scores: Pathway scores.
    """
    tinit = time.time()  
    
    #Gets corrected gene mutation scores. Subseted to cancer-associated genes.
    mutrate = corgenelength(mutdata,cangenes,exonsize,normMut,subcangenes)
    #Gets adjacency matrix of genes and pathways. Subseted to cancer-associated genes.
    sign_matrix = convertgmt(gmtfile,cangenes,gmtMSigDB,subcangenes)

    #Match genes in both lists
    genes_ms = sorted(list(set(mutrate.columns).intersection(set(sign_matrix.columns))))
    ## Computes the number of genes in a pathway (before matching with mutation data).
    sm = sign_matrix.sum(1)
    
    #Matches genes in mutation data with genes in pathway data
    mutrate = mutrate[genes_ms]
    sign_matrix = sign_matrix[genes_ms]
    mt = mutrate
    #Corrects mutation rates by frequency in a pathway
    genefreq = list(sign_matrix.sum())
    mutrate = mutrate.div(genefreq).transpose()
    
    #This should be done here, otherwise the matching before checkpoint 2 might generate rows with only 0. SHOULD NOT BE COMMENTED
    sign_matrix = sign_matrix.drop(list(sign_matrix[(sign_matrix == 0).all(1)].index))
    sign_matrix = sign_matrix.transpose()
    sign_matrix = sign_matrix.drop(list(sign_matrix[(sign_matrix == 0).all(1)].index))
    sign_matrix = sign_matrix.transpose()
    
    #Matrix multiply dataframes to get the pathway mutation scores
    pathway_scores = sign_matrix.dot(mutrate).div(sm,axis=0)
    pathway_scores = pathway_scores.dropna() #Removes pathways with NaN values.
    #Remove rows and columns with only 0.
    pathway_scores = pathway_scores.drop(list(pathway_scores[(pathway_scores == 0).all(1)].index))
    pathway_scores = pathway_scores.transpose()
    pathway_scores = pathway_scores.drop(list(pathway_scores[(pathway_scores == 0).all(1)].index))
    pathway_scores = pathway_scores.transpose()
    
    tfin = time.time()
    print("Sambar runtime: ",tfin-tinit)
    
    return mt,pathway_scores

def binomial_dist(u,v):
    """
    Description:
        Implementation of the binomial dissimilarity funcion or Millar distance from the vegan:vegdist package in R.

    Inputs:
        u : First vector to compare distance from.
        v : Second vector to compare distance to.

    Outputs:
        bd: binomial dissimilarity.
    """
    np.seterr(divide='ignore', invalid='ignore') #The following steps raise divide by zero errors. This avoids the error to show up as it is not a problem since it is processed afterwards.
    x,y= u,v
    nk = x+y
    lognk = np.nan_to_num(np.log(nk))
    t1 = np.nan_to_num((x*(np.log(x)-lognk)))
    t2 = np.nan_to_num((y*(np.log(y)-lognk)))
    dn = t1+t2+nk*np.log(2)
    d  = sum(np.nan_to_num((dn/nk)))
    bd = np.where(d<0, 0, d)
    return bd

def clustering(pt, kmin,kmax,distance,linkagem):
    """
    Description:
        Computes the clustering for the pathways matrix and returns a dataframe with the groups with k clusters from kmin to kmax.

    Inputs:
        pt      : Pathway scores.
        kmin    : Min number of groups in the clustering.
        kmax    : Max number of groups in the clustering.
        distance: Similarity metric for the clustering. Default is binomial distance but any distance from scipy.spatial.distance can be used.
        linkagem: Linkage method. Default is complete.

    Outputs:
        df      : Cluster assignement dataframe.    
    """
    tinit = time.time()
    if(distance=="binomial"):
        Y = pdist(pt.transpose(),binomial_dist) # Computes the distance matrix. pt is transposed because the pdist function takes rows as input and what we want to cluster are the samples.
    else:
        Y = pdist(pt.transpose(),distance)
    Z = linkage(Y,linkagem) #Linkage of the clusters using the distance matrix and the complete method.
    np.savetxt("dist_matrix.csv",squareform(Y),delimiter=",")#Saves the distance matrix. Note that the output of pdist is a condensed matrix!!

    #Building the output dataframe.
    df = pd.DataFrame()
    for k in range(kmin,kmax+1):
        R = cut_tree(Z,k) #Cuts the tree at k groups.
        u = [item for sublist in R for item in sublist]
        df["X"+str(k)] =u
        df.index = pt.columns
    tfin = time.time()
    print("Clustering runtime: ",tfin-tinit)
    return df.transpose()

def sambar(mut_file=mut,esize_file=esize,genes_file=genes,gmtfile=sign,normPatient=True,kmin=2,kmax=4,gmtMSigDB=True,subcangenes=True,distance="binomial",linkagem="complete",cluster=True):
    """
    Description:
        Runs SAMBAR and outputs the pt matrix, the mt matrix and the clustering matrix.

    Inputs:
        mutfile    : matrix of mutations with genes as columns and samples as rows. Format CSV.
        esize_file : file with genes and their length.
        genes_file : file with list of cancer-associated genes.
        gmt_file   : genelist by pathway, format from MSigDB.
        normPatien : Normalize mutation data by number of mutations in a sample.
        kmin,kmax  : Number of groups in the clustering.
        gmtMSigDB  : Whether the signature file comes from MSigDB or not. (Important for processing the file).
        subcangenes: Makes optional subsetting to cancer-associated genes.
        distance   : Similarity metric for the clustering. Default is binomial distance but any distance from scipy.spatial.distance can be used.
        linkagem   : Linkage method. Default is complete.
        cluster    : Whether the clustering has to be compute or the output will just be the pathway mutation scores.
    
    Outputs:
        pt               : the pathway matrix.
        groups           : the groups dataframe as python objects.
        mt_out.csv       : processed gene mutation scores.
        pt_out.csv       : pathway mutation scores
        clustergroups.csv: matrix of pertinence to a group in the clustering.
        dist_matrix.csv  : Computation of the distance matrix is resource-consuming so the matrix is writen so it doesn't have to be computed again.

    Reference:
        Kuijjer, Marieke Lydia, et al. "Cancer subtype identification using somatic mutation data." British journal of cancer 118.11 (2018): 1492-1501.
    """
    mut = pd.read_csv(mut_file,index_col=0)
    esize = pd.read_csv(esize_file,index_col=0)
    file = open(genes_file) # Loads the cancer associated gene list file
    gene_line = file.readline()
    genes = set(gene_line.split("\t"))    
    mt,pt = desparsify(mut,esize,gmtfile,genes,normPatient,gmtMSigDB,subcangenes)
    
    pt.to_csv("pt_out.csv")
    mt.to_csv("mt_out.csv")
    
    if (not cluster): return pt
    
    groups = clustering(pt,kmin,kmax,distance,linkagem)
    groups.to_csv("clustergroups.csv")
    return pt,groups