"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Copyright (c) 2024 [Mengtian Li]
[Capital University of Economics and Business]. All rights reserved.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
from scipy.sparse import linalg,csr_matrix,csc_matrix,diags
import networkx as nx
from FindChimericRead import *
import pandas as pd

def importLongReads(filename,sep='.'):
    ReadSeq={}    
    for seq_record in SeqIO.parse(filename, 'fasta'):
        read_id = seq_record.id.split(sep)[1]
        ReadSeq[read_id+"_F"] = str(seq_record.seq).upper()
        ReadSeq[read_id+"_R"] = str(Seq(seq_record.seq).reverse_complement()).upper()        
    return ReadSeq

def writeSequence(f, sequence):
    length = len(sequence)
    line=int(length/100)
    for k in range(line):
        f.write(sequence[100*k:100*(k+1)]+'\n')
    if length>line*100:
        f.write(sequence[100*line:]+'\n')
        
def ImportAlignmentsFromMinimap(filename,HangingOut,filter_chimera=False, sep='.'):
    candidate_chimerias = []
    chimeria_cal = True
    G = nx.Graph()
    edge_dict = {}
    direct_dict = {}
    ReadLength = dict()
    with open(filename,'r') as fin:
        for line in fin:
            record = line.strip().split()
            (A_length,A_start,A_end) = list(map(int,record[1:4]))
            (B_length,B_start,B_end,align_bases) = list(map(int,record[6:10]))
            A_name = int(record[0].split(sep)[1])
            B_name = int(record[5].split(sep)[1])
            A_Orientation = 0
            B_Orientation = 0 if record[4]=='+' else 1
            align_bases = round(np.sqrt(align_bases)) #int(record[9])#
            # align_block_len = int(record[10])
                
            if filter_chimera: #filter candidate chimerias
                if chimeria_cal: ##只计算一次
                    candidate_chimerias_list,ReadLength = FindChimeria_minimap2(filename, sep=sep, Hangingout=HangingOut)
                    chimeria_cal = False
                    candidate_chimerias = collections.defaultdict(bool)
                    candidate_chimerias.update((x,True) for x in candidate_chimerias_list)
                    print("Filter Chimeric reads:", len(candidate_chimerias))
                if str(A_name) in candidate_chimerias or str(B_name) in candidate_chimerias:
                    continue
            else:
                ReadLength[A_name] = A_length
                ReadLength[B_name] = B_length
            if A_name==B_name:
                continue
                
            else: ## target start on original strand
                if B_Orientation==1:
                    t_B = B_start
                    B_start = B_length-B_end
                    B_end = B_length-t_B
                overlap={'A_name':A_name,'B_name':B_name,'A_start':A_start,'A_end':A_end,'A_length':A_length,'B_Orientation':B_Orientation,'B_start':B_start,'B_end':B_end,'B_length':B_length,'align_bases':align_bases}
                key = str(min(A_name,B_name))+"-"+str(max(A_name,B_name)) ##direct_dict不需要分AB的顺序
                e_key = str(A_name)+"-"+str(B_name)
                hanging = max(min(A_start,B_start),min(A_length-A_end,B_length-B_end))
                # print("HangingOut",hanging)
                if A_Orientation==B_Orientation:
                    ##权重绝对值越大可靠性越高
                    direct = align_bases
                else:
                    direct = -align_bases
                if hanging<HangingOut: #min(A_start,B_start)<HangingOut and min(A_length-A_end,B_length-B_end)<HangingOut:
                    if key in direct_dict:
                        direct_dict[key] += [direct]
                    else:
                        direct_dict[key] = [direct]
                    G.add_edge(A_name,B_name,weight = direct)
                    if e_key in edge_dict:
                        edge_dict[e_key].append(overlap)
                    else:
                        edge_dict[e_key]=[overlap]
    fin.close()
    return (edge_dict,direct_dict,G,ReadLength)

def setupRegressionModel_direct(reads,D, edge_dict):    
    row=0
    indptr=[0]
    indices=[]
    data=[]
    response=[]

    size=len(reads)
    response= []
    inconsist_link = 0
    for i in range(size-1):
        for j in range(i+1,size):
            A = reads[i]
            B = reads[j]
            key = str(A)+"-"+str(B)
            if key in edge_dict:
                for overlap in edge_dict[key]:##总是reads[j]['start']-reads[i]['start']，且以i的方向为正方向
                    if overlap['B_Orientation']==1 and D[i]!=D[j]: #异向overlap
                        if D[i]==1: #reads[i] --> reads[j] <--
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']
                        elif D[i]==-1: #reads[i] <-- reads[j] --> 
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']#overlap['B_length']-overlap['B_end']-overlap['A_length']+overlap['A_end']
                            diff = -diff
                        else:
                            print("Error!")
                    elif overlap['B_Orientation']==0 and D[i]==D[j]: #同向overlap                        
                        if D[i]==1: #reads[i] --> reads[j] -->    
                            diff = overlap['A_start']-overlap['B_start']
                        elif D[i]==-1:#reads[i] <-- reads[j] <--
                            diff = overlap['A_start']-overlap['B_start']#overlap['A_end']-overlap['A_length']+overlap['B_end']
                            diff = -diff
                        else:
                            print("Error!")
                        # print(A,B,D[i],D[j],diff)
                    else:#else的情况说明该overlap的方向与估计的方向关系不符，应该被滤除
                        inconsist_link +=1
                        # print("Inconsistent-direction overlap",overlap,D[i],D[j])
                        continue
                    row+=1
                    indices.append(j) 
                    data.append(1)
                    indices.append(i)
                    data.append(-1)                   
                    indptr.append(indptr[-1]+2)
                    response+=[diff]                             
    print("Inconsistent-direction overlaps:",inconsist_link)
    design=csr_matrix((data, indices, indptr), shape=(row, size))
    response=np.matrix(response).T

    return(design, response)

def setupRegressionModel_W(reads,D, edge_dict):    
    row=0
    indptr=[0]
    indices=[]
    data=[]
    response=[]
    W = []

    size=len(reads)
    response= []
    for i in range(size-1):
        for j in range(i+1,size):
            A = reads[i]
            B = reads[j]
            key = str(A)+"-"+str(B)
            if key in edge_dict:
                for overlap in edge_dict[key]:##总是reads[j]['start']-reads[i]['start']，且以i的方向为正方向
                    if overlap['B_Orientation']==1 and D[i]!=D[j]: #异向overlap
                        if D[i]==1: #reads[i] --> reads[j] <--
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']
                        elif D[i]==-1: #reads[i] <-- reads[j] --> 
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']#overlap['B_length']-overlap['B_end']-overlap['A_length']+overlap['A_end']
                            diff = -diff
                        else:
                            print("Error!")
                    elif overlap['B_Orientation']==0 and D[i]==D[j]: #同向overlap                        
                        if D[i]==1: #reads[i] --> reads[j] -->    
                            diff = overlap['A_start']-overlap['B_start']
                        elif D[i]==-1:#reads[i] <-- reads[j] <--
                            diff = overlap['A_start']-overlap['B_start']#overlap['A_end']-overlap['A_length']+overlap['B_end']
                            diff = -diff
                        else:
                            print("Error!")
                        # print(A,B,D[i],D[j],diff)
                    else:#else的情况说明该overlap的方向与估计的方向关系不符，应该被滤除
                        # print("Inconsistent-direction overlap",overlap,D[i],D[j])
                        continue
                    row+=1
                    indices.append(j) 
                    data.append(1)
                    indices.append(i)
                    data.append(-1)                   
                    indptr.append(indptr[-1]+2)
                    response+=[diff]
                    W+=[np.sqrt(overlap['align_bases'])]                  
                                        
    design=csr_matrix((data, indices, indptr), shape=(row, size))
    response=np.matrix(response).T

    return(design, response, W)

def setupRegressionModel_indel(reads,D, edge_dict,indel_ratio = 0.006875):    
    row=0
    indptr=[0]
    indices=[]
    data=[]
    response=[]

    size=len(reads)
    response= []
    inconsist_link = 0
    for i in range(size-1):
        for j in range(i+1,size):
            A = reads[i]
            B = reads[j]
            key = str(A)+"-"+str(B)
            if key in edge_dict:
                for overlap in edge_dict[key]:##总是reads[j]['start']-reads[i]['start']，且以i的方向为正方向
                    if overlap['B_Orientation']==1 and D[i]!=D[j]: #异向overlap
                        if D[i]==1: #reads[i] --> reads[j] <--
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']
                        elif D[i]==-1: #reads[i] <-- reads[j] --> 
                            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']#overlap['B_length']-overlap['B_end']-overlap['A_length']+overlap['A_end']
                            diff = -diff
                        else:
                            print("Error!")
                    elif overlap['B_Orientation']==0 and D[i]==D[j]: #同向overlap                        
                        if D[i]==1: #reads[i] --> reads[j] -->    
                            diff = overlap['A_start']-overlap['B_start']
                        elif D[i]==-1:#reads[i] <-- reads[j] <--
                            diff = overlap['A_start']-overlap['B_start']#overlap['A_end']-overlap['A_length']+overlap['B_end']
                            diff = -diff
                        else:
                            print("Error!")
                        # print(A,B,D[i],D[j],diff)
                    else:#else的情况说明该overlap的方向与估计的方向关系不符，应该被滤除
                        inconsist_link +=1
                        # print("Inconsistent-direction overlap",overlap,D[i],D[j])
                        continue
                    row+=1
                    indices.append(j) 
                    data.append(1)
                    indices.append(i)
                    data.append(-1)                   
                    indptr.append(indptr[-1]+2)
                    response+=[int(diff*(1-indel_ratio))]
    print("Inconsistent-direction overlaps:",inconsist_link)
    design=csr_matrix((data, indices, indptr), shape=(row, size))
    response=np.matrix(response).T

    return(design, response)

## 将一对contig pair之间的link集成后，再取其代表值（中位数）作为集成观测进行回归
def setupRegressionModel_ContigMedian(CD, split_contigs_reads, read_in_contig, D, ReadLength_dict, c_edge_dict,truc=3): 
    row=0
    indptr=[0]
    indices=[]
    data=[]
    c_response=[]    
    size = len(split_contigs_reads) ## contig count
    for ci,A_contig in enumerate(CD):
        for cj,B_contig in enumerate(CD):
            s_key = str(A_contig)+"-"+str(B_contig)
            if s_key in c_edge_dict:
                diff_l = []
                for olp in c_edge_dict[s_key]:
                    contig_diff = CD[A_contig]*CalContigDiff_2D(olp,D,read_in_contig,ReadLength_dict,CD[A_contig]*CD[B_contig])
                    if contig_diff*CD[A_contig]==-1:
                        continue
                    # print("contig index", A_contig, B_contig, contig_diff)
                    diff_l.append(contig_diff)                    
                              
                if len(diff_l)<truc :#or np.std(diff_l)>2000: # link<3不采用,方差大也不采用 其实应该
                    contig_diff = -1
                    continue
                else:
                    contig_diff = np.median(diff_l) #中位数
                row+=1
                indices.append(cj) 
                data.append(1)
                indices.append(ci)
                data.append(-1)                   
                indptr.append(indptr[-1]+2)
                c_response+=[contig_diff]   
        
    c_design=csr_matrix((data, indices, indptr), shape=(row, size))
    c_response=np.matrix(c_response).T

    return c_design, c_response

def setupRegressionModel_layout(split_contigs_reads, read_in_contig, D_dict, edge_dict, ReadLength_dict): #2024.4.10未经测试
    row=0
    indptr=[0]
    indices=[]
    data=[]
    size = len(split_contigs_reads) ## contig count
    c_response=[]
    for key,overlap_l in edge_dict.items():
        A_name,B_name = map(int,key.split('-'))
        if A_name in read_in_contig and B_name in read_in_contig:
            A_contig = read_in_contig[A_name][0]
            B_contig = read_in_contig[B_name][0]
            if A_contig==B_contig:
                continue
            else:

                for overlap in overlap_l:
                    contig_diff = CalContigDiff(overlap,D_dict[A_name],D_dict[B_name],read_in_contig,ReadLength_dict)
                    if contig_diff==-1:
                        continue
                    # if abs(contig_diff)>100000:
                    #     print("contig index", A_contig, B_contig, A_name, B_name, contig_diff)
                    row+=1
                    indices.append(B_contig) 
                    data.append(1)
                    indices.append(A_contig)
                    data.append(-1)                   
                    indptr.append(indptr[-1]+2)
                    c_response+=[contig_diff]   

    c_design=csr_matrix((data, indices, indptr), shape=(row, size))
    c_response=np.matrix(c_response).T
    return c_design, c_response

## 将一对contig pair之间的link集成后，再取其代表值（中位数）作为集成观测进行回归
def setupRegressionModel_layoutMedian(split_contigs_reads, read_in_contig, D_dict, edge_dict, ReadLength_dict): 
    layout_edge_dict = {}    
    for key,overlap_l in edge_dict.items():
        A_name,B_name = map(int,key.split('-'))
        if A_name in read_in_contig and B_name in read_in_contig:
            A_contig = read_in_contig[A_name][0]
            B_contig = read_in_contig[B_name][0]
            if A_contig==B_contig:
                continue
            else:
                for overlap in overlap_l:
                    contig_diff = CalContigDiff(overlap,D_dict[A_name],D_dict[B_name],read_in_contig,ReadLength_dict)
                    if contig_diff==-1:
                        continue
                    # print("contig index", A_contig, B_contig, contig_diff)
                    e_key = str(A_contig)+'-'+str(B_contig)
                    if e_key in layout_edge_dict:
                        layout_edge_dict[e_key].append(contig_diff)
                    else:
                        layout_edge_dict[e_key] = [contig_diff]
    row=0
    indptr=[0]
    indices=[]
    data=[]
    size = len(split_contigs_reads) ## contig count
    c_response=[]                    
    for e_key, e_overlaps in layout_edge_dict.items():                
        A_contig, B_contig = map(int,e_key.split('-'))
        contig_diff = np.median(e_overlaps) #中位数
        row+=1
        indices.append(B_contig) 
        data.append(1)
        indices.append(A_contig)
        data.append(-1)                   
        indptr.append(indptr[-1]+2)
        c_response+=[contig_diff]   
        
    c_design=csr_matrix((data, indices, indptr), shape=(row, size))
    c_response=np.matrix(c_response).T

    return c_design, c_response

def deleteRowsCsr(mat, indices):
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]


def IRLS_Huber(X,Y,reads,thr1=100,thr2=200,iter=100):
    if Y.shape[0]==0:
        return([], [], [])

    residual_list = []
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte
    residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    residual_list.append(residual)
    max_residual=max(residual)
    print("Initial max_residual:", round(max_residual,5),"Regression Size:",Y.shape[0])
    # print("Quantile 50 55 60 ... 95:",[np.percentile(residual,5*x+50) for x in range(10)],np.percentile(residual,99))
    threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] #threshold of Huber's weight function
    old_estimate=estimate
    n=0
    diff = 100
    res_99 = thr2+10
    
    while n<iter and res_99>thr2:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        print("Quantile 50 55 60 ... 95,99:",[np.percentile(residual,50+5*x) for x in range(1,9)],np.percentile(residual,99))
    
        # thr1 = np.percentile(residual,perc) #set thr1 to 90/95 percentile of all residuals
        threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] 
        reweight[index]=threshold[index]/abs(residual[index]) # update weights for alignments with residuals greater than thr1
        reweight=diags(reweight)
        t=X.T
        A=t.dot(reweight).dot(X)
        y=csr_matrix(Y)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        residual_list.append(residual)
        max_residual=max(residual)
        res_99 = np.percentile(residual,99)
        #print("max_residual:", round(max_residual,5))
        diff=np.percentile(abs(estimate-old_estimate),95)#max(abs(estimate-old_estimate))
        if diff<10: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1
    print("IRLS Stop at the %d rounds with max diff %d, max residual %f"%(n,diff,round(max_residual,5)))
    
    # outlier_overlap=[]
    # for i in range(len(residual)):        
    #     n1=X.indices[2*i]
    #     n2=X.indices[2*i+1]
    #     if residual[i]>thr2:
    #         outlier_overlap.append((reads[n1],reads[n2],residual[i]))
    
    index=np.where(residual>thr2)[0]
    X=deleteRowsCsr(X,index)
    Y=np.delete(Y,index,0)
    
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05)
    # residual_2=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0] ## residuals after step 2
    #max_residual=max(residual)

    #remove alignments with residuals greater than thr2
    G = nx.Graph()
    for i in range(X.shape[0]):
        n1=X.indices[2*i] #n1, n2是列序号
        n2=X.indices[2*i+1]
        G.add_edge(n1,n2)

#divide reads into separate connected components
    reads_list=[]
    estimates_list=[]
    index_list = []
    for c in nx.connected_components(G):
        # sub_index=list(G.subgraph(c).nodes)
        # sub_estimates=[estimate[i] for i in sub_index]
        # sub_reads=[reads[i] for i in sub_index]
        # estimates_list.append(list(map(int,np.round(sub_estimates))))
        # reads_list.append(sub_reads)
        sub_index=list(G.subgraph(c).nodes)
        sub_estimates=[estimate[i] for i in sub_index]
        sub_reads=[reads[i] for i in sub_index]
        sortIndex=np.argsort(sub_estimates)
        sub_index=[sub_index[k] for k in sortIndex]
        sub_reads=[sub_reads[k] for k in sortIndex]
        sub_estimates=[sub_estimates[k] for k in sortIndex]
        estimates_list.append(list(map(int,np.round(sub_estimates))))
        reads_list.append(sub_reads)
        index_list.append(sub_index)

    # return(estimates_list, reads_list, index_list, outlier_overlap, residual_list, residual_2)
    return(estimates_list, reads_list, index_list)


def WIRLS_Huber(X,Y,W,reads,thr1=100,thr2=200,iter=100):
    if Y.shape[0]==0:
        return([], [], [])

    residual_list = []
    WX = diags(W).dot(X)
    WY = diags(W).dot(Y)
    t = WX.T
    A = t.dot(WX)
    y = csr_matrix(WY)
    b = t.dot(y).todense()
    W_mean = np.mean(W)
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte
    residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    residual_list.append(residual)
    max_residual=max(residual)
    print("Initial max_residual:", round(max_residual,5),"Regression Size:",Y.shape[0])
    # print("Quantile 50 55 60 ... 95:",[np.percentile(residual,5*x+50) for x in range(10)],np.percentile(residual,99))
    threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] #threshold of Huber's weight function
    old_estimate=estimate
    n=0
    diff = 100
    res_99 = thr2+10
    
    while n<iter and res_99>thr2:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        print("Quantile 50 55 60 ... 95,99:",[np.percentile(residual,50+5*x) for x in range(1,9)],np.percentile(residual,99))
    
        # thr1 = np.percentile(residual,perc) #set thr1 to 90/95 percentile of all residuals
        threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] 
        reweight[index]=threshold[index]/abs(residual[index]) # update weights for alignments with residuals greater than thr1
        reweight=diags(reweight)
        A=t.dot(reweight).dot(WX)
        # y=csr_matrix(WY)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((WX.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        residual_list.append(residual)
        max_residual=max(residual)
        res_99 = np.percentile(residual,99)/W_mean
        #print("max_residual:", round(max_residual,5))
        diff=np.percentile(abs(estimate-old_estimate),95)#max(abs(estimate-old_estimate))
        if diff<10: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1
    print("IRLS Stop at the %d rounds with max diff %d, max residual %f"%(n,diff,round(max_residual,5)))
    
    # outlier_overlap=[]
    # for i in range(len(residual)):        
    #     n1=X.indices[2*i]
    #     n2=X.indices[2*i+1]
    #     if residual[i]>thr2:
    #         outlier_overlap.append((reads[n1],reads[n2],residual[i]))
    
    index=np.where(residual>thr2)[0]
    X=deleteRowsCsr(X,index)
    Y=np.delete(Y,index,0)
    
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05)
    # residual_2=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0] ## residuals after step 2
    #max_residual=max(residual)

    #remove alignments with residuals greater than thr2
    G = nx.Graph()
    for i in range(X.shape[0]):
        n1=X.indices[2*i] #n1, n2是列序号
        n2=X.indices[2*i+1]
        G.add_edge(n1,n2)

#divide reads into separate connected components
    reads_list=[]
    estimates_list=[]
    index_list = []
    for c in nx.connected_components(G):
        # sub_index=list(G.subgraph(c).nodes)
        # sub_estimates=[estimate[i] for i in sub_index]
        # sub_reads=[reads[i] for i in sub_index]
        # estimates_list.append(list(map(int,np.round(sub_estimates))))
        # reads_list.append(sub_reads)
        sub_index=list(G.subgraph(c).nodes)
        sub_reads=[reads[i] for i in sub_index]
        sub_estimates=[estimate[i] for i in sub_index]
        sortIndex=np.argsort(sub_estimates)
        sub_index=[sub_index[k] for k in sortIndex]
        sub_reads=[sub_reads[k] for k in sortIndex]
        sub_estimates=[sub_estimates[k] for k in sortIndex]
        estimates_list.append(list(map(int,np.round(sub_estimates))))
        reads_list.append(sub_reads)
        index_list.append(sub_index)

    # return(estimates_list, reads_list, index_list, outlier_overlap, residual_list, residual_2)
    return(estimates_list, reads_list, index_list)


def SortLayout(estimate,reads, reads_index, D_dict,ReadLength_dict):
    # 根据估出来的方向，将反向read调整坐标后，按起点排序
    reversed_estimate = [estimate[i] if D_dict[r]==1 else estimate[i]-ReadLength_dict[r] for i,r in enumerate(reads)]
    offset = min(reversed_estimate)
    lay_est = [x-offset for x in reversed_estimate]
    ## 按起点坐标重新排序
    resort_read_id = np.argsort(np.array(lay_est))
    resort_ests = sorted(lay_est)
    resort_reads = np.array(reads)[resort_read_id]
    resort_reads_index = np.array(reads_index)[resort_read_id]
    return resort_ests,resort_reads,resort_reads_index



def SortLayout_sim(estimate,reads, D_dict,ReadLength_dict): #不返回index版本
    # 根据估出来的方向，将反向read调整坐标后，按起点排序
    reversed_estimate = [estimate[i] if D_dict[r]==1 else estimate[i]-ReadLength_dict[r] for i,r in enumerate(reads)]
    offset = min(reversed_estimate)
    lay_est = [x-offset for x in reversed_estimate]
    ## 按起点坐标重新排序
    resort_read_id = np.argsort(np.array(lay_est))
    resort_ests = sorted(lay_est)
    resort_reads = np.array(reads)[resort_read_id]
    return resort_ests,resort_reads

## 注意Merge之后的reads，反向的estimate在Merge时已经换成起点了，只需要排序
def SortLayout_aftermerge(estimate,reads): 
    lay_est = [x-min(estimate) for x in estimate]
    ## 按起点坐标重新排序
    resort_read_id = np.argsort(np.array(lay_est))
    resort_ests = sorted(lay_est)
    resort_reads = np.array(reads)[resort_read_id]
    return resort_ests,resort_reads

# def CalReadDiff(overlap, D_i, D_j):
#     diff = 0
#     if overlap['B_Orientation']==1 and D_i!=D_j: #异向overlap
#         print("D_i != D_j and overlap is reverse.")
#         if D_i==1: #reads[i] --> reads[j] <--
#             diff = overlap['A_end']+overlap['B_length']-overlap['B_end']
#         elif D_i==-1: #reads[i] <-- reads[j] --> 
#             diff = overlap['A_end']+overlap['B_length']-overlap['B_end']#overlap['B_length']-overlap['B_end']-overlap['A_length']+overlap['A_end']
#             diff = -diff
#         else:
#             print("Error!")
#     elif overlap['B_Orientation']==0 and D_i==D_j: #同向overlap                        
#         if D_i==1: #reads[i] --> reads[j] -->    
#             diff = overlap['A_start']-overlap['B_start']
#         elif D_i==-1:#reads[i] <-- reads[j] <--
#             diff = overlap['A_start']-overlap['B_start']#overlap['A_end']-overlap['A_length']+overlap['B_end']
#             diff = -diff
#         else:
#             print("Error!")
#     else:
#         print('Two layouts are from different strains')
#     return diff


# 根据overlap的信息计算两条read之间的位移差
def CalReadDiff(overlap, D_i, D_j):
    diff = -1
    if overlap['B_Orientation']==1 and D_i!=D_j: #异向overlap
        if D_i==1: #reads[i] --> reads[j] <--
            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']
        elif D_i==-1: #reads[i] <-- reads[j] --> 
            diff = overlap['A_end']+overlap['B_length']-overlap['B_end']#overlap['B_length']-overlap['B_end']-overlap['A_length']+overlap['A_end']
            diff = -diff
        else:
            print("Error!")
    elif overlap['B_Orientation']==0 and D_i==D_j: #同向overlap
        if D_i==1: #reads[i] --> reads[j] -->    
            diff = overlap['A_start']-overlap['B_start']
        elif D_i==-1:#reads[i] <-- reads[j] <--
            diff = overlap['A_start']-overlap['B_start']#overlap['A_end']-overlap['A_length']+overlap['B_end']
            diff = -diff
        else:
            print("Error!")
        # print(A,B,D_i,D_j,diff)
    # else:#else的情况说明该overlap的方向与估计的方向关系不符，应该被滤除
    #     print("Inconsistent-direction overlap",overlap,D_i,D_j)
    return diff

#read_in_contig is dict, read_name:[contig_id,r_start_estimate,r_index] 反向后的起点
#若overlap指示的两条layout方向相反，说明与估计方向不一致，弃用
def CalContigDiff(overlap,D_i,D_j,read_in_contig,ReadLength_dict):
    diff = -1
    r1 = overlap['A_name']
    r2 = overlap['B_name']
    pos_1 = read_in_contig[r1][1]
    pos_2 = read_in_contig[r2][1]
    read_diff = CalReadDiff(overlap,D_i,D_j)
    if read_diff==-1:
        return diff
    if D_i==1 and D_j==1 and overlap['B_Orientation']==0:
        diff = read_diff+pos_1-pos_2
    elif D_i==1 and D_j==-1 and overlap['B_Orientation']==1:
        diff = read_diff+pos_1-pos_2-ReadLength_dict[r2]
    elif D_i==-1 and D_j==1 and overlap['B_Orientation']==1:
        diff = read_diff+pos_1+ReadLength_dict[r1]-pos_2 #-read_diff+pos_1+ReadLength_dict[r1]-pos_2
    elif D_i==-1 and D_j==-1 and overlap['B_Orientation']==0:
        diff = read_diff+pos_1+ReadLength_dict[r1]-pos_2-ReadLength_dict[r2] #-read_diff+pos_1+ReadLength_dict[r1]-pos_2-ReadLength_dict[r2]
    return diff

#read_in_contig is dict, read_name:[contig_id,r_start_estimate,r_index_in_contig, r_backwards_index]
def CalContigDiff_2D(overlap,D,read_in_contig,ReadLength_dict,CD_ij):
    diff = -1    
    r1 = int(overlap['A_name'])
    r2 = int(overlap['B_name'])
    D_i = D[r1]
    D_j = D[r2]
    pos_1 = read_in_contig[r1][1]
    pos_2 = read_in_contig[r2][1]
    read_diff = CalReadDiff(overlap,D_i,D_j*CD_ij) ##注意D_i 是read_i 在contig_i中的方向，所以计算read_diff要放进统一坐标系，得乘CD_ij
    if read_diff==-1:
        return diff
    # if CD_ij==1: ## 按RegScaf中的公式
    #     if D_i==1 and D_j==1 and overlap['B_Orientation']==0:
    #         diff = read_diff+pos_1-pos_2
    #     elif D_i==1 and D_j==-1 and overlap['B_Orientation']==1:
    #         diff = read_diff+pos_1+pos_2+ReadLength_dict[r2]
    #     elif D_i==-1 and D_j==1 and overlap['B_Orientation']==1:
    #         diff = -read_diff+pos_1+ReadLength_dict[r1]+pos_2 #-read_diff+pos_1+ReadLength_dict[r1]-pos_2
    #     elif D_i==-1 and D_j==-1 and overlap['B_Orientation']==0:
    #         diff = -read_diff+pos_1+ReadLength_dict[r1]-pos_2-ReadLength_dict[r2] #-read_diff+pos_1+ReadLength_dict[r1]-pos_2-ReadLength_dict[r2]
    # else:
    #     if D_i==1 and D_j==1 and overlap['B_Orientation']==1:
    #         diff = read_diff+pos_1+pos_2
    #     elif D_i==1 and D_j==-1 and overlap['B_Orientation']==0:
    #         diff = read_diff+pos_1-pos_2-ReadLength_dict[r2]
    #     elif D_i==-1 and D_j==1 and overlap['B_Orientation']==0:
    #         diff = -read_diff+pos_1+ReadLength_dict[r1]-pos_2 #-read_diff+pos_1+ReadLength_dict[r1]-pos_2
    #     elif D_i==-1 and D_j==-1 and overlap['B_Orientation']==1:
    #         diff = -read_diff+pos_1+ReadLength_dict[r1]+pos_2+ReadLength_dict[r2]
    if CD_ij==1:
        if D_i==1 and D_j==1 and overlap['B_Orientation']==0:
            diff = read_diff+pos_1-pos_2
        elif D_i==1 and D_j==-1 and overlap['B_Orientation']==1:
            diff = read_diff+pos_1-pos_2-ReadLength_dict[r2]
        elif D_i==-1 and D_j==1 and overlap['B_Orientation']==1:
            diff = read_diff+pos_1+ReadLength_dict[r1]-pos_2 #-read_diff+pos_1+ReadLength_dict[r1]-pos_2
        elif D_i==-1 and D_j==-1 and overlap['B_Orientation']==0:
            diff = read_diff+pos_1+ReadLength_dict[r1]-pos_2-ReadLength_dict[r2] #-read_diff+pos_1+ReadLength_dict[r1]-pos_2-ReadLength_dict[r2]
    else:
        if D_i==1 and D_j==1 and overlap['B_Orientation']==1:
            diff = read_diff+pos_1+pos_2
        elif D_i==1 and D_j==-1 and overlap['B_Orientation']==0:
            diff = read_diff+pos_1+pos_2+ReadLength_dict[r2]
        elif D_i==-1 and D_j==1 and overlap['B_Orientation']==0:
            diff = -read_diff+pos_1+ReadLength_dict[r1]+pos_2 #-read_diff+pos_1+ReadLength_dict[r1]-pos_2
        elif D_i==-1 and D_j==-1 and overlap['B_Orientation']==1:
            diff = -read_diff+pos_1+ReadLength_dict[r1]+pos_2+ReadLength_dict[r2]
    return diff

def CalStartDiff(overlap,D_i,D_j): #若D_i
    diff = CalReadDiff(overlap,D_i,D_j)
    if D_i==-1:
        diff += overlap['A_length']
    if D_j==-1:
        diff += -overlap['B_length']
    return diff

# 修改后的split函数，之前应该没弄明白i,next_i
# 下面尝试SplitRegression
# def SpiltRegressionLayout(sorted_Reads,Starts,ReadLength_dict,edge_dict,D_dict,direct_dict,sd=600):
#     remained_index = [x for x in range(len(sorted_Reads))]
#     if len(sorted_Reads)<2:
#         return [remained_index]
#     Ends = [Starts[i]+ReadLength_dict[r_i] for i,r_i in enumerate(sorted_Reads)]
#     Split_Layouts = []
#     k = 0
#     # print('Now:',len(remained_index))
#     while len(remained_index)>0:
#         # i,ni是指针序号,c_i,next_i是read在sorted_reads中的原始序号
#         i = 0 #remained_index[0]
#         c_i = remained_index[i]
#         layout = [c_i]
#         if len(remained_index)>1:
#             ni = i+1
#             next_i = remained_index[ni]
#             while ni<len(remained_index)-1 and Starts[next_i]<Ends[c_i]+sd:
#                 current_r = sorted_Reads[c_i]
#                 next_r = sorted_Reads[next_i]
#                 est_diff = Starts[next_i]-Starts[c_i]
#                 est_direct = 1 if D_dict[current_r]==D_dict[next_r] else -1
#                 link_diff = -1
#                 link_direct = 0
#                 key = str(current_r)+'-'+str(next_r)
#                 key_1 = str(next_r)+'-'+str(current_r)                      
#                 d_key = str(min(current_r,next_r))+"-"+str(max(current_r,next_r))
#                 if d_key in direct_dict:
#                     link_direct = np.sign(direct_dict[d_key][0])
#                 if key in edge_dict:
#                     link_diff = CalStartDiff(edge_dict[key][0],D_dict[current_r],D_dict[next_r])                      
#                 elif key_1 in edge_dict:
#                     link_diff = -CalStartDiff(edge_dict[key_1][0],D_dict[next_r],D_dict[current_r])
#                 # print(c_i,next_i,key,D_dict[current_r],D_dict[next_r],link_diff,est_diff)
#                 if (key in edge_dict or key_1 in edge_dict) and abs(link_diff-est_diff)<sd and est_direct==link_direct:
#                     layout.append(next_i)
#                     i = ni
#                     c_i = next_i
#                     ni +=1
#                     next_i = remained_index[ni]
#                 elif ni<len(remained_index)-1:
#                     ni +=1
#                     next_i = remained_index[ni]
#                 else: # ni是最末尾的read，这一条layout应该到这结束
#                     break
#         for x in layout:   
#             remained_index.remove(x)
#         # print('Now:',len(remained_index))
#         if len(layout)>1: # drop layout with only one read (Right?)
#             Split_Layouts.append(layout)
#             k+=1
#             # print("Join the %d-th layout:"%k,len(layout))
#     # print('Split into %d layouts'%(len(Split_Layouts)),[len(l) for l in Split_Layouts])
#     return Split_Layouts

def IfValidlink(r1,r2,D_dict,edge_dict,direct_dict):
    link_diff = -1
    link_direct = 0
    key = str(r1)+'-'+str(r2)
    key_1 = str(r2)+'-'+str(r1)                      
    d_key = str(min(r1,r2))+"-"+str(max(r1,r2))
    if d_key in direct_dict:
        link_direct = np.sign(direct_dict[d_key][0])
    if key in edge_dict:
        link_diff = CalStartDiff(edge_dict[key][0],D_dict[r1],D_dict[r2])                      
    elif key_1 in edge_dict:
        link_diff = -CalStartDiff(edge_dict[key_1][0],D_dict[r2],D_dict[r1])
    return link_direct,link_diff
    

# SplitRegression 加入next_r时往前多找r几条read进行判断
def SpiltRegressionLayout(sorted_Reads,Starts,ReadLength_dict,edge_dict,D_dict,direct_dict,sd=600,back_w=4):
    remained_index = [x for x in range(len(sorted_Reads))]
    if len(sorted_Reads)<2:
        return [remained_index]
    Ends = [Starts[i]+ReadLength_dict[r_i] for i,r_i in enumerate(sorted_Reads)]
    Split_Layouts = []
    k = 0
    # print('Now:',len(remained_index))
    while len(remained_index)>0:
        # i,ni是指针序号,c_i,next_i是read在sorted_reads中的原始序号
        i = 0 #remained_index[0]
        c_i = remained_index[i]
        layout = [c_i]
        if len(remained_index)>1:
            ni = i+1
            next_i = remained_index[ni]
            while ni<len(remained_index)-1 and Starts[next_i]<Ends[c_i]+sd:
                current_r = sorted_Reads[c_i]
                next_r = sorted_Reads[next_i]
                est_diff = Starts[next_i]-Starts[c_i]
                est_direct = 1 if D_dict[current_r]==D_dict[next_r] else -1
                link_direct,link_diff = IfValidlink(current_r,next_r,D_dict,edge_dict,direct_dict)    
                wi=1
                while link_direct==0 and wi<=min(back_w,len(layout)-1):
                    wi+=1
                    last_r = sorted_Reads[layout[-wi]]
                    est_diff = Starts[next_i]-Starts[layout[-wi]]                   
                    est_direct = 1 if D_dict[last_r]==D_dict[next_r] else -1
                    link_direct,link_diff = IfValidlink(last_r,next_r,D_dict,edge_dict,direct_dict)    
                # print(c_i,next_i,key,D_dict[current_r],D_dict[next_r],link_diff,est_diff)
                if abs(link_diff-est_diff)<sd and est_direct==link_direct:
                    layout.append(next_i)
                    i = ni
                    c_i = next_i
                    ni +=1
                    next_i = remained_index[ni]
                elif ni<len(remained_index)-1:
                    ni +=1
                    next_i = remained_index[ni]
                else: # ni是最末尾的read，这一条layout应该到这结束
                    break
        for x in layout:   
            remained_index.remove(x)
        # print('Now:',len(remained_index))
        # if len(layout)>1: # drop layout with only one read (Right?)
        Split_Layouts.append(layout)
        k+=1
        # print("Join the %d-th layout:"%k,len(layout))
    # print('Split into %d layouts'%(len(Split_Layouts)),[len(l) for l in Split_Layouts])
    return Split_Layouts    

def ContigIRLS_Huber(X,Y,contigs,thr1=300,thr2=1000,iter=100):
    if Y.shape[0]==0:
        return([], [], [])

    residual_list = []
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute ordinary least squares estiamte
    residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
    residual_list.append(residual)
    max_residual=max(residual)
    res_perc95 = np.percentile(residual,95)
    print("Initial max_residual:", round(max_residual,5),"Regression Size:",Y.shape[0])
    print("Quantile 50 55 60 ... 95:",[np.percentile(residual,5*x+50) for x in range(10)],np.percentile(residual,99))
    threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] #threshold of Huber's weight function
    old_estimate=estimate
    n=0
    diff = 100
    
    while n<iter and res_perc95>thr2:
        index=np.where(residual>threshold)[0]
        reweight=np.ones(len(residual))
        print("Quantile 50 55 60 ... 95,99:",[np.percentile(residual,50+5*x) for x in range(1,9)],np.percentile(residual,99))
        res_perc95 = np.percentile(residual,95)
        # thr1 = np.percentile(residual,perc) #set thr1 to 90/95 percentile of all residuals
        threshold=csr_matrix(np.ones(len(residual))*thr1).toarray()[0] 
        reweight[index]=threshold[index]/abs(residual[index]) # update weights for alignments with residuals greater than thr1
        reweight=diags(reweight)
        t=X.T
        A=t.dot(reweight).dot(X)
        y=csr_matrix(Y)
        b=t.dot(reweight).dot(y).todense()
        estimate, exitCode = linalg.lgmres(A, b, atol=1e-05) #compute weighted least squares estimate
        residual=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0]
        residual_list.append(residual)
        max_residual=max(residual)
        #print("max_residual:", round(max_residual,5))
        diff=np.percentile(abs(estimate-old_estimate),95)#max(abs(estimate-old_estimate))
        if diff<10: # convergence condition of estimates
            break
        else:
            old_estimate=estimate
            n+=1
    print("IRLS Stop at the %d rounds with max diff %d, max residual %f"%(n,diff,round(max_residual,5)))
        
    index=np.where(residual>thr2)[0]
    X=deleteRowsCsr(X,index)
    Y=np.delete(Y,index,0)
    
    t=X.T
    A=t.dot(X)
    y=csr_matrix(Y)
    b=t.dot(y).todense()
    estimate, exitCode = linalg.lgmres(A, b, atol=1e-05)
    # residual_2=abs((X.dot(csr_matrix(estimate).T)-y).todense()).T.getA()[0] ## residuals after step 2
    #max_residual=max(residual)

    #remove alignments with residuals greater than thr2
    G = nx.Graph()
    G.add_nodes_from(range(len(contigs))) ##in case an island node is dropped
    for i in range(X.shape[0]):
        n1=X.indices[2*i]
        n2=X.indices[2*i+1]
        G.add_edge(n1,n2)

    #divide contigs into separate connected components
    reads_list=[]
    estimates_list=[]
    index_list = []
    for c in nx.connected_components(G):
        sub_index=list(G.subgraph(c).nodes)
        sub_contigs=[contigs[i] for i in sub_index]
        sub_estimates=[estimate[i] for i in sub_index]
        sortIndex=np.argsort(sub_estimates)
        sub_index=[sub_index[k] for k in sortIndex]
        sub_contigs=[sub_contigs[k] for k in sortIndex]
        sub_estimates=[sub_estimates[k] for k in sortIndex]
        estimates_list.append(list(map(int,np.round(sub_estimates))))
        reads_list.append(sub_contigs)
        index_list.append(sub_index)

    # return(estimates_list, reads_list, index_list, outlier_overlap, residual_list, residual_2)
    return(estimates_list, reads_list, index_list)

## 将再次回归得到的contigs合并
def MergeContigs(c_estimates, c_contigs, split_contigs_reads,split_contigs_estimates):
    c_0 = c_contigs[0]
    m_estimates = split_contigs_estimates[c_0]
    m_reads = split_contigs_reads[c_0]
    
    for i in range(1,len(c_contigs)):
        l_diff = c_estimates[i]-c_estimates[0]
        m_reads = np.concatenate((m_reads,split_contigs_reads[c_contigs[i]]), axis = 0) #追加到列表末尾
        c_ests = split_contigs_estimates[c_contigs[i]]
        m_estimates = np.concatenate((m_estimates, np.array([x+l_diff-c_ests[0] for x in c_ests])), axis = 0) #每条contig 的reads坐标应该从0开始，每条read平移l_diff
        # m_reads += split_contigs_reads[c_contigs[i]]
        # m_estimates += [x+l_diff for x in split_contigs_estimates[c_contigs[i]]]
        
    return m_reads, m_estimates

## 将再次回归得到的contigs合并,考虑方向
def MergeContigs_D(c_estimates,c_contigs, CD, split_contigs_reads,split_contigs_estimates,D,D_update, ReadLength_dict):
    # c_contigs = list(CD.keys())
    c_0 = c_contigs[0]    
    if CD[c_contigs[0]]==1:
        m_estimates = split_contigs_estimates[c_0]
        m_reads = split_contigs_reads[c_0]
    else: #负向contig
        c_ests = split_contigs_estimates[c_0]
        m_reads = split_contigs_reads[c_0]
        s_contigs_ends = [c_ests[i]+ReadLength_dict[r_i] for i,r_i in enumerate(m_reads)]
        m_estimates = (-1)*np.array(s_contigs_ends)
        m_estimates = [x-m_estimates[0]-ReadLength_dict[m_reads[0]] for x in m_estimates]
        for sc in m_reads:
            D_update[sc]*=-1
    if len(c_estimates)<2:
        return m_reads, m_estimates, D_update
    
    for i in range(1,len(c_contigs)):
        c_i = c_contigs[i]
        l_diff = c_estimates[i]-c_estimates[0]
        c_ests = split_contigs_estimates[c_i]
        c_reads = split_contigs_reads[c_i]
        ini_start = c_ests[0]        
        if CD[c_contigs[i]]==-1: ##负向contig
            s_contigs_reads = split_contigs_reads[c_i]
            s_contig_starts = split_contigs_estimates[c_i]
            s_contigs_ends = [s_contig_starts[r]+ReadLength_dict[r_i] for r,r_i in enumerate(s_contigs_reads)] ##整条contig反向，把起点换成终点
            c_ests = (-1)*np.array(s_contigs_ends)
            c_ests = [x-min(c_ests) for x in c_ests]
            ini_start = c_ests[0]+ReadLength_dict[c_reads[0]]
            for sc in c_reads:
                D_update[sc]*=-1
        #与基准（即当前分支的第一条contig）方向相同，直接平移追加
        m_reads = np.concatenate((m_reads,c_reads), axis = 0) #追加到列表末尾        
        m_estimates = np.concatenate((m_estimates, np.array([x+l_diff-ini_start for x in c_ests])), axis = 0) #每条contig 的reads坐标应该从0开始，每条read平移l_diff            
        # print(len(m_reads),len(m_estimates))
    return m_reads, m_estimates, D_update

## 对每个子图基于read diff回归一次，并sort,split,merge
def regAssem(V,D,edge_dict,ReadLength_dict,direct_dict):
    D_list = [D[v] for v in V]
    design, response = setupRegressionModel_direct(V,D_list,edge_dict)
    estimates_list, reads_list, index_list = IRLS_Huber(design, response, V, thr1=200, thr2=600)

    sorted_estimates_list = []
    sorted_reads_list = []
    sorted_reads_index = []
    for i in range(len(estimates_list)):
        resort_ests,resort_reads,resort_reads_index = SortLayout(estimates_list[i],reads_list[i],index_list[i],D,ReadLength_dict)
        sorted_estimates_list.append(resort_ests)
        sorted_reads_list.append(resort_reads)
        sorted_reads_index.append(resort_reads_index)

    read_in_contig = dict()
    split_contigs_reads = []
    split_contigs_estimates = []
    contig_num = 0
    for i in range(len(sorted_estimates_list)):
        split_contigs = SpiltRegressionLayout(sorted_reads_list[i],sorted_estimates_list[i],ReadLength_dict,edge_dict,D,direct_dict,sd=1000)
        for sl in split_contigs:
            contig_reads = np.array(sorted_reads_list[i])[sl]
            contig_estimates = np.array(sorted_estimates_list[i])[sl]
            split_contigs_reads.append(contig_reads)
            split_contigs_estimates.append(contig_estimates)
            print("Split Layout length:", max(contig_estimates)-min(contig_estimates)+ReadLength_dict[contig_reads[-1]])
            for r_i,r_name in enumerate(contig_reads):
                read_in_contig[r_name]=[contig_num,contig_estimates[r_i]]
            contig_num +=1
            
    # 对split得到的layout二次回归,但是目前是没有重新定向,所以循环放在V之内
    c_design, c_response = setupRegressionModel_layout(split_contigs_reads, read_in_contig, D, edge_dict, ReadLength_dict)
    sub_contigs = np.array(range(contig_num))
    c_estimates_list, c_reads_list, c_index_list = ContigIRLS_Huber(c_design, c_response, sub_contigs, thr1=200, thr2=600)

    Merged_reads_list = []
    Merged_estimates_list = []
    for i, c_contigs in enumerate(c_reads_list):
        c_estimates = c_estimates_list[i]
        m_reads, m_estimates = MergeContigs(c_estimates, c_contigs, split_contigs_reads, split_contigs_estimates)
        Merged_reads_list.append(m_reads)
        Merged_estimates_list.append(m_estimates)
        print("Merged Layout length:", max(m_estimates)-min(m_estimates)+ReadLength_dict[m_reads[-1]])
    return(Merged_reads_list,Merged_estimates_list)

## iterate the regression at the contig level
## 从SortLayout, SplitRegression 开始
def regAssem_contigIter(reads_list_l, estimates_l,D,ReadLength_dict,edge_dict,direct_dict):
    sorted_estimates_list = []
    sorted_reads_list = []
    for i in range(len(estimates_l)):
        resort_ests,resort_reads = SortLayout_sim(estimates_l[i],reads_list_l[i],D,ReadLength_dict)
        sorted_estimates_list.append(resort_ests)
        sorted_reads_list.append(resort_reads)
        
    read_in_contig = dict()
    split_contigs_reads = []
    split_contigs_estimates = []
    contig_num = 0
    for i in range(len(sorted_estimates_list)):
        split_contigs = SpiltRegressionLayout(sorted_reads_list[i],sorted_estimates_list[i],ReadLength_dict,edge_dict,D,direct_dict,sd=1000)
        for sl in split_contigs:
            contig_reads = np.array(sorted_reads_list[i])[sl]
            contig_estimates = np.array(sorted_estimates_list[i])[sl]
            split_contigs_reads.append(contig_reads)
            split_contigs_estimates.append(contig_estimates)
            for r_i,r_name in enumerate(contig_reads):
                read_in_contig[r_name]=[contig_num,contig_estimates[r_i]]
            contig_num +=1
            
    c_design, c_response = setupRegressionModel_layout(split_contigs_reads, read_in_contig, D, edge_dict, ReadLength_dict)
    sub_contigs = np.array(range(contig_num))
    c_estimates_list, c_reads_list, c_index_list = ContigIRLS_Huber(c_design, c_response, sub_contigs, thr1=200, thr2=1000)

    Merged_reads_list = []
    Merged_estimates_list = []
    for i, c_contigs in enumerate(c_reads_list):
        c_estimates = c_estimates_list[i]
        m_reads, m_estimates = MergeContigs(c_estimates, c_contigs, split_contigs_reads, split_contigs_estimates)
        Merged_reads_list.append(m_reads)
        Merged_estimates_list.append(m_estimates)
        print("Layout length:", max(m_estimates)-min(m_estimates)+ReadLength_dict[m_reads[-1]])
        
    return(Merged_reads_list,Merged_estimates_list)

def WriteLayoutResult(filename,ctg_id, n, resort_reads, resort_ests,resort_est_ends,lay_D,ReadSeq, block=2048, block_overlap=1000):
    file_in = open(filename,'a+')    
    block_start = 0
    block_end = block
    lay_length = max(resort_est_ends)-min(resort_ests)
    print('>ctg%d\tnodes=%d\tlen=%d'%(ctg_id,int(lay_length/(block-block_overlap))+1, lay_length),file=file_in)

    while block_start<lay_length : #and n<10
        print('E\t%d\tN%d\t+\tN%d\t+'%(block_start, n, n+1),file=file_in)
        for r_i in range(len(resort_ests)):
            if resort_est_ends[r_i]<=block_start:
                continue
            if resort_ests[r_i]>=block_end:
                break
            read_name = resort_reads[r_i]
            read_strand = lay_D[read_name]
            read_block_start = read_block_length = 0
            if resort_ests[r_i]<=block_start and resort_est_ends[r_i]>=block_end:
                read_block_start = block_start-resort_ests[r_i]
                read_block_length = block
            if resort_ests[r_i]<=block_start and resort_est_ends[r_i]<=block_end:
                read_block_start = block_start-resort_ests[r_i]
                read_block_length = resort_est_ends[r_i]-block_start
            if resort_ests[r_i]>=block_start and resort_est_ends[r_i]>=block_end:
                read_block_start = 0
                read_block_length = block_end - resort_ests[r_i]
            if resort_ests[r_i]>=block_start and resort_est_ends[r_i]<=block_end:
                read_block_start = 0
                read_block_length = resort_est_ends[r_i]-resort_ests[r_i]
            if read_strand=='+':
                read_seq = ReadSeq[str(read_name)+'_F'][read_block_start:(read_block_start+read_block_length)]
            else:
                read_seq = Seq(ReadSeq[str(read_name)+'_F']).reverse_complement().upper()
                read_seq = read_seq[read_block_start:(read_block_start+read_block_length)]
            if read_block_length>1:
                # print(read_block_length)
                print('S\tSRR10971019.%s\t%s\t%d\t%d\t%s'%(read_name, read_strand, read_block_start, read_block_length, read_seq),file=file_in)
                
        block_start = block_end-block_overlap
        block_end = block_start+block
        n+=1
    return n
