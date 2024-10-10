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

import networkx as nx
from FindChimericRead import *
from PreFunctions import *
from collections import defaultdict

from Direction import *
#from Evaluation import *
import argparse
from multiprocessing import Pool, Manager
from Bio import Align
from Bio.Seq import Seq
from spoa import poa
import datetime

aligner = Align.PairwiseAligner(mode = 'local', match_score=1, mismatch_score=-2, open_gap_score=-5, extend_gap_score=-1)

descript="This script is designed for de novo assembly based on a robust regression model in RegAssembler.\n"

parser = argparse.ArgumentParser(description=descript)
parser.add_argument('-r', required=True, help='Input long reads, should be fasta file')
parser.add_argument('-m', required=True, help='Input self-mapping file, minimap2.paf')
parser.add_argument('-o', required=False, default='output.fa', help='Output contig file name [output.fa]')
parser.add_argument('-p', type=str, default='_', help='The delimiter of read name in input file [_]')
parser.add_argument('-t', type=int, default=1, help='number of threads [1]')
parser.add_argument('-HO', type=int, default=80, help="The threshold for hangingout in filter overlaps [80]")
parser.add_argument('-c', action="store_true", help="Determine if drop overlaps with a cut read, for circular genome.")
parser.add_argument('-C', action="store_true", help="Determine if drop suspicious chimeric reads.")
parser.add_argument('-i', type=int, default=3, help='Rounds to Iterating on contig level regression')
parser.add_argument('-sd', type=int, default=1000, help='The variation of position difference allowed in splitting regression layouts')
parser.add_argument('-ms', type=int, default=20000, help='The max size of a connected subgraph allowed in BFS, for computation efficiency')
# parser.add_argument('-log', type=str, default='regAssembler.log', help='Log file')
parser.add_argument('-bw', default=12000, type=int, help='The block width in block Consensus')

args = parser.parse_args()

## 类似要注意，有hanging out的情形，所以最好是从overlap区域开始更新consensus
def extendSeq(cns_seq,next_seq,truc=1000,minscore=40):    
    seq_L = cns_seq[-truc:]
    seq_R = next_seq[:truc]
    alignments = aligner.align(seq_L,seq_R) ##对long read要截取比对
    alignment = alignments[0]
    extendSeq = ''
    print('Alignment.score', alignment.score)
    if alignment.score<minscore:
        print("Abnormal!")
        return ''
    else:
        a1, _ = alignment.aligned[0][0]
        a2, _ = alignment.aligned[1][0]
        _,b1 = alignment.aligned[0][-1]
        _,b2 = alignment.aligned[1][-1]
        extendSeq = seq_L[:b1]+seq_R[b2:]+next_seq[truc:] ## overlap
    return extendSeq

def POACNS(k,block_seqs,cns_dict):
    consensus, _ = poa(block_seqs)   
    cns_dict[k] = consensus
        
## 把BlockConsensu这一步改成多线程
def BlockConsensus_MultiThread(layout_num,contig_reads,contig_estimates,ReadLength_dict,D,consensus_dict,bw=3000,block_overlap=600,t=40):
    est_ends = [x+ReadLength_dict[contig_reads[r_i]] for r_i,x in enumerate(contig_estimates)]
    contigs_l = []
    block_start = 0
    block_end = bw
    lay_length = max(est_ends)
    
    # cns_dict=dict()
    m=Manager()
    pool = Pool(t)
    cns_dict = m.dict()
    k = 0
    while block_start<lay_length :
        block_seqs = []
        for r_i,read_name in enumerate(contig_reads):
            if est_ends[r_i]<=block_start:
                continue
            if contig_estimates[r_i]>=block_end:
                break          
            read_block_start = read_block_length = 0
            if contig_estimates[r_i]<=block_start and est_ends[r_i]>=block_end:
                read_block_start = block_start-contig_estimates[r_i]
                read_block_length = bw
            elif contig_estimates[r_i]<=block_start and est_ends[r_i]<=block_end:
                read_block_start = block_start-contig_estimates[r_i]
                read_block_length = est_ends[r_i]-block_start
            elif contig_estimates[r_i]>=block_start and est_ends[r_i]>=block_end:
                read_block_start = 0
                read_block_length = block_end - contig_estimates[r_i]
            elif contig_estimates[r_i]>=block_start and est_ends[r_i]<=block_end:
                read_block_start = 0
                read_block_length = est_ends[r_i]-contig_estimates[r_i]
            
            if read_block_length>1:
                if D[read_name]==1:
                    read_seq = str(ReadSeq[str(read_name)+'_F'][read_block_start:(read_block_start+read_block_length)])
                else:
                    read_seq = Seq(ReadSeq[str(read_name)+'_F']).reverse_complement().upper()
                    read_seq = str(read_seq[read_block_start:(read_block_start+read_block_length)])
                block_seqs.append(read_seq)        
        block_start = block_end-block_overlap
        block_end = block_start+bw
        k+=1
        # POACNS(k,block_seqs,cns_dict)
        pool.apply_async(POACNS,args=(k,block_seqs,cns_dict))        
    
    pool.close()
    pool.join()
    
    contig_seq = last_seq = ReadSeq[str(contig_reads[0])+'_F'] if D[contig_reads[0]]==1 else ReadSeq[str(contig_reads[0])+'_R'][:block_overlap]
    for i_k in range(1,k+1):
        consensus = cns_dict[i_k]
        bw_truc = block_overlap+max(0,len(consensus)-bw)
        extend_seq = extendSeq(last_seq,consensus,truc=bw_truc)
        if len(extend_seq)>0:
            contig_seq =contig_seq[:-bw_truc]+extend_seq
            print("Current len", len(contig_seq))
        else: ##若不能前后比对上，只能断开，开始新的一条contig
            print("Contig stops at len %d"%(len(contig_seq)))
            contigs_l.append(contig_seq)
            contig_seq = consensus
        last_seq = consensus   
        
    print("Contig stops at len %d with %d reads"%(len(contig_seq),len(contig_reads)))
    contigs_l.append(contig_seq) ## 最后一条不能忘了
    consensus_dict[layout_num] = contigs_l     
    return contig_seq

if __name__=='__main__':
    # log_file = open(args.log, 'w')
    split_sd = args.sd
    HangingOut = args.HO

    now = datetime.datetime.now()
    print("RegAssembler starts at time", now.strftime("%Y-%m-%d %H:%M:%S"))

    print("Loading self-mapping results from %s ..."%args.m)
    edge_dict,direct_dict,G_ini,ReadLength_dict = ImportAlignmentsFromMinimap(args.m,HangingOut,filter_chimera=args.C,sep=args.p) 
    now = datetime.datetime.now()
    print("Loaded successfully at time", now.strftime("%Y-%m-%d %H:%M:%S"))
    
    
    # 拷贝一份
    G = G_ini.copy()    
    all_nodes = G.nodes()
    
    ## 定向分子图   G = G_ini.copy()   
    all_nodes = G.nodes()

    ## 定向分子图   
    V_list = []
    D = dict()
    suspious_V = dict()
    max_size = 20000
    single_node = set()
    Visited = {key: 0 for key in all_nodes}
    while len(G)>0:
        S = list(G.nodes())
        V = [S.pop()]
        if len(G[V[0]])==0: # V[0]没有邻居
            G.remove_nodes_from(V)
            single_node.add(V[0])
            continue
        V_neighbors = dict(G[V[0]])
        # print(len(V_neighbors))
        
        G.nodes[V[0]]['direct'] = 1
        for n1 in V_neighbors:
            G.nodes[n1]['direct'] = np.sign(V_neighbors[n1]['weight'])
        node = max(V_neighbors, key=lambda k:abs(V_neighbors[k]['weight']))
        # G.nodes[node]['direct'] = np.sign(V_neighbors[node]['weight'])
        Visited[V[0]]=1

        while V_neighbors and len(V)<max_size:    # 限制V的大小
            node = max(V_neighbors, key=lambda k:abs(V_neighbors[k]['weight']))
            # print('Add ', node, 'V:', len(V), 'V_neighbors:', len(V_neighbors), 'Weight', V_neighbors[node]['weight'])
            V.append(node)
            Visited[node] = 1
            del V_neighbors[node]
            nei = G[node]
            # print('Neighbors of node:', nei)
            for n2 in nei:
                if Visited[n2]==1 or n2 in suspious_V: ## visited or suspious
                    continue
                if n2 in V_neighbors: ##loop, n2 already in V_nei, see if contradict            
                    if np.sign(G.nodes[node]['direct']*nei[n2]['weight'])==np.sign(V_neighbors[n2]['weight']): ##Updata n2 by adding weight
                        V_neighbors[n2]['weight']+=nei[n2]['weight']*G.nodes[node]['direct']                    
                    else: ##n2的方向矛盾，不应该加入V，之后也不应该再考虑
                        # print('Delelte',n2)
                        del V_neighbors[n2]
                        suspious_V[n2]=1
                        # G.remove_node(n2) 
                else:
                    V_neighbors[n2]={'weight': nei[n2]['weight']*G.nodes[node]['direct']}  #add n2     
                    G.nodes[n2]['direct'] = G.nodes[node]['direct']*np.sign(nei[n2]['weight'])
                    # print("New nei",n2,nei[n2]['weight'])
        print("New V size", len(V))
        V_list.append(V)
        D.update(dict(zip(V,[G.nodes[v]['direct'] for v in V])))
        G.remove_nodes_from(V)
        # if len(V)>10000:
        #     break    
        
    now = datetime.datetime.now()
    print("Successfully obtained %d connected components at time %s"%(len(V_list), now.strftime("%Y-%m-%d %H:%M:%S")))
        
    print("Total reads:",sum([len(V) for V in V_list]))
    
    contig_num = 0
    read_in_contig = dict()
    split_contigs_reads = []
    split_contigs_estimates = []
    for V in V_list:        
        # if len(V)<10:
        #     continue
        D_list = [D[v] for v in V]
        design, response= setupRegressionModel_direct(V,D_list,edge_dict)
        estimates_list, reads_list, index_list = IRLS_Huber(design, response, V, thr1=200, thr2=800, iter=50)
            
        for i in range(len(estimates_list)):
            resort_ests,resort_reads,resort_reads_index = SortLayout(estimates_list[i],reads_list[i],index_list[i],D,ReadLength_dict)
            split_contigs = SpiltRegressionLayout(resort_reads,resort_ests,ReadLength_dict,edge_dict,D,direct_dict,sd=split_sd)
            for sl in split_contigs:
                contig_reads = np.array(resort_reads)[sl]
                contig_estimates = np.array(resort_ests)[sl]
                contig_estimates = [x-contig_estimates[0] for x in contig_estimates]
                # if len(contig_reads)<10: ##筛掉reads很少的contig
                #     continue
                split_contigs_reads.append(contig_reads)
                split_contigs_estimates.append(contig_estimates)
                print(i, "Split Layout length:", max(contig_estimates)-min(contig_estimates)+ReadLength_dict[contig_reads[-1]])
                for r_i,r_name in enumerate(contig_reads):
                    read_in_contig[r_name]=[contig_num,contig_estimates[r_i],r_i,len(contig_reads)-r_i]
                contig_num +=1
    
    print("Total reads:",sum([len(x) for x in split_contigs_reads]))
    iter = 1
    while iter<args.i:
        now = datetime.datetime.now()
        print("Start the %d round of contig level regression at time %s......."%(iter,now.strftime("%Y-%m-%d %H:%M:%S")))
        contig_direct_dict = defaultdict(list)
        c_edge_dict = defaultdict(list)
        end_w = 15
        for key,overlap_l in edge_dict.items():
            A_name,B_name = map(int,key.split('-'))
            if A_name in read_in_contig and B_name in read_in_contig:
                A_contig = read_in_contig[A_name][0]
                B_contig = read_in_contig[B_name][0]
                if A_contig==B_contig or min(read_in_contig[A_name][2:])>end_w or min(read_in_contig[B_name][2:])>end_w: #read在contig的中间位置，筛去
                    continue
                s_key = str(min(A_contig,B_contig))+"-"+str(max(A_contig,B_contig))
                c_key = str(A_contig)+"-"+str(B_contig)
                # r_key = str(min(A_name,B_name))+"-"+str(max(A_name,B_name))
                for overlap in overlap_l:
                    c_edge_dict[c_key].append(overlap)
                    contig_direct_dict[s_key].append((1 if overlap['B_Orientation']==0 else -1)*overlap['align_bases']*D[A_name]*D[B_name]) ##contig A B的方向是否一致
                    
        (C_D,J_mat) = DirectMST_contig(contig_direct_dict,range(contig_num),truc=10)

        Merged_reads_list = []
        Merged_estimates_list = []
        D_update = D.copy()
        for CD in C_D:
            contig_l = list(CD.keys())
            c_design, c_response = setupRegressionModel_ContigMedian(CD, contig_l, read_in_contig, D, ReadLength_dict, c_edge_dict)
            # c_estimates_list, c_reads_list, c_index_list = ContigIRLS_Huber(c_design, c_response, contig_l, thr1=200, thr2=1200)
            c_estimates_list, c_reads_list, c_index_list = ContigIRLS_Huber(c_design, c_response, contig_l, thr1=500, thr2=2000) ##试一下thr2=2000，本来是1200
            
            for i, c_contigs in enumerate(c_reads_list):
                c_estimates = c_estimates_list[i]
                m_reads, m_estimates, D_update = MergeContigs_D(c_estimates, c_contigs, CD, split_contigs_reads, split_contigs_estimates,D,D_update,ReadLength_dict)
                Merged_reads_list.append(m_reads)
                Merged_estimates_list.append(m_estimates)
                max_i = np.argmax(m_estimates)
                print("Layout length:", max(m_estimates)-min(m_estimates)+ReadLength_dict[m_reads[max_i]])
        
        split_contigs_reads = []
        split_contigs_estimates = []
        read_in_contig = dict()
        contig_num = 0
        for i in range(len(Merged_estimates_list)):
            resort_ests,resort_reads = SortLayout_aftermerge(Merged_estimates_list[i],Merged_reads_list[i])
            split_contigs = SpiltRegressionLayout(resort_reads,resort_ests,ReadLength_dict,edge_dict,D_update,direct_dict,sd=split_sd*2)
            for sl in split_contigs:
                contig_reads = np.array(resort_reads)[sl]
                contig_estimates = np.array(resort_ests)[sl]
                contig_estimates = [x-contig_estimates[0] for x in contig_estimates]
                split_contigs_reads.append(contig_reads)
                split_contigs_estimates.append(contig_estimates)
                print(i, "Split Layout length:", max(contig_estimates)-min(contig_estimates)+ReadLength_dict[contig_reads[-1]])                
                for r_i,r_name in enumerate(contig_reads):
                    read_in_contig[r_name]=[contig_num,contig_estimates[r_i],r_i,len(contig_reads)-r_i]
                contig_num +=1
        
        print("Total reads:",sum([len(x) for x in split_contigs_reads]))
        iter +=1
        D = D_update.copy()        
    
    print("Loading long reads from %s ..."%args.r)    
    ReadSeq = importLongReads(args.r,sep=args.p)
    now = datetime.datetime.now()
    print("Loaded successfully at time", now.strftime("%Y-%m-%d %H:%M:%S"))

    now = datetime.datetime.now()
    print("Start the SplitRegression Consensus at time ", now.strftime("%Y-%m-%d %H:%M:%S"))
    bw = args.bw
    consensus_dict = dict()
    for layout_num in range(len(split_contigs_reads)):  
        m_reads = split_contigs_reads[layout_num]
        if len(m_reads)<3:
            continue
        m_estimates = split_contigs_estimates[layout_num]
        resort_ests,resort_reads = SortLayout_aftermerge(m_estimates,m_reads)
        BlockConsensus_MultiThread(layout_num,resort_reads,resort_ests,ReadLength_dict,D,consensus_dict,bw,600,t=args.t)
    
    contig_file = open(args.o, 'w')
    
    contig_num = 1
    for k,cns_l in consensus_dict.items():
        for seq_str in cns_l:
            contig_file.write('>contig%d\n'%contig_num)
            writeSequence(contig_file, seq_str)
            contig_num+=1
    contig_file.close()
