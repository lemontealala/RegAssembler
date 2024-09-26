import numpy as np
import networkx as nx

##最小生成树定方向
def STreeDFS(T,node,DFlag,Direction):
	sorted_nei = sorted(T.neighbors(node),key=lambda x:abs(T.edges[x,node]['weight']),reverse=True)
	for nei in sorted_nei:
		if Direction[nei]!=0:
			continue
		Direction[nei]=int(Direction[node]*DFlag[nei,node])
		# print(Direction[node],Direction[nei])
		STreeDFS(T,nei,DFlag,Direction)
	return Direction

##循环版本
def STreeDFS_stack(T, start_node, DFlag):
    Direction = {start_node: 1}
    stack = [start_node]
    while stack:
        node = stack.pop()
        sorted_nei = sorted(T.neighbors(node), key=lambda x: abs(T.edges[x, node]['weight']), reverse=True)
        for nei in sorted_nei:
            if nei in Direction:
                continue
            Direction[nei] = Direction[node] * DFlag[node, nei]
            stack.append(nei)
    return Direction

def TreeDFS(T,node,DFlag,Direction):
	for nei in T.neighbors(node):
		if Direction[nei]!=0:
			continue
		Direction[nei]=int(Direction[node]*DFlag[nei,node])
		#print(Direction[node],Direction[nei])
		TreeDFS(T,nei,DFlag,Direction)
	return Direction
	
def sigmoid(x):
	return 1/(1+np.exp(-x))

def DirectMST(direct_dict,reads,s_node = 0):
	DG = nx.Graph()
	ReadsCount = len(reads)
	# print(ReadsCount)
	ReadsDirection = [0 for i in range(ReadsCount)]
	DFlag = np.zeros((ReadsCount,ReadsCount))
	J_mat = np.zeros((ReadsCount,ReadsCount))
	for n1 in range(ReadsCount):
		for n2 in range(n1,ReadsCount):
			c1 = reads[n1]
			c2 = reads[n2]
			# if c1 in trouble_reads or c2 in trouble_reads:
			# 	continue
			key = str(min(c1,c2))+'-'+str(max(c1,c2))
			d_list = []
			if key in direct_dict:	
				d_list = direct_dict[key]
				# print(key,d_list)				
			if len(d_list)<1:
				continue
			d_sum = sum(d_list)
			d_abssum = sum(np.abs(d_list))
			posi = (d_sum+d_abssum)/2
			nega = (d_abssum-d_sum)/2
			J_mat[n1,n2] = d_sum #np.tanh(d_sum)#sigmoid(d_sum)#np.sign(d_sum)*np.sqrt(abs(d_sum))
			J_mat[n2,n1] = d_sum #np.tanh(d_sum)#sigmoid(d_sum)#np.sign(d_sum)*np.sqrt(abs(d_sum))#d_sum
			if max(posi,nega)>0:
				DG.add_edge(n1,n2,weight=-max(posi,nega)) #np.tanh
				DFlag[n1,n2] = np.sign(d_sum)
				DFlag[n2,n1] = DFlag[n1,n2]
				# print(c1,c2,DFlag[n1,n2],"weight",-max(posi,nega))
			else:
				print("Why?",key,d_list,d_sum,d_abssum,posi,nega)
	MST = nx.minimum_spanning_tree(DG)
	#sorted(G.degree(),key = lambda x:x[1],reverse=True)[0][0]
	Direction = STreeDFS_stack(MST,s_node,DFlag)
 	# ReadsDirection[s_node]=1
	# Direction = STreeDFS(MST,s_node,DFlag,ReadsDirection)
	return (Direction,J_mat)

# 修改，direct_dict权重是align_bases，最小不得小于10k
def DirectMST_contig(direct_dict, reads, truc=10000):
    DG = nx.Graph()
    ReadsCount = len(reads)
    # print(ReadsCount)
    DFlag = np.zeros((ReadsCount,ReadsCount))
    J_mat = np.zeros((ReadsCount,ReadsCount))
    for n1 in range(ReadsCount):
        for n2 in range(n1,ReadsCount):
            c1 = reads[n1]
            c2 = reads[n2]
            key = str(min(c1,c2))+'-'+str(max(c1,c2))
            d_list = []
            if key in direct_dict:	
                d_list = direct_dict[key]
                # print(key,d_list)				
            if len(d_list)<1:
                continue
            d_sum = sum(d_list)
            d_abssum = sum(np.abs(d_list))
            posi = (d_sum+d_abssum)/2
            nega = (d_abssum-d_sum)/2
            J_mat[n1,n2] = d_sum #np.tanh(d_sum)#sigmoid(d_sum)#np.sign(d_sum)*np.sqrt(abs(d_sum))
            J_mat[n2,n1] = d_sum #np.tanh(d_sum)#sigmoid(d_sum)#np.sign(d_sum)*np.sqrt(abs(d_sum))#d_sum
            if max(posi,nega)>truc:
                DG.add_edge(n1,n2,weight=-max(posi,nega)) #np.tanh
                DFlag[n1,n2] = np.sign(d_sum)
                DFlag[n2,n1] = DFlag[n1,n2]
                # print(c1,c2,DFlag[n1,n2],"weight",-max(posi,nega))
            # else:
            #     print("Why?",key,d_list,d_sum,d_abssum,posi,nega)
    Direction = []  #可能有多个连通分支，分别用MST定向
    for subG in nx.connected_components(DG):
        MST = nx.minimum_spanning_tree(DG.subgraph(subG))
        print(len(subG),MST.nodes())
        Direction.append(STreeDFS_stack(MST,list(MST.nodes())[0],DFlag))
    # ReadsDirection[s_node]=1
    # Direction = STreeDFS(MST,s_node,DFlag,ReadsDirection)
    return (Direction,J_mat)
