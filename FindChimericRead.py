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
# 想要根据chimeric read的特征，将每条read的overlap分为比对到前半段和后半段的，
# 判别这两部分是否有交叉。
# 比对到前半段的，即overlap在该read的起点处hangingout很短，计在read_F中，
# 后半段，即overlap在该read的终点处hangingout很短，计到read_R中。
# 记录overlap远端到这一端的距离
import collections
import numpy as np

def DCountOverlap(file):
    readLength = {}
    Hangingout = 50
    # suspious_reads = []
    with open(file, 'r') as blasrout:
        read_olp_count = dict()
        for line in blasrout:
            record = line.strip().split()
            (A_Orientation, B_Orientation, score)=list(map(int,record[2:5]))
            A_name,B_name = map(str,[record[0].split('.')[1],record[1].split('.')[1]])
            if A_name==B_name:
                continue
            (B_start, B_end, B_length, A_start, A_end, A_length) = list(map(int,record[6:12]))
            readLength[A_name] = A_length
            readLength[B_name] = B_length
            if B_Orientation==0 and B_length-B_end<Hangingout and A_start<Hangingout:
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_R'] = [B_length-B_start]
            elif B_Orientation==0 and A_length-A_end<Hangingout and B_start<Hangingout:
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_end]
                else:
                    read_olp_count[B_name+'_F'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            elif B_Orientation==1 and B_length-B_end<Hangingout and A_start<Hangingout:
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_F'] = [B_length-B_start]
            elif B_Orientation==1 and A_length-A_end<Hangingout and B_start<Hangingout:
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_end]
                else:
                    read_olp_count[B_name+'_R'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            # else:
                # print("Hangingout!",record)
                # suspious_reads+=[A_name,B_name]
    return read_olp_count,readLength

def FindChimeria(file_name,N):
    Read_olp_count,readLength = DCountOverlap(file_name)
    candidate_chimera = []
    for i in range(1,N+1):
        x = str(i)
        if x+'_F' in Read_olp_count and x+'_R' in Read_olp_count:
            LF = max(Read_olp_count[x+'_F'])+max(Read_olp_count[x+'_R'])
            # print(x,LF,readLength[x], LF>readLength[x])
            if LF<readLength[x]:
                candidate_chimera.append(x)
        else:
            candidate_chimera.append(x)
    return candidate_chimera

def DCountOverlap_MHAP(file):
    readLength = {}
    Hangingout = 50
    # suspious_reads = []
    with open(file, 'r') as blasrout:
        read_olp_count = dict()
        for line in blasrout:
            record = line.strip().split()
            (A_Orientation,A_start,A_end,A_length,B_Orientation,B_start,B_end,B_length) = list(map(int,record[4:]))            
            A_name,B_name = map(str,record[:2])
            if B_Orientation==1:
                T = B_end 
                B_end = B_length-B_start
                B_start = B_length-T
            if A_name==B_name:
                continue
            readLength[A_name] = A_length
            readLength[B_name] = B_length
            if B_Orientation==0 and B_length-B_end<Hangingout and A_start<Hangingout:
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_R'] = [B_length-B_start]
            elif B_Orientation==0 and A_length-A_end<Hangingout and B_start<Hangingout:
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_end]
                else:
                    read_olp_count[B_name+'_F'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            elif B_Orientation==1 and B_length-B_end<Hangingout and A_start<Hangingout:
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_F'] = [B_length-B_start]
            elif B_Orientation==1 and A_length-A_end<Hangingout and B_start<Hangingout:
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_end]
                else:
                    read_olp_count[B_name+'_R'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            # else:
                # print("Hangingout!",record)
                # suspious_reads+=[A_name,B_name]
    return read_olp_count,readLength

def FindChimeria_MHAP(file_name,N):
    Read_olp_count,readLength = DCountOverlap_MHAP(file_name)
    candidate_chimera = []
    for i in range(1,N+1):
        x = str(i)
        if x+'_F' in Read_olp_count and x+'_R' in Read_olp_count:
            LF = max(Read_olp_count[x+'_F'])+max(Read_olp_count[x+'_R'])
            # print(x,LF,readLength[x], LF>readLength[x])
            if LF<readLength[x]:
                candidate_chimera.append(x)
        else:
            candidate_chimera.append(x)
    return candidate_chimera


def DCountOverlap_minimap2(file, sep='.', Hangingout = 150):
    readLength = {}
    with open(file, 'r') as minimappaf:
        read_olp_count = collections.defaultdict(list)
        for line in minimappaf:
            record = line.strip().split()
            (A_length,A_start,A_end) = list(map(int,record[1:4]))
            (B_length,B_start,B_end) = list(map(int,record[6:9]))
            A_name = str(record[0].split(sep)[1])
            B_name = str(record[5].split(sep)[1])
            if A_name==B_name :
                continue
            B_Orientation = record[4]
            readLength[int(A_name)] = A_length
            readLength[int(B_name)] = B_length
            #Case 1:
            if B_Orientation=='+' and B_length-B_end<Hangingout and A_start<Hangingout:
                # print("Case 1:", A_name+'_F',A_end, B_name+'_R', B_length-B_start)
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_R'] = [B_length-B_start]
            elif B_Orientation=='+' and A_length-A_end<Hangingout and B_start<Hangingout:
                # print("Case 2:", B_name+'_F', B_end, A_name+'_R', A_length-A_start)
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_end]
                else:
                    read_olp_count[B_name+'_F'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            elif B_Orientation=='-' and B_start<Hangingout and A_start<Hangingout:
                # print("Case 3:", A_name+'_F',A_end, B_name+'_F', B_start)
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_end]
                else:
                    read_olp_count[B_name+'_F'] = [B_end]
            elif B_Orientation=='-' and A_length-A_end<Hangingout and B_length-B_end<Hangingout:
                # print("Case 4:", A_name+'_R',A_length-A_start, B_name+'_R', B_length-B_start)
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_R'] = [B_length-B_start]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            else:
                # print("Hangingout!")
                continue
    return read_olp_count,readLength

# 在大基因组上跑的实在太慢 N是read count,计算复杂度是O（N+M)，M是minimap2 paf文件行数
def FindChimeria_minimap2(file_name,sep='.',Hangingout=150):
    print("Counting overlap ...")
    Read_olp_count,readLength= DCountOverlap_minimap2(file_name,sep=sep,Hangingout=Hangingout)
    
    print("Filter by alignment lowcoverage ...")
    candidate_chimera = []    
    for x,l in readLength.items():
        LF = max(Read_olp_count[str(x)+'_F'],default=0)+max(Read_olp_count[str(x)+'_R'],default=0)
        # print(x,LF,len(Read_olp_count[x+'_F']), len(Read_olp_count[x+'_R']), readLength[x], LF>readLength[x])
        if LF<l+300:
            candidate_chimera.append(x)
    return candidate_chimera,readLength

## 可以在去除chimeria的同时记录align_bases, align_ratio低的个数 短而碎的比对特别多的reads

def DCountOverlap_minimap2_low(file, sep='.', Hangingout = 150):
    readLength = {}
    # suspious_reads = []
    with open(file, 'r') as minimappaf:
        read_olp_count = collections.defaultdict(list)
        read_lowolp_count = collections.defaultdict(int)
        for line in minimappaf:
            record = line.strip().split()
            (A_length,A_start,A_end) = list(map(int,record[1:4]))
            (B_length,B_start,B_end) = list(map(int,record[6:9]))
            A_name = str(record[0].split(sep)[1])
            B_name = str(record[5].split(sep)[1])
            align_bases = int(record[9])
            align_ratio = int(record[9])/int(record[10])
            if align_bases<300 or align_ratio<0.6:
                read_lowolp_count[A_name] +=1
                read_lowolp_count[B_name] +=1                
            if A_name==B_name :
                continue
            B_Orientation = record[4]
            readLength[A_name] = A_length
            readLength[B_name] = B_length
            #Case 1:
            if B_Orientation=='+' and B_length-B_end<Hangingout and A_start<Hangingout:
                # print("Case 1:", A_name+'_F',A_end, B_name+'_R', B_length-B_start)
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_R'] = [B_length-B_start]
            elif B_Orientation=='+' and A_length-A_end<Hangingout and B_start<Hangingout:
                # print("Case 2:", B_name+'_F', B_end, A_name+'_R', A_length-A_start)
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_end]
                else:
                    read_olp_count[B_name+'_F'] = [B_end]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            elif B_Orientation=='-' and B_start<Hangingout and A_start<Hangingout:
                # print("Case 3:", A_name+'_F',A_end, B_name+'_F', B_start)
                if A_name+'_F' in read_olp_count:
                    read_olp_count[A_name+'_F'] += [A_end]
                else:
                    read_olp_count[A_name+'_F'] = [A_end]
                if B_name+'_F' in read_olp_count:
                    read_olp_count[B_name+'_F'] += [B_end]
                else:
                    read_olp_count[B_name+'_F'] = [B_end]
            elif B_Orientation=='-' and A_length-A_end<Hangingout and B_length-B_end<Hangingout:
                # print("Case 4:", A_name+'_R',A_length-A_start, B_name+'_R', B_length-B_start)
                if B_name+'_R' in read_olp_count:
                    read_olp_count[B_name+'_R'] += [B_length-B_start]
                else:
                    read_olp_count[B_name+'_R'] = [B_length-B_start]
                if A_name+'_R' in read_olp_count:
                    read_olp_count[A_name+'_R'] += [A_length-A_start]
                else:
                    read_olp_count[A_name+'_R'] = [A_length-A_start]
            else:
                # print("Hangingout!")
                continue
    return read_olp_count,readLength,read_lowolp_count
