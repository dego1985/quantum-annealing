import os
from dwave.cloud import Client
from dwave.system.samplers import DWaveSampler
from dimod import SimulatedAnnealingSampler
import neal
import dimod
import networkx as nx
import dwave_networkx as dnx
from dwave.system.composites import FixedEmbeddingComposite, LazyFixedEmbeddingComposite
import matplotlib.pyplot as plt
from dwave.system import DWaveSampler, EmbeddingComposite
import numpy as np



def v1():
    h = {0:1.0, 1:1.0, 2:0.0, 3:1.0, 
         4:1.0, 5:1.0, 6:0.0, 7:1.0}
    J = {(0,4):+0.5, (0,5):+1.0, (0,6):-0.5, (0,7):+0.5, 
         (1,4):+1.0, (1,5):-1.0, (1,6): 0.0, (1,7):+1.0, 
         (2,4):-0.5, (2,5): 0.0, (2,6): 0.0, (2,7):+0.5, 
         (3,4):+0.5, (3,5):+1.0, (3,6):+0.5, (3,7):+0.5}
    return h, J

def V(unit_size=16,i0=0, j0=0, machine_size=16):
    h1, J1 = v1()
    h = {}
    J = {}

    def offset(i,j):
        return (i+i0  ) * machine_size * 8 + (j+j0  ) * 8
    def offset_right(i,j):
        return (i+i0  ) * machine_size * 8 + (j+j0+1) * 8
    def offset_bottom(i,j):
        return (i+i0+1) * machine_size * 8 + (j+j0  ) * 8

    for i in range(unit_size):
        for j in range(unit_size):
            os = offset(i,j)

            for x in h1.keys():
                h[os + x] = h1[x]

            for x in J1.keys():
                J[(x[0]+os, x[1]+os)] = J1[x]

    # unit to unit
    for i in range(unit_size):
        for j in range(unit_size):
            up = (i+i0+j+j0)%2 == 0
            os  = offset(i,j)
            osr = offset_right(i,j)
            osb = offset_bottom(i,j)

            if up:
                if i < unit_size - 1:
                    J[(3+os, 3+osb)] = -0.5
                if j < unit_size - 1:
                    J[(7+os, 7+osr)] = -0.5
            else:
                if i < unit_size - 1:
                    J[(0+os, 0+osb)] = -0.5
                if j < unit_size - 1:
                    J[(4+os, 4+osr)] = -0.5

    # boarder
    lmap = [[0,1,2,3,4,5,6,7], [3,2,1,0,7,6,5,4]]
    for i in range(unit_size):

        # left
        j = 0
        up = (i+i0+j+j0)%2 
        os  = offset(i,j)
        h[lmap[up][4]+os] += 1.0

        # right
        j = unit_size-1
        up = (i+i0+j+j0)%2 
        os  = offset(i,j)
        h[lmap[up][7]+os] += 1.0
    
    for j in range(unit_size):

        # top
        i = 0
        up = (i+i0+j+j0)%2 
        os  = offset(i,j)
        h[lmap[up][0]+os] += 1.0

        # bottom
        i = unit_size-1
        up = (i+i0+j+j0)%2 
        os  = offset(i,j)
        h[lmap[up][3]+os] += 1.0


    return h, J


def print_chimera_simple(spins, machine_size = 16):

    out = np.full((machine_size, machine_size, 3, 3), '. ')

    # set char
    for i in range(machine_size):
        for j in range(machine_size):

            # offset
            offset = i * machine_size * 8 + j * 8
            out[i,j,1,1] = format(offset//8, '02x')

            # up
            up = (i+j)%2 == 0
            if up:
                lmap = [0,1,2,3,4,5,6,7]
            else:
                lmap = [3,2,1,0,7,6,5,4]

            # nodes
            idx = offset+lmap[0]
            out[i,j,0,1] = '┃' if idx in list(spins.keys()) and spins[idx] > 0 else '  '
            idx = offset+lmap[3]
            out[i,j,2,1] = '┃' if idx in list(spins.keys()) and spins[idx] > 0 else '  '
            idx = offset+lmap[4]
            out[i,j,1,0] = '━' if idx in list(spins.keys()) and spins[idx] > 0 else '  '
            idx = offset+lmap[7]
            out[i,j,1,2] = '━' if idx in list(spins.keys()) and spins[idx] > 0 else '  '

    # print
    out_t = out.transpose(0,2,1,3)
    for o1 in out_t:
        for o2 in o1:
            print(''.join(o2.flatten()))

def find_loop(unit_size=8, i0=0, j0=0, use_qpu=True):
    '''
    dwaveの各セルを頂点として、頂点をつなげるループを
    見つけようとするモデルを実行する。 

    結果は、エネルギー最小値にならないため
    期待した結果にはならない。 

    今後のマシンの精度向上に期待。 
    または、どうにかして頑健なグラフ構造にすべきか。 

    ※注意：マシンは個体ごとにところどころ壊れているので、使用するセルを以下の引数で指定する。 

    Args:
        unit_size:
            セルの正方形の辺の長さ
        i0:
            セルの縦方向のオフセット
        j0:
            セルの横方向のオフセット
        use_qpu:
            qpuを使うか指定する
    '''
    qpu_unit_size = 16
    C16 = dnx.chimera_graph(qpu_unit_size)

    h, J = V(unit_size=unit_size, i0=i0, j0=j0)
    #print(h)
    #print(J)
    if use_qpu:
        sampler = DWaveSampler()
        samples = sampler.sample_ising(h, J, num_reads=10, annealing_time=20)

        if 0:
            count = 0
            for s,e,o in samples.data(['sample', 'energy', 'num_occurrences']):     
                print(e, o)
                print_chimera_simple(s)


                samples = sampler.sample_ising(h, J, 
                        num_reads=100,
                        anneal_schedule = [[0.0,1.0],[20.0,0.0],[30.0,0.3],[35.0,0.3],[40.0,1.0]],
                        reinitialize_state=False,
                        initial_state=s)

                count += 1
                if count >= 1 :
                    break

    else:
        classical_sampler = neal.SimulatedAnnealingSampler()
        sampler = dimod.StructureComposite(classical_sampler, C16.nodes, C16.edges)    
        samples = sampler.sample_ising(h, J, num_reads=10)

    count = 0
    for s,e,o in samples.data(['sample', 'energy', 'num_occurrences']):     
        print(e, o)
        print_chimera_simple(s)
        count += 1
        if count >= 10 :
            break

if __name__ == "__main__":

    find_loop(8,2,3,use_qpu=True)
