import os
from dwave.system.samplers import DWaveSampler
import numpy as np

def print_chimera_active(nodes, edges, unit_size = 16):
    '''
    dwaveの使用できるノードとエッジを示す。

    - ○はノード。
    - ◇はChimeraグラフの1セル内の２部グラフのエッジ。
    - □はセル間のエッジ。
    - ●、◆、■は壊れて使用できないノードとエッジを表す。

    左上のノードの番号は

        4
        5
    0 1   2 3
        6
        7
    
    エッジの状態は、ノードを通る交点に書いた。

    中央の数字はノード番号のオフセット。
    
    ※注意全角文字は、全角文字の幅で表示しないと崩れる。
    '''
    out = np.full((unit_size, unit_size, 6, 6), '. ')

    # set char
    for i in range(unit_size):
        for j in range(unit_size):
            # offset
            offset = i * unit_size * 8 + j * 8
            out[i,j,2,2] = format(offset//8, '02x')

            # node
            out[i,j,2,0] = '○' if offset+0 in nodes else '●'
            out[i,j,2,1] = '○' if offset+1 in nodes else '●'
            out[i,j,2,3] = '○' if offset+2 in nodes else '●'
            out[i,j,2,4] = '○' if offset+3 in nodes else '●'
            out[i,j,0,2] = '○' if offset+4 in nodes else '●'
            out[i,j,1,2] = '○' if offset+5 in nodes else '●'
            out[i,j,3,2] = '○' if offset+6 in nodes else '●'
            out[i,j,4,2] = '○' if offset+7 in nodes else '●'

            # connection inner
            out[i,j,0,0] = '◇' if (offset+0, offset+4) in edges else '◆'
            out[i,j,1,0] = '◇' if (offset+0, offset+5) in edges else '◆'
            out[i,j,3,0] = '◇' if (offset+0, offset+6) in edges else '◆'
            out[i,j,4,0] = '◇' if (offset+0, offset+7) in edges else '◆'

            out[i,j,0,1] = '◇' if (offset+1, offset+4) in edges else '◆'
            out[i,j,1,1] = '◇' if (offset+1, offset+5) in edges else '◆'
            out[i,j,3,1] = '◇' if (offset+1, offset+6) in edges else '◆'
            out[i,j,4,1] = '◇' if (offset+1, offset+7) in edges else '◆'

            out[i,j,0,3] = '◇' if (offset+2, offset+4) in edges else '◆'
            out[i,j,1,3] = '◇' if (offset+2, offset+5) in edges else '◆'
            out[i,j,3,3] = '◇' if (offset+2, offset+6) in edges else '◆'
            out[i,j,4,3] = '◇' if (offset+2, offset+7) in edges else '◆'

            out[i,j,0,4] = '◇' if (offset+3, offset+4) in edges else '◆'
            out[i,j,1,4] = '◇' if (offset+3, offset+5) in edges else '◆'
            out[i,j,3,4] = '◇' if (offset+3, offset+6) in edges else '◆'
            out[i,j,4,4] = '◇' if (offset+3, offset+7) in edges else '◆'

            # connection outer to right
            if j < unit_size - 1:
                offset_R = offset + 8
                out[i,j,0,5] = '□' if (offset+4, offset_R+4) in edges else '■'
                out[i,j,1,5] = '□' if (offset+5, offset_R+5) in edges else '■'
                out[i,j,3,5] = '□' if (offset+6, offset_R+6) in edges else '■'
                out[i,j,4,5] = '□' if (offset+7, offset_R+7) in edges else '■'

            # connection outer to buttom
            if i < unit_size - 1:
                offset_B = offset + 8 * unit_size
                out[i,j,5,0] = '□' if (offset+0, offset_B+0) in edges else '■'
                out[i,j,5,1] = '□' if (offset+1, offset_B+1) in edges else '■'
                out[i,j,5,3] = '□' if (offset+2, offset_B+2) in edges else '■'
                out[i,j,5,4] = '□' if (offset+3, offset_B+3) in edges else '■'

    # print
    out_t = out.transpose(0,2,1,3)
    for o1 in out_t:
        for o2 in o1:
            print(''.join(o2.flatten()))

if __name__ == '__main__':
    sampler = DWaveSampler()
    print_chimera_active(sampler.nodelist, sampler.edgelist)

