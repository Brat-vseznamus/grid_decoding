import itertools
import galois
import numpy as np

from typing import Dict, List, Tuple
 
gf = galois.GF(2)

def read_matrix():
    with open("./matrix.txt", 'r') as f:
        G = []
        for line in f.readlines():
            G.append(list(map(int, line.split(' '))))
        return gf(G)

G = read_matrix()

(k, n) = G.shape

def active_range(row):
    idx = [i for i in range(row.shape[0]) if row[i] == 1]
    return (idx[0], idx[len(idx) - 1] - 1)


def make_difference():
    while(True):
        rngs = [active_range(G[i]) for i in range(k)]

        diff_starts = set([r[0] for r in rngs])

        if len(diff_starts) != k:
            for start in diff_starts:
                indx = [i for i in range(k) if start == rngs[i][0]]
                if len(indx) == 1:
                    continue
                
                indx_with_ends = [(i, rngs[i][1]) for i in indx]
                indx_with_ends.sort(key=lambda x: x[1])

                index_with_smallest_end = indx_with_ends[0][0]

                for i in indx:
                    if i != index_with_smallest_end:
                        G[i] += G[index_with_smallest_end]
        
        diff_ends = set([r[1] for r in rngs])

        if len(diff_ends) != k:
            for end in diff_ends:
                indx = [i for i in range(k) if end == rngs[i][1]]
                if len(indx) == 1:
                    continue
                
                indx_with_starts = [(i, rngs[i][0]) for i in indx]
                indx_with_starts.sort(key=lambda x: x[1])
                indx_with_starts.reverse()

                index_with_biggest_end = indx_with_starts[0][0]

                for i in indx:
                    if i != index_with_biggest_end:
                        G[i] += G[index_with_biggest_end]
        
        if len(diff_ends) == len(diff_starts) and len(diff_ends) == k:
            break


make_difference()

print('МФС:')
print(G)

def active_rows(indx, rngs):
    return [i for i in range(len(rngs)) if rngs[i][0] <= indx and indx <= rngs[i][1]]

rngs = [active_range(G[i]) for i in range(k)]

class Vertex:
    vertex_cnt: int = 0
    def __init__(self, values:Dict[int, int]=None):
        self.values = {} if values is None else values
        self.children: List[Tuple[int, Vertex]] = []
        self.w: int = None
        self.num = Vertex.vertex_cnt
        Vertex.vertex_cnt += 1
    
    def __str__(self):
        return str(self.values)

    def __repr__(self):
        return self.__str__()


layers: List[List[Vertex]] = [[Vertex()]] # стартовый слой

for i in range(0, n):
    rows = active_rows(i, rngs)

    layer = []

    for v in list(itertools.product([0, 1], repeat=len(rows))):
        d = {}
        for i in range(len(rows)):
            d[rows[i]] = v[i]
        
        layer.append(Vertex(d))
    
    layers.append(layer)


def connect_layers(layer_id):
    current_layer = layers[layer_id]
    prev_layer = layers[layer_id - 1]

    curr_active_rows = active_rows(layer_id - 1, rngs)
    prev_active_rows = active_rows(layer_id - 2, rngs)

    common_rows = set(curr_active_rows).intersection(set(prev_active_rows))

    diff_curr = set(curr_active_rows).difference(common_rows)

    for pv in prev_layer:
        for cv in current_layer:
            # print(pv, cv)
            if len(list(filter(lambda i: pv.values[i] == cv.values[i], common_rows))) == len(common_rows):
                values = gf([0 for _ in range(k)])
                for i in prev_active_rows:
                    values[i] = pv.values[i]
                for i in diff_curr:
                    values[i] = cv.values[i]
                
                column = gf([G[i][layer_id - 1] for i in range(k)])

                res = values.dot(column.T)
                pv.children.append((res, cv))


def read_vector():
    with open('./y.txt', 'r') as f:
        return gf(list(map(int, f.readline().split(' '))))

for i in range(1, len(layers)):
    connect_layers(i)


Y = read_vector()
layers[0][0].w = 0

for i in range(1, len(layers)):
    prev_layer = layers[i - 1]
    for v in prev_layer:
        for edge in v.children:
            value = edge[0]
            vc = edge[1]

            dw = 0 if value == Y[i - 1] else 1

            if vc.w == None:
                vc.w = v.w + dw
            else:
                vc.w = min(v.w + dw, vc.w)


print()

for i in range(len(layers)):
    print(f'log(layer[{i}]) =', int(np.log2(len(layers[i]))))

print()

with open('graph.txt', 'w') as f:
    for layer in layers:
        for v in layer:
            text = ''.join(map(str, list(v.values.values())))
            # text = ''.join(map(str, list(v.values.values()))) + ",w:" + str(v.w)
            f.write(f"v{v.num}[label=\"{text}\"]\n")
    
    for layer in layers:
        for v in layer:
            for vc in v.children:
                f.write(f"v{v.num} -> v{vc[1].num} [color={'red' if vc[0] == 1 else 'blue' }]\n")


Y_ = []

end = layers[len(layers) - 1][0]
curr_layer_id = len(layers) - 1

while(curr_layer_id > 0):
    prev_layer = layers[curr_layer_id - 1]

    candidates: List[Tuple[Vertex, int]]= []

    for v in prev_layer:
        for edge in v.children:
            if edge[1] == end:
                candidates.append((v, edge[0]))
    
    for c in candidates:
        dw = 0 if c[1] == Y[curr_layer_id - 1] else 1
        if c[0].w + dw == end.w:
            end = c[0]
            Y_.append(0 if 0 == c[1] else 1)
            break
    
    curr_layer_id -= 1

print("original:", Y)

Y_.reverse()

print("fixed:   ", Y_)


