import os
import sys
import gzip

def calc_jacard_index(setA, setB):
    if len(setA) == 0 or len(setB) == 0:
        return 1
    else:
        return len(setA & setB) / len(setA | setB)

id2symbol_file, cluster_file = sys.argv[1:]

id2symbol = {}
for line in gzip.open(id2symbol_file, "rt"):
    if line.startswith("#"):
        continue
    
    prot_id, symbol, _, _ = line.strip().split("\t")
    id2symbol[prot_id] = symbol

cluster = {}
for line in gzip.open(cluster_file, "rt"):
    if line.startswith("#"):
        continue

    _, cluster_id, prot_id = line.strip().split()
    if not cluster_id in cluster:
        cluster[cluster_id] = set()
    cluster[cluster_id].add(id2symbol[prot_id])

clusters_sorted = sorted([ [cluster_id, len(cluster[cluster_id])] for cluster_id in cluster], key=lambda x:x[1], reverse=True)
select_clusters = [clusters_sorted[0][0]]
for cluster_id, _ in clusters_sorted:
    ## 计算与已经选中的cluster的Jacard index
    jacard_index = [ calc_jacard_index(cluster[c], cluster[cluster_id]) for c in select_clusters ]
    if sum([ i<0.3 for i in jacard_index]) == len(select_clusters):
        select_clusters.append(cluster_id)
    
    if len(select_clusters) == 20:
        break
tmp = set()
for cluster_id in select_clusters:
    print(cluster_id, ",".join(cluster[cluster_id]), sep="\t")
    tmp = tmp | cluster[cluster_id]
sys.stderr.write(f"{len(tmp)}, {len(tmp)/21895}\n")

