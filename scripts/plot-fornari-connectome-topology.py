#!/usr/bin/env python3
import argparse, csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def read_data(nodes_file,edges_file):
    with nodes_file.open(newline="") as stream: nodes=list(csv.DictReader(stream))
    with edges_file.open(newline="") as stream: edges=list(csv.DictReader(stream))
    return nodes,edges

def scale_sizes(values):
    values=np.asarray(values,float)
    return 10.+55.*(values-values.min())/(values.max()-values.min())

def draw_graph(axis,xyz,edges,values,title,label):
    for edge in edges:
        source,target=int(edge["source"]),int(edge["target"])
        axis.plot(xyz[[source,target],0],xyz[[source,target],1],xyz[[source,target],2],
                  color=".55",linewidth=.25,alpha=.13,zorder=1)
    scatter=axis.scatter(xyz[:,0],xyz[:,1],xyz[:,2],c=values,s=scale_sizes(values),
                         cmap="plasma",edgecolors="white",linewidths=.25,
                         depthshade=False,zorder=5)
    axis.set_title(title,fontsize=10)
    axis.view_init(elev=72,azim=-90)
    axis.set_axis_off()
    colorbar=axis.figure.colorbar(scatter,ax=axis,shrink=.58,pad=.01)
    colorbar.set_label(label,fontsize=8)
    colorbar.ax.tick_params(labelsize=7)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--nodes",type=Path,required=True)
    parser.add_argument("--edges",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    nodes,edges=read_data(args.nodes,args.edges)
    n=len(nodes)
    xyz=np.array([[float(node[k]) for k in ("x","y","z")] for node in nodes])
    adjacency=np.zeros((n,n)); degree=np.zeros(n)
    for edge in edges:
        source,target=int(edge["source"]),int(edge["target"])
        weight=float(edge["connectivity_weight"])
        adjacency[source,target]+=weight; adjacency[target,source]+=weight
        degree[source]+=1; degree[target]+=1
    weighted_degree=adjacency.sum(axis=1)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/"connectivity_diagnostics.csv").open("w",newline="") as stream:
        writer=csv.writer(stream)
        writer.writerow(["node_id","name","unweighted_degree","weighted_degree"])
        for node,unweighted,weighted in zip(nodes,degree,weighted_degree):
            writer.writerow([node["node_id"],node["name"],int(unweighted),weighted])
    figure=plt.figure(figsize=(11.2,3.8),constrained_layout=True)
    draw_graph(figure.add_subplot(1,3,1,projection="3d"),xyz,edges,degree,
               "(a) Non-weighted degree",r"$D_{II}$")
    draw_graph(figure.add_subplot(1,3,2,projection="3d"),xyz,edges,weighted_degree,
               "(b) Connectivity-weighted degree",r"$D_{II}$")
    axis=figure.add_subplot(1,3,3)
    image=axis.imshow(adjacency,cmap="magma",origin="lower",interpolation="nearest")
    axis.set(title=r"(c) Adjacency $A_{IJ}$",xlabel="region index $J$",ylabel="region index $I$")
    colorbar=figure.colorbar(image,ax=axis,shrink=.73); colorbar.set_label(r"$A_{IJ}$")
    output=args.output_dir/"connectome_topology.pdf"
    figure.savefig(output,dpi=300); figure.savefig(output.with_suffix(".png"),dpi=300)
    print(f"Saved {output}")
    print(f"nodes={n}, edges={len(edges)}, degree=[{degree.min():g},{degree.max():g}], "
          f"weighted degree=[{weighted_degree.min():.6g},{weighted_degree.max():.6g}]")
if __name__=="__main__": main()
