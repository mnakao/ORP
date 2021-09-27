#!/usr/bin/python3
from __future__ import print_function
from __future__ import unicode_literals
from itertools import *
import argparse
import os, sys, tempfile
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import subprocess
host_size    = 30
switch_size  = 100
host_color   = "blue"
switch_color = "red"
edge_color   = '#222222'
alpha        = 0.7

argumentparser = argparse.ArgumentParser()
argumentparser.add_argument("edges_path", help="Input edgelist file path")
argumentparser.add_argument("-n", action="store_true", help="Add number label")
argumentparser.add_argument("-s", help="Image size", type=int, default=4)
args = argumentparser.parse_args()

def main():
	in_edges_path = os.path.normpath(args.edges_path)
	assert os.path.isfile(in_edges_path)
	assert args.s > 0

	fp_edge    = open(in_edges_path, "r")
	first_line = fp_edge.readline().split(" ")
	hosts      = int(first_line[0])
	switches   = int(first_line[1])

	fp_temp = tempfile.NamedTemporaryFile()
	while True:
		line = fp_edge.readline()
		if line:
			s = line.split(" ")[0] + " " + line.split(" ")[1] + "\n"
			fp_temp.write(s.encode('utf-8'))
		else:
			break
	fp_temp.seek(0)
	g = nx.read_edgelist(fp_temp.name, nodetype=int, data=False)

	if in_edges_path.endswith(".gz"):
		in_edges_path = os.path.splitext(in_edges_path)[0]

	out_image_path = os.path.splitext(in_edges_path)[0] + ".png"
	out_image_path = os.path.basename(out_image_path)

	node_size = []
	node_color = []
	for i in g.nodes:
		if i < hosts:
			node_size.append(host_size)
			node_color.append(host_color)
		else:
			node_size.append(switch_size)
			node_color.append(switch_color)

	layout_host   = nx.circular_layout([i for i in range(hosts)])
	layout_switch = nx.circular_layout([i + hosts for i in range(switches)])
	for key, value in layout_switch.items():
		layout_switch[key] = value * 0.8

	plt.figure(figsize=(args.s,args.s))
	nx.draw(g, with_labels=args.n, node_size=node_size, linewidths=0, alpha=alpha, node_color=node_color,
			edge_color=edge_color, pos={**layout_host, **layout_switch}, width=2)
	plt.draw()
	plt.savefig(out_image_path)
	cmd = "open " + out_image_path
	subprocess.call(cmd, shell=True)
	
	fp_edge.close()
	fp_temp.close()
	return

if __name__ == '__main__':
	main()
