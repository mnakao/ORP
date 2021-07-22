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

argumentparser = argparse.ArgumentParser()
argumentparser.add_argument("edges_path", help="Input edgelist file path")
argumentparser.add_argument("-n", action="store_true", help="Add number label")
argumentparser.add_argument("-s", help="Image size", type=int, default=4)
args = argumentparser.parse_args()

def main():
	in_edges_path = os.path.normpath(args.edges_path)
	assert os.path.isfile(in_edges_path)
	assert args.s > 0
	fp = open(in_edges_path, "r")
	line = fp.readline().split(" ")
	hosts = line[0]
	switches = line[1]
	radix = line[2]
	print("Hosts: " + hosts + " Switches: " + switches + " Radix: " + radix)
	row_no = 0
	fp_tmp = tempfile.NamedTemporaryFile()
	hosts = int(hosts)
	while True:
		line = fp.readline()
		if line:
			n0 = int(line.split(" ")[0])
			n1 = int(line.split(" ")[1])
			if(n0 >= hosts and n1 >= hosts):
				s = str(n0) + " " + str(n1) + "\n"
				fp_tmp.write(s.encode('utf-8'))
		else:
			break
	fp_tmp.seek(0)
	g = nx.read_edgelist(fp_tmp.name, nodetype=int, data=False)
	
	node_size = 50
	if g.order() > 400:
		node_size = 5
	if in_edges_path.endswith(".gz"):
		in_edges_path = os.path.splitext(in_edges_path)[0]

	out_image_path = os.path.splitext(in_edges_path)[0] + ".png"
	out_image_path = os.path.basename(out_image_path)
	layout = nx.circular_layout(sorted(g))
	plt.figure(figsize=(args.s,args.s))
	nx.draw(g, with_labels=args.n, node_size=node_size, linewidths=0, alpha=0.5, node_color='#3399ff', edge_color='#666666', pos=layout, width=2)
	plt.draw()
	plt.savefig(out_image_path)
	cmd = "open " + out_image_path
	subprocess.call(cmd, shell=True)
	fp_tmp.close()
	fp.close()
	return

def grid_layout(g):
	layout = {}
	xmin = float("inf")
	ymin = float("inf")
	xmax = 0
	ymax = 0
	for n in g.nodes():
		x, y = map(int, n.split(','))
		if x < xmin:
			xmin = x
		if y < ymin:
			ymin = y
		if xmax < x:
			xmax = x
		if ymax < y:
			ymax = y
	x0 = float(xmin + xmax) / 2
	y0 = float(ymin + ymax) / 2
	r = float(max(xmax - xmin, ymax - ymin)) / 2
# 	print("min=({}, {}), max=({}, {}), center=({}, {}), radius={}".format(xmin, xmax, ymin, ymax, x0, y0, r))
	for n in g.nodes():
		x, y = map(int, n.split(','))
		layout[n] = [float(x - x0) / r, float(y - y0) / r]
	return layout

if __name__ == '__main__':
	main()
