#!/usr/bin/python3
from __future__ import print_function
from __future__ import unicode_literals
from itertools import *
import argparse
import os
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
	line = fp.readline().split(",")
	if len(line) == 1:
		is_general = True
	else:
		is_general = False
	fp.close()

	try:
		g = nx.read_edgelist(in_edges_path, nodetype=int, data=False)
#		g = nx.convert_node_labels_to_integers(g, ordering='sorted')
	except TypeError:
		g = nx.read_edgelist(in_edges_path, data=False)

	node_size = 50
	if g.order() > 400:
		node_size = 5
	if in_edges_path.endswith(".gz"):
		in_edges_path = os.path.splitext(in_edges_path)[0]

	out_image_path = os.path.splitext(in_edges_path)[0] + ".png"
	out_image_path = os.path.basename(out_image_path)
	if is_general:
		layout = nx.circular_layout(sorted(g))
	else:
		layout = grid_layout(g)
	plt.figure(figsize=(args.s,args.s))
#	nx.draw(g, with_labels=args.n, node_size=node_size, linewidths=0, alpha=1, node_color='#000000', edge_color='#000000', pos=layout)
	nx.draw(g, with_labels=args.n, node_size=node_size, linewidths=0, alpha=0.5, node_color='#3399ff', edge_color='#666666', pos=layout, width=2)
#	plt.xlim(-1.03, 1.03)
#	plt.ylim(-1.03, 1.03)
	plt.draw()
	plt.savefig(out_image_path)
	cmd = "open " + out_image_path
	subprocess.call(cmd, shell=True)
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
