## Overview
This library is for Order/Radix Problem with a noweight undirected graph.
You can use following.
* liborp.a : Serial version
* liborp_threads.a : Threads version

## Algotherm for hASPL (hosts Average Shortest Path Length)
It uses an algorithm that is an extension of the following paper (Open Access). _When you write a paper using this library, please refer to the paper._
* https://dl.acm.org/doi/10.1145/3368474.3368478

@inproceedings{10.1145/3368474.3368478,
  author = {Nakao, Masahiro and Murai, Hitoshi and Sato, Mitsuhisa},
  title = {Parallelization of All-Pairs-Shortest-Path Algorithms in Unweighted Graph},
  year = {2020},
  isbn = {9781450372367},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3368474.3368478},
  doi = {10.1145/3368474.3368478},
  booktitle = {Proceedings of the International Conference on High Performance Computing in Asia-Pacific Region},
  pages = {63–72},
  numpages = {10},
  keywords = {network, hybrid parallelization, graph theory, GPU},
  location = {Fukuoka, Japan},
  series = {HPCAsia2020}
}

## Quick start
### Serial version
```
$ git clone https://github.com/mnakao/ORP.git
$ cd ORP
$ make
$ cd ./sample
$ make
$ ./simple.x ./graph/h32s26r4.edges
Hosts = 32, Switches 26, Radix = 4
Diameter     = 11
Diameter Gap = 6 (11 - 5)
ASPL         = 5.9677419355 (2960/496)
ASPL Gap     = 1.7741935484 (5.9677419355 - 4.1935483871)
```

### Threads version
```
$ git clone https://github.com/mnakao/ORP.git
$ cd ORP
$ make threads
$ cd ./sample
$ make threads
$ ./threads_simple.x ./graph/h32s26r4.edges
Hosts = 32, Switches 26, Radix = 4
Diameter     = 11
Diameter Gap = 6 (11 - 5)
ASPL         = 5.9677419355 (2960/496)
ASPL Gap     = 1.7741935484 (5.9677419355 - 4.1935483871)
```

## How to create libraries
```
$ cd ORP
$ make [serial|threads|all]
```

## How to use libraries
In order to create the executable, you need to include the header file `./include/orp.h` and link it with one of the libraries.

## How to compile sample programs
```
$ cd ORP/sample
$ make [serial|threads|all]
```

## File format for graph
* The first line describes the number of hosts (h), the number of switches (s), and the radix of the switches (r).
* The hosts must be integers 0, 1, 2, ..., h-2, h-1.
* The switches must be integers h, h+1, h+2, ..., h+s-2, h+s-1.
* The switch can have up to r edges.
* From the second line, describe the vertex numbers (from 0 to h+s-1) of the start and end points of the edge.
* The switch connects to the host or switch.
* The host cannot connect to the host.
* For example, see "sample/graph/".

Following figure and its edge file are examples of a graph (h = 15, s = 4, r = 6), which is referred to [1].

![](misc/img/sample.png)

```
15 4 6
0 15
5 15
6 15
7 15
1 16
8 16
9 16
2 17
10 17
11 17
12 17
3 18
4 18
13 18
14 18
15 16
16 17
16 18
17 18
```

## Environment variable
### ORP_ASPL=BFS

This library provides two algorithms for ASPL.
The default setting calculates the APSL using a matrix product method.
Apart from that, it also provides a method using breadth-first search (BFS).
In most cases, the matrix product method is faster, but the BFS method consumes less memory.
In default setting, the matrix product method is used.
In case of ORP_ASPL=BFS, BFS is used.

### ORP_PROFILE=1

Output the performance profile for ORP_Set_aspl(). 
This profile is output when ORP_Finalize_aspl() is executed.
```
$ ORP_ASPL=BFS ORP_PROFILE=1 ./simple.x ./graph/h32s26r4.edges
------ Profile for SET_ASPL ------
Date            = Thu May  6 17:49:08 2021
Hostname        = crimson.r-ccs27.riken.jp
Number of Times = 1
Total Time      = 0.000026 sec.
Average Time    = 0.000026 sec.
Algorithm       = BFS (SERIAL)
Memory Usage    = 0.000 MB
Num of Threads  = 1
--------- End of Profile ---------
Hosts = 32, Switches 26, Radix = 4
Diameter     = 11
Diameter Gap = 6 (11 - 5)
ASPL         = 5.9677419355 (2960/496)
ASPL Gap     = 1.7741935484 (5.9677419355 - 4.1935483871)
```

The meaning of each item in the profile is as follows.
* Date : The time when the profile was output.
* Hostname : Name of the machine on which the program ran.
* Number of Times : Number of times ORP_Set_aspl() was executed.
* Total Time : Total execution time of ORP_Set_aspl().
* Average Time : Average execution time of ORP_Set_aspl().
* Algorithm : MATRIX or BFS. The parentheses are the types of libraries used. That is, SERIAL or THREADS.
* Memory Usage : Estimated ammount of memory used in the library.
* Num of Threads : Number of threads used in the library.

## Basic Function
### Overview
Basic function calculates ASPL and diameter of a graph.
```
#include "orp.h"

int main()
{
  ...
  ORP_Init_aspl(...);
  for(int i=0;i<ITERATIONS;i++){
    /* Optimization */
    ORP_Set_aspl(...);
  }
  ORP_Finalize_aspl();
  ...
  return 0;
}
```
ORP_Set_aspl() calculates APSL and diameter of a graph.
While ORP_Init_aspl() initializes for ORP_Set_aspl(), ORP_Finalize_aspl() finalizes for ORP_Set_aspl().
Thus, ORP_Set_aspl() should be called between ORP_Init_aspl() and ORP_Finalize_aspl().
Note that ORP_Init_aspl() and ORP_Finalize_aspl() will basically be called only once each in a program.

### Initialize
Perform the initialization process before executing ORP_Set_aspl().
```
void ORP_Init_aspl(int hosts, int switches, int radix)
```
* [IN] hosts : Number of hosts.
* [IN] switches: Number of switches.
* [IN] radix : Radix of the switch.

### Set diameter, sum, and ASPL
Set diameter, sum, and ASPL. 
In the case of an unconnected graph, INT_MAX, LONG_MAX, and DBL_MAX are assigned to the values of diameter, sum, and ASPL, respectively.
```
void ORP_Set_aspl(int h_degree[switches], int s_degree[switches], int adjacency[switches][radix], int *diameter, long *sum, double *ASPL)
```
* [IN] h_degree : Number of hosts connected to each switch.
* [IN] s_degree : Number of switches connected to each switch.
* [IN] adjacency : Adjacency matrix.
* [OUT] diameter : Diameter of a graph.
* [OUT] sum : Total value of the distances between each host in a graph.
* [OUT] ASPL : Average shortest path length of a graph (sum = ASPL*(hosts*(hosts-1)/2)).

In the case of the above figure, h_degree, s_degree and adjacency are as follows.
* h_degree[] = {4, 3, 4, 4}
* s_degree[] = {1, 3, 2, 2}
* adjacency[][] = {{1,-1,-1,-1},{0,2,3,-1},{1,3,-1,-1},{1,2,-1,-1}}

"h_degree[0]=4" means that s0 connects to four hosts.
"s_degree[1]=3" means that s1 connects to three switches.
"adjacency[2][]={1,3,-1,-1}" means that s2 connects to s1 and s3 (-1 means no value).

### Finalize
Release the resources allocated in ORP_Init_aspl().
```
void ORP_Finalize_aspl()
```

## Utility
### Read a property of a file
```
void ORP_Read_property(char* fname, int* host, int* switches, int* radix, int *lines)
```
* [IN] fname : File name of a graph.
* [OUT] hosts : Number of hosts.
* [OUT] switches : Number of switches.
* [OUT] radix : Radix of the switch.
* [OUT] lines : Number of lines in an edge list.

### Read an edge from a file
```
void ORP_Read_edge(char* fname, int edge[lines][2])
```
* [IN] fname : File name of a graph.
* [OUT] edge : Edge list.

In the case of the above figure, lines and edge are as follows.
* lines = 19
* edge[][] = {{0,15},{5,15},{6,15},{7,15},{1,16},{8,16},{9,16},{2,17},{10,17},{11,17},{12,17},{3,18},{4,18},{13,18},{14,18},{15,16},{16,17},{16,18},{17,18}}

That is, lines is the number of lines in the file excluding the first line, and edge is the value after the second line.

### Write an edge to a file
```
void ORP_Write_edge(int hosts, int switches, int radix, int lines, void *edge, char *fname)
```
* [IN] hosts : Number of hosts.
* [IN] switches : Number of switches.
* [IN] radix : Radix of the switch.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list.
* [OUT] fname : File name of a graph.

### Print an adjacency matrix
```
void ORP_Print_adjacency(int hosts, int switches, int radix, int s_degree[switches], int adjacency[switches][radix])
```
* [IN] hosts : Number of hosts.
* [IN] switches : Number of switches.
* [IN] radix : Radix of the switch.
* [IN] s_degree : Number of switches connected to each switch.
* [IN] adjacency : Adjacency matrix.

### Print an edge list
```
void ORP_Print_edge(int lines, int edge[lines][2])
```
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list.

### Print an switch list
```
void ORP_Print_switch(int hosts, int lines, int edge[lines][2])
```
* [IN] hosts : Number of hosts.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list.

### Print an degree
```
void ORP_Print_degree(int hosts, int switches, int degree[switches])
```
* [IN] hosts : Number of hosts.
* [IN] switches : Number of switches.
* [IN] degree : Degree of host or switch (h_degree or s_degree).

### Convert an edge list to an adjacency matrix
```
void ORP_Conv_edge2adjacency(int hosts, int switches, int radix, int lines, int edge[lines][2], int adjacency[switches][radix])
```
* [IN] hosts : Number of hosts.
* [IN] switches : Number of switches.
* [IN] radix : Radix of the switch.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list.
* [OUT] adjacency : Adjacency matrix.

### Convert an adjacency matrix to an edge list
```
void ORP_Conv_adjacency2edge(int hosts, int switches, int radix, int h_degree[switches], int s_degree[switches], int adjacency[switches][radix], int edge[lines][2])
```
* [IN] hosts : Number of hosts.
* [IN] switches : Number of switches.
* [IN] radix : Radix of the switch.
* [IN] h_degree : Number of hosts connected to each switch.
* [IN] s_degree : Number of switches connected to each switch.
* [IN] adjacency : Adjacency matrix.
* [OUT] edge : Edge list.

### Set theoretical lower bounds
```
void ORP_Set_lbounds(int hosts, int radix, int *low_diameter, double *low_ASPL)
```
* [IN] hosts : Number of hosts.
* [IN] radix : Radix of the switch.
* [OUT] low_diameter : Theoretical lower bound of diameter.
* [OUT] low_ASPL : Theoretical lower bound of ASPL.

### Set degrees for a graph
ORP_Set_host_degree() and ORP_Set_switch_degree() set h_degree and s_degree, respectively.
ORP_Set_degrees() sets both h_degree and s_degree.
```
void ORP_Set_host_degree  (int hosts, int switches, int lines, int edge[lines][2], int h_degree[switches])
void ORP_Set_switch_degree(int hosts, int switches, int lines, int edge[lines][2], int s_degree[switches])
void ORP_Set_degrees      (int hosts, int switches, int lines, int edge[lines][2], int h_degree[switches], int s_degree[switches])
```
* [IN] hosts : Number of hosts.
* [IN] switches : Number of switches.
* [IN] lines : Number of lines in an edge list.
* [IN] edge : Edge list.
* [OUT] h_degree : Number of hosts connected to each switch.
* [OUT] s_degree : Number of switches connected to each switch.

### Seed for a random number
A random number is used in ORP_Generate_random(), ORP_Swap_adjacency(), and ORP_Swing_adjacency().
Thus, ORP_Srand() must be executed before these functions.
```
void ORP_Srand(unsigned int seed)
```
* [IN] seed : Seed for random.

### Generate a random graph
Generate a graph with randomly connected vertices. Note that the graph may contain multiple edges and loops.
When assigned_evenly=true, it assigns evenly hosts to switches.
```
void* ORP_Generate_random(int hosts, int switches, int radix, bool assign_evenly, int *lines, int *h_degree, int *s_degree)
```
* [IN] hosts : Number of hosts.
* [IN] switches : Number of switches.
* [IN] radix : Radix of the switch.
* [IN] assign_evenly : Whether to assign evenly hosts to switches.
* [OUT] lines : Number of lines in an edge list.
* [OUT] h_degree : Number of hosts connected to each switch.
* [OUT] s_degree : Number of switches connected to each switch.
* [RETURN] edge : Edge list.

### Mutate an adjacency matrix
Mutate an adjacency matrix slightly.
ORP_Swap_adjacency() and ORP_Swing_adjacency() perform the SWAP and SWING operations described in ref. [1], respectively.
The changes are stored in the ORP_Resotre structure variable.
```
void ORP_Swap_adjacency(int switches, int radix, int s_degree[switches], ORP_Restore *restore, int adjacency[switches][radix])
```
* [IN] switches : Number of switches.
* [IN] radix : Radix of the switch.
* [IN] s_degree : Number of switches connected to each switch.
* [OUT] restore : Changes.
* [OUT] adjacency : Adjacency matrix.

```
void ORP_Swing_adjacency(int switches, int radix, int h_degree[switches], int s_degree[switches], ORP_Restore *restore, int adjacency[switches][radix])
```
* [IN] switches : Number of switches.
* [IN] radix : Radix of the switch.
* [OUT] h_degree : Number of hosts connected to each switch.
* [OUT] s_degree : Number of switches connected to each switch.
* [OUT] restore : Changes.
* [OUT] adjacency : Adjacency matrix.

### Restore an adjacency matrix
Changes are undone by the ORP_Resotre structure variable.
```
void ORP_Restore_adjacency(ORP_Restore r, int radix, int h_degree[switches], int s_degree[switches], int adjacency[switches][radix])
```
* [IN] restore : Changes.
* [IN] radix : Radix of the switch.
* [OUT] h_degree : Number of hosts connected to each switch.
* [OUT] s_degree : Number of switches connected to each switch.
* [OUT] adjacency : Adjacency matrix.

## Performance
* On Cygnus system in University of Tsukuba, Japan
* Xeon Gold 6126 2.6GHz (12cores) x 2, gcc/8.3.1
* Graph with (hosts, switches, radix) = (65536, 3000, 64)
* http://research.nii.ac.jp/graphgolf/solutions/h65536r64k5l453.20210426-ftr2ae.edges.gz
  * libapsp.a with BFS : 602.2 msec.
  * libapsp.a with matrix product method : 52.6 msec.
  * libapsp_threads.a with matrix product method : 9.9 msec. (12 threads, 1 CPU), 5.4 msec. (24 threads, 2 CPUs)

## Reference
[1] Ryota Yasudo et al. ``Designing High-Performance Interconnection Networks with Host-Switch Graphs``, in IEEE Transactions on Parallel and Distributed Systems, vol. 30, no. 2, pp. 315-330, 1 Feb. 2019, doi: 10.1109/TPDS.2018.2864286.

