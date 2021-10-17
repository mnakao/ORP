#include <stdio.h>
#include <stdlib.h>
#include "orp.h"
#define NOT_VISITED -1
#define VISITED      1
int aa = 0;

static int top_down_step(const int level, const int switches, const int num_frontier, const int radix,
                         int *s_degree, int (*adjacency)[radix], const int* restrict frontier,
                         int* restrict next, int* restrict distance)
{
  int count = 0;
  for(int i=0;i<num_frontier;i++){
    int v = frontier[i];
    for(int j=0;j<s_degree[v];j++){
      int n = adjacency[v][j];
      if(distance[n] == NOT_VISITED){
        distance[n] = level;
        next[count++] = n;
      }
    }
  }
  
  return count;
}

static void simple_bfs(const int hosts, const int switches, const int radix, int *h_degree, int* s_degree, int (*adjacency)[radix])
{
  int *frontier = malloc(sizeof(int)  * switches);
  int *next     = malloc(sizeof(int)  * switches);
  int *distance = malloc(sizeof(int)  * switches);

  for(int s=0;s<switches;s++){
    if(h_degree[s] == 0)
      continue;
    
    int num_frontier = 1, level = 1;
    frontier[0] = s;
    for(int i=0;i<switches;i++) distance[i] = NOT_VISITED;
    distance[s] = VISITED;

    while(1){
      num_frontier = top_down_step(level++, switches, num_frontier, radix, s_degree,
                                   adjacency, frontier, next, distance);
      if(num_frontier == 0) break;
      
      int *tmp = frontier;
      frontier = next;
      next     = tmp;
    }

    for(int i=0;i<switches;i++)
      if(distance[i] == 3 && h_degree[i] != 0)
        aa++;
  }
  printf("aa = %d\n", aa);
  free(frontier);
  free(next);
  free(distance);
}

int main(int argc, char **argv)
{
  int hosts, switches, radix, lines, diameter = -1, low_diameter = -1;
  long sum = -1;
  double ASPL = -1, low_ASPL = -1;

  if(argc != 2){
    fprintf(stderr, "Please add an edge file. \"%s xxx.edges\"\n", argv[0]);
    exit(1);
  }
    
  ORP_Read_property(argv[1], &hosts, &switches, &radix, &lines);
  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ORP_Read_edge(argv[1], edge);
  
  int *s_degree = malloc(sizeof(int) * switches);
  int *h_degree = malloc(sizeof(int) * switches);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix); // int adjacency[switches][radix];
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);

  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  ORP_Finalize_aspl();

  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  printf("Hosts = %d, Switches %d, Radix = %d\n", hosts, switches, radix);
  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);

  simple_bfs(hosts, switches, radix, h_degree, s_degree, adjacency);
    
  free(adjacency);
  free(h_degree);
  free(s_degree);
  free(edge);
  return 0;
}
