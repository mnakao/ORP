#include <stdio.h>
#include <stdlib.h>
#include "orp.h"

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
  int (*edge)[2] = malloc(sizeof(int)*lines*2);
  ORP_Read_edge(argv[1], edge);
  
  for(int i=0;i<lines;i++)
    if(edge[i][0] >= hosts && edge[i][1] >= hosts)
      printf("%d %d\n", edge[i][0]-hosts, edge[i][1]-hosts);

  return 0;
}
