#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "orp.h"

static int search_max(const int switches, const int *s_degree)
{
  int max = s_degree[0];
  for(int i=1;i<switches;i++){
    if(s_degree[i] > max)
      max = s_degree[i];
  }
  return max;
}

static int search_min(const int switches, const int *s_degree)
{
  int min = s_degree[0];
  for(int i=1;i<switches;i++){
    if(s_degree[i] < min)
      min = s_degree[i];
  }
  return min;
}

int main(int argc, char *argv[])
{
  if(argc != 2){
    fprintf(stderr, "%s <edgefile>\n", argv[0]);
    exit(1);
  }

  int hosts, switches, radix, lines;
  ORP_Read_property(argv[1], &hosts, &switches, &radix, &lines);
  printf("(h, r, s) = (%d, %d, %d)\n", hosts, radix, switches);

  int *h_degree  = malloc(sizeof(int) * switches);
  int *s_degree  = malloc(sizeof(int) * switches);
  int (*edge)[2] = malloc(sizeof(int) * lines * 2);
  ORP_Read_edge(argv[1], edge);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);

  int max      = search_max(switches, s_degree);
  int min      = search_min(switches, s_degree);
  int elements = max - min + 1;
  int *hist    = malloc(elements * sizeof(int));
  for(int i=0;i<elements;i++) hist[i] = 0;

  for(int i=0;i<switches;i++)
    hist[s_degree[i]-min]++;

  for(int i=min;i<max+1;i++)
    printf("%d : %d\n", i, hist[i-min]);
  
  return 0;
}
