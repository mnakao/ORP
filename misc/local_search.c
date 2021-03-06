#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "orp.h"
int ORP_top_down_step(const int level, const int num_frontier, const int* restrict adjacency,
                      const int switches, const int radix, const int* restrict ports,
                      int* restrict frontier, int* restrict next, int* restrict distance);
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX_FILENAME_LENGTH 256
#define NOT_DEFINED -1

static void print_help(char *argv)
{
  fprintf(stderr, "%s -f input -o output [-G symmetries]\n", argv);
  fprintf(stderr, "  -f : Input file\n");
  fprintf(stderr, "  -o : Output file\n");
  fprintf(stderr, "  -G : Numer of symmetries (Default: 1)\n");
  exit(1);
}

static void set_args(const int argc, char **argv, char **infname, char **outfname, int *symmetries)
{
  int result;
  while((result = getopt(argc,argv,"f:o:G:"))!=-1){
    switch(result){
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s).\n", optarg);
      *infname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*infname, optarg);
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s).\n", optarg);
      *outfname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*outfname, optarg);
      break;
    case 'G':
      *symmetries = atoi(optarg);
      if(*symmetries < 1)
        ERROR("-G value >= 1\n");
      break;
    default:
      print_help(argv[0]);
    }
  }
}

static int search_index(const int v, const int target, const int exclusion, const int *s_degree,
                        const int radix, const int (*adjacency)[radix])
{
  if(v == target){
    for(int i=0;i<s_degree[v];i++){
      if(adjacency[v][i] == target && i != exclusion){
        return i;
      }
    }
  }
  else{
    for(int i=0;i<s_degree[v];i++){
      if(adjacency[v][i] == target){
        return i;
      }
    }
  }

  ERROR("Something Wrong (id=2)\n");
  return -1; // dummy
}

static void set_distance(const int* restrict s_degree, const int switches,
                         const int radix, const void* adjacency, int (*distance)[switches])
{
  int *frontier     = malloc(sizeof(int) * switches);
  int *tmp_distance = malloc(sizeof(int) * switches);
  int *next         = malloc(sizeof(int) * switches);

  for(int i=0;i<switches;i++)
    for(int j=0;j<switches;j++)
      distance[i][j] = NOT_DEFINED;
  
  for(int s=0;s<switches;s++){
    int num_frontier = 1, level = 1;
    for(int i=0;i<switches;i++)
      tmp_distance[i] = NOT_DEFINED;

    frontier[0] = s;
    tmp_distance[s] = level;

    while(1){
      num_frontier = ORP_top_down_step(level++, num_frontier, adjacency, switches, radix, s_degree,
                                       frontier, next, tmp_distance);
      if(num_frontier == 0) break;

      int *tmp = frontier;
      frontier = next;
      next     = tmp;
    }
    memcpy(distance[s], tmp_distance, sizeof(int) * switches);
  }
  
  free(frontier);
  free(tmp_distance);
  free(next);
}

static double calc_aspl(const int hosts, const int switches, const int h_degree[switches], const int distance[switches][switches])
{
  long sum = 0;
  for(int i=0;i<switches;i++){
    if(h_degree[i] == 0) continue;
    for(int j=i+1;j<switches;j++){
      if(h_degree[j] != 0)
        sum += (long)(distance[i][j] + 2) * h_degree[i] * h_degree[j];
    }
  }
  for(int i=0;i<switches;i++)
    sum += (long)h_degree[i] * (h_degree[i] - 1);

  return sum/(double)(((long)hosts*(hosts-1))/2);
}

int main(int argc, char **argv)
{
  char *infname = NULL, *outfname = NULL;
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, lines = NOT_DEFINED;
  int diameter = NOT_DEFINED, low_diameter = NOT_DEFINED, new_diameter = NOT_DEFINED, tmp_diameter = NOT_DEFINED;
  int symmetries = 1;
  long sum = NOT_DEFINED, new_sum = NOT_DEFINED, tmp_sum = NOT_DEFINED;
  double ASPL = NOT_DEFINED, low_ASPL = NOT_DEFINED, new_ASPL = NOT_DEFINED, tmp_ASPL = NOT_DEFINED;
  bool flag = false;

  set_args(argc, argv, &infname, &outfname, &symmetries);
  if(infname == NULL || outfname == NULL)
    ERROR("./%s -f input -o output\n", argv[0]);
  
  ORP_Read_property(infname, &hosts, &switches, &radix, &lines);
  printf("Hosts = %d, Switches = %d, Radix = %d\n", hosts, switches, radix);

  int (*edge)[2] = malloc(sizeof(int)*lines*2); // int edge[lines][2];
  ORP_Read_edge(infname, edge);
  
  int *s_degree = malloc(sizeof(int) * switches);
  int *h_degree = malloc(sizeof(int) * switches);
  int *best_s_degree = malloc(sizeof(int) * switches);
  int *best_h_degree = malloc(sizeof(int) * switches);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
  
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix); // int adjacency[switches][radix];
  int (*best_adjacency)[radix] = malloc(sizeof(int) * switches * radix);
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adjacency);

  ORP_Init_aspl(hosts, switches, radix);
  ORP_Set_aspl(h_degree, s_degree, adjacency, &diameter, &sum, &ASPL);
  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);

  printf("Diameter     = %d\n", diameter);
  printf("Diameter Gap = %d (%d - %d)\n", diameter - low_diameter, diameter, low_diameter);
  printf("ASPL         = %.10f (%ld/%.0f)\n", ASPL, sum, (double)hosts*(hosts-1)/2);
  printf("ASPL Gap     = %.10f (%.10f - %.10f)\n", ASPL - low_ASPL, ASPL, low_ASPL);
  printf("--\n");
  
  double best_ASPL = ASPL;
  int (*distance)[switches] = malloc(sizeof(int) * switches * switches);
  for(int u=0;u<switches/symmetries;u++){
    if(s_degree[u] == 1) continue;
    for(int u_d=0;u_d<s_degree[u];u_d++){
      int v = adjacency[u][u_d];
      if(s_degree[v] == 1) continue; // (u == v && s_degree[u] == 2) cannot be happen.
      int v_d = search_index(v, u, u_d, s_degree, radix, adjacency);
      adjacency[u][u_d] = adjacency[u][s_degree[u]-1]; s_degree[u]--;
      adjacency[v][v_d] = adjacency[v][s_degree[v]-1]; s_degree[v]--;
      ORP_Set_aspl(h_degree, s_degree, adjacency, &new_diameter, &new_sum, &new_ASPL);
      if(new_ASPL == ASPL && new_diameter == diameter){
        printf("Delete : %d -- %d\n", u+hosts, v+hosts);
        set_distance(s_degree, switches, radix, adjacency, distance);
        for(int i=0;i<switches;i++){
          if(h_degree[i] == 0) continue;
          for(int j=i;j<switches;j++){
            if(h_degree[j] == 0 || (i == j && h_degree[i] == 1)) continue;
            h_degree[i]--; h_degree[j]--; h_degree[u]++; h_degree[v]++;
            double tmp_ASPL = calc_aspl(hosts, switches, h_degree, distance);
            if(tmp_ASPL < best_ASPL){
              flag = true;
              best_ASPL = tmp_ASPL;
              printf("Found a graph with smaller h-ASPL (%f)\n", tmp_ASPL);
              memcpy(best_h_degree, h_degree, sizeof(int) * switches);
              memcpy(best_s_degree, s_degree, sizeof(int) * switches);
              memcpy(best_adjacency, adjacency, sizeof(int) * switches * radix);
            }
            h_degree[i]++; h_degree[j]++; h_degree[u]--; h_degree[v]--;
          }
        }
      }
      adjacency[u][u_d] = v; s_degree[u]++;
      adjacency[v][v_d] = u; s_degree[v]++;
    }
  }
  
  ORP_Finalize_aspl();
  
  if(!flag)
    printf("Not found.\n");
  else{
    int (*best_edge)[2] = malloc(sizeof(int)*(lines-1)*2);
    ORP_Conv_adjacency2edge(hosts, switches, radix, best_h_degree, best_s_degree, best_adjacency, best_edge);
    ORP_Write_edge(hosts, switches, radix, lines-1, best_edge, outfname);
    free(best_edge);
    printf("Created a new graph to %s\n", outfname);
  }
  
  free(adjacency); free(best_adjacency); free(h_degree); free(best_h_degree);
  free(s_degree);  free(best_s_degree);  free(edge);     free(distance);

  return 0;
}
