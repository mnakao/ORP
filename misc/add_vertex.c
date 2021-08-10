#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "orp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX_FILENAME_LENGTH 256
#define NOT_DEFINED -1

static void print_help(char *argv)
{
  fprintf(stderr, "%s -f input -a added_vertices [-s seed] [-n calcs]\n", argv);
  fprintf(stderr, "  -f : Input file\n");
  fprintf(stderr, "  -a : Added vertices\n");
  fprintf(stderr, "  -s : Random seed\n");
  fprintf(stderr, "  -n : Number of calculations\n");
  exit(1);
}

static void set_args(const int argc, char **argv, char **infname, int *added_vertices, int *seed, int *ncalcs)
{
  int result;
  while((result = getopt(argc,argv,"f:a:s:n"))!=-1){
    switch(result){
    case 'f':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Input filename is long (%s).\n", optarg);
      *infname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*infname, optarg);
      break;
    case 'a':
      *added_vertices = atoi(optarg);
      if(*added_vertices <= 0)
        ERROR("-a value > 0\n");
      break;
    case 's':
      *seed = atoi(optarg);
      if(*seed < 0)
        ERROR("-s value >= 0\n");
      break;
    case 'n':
      *ncalcs = atoi(optarg);
      if(*ncalcs < 0)
        ERROR("-n value >= 0\n");
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

static int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

int main(int argc, char *argv[])
{
  char *infname = NULL;
  int added_vertices = NOT_DEFINED, seed = 0, ncalcs = 100000;  
  int hosts, switches, radix, lines, diameter, low_diameter;
  long sum;
  double ASPL, low_ASPL;

  set_args(argc, argv, &infname, &added_vertices, &seed, &ncalcs);
  if(infname == NULL || added_vertices == NOT_DEFINED) print_help(argv[0]);
  else if(added_vertices%2 != 0) ERROR("Added vertices must be even number.\n");

  ORP_Srand(seed);
  ORP_Read_property(infname, &hosts, &switches, &radix, &lines);
  int *h_degree  = malloc(sizeof(int) * switches);
  int *s_degree  = malloc(sizeof(int) * switches);
  int (*edge)[2] = malloc(sizeof(int) * lines * 2);
  int (*adjacency)[radix] = malloc(sizeof(int) * switches * radix);
  
  ORP_Read_edge(infname, edge);
  ORP_Set_degrees(hosts, switches, lines, edge, h_degree, s_degree);
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
  printf("\n");

  // Delete edges
  for(int i=0;i<added_vertices/2;i++){
    int u   = get_random(switches);
    int u_d = get_random(s_degree[u]);
    int v   = adjacency[u][u_d];
    int v_d = search_index(v, u, u_d, s_degree, radix, adjacency);
    adjacency[u][u_d]           = adjacency[u][s_degree[u]-1];
    adjacency[u][s_degree[u]-1] = NOT_DEFINED;
    adjacency[v][v_d]           = adjacency[v][s_degree[v]-1];
    adjacency[v][s_degree[v]-1] = NOT_DEFINED;
    s_degree[u]--; h_degree[u]++;
    s_degree[v]--; h_degree[v]++;
  }

  int new_lines = lines + added_vertices/2;
  int (*new_edge)[2] = malloc(sizeof(int) * new_lines * 2);
  ORP_Conv_adjacency2edge(hosts+added_vertices, switches, radix, h_degree, s_degree, adjacency, new_edge);

  printf("%d %d %d\n", hosts+added_vertices, switches, radix);
  for(int i=0;i<new_lines;i++)
    printf("%d %d\n", new_edge[i][0], new_edge[i][1]);
  
  return 0;
}
