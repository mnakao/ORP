#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include "orp.h"
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX_FILENAME_LENGTH 256
#define NOT_DEFINED -1

static void print_help(char *argv)
{
  fprintf(stderr, "%s -H hosts -S switches -R radix [-o output]\n", argv);
  fprintf(stderr, "  -H : Number of hosts\n");
  fprintf(stderr, "  -S : Number of switches\n");
  fprintf(stderr, "  -R : Radix\n");
  fprintf(stderr, "  -o : Output file\n");
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix, char **outfname)
{
  int result;
  while((result = getopt(argc,argv,"H:S:R:o:"))!=-1){
    switch(result){
    case 'H':
      *hosts = atoi(optarg);
      if(*hosts <= 3)
        ERROR("-H value > 3\n");
      break;
    case 'S':
      *switches = atoi(optarg);
      if(*switches <= 2)
        ERROR("-S value > 2\n");
      break;
    case 'R':
      *radix = atoi(optarg);
      if(*radix <= 3)
        ERROR("-R value > 3\n");
      break;
    case 'o':
      if(strlen(optarg) > MAX_FILENAME_LENGTH)
        ERROR("Output filename is long (%s).\n", optarg);
      *outfname = malloc(MAX_FILENAME_LENGTH);
      strcpy(*outfname, optarg);
      break;
    default:
      print_help(argv[0]);
    }
  }
}

static void generate_graph(const int hosts, const int switches, const int radix, const int lines,
                           int h_degree[switches], int s_degree[switches], int edge[lines][2], int adj[switches][radix])
{
  for(int i=0;i<switches;i++)
    h_degree[i] = hosts/switches;
  for(int i=0;i<hosts%switches;i++)
    h_degree[i]++;

  int tmp_lines = 0;
  for(int i=0;i<lines;i++)
    edge[i][0] = edge[i][1] = NOT_DEFINED;
  for(int i=0;i<switches;i++)
    s_degree[i] = 0;
  for(int i=0;i<switches;i++){
    for(int j=i+1;j<switches;j++){
      edge[tmp_lines][0] = i+hosts;
      edge[tmp_lines][1] = j+hosts;
      s_degree[i]++;
      s_degree[j]++;
      tmp_lines++;
    }
  }
  
  ORP_Conv_edge2adjacency(hosts, switches, radix, lines, edge, adj);
  ORP_Conv_adjacency2edge(hosts, switches, radix, h_degree, s_degree, adj, edge);
}

static void set_aspl(const int hosts, const int switches, const int h_degree[switches],
                     int *diameter, long *sum, double *ASPL)
{
  *sum = 0;
  for(int i=0;i<switches;i++){
    *sum += (long)h_degree[i] * (h_degree[i] - 1);
    for(int j=i+1;j<switches;j++)
      *sum += (long)3 * h_degree[i] * h_degree[j];
  }

  *ASPL = *sum / (double)(((long)hosts*(hosts-1))/2);
  *diameter = 3;
}

int main(int argc, char *argv[])
{
  char *outfname = NULL;
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED;
  int diameter, best_diameter, low_diameter;
  long sum, best_sum;
  double ASPL, best_ASPL, low_ASPL;

  set_args(argc, argv, &hosts, &switches, &radix, &outfname);
  
  if(hosts == NOT_DEFINED || radix == NOT_DEFINED || switches == NOT_DEFINED)
    print_help(argv[0]);

  int lines = hosts + (switches * (switches-1))/2;
  if(lines > (switches * radix - hosts)/2 + hosts)
    ERROR("This is not complete graph.\n");
  
  int *h_degree     = malloc(sizeof(int) * switches);
  int *s_degree     = malloc(sizeof(int) * switches);
  int (*edge)[2]    = malloc(sizeof(int) * lines * 2);
  int (*adj)[radix] = malloc(sizeof(int) * switches * radix);

  generate_graph(hosts, switches, radix, lines, h_degree, s_degree, edge, adj);

  ORP_Set_lbounds(hosts, radix, &low_diameter, &low_ASPL);
  int max = radix - s_degree[switches-1] + 1;
  for(h_degree[0]=0;h_degree[0]<max;h_degree[0]++){
    for(h_degree[1]=h_degree[0];h_degree[1]<max;h_degree[1]++){
      for(h_degree[2]=h_degree[1];h_degree[2]<max;h_degree[2]++){
        for(h_degree[3]=h_degree[2];h_degree[3]<max;h_degree[3]++){
          for(h_degree[4]=h_degree[3];h_degree[4]<max;h_degree[4]++){
            for(h_degree[5]=h_degree[4];h_degree[5]<max;h_degree[5]++){
              for(h_degree[6]=h_degree[5];h_degree[6]<max;h_degree[6]++){
                for(h_degree[7]=h_degree[6];h_degree[7]<max;h_degree[7]++){
                  for(h_degree[8]=h_degree[7];h_degree[8]<max;h_degree[8]++){
                  int s = 0;
                  for(int i=0;i<9;i++) s += h_degree[i];
                  if(s == hosts){
                    printf("%d %d %d %d %d %d %d %d\n",
                           h_degree[0], h_degree[1], h_degree[2], h_degree[3], h_degree[4], h_degree[5], h_degree[6], h_degree[7]);
                    set_aspl(hosts, switches, h_degree, &best_diameter, &best_sum, &best_ASPL);
                    printf("ASPL Gap        = %.10f (%.10f - %.10f)\n", best_ASPL - low_ASPL, best_ASPL, low_ASPL);}
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  if(outfname)
    ORP_Write_edge(hosts, switches, radix, lines, edge, outfname);

  return 0;
}
