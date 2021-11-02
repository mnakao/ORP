#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#define PERFORMANCE "500Gf"
#define BANDWIDTH   "100Gbps"
#define NIC_LATENCY "600ns"
#define SW_LATENCY  "100ns"
#define MAX_FILENAME_LENGTH (256)
#define NOT_DEFINED         (-1)
#define INF                 (1000000000)
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
static int _n = 0;

void print_header()
{
  printf("<?xml version='1.0'?>\n");
  printf("<!DOCTYPE platform SYSTEM \"https://simgrid.org/simgrid.dtd\">\n");
  printf("<platform version=\"4.1\">\n");
  printf("  <zone id=\"AS0\" routing=\"Full\">\n");
}

void print_footer()
{
  printf("  </zone>\n");
  printf("</platform>\n");
}

void print_hosts(const int hosts)
{
  for(int i=0;i<hosts;i++)
    printf("    <host id=\"host%d\" speed=\"%s\"/>\n", i, PERFORMANCE);
}

void print_lines(const int hosts, const int lines, const int edge[lines][2])
{
  for(int i=0;i<lines;i++){
    int n = edge[i][0];
    int m = edge[i][1];
    if(n > m){
      if(n >= hosts && m >= hosts)
        printf("    <link bandwidth=\"%s\" latency=\"%s\" id=\"link%dto%d\"/>\n", BANDWIDTH, SW_LATENCY, m, n);
      else
        printf("    <link bandwidth=\"%s\" latency=\"%s\" id=\"link%dto%d\"/>\n", BANDWIDTH, NIC_LATENCY, m, n);
    }
    else{
      if(n >= hosts && m >= hosts)
        printf("    <link bandwidth=\"%s\" latency=\"%s\" id=\"link%dto%d\"/>\n", BANDWIDTH, SW_LATENCY, n, m);
      else
        printf("    <link bandwidth=\"%s\" latency=\"%s\" id=\"link%dto%d\"/>\n", BANDWIDTH, NIC_LATENCY, n, m);
    }
  }
}

// https://algo-logic.info/warshall-floyd/#toc_id_2_1
void print_routes(const int hosts, const int switches, const int lines, const int edge[lines][2])
{
  int vertices = hosts + switches;
  int (*dist)[vertices] = malloc(sizeof(int) * vertices * vertices);
  int (*prev)[vertices] = malloc(sizeof(int) * vertices * vertices);

  for(int i=0;i<vertices;i++)
    for(int j=0;j<vertices;j++)
      dist[i][j] = prev[i][j] = INF;

  for(int i=0;i<lines;i++){
    dist[edge[i][0]][edge[i][1]] = dist[edge[i][1]][edge[i][0]] = 1;
    prev[edge[i][0]][edge[i][1]] = edge[i][0];
    prev[edge[i][1]][edge[i][0]] = edge[i][1];
  }

  for(int k=0;k<vertices;k++){
    for(int i=0;i<vertices;i++){
      for(int j=0;j<vertices;j++){
        if(dist[i][j] > dist[i][k] + dist[k][j]){
          dist[i][j] = dist[i][k] + dist[k][j];
          prev[i][j] = prev[k][j];
        }
      }
    } 
  }
  
  for(int src=0;src<hosts;src++){
    for(int dst=src+1;dst<hosts;dst++){
      printf("    <route src=\"host%d\" dst=\"host%d\">", dst, src);
      for(int cur=dst; cur!=src; cur=prev[src][cur]){
        int s = MIN(cur, prev[src][cur]);
        int e = MAX(cur, prev[src][cur]);
        printf("<link_ctn id=\"link%dto%d\"/>", s, e);
      }
      printf("</route>\n");
    }
  }

  free(dist); free(prev);
}

void read_property(const char* fname, int* hosts, int* switches, int* radix, int *lines)
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  fscanf(fp, "%d %d %d", hosts, switches, radix);

  *lines = 0;
  int n1, n2;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF)
    (*lines)++;

  fclose(fp);
}

void read_edge(const char* fname, int (*edge)[2])
{
  FILE *fp;
  if((fp = fopen(fname, "r")) == NULL)
    ERROR("File not found\n");

  int hosts, switches, radix;
  fscanf(fp, "%d %d %d", &hosts, &switches, &radix); // Skip the first line

  int n1, n2, i = 0;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF){
    edge[i][0] = n1;
    edge[i][1] = n2;
    i++;
  }

  fclose(fp);
}

void check_args(const int argc, char* argv[])
{
  if(argc != 2)
    ERROR("%s filename\n", argv[0]);
}

void dfs(const int src, int *s_degree, const int radix, int (*adj)[radix], const int hosts,
         int *check, int *s_output)
{
  s_output[_n++] = src + hosts;
  check[src] = 1;
  for(int i=0;i<s_degree[src];i++)
    if(check[adj[src][i]] == NOT_DEFINED)
      dfs(adj[src][i], s_degree, radix, adj, hosts, check, s_output);
}

// hosts are sequentially connected to switches in depth-first order by backtracking.
void exchange_switch(const int hosts, const int switches, const int radix, const int lines, int (*edge)[2])
{
  int *check        = malloc(sizeof(int) * switches);
  int *s_output     = malloc(sizeof(int) * switches);
  int *h_output     = malloc(sizeof(int) * hosts);
  int *s_degree     = malloc(sizeof(int) * switches);
  int *h_degree     = malloc(sizeof(int) * switches);
  int (*adj)[radix] = malloc(sizeof(int) * switches * radix);
  
  for(int i=0;i<switches;i++){
    check[i] = s_output[i] = NOT_DEFINED;
    s_degree[i] = h_degree[i] = 0;
    for(int j=0;j<radix;j++)
      adj[i][j] = NOT_DEFINED;
  }
  for(int i=0;i<hosts;i++)
    h_output[i] = NOT_DEFINED;

  for(int i=0;i<lines;i++){
    int n1 = edge[i][0];
    int	n2 = edge[i][1];
    if(n1 >= hosts && n2 >= hosts){
      adj[n1-hosts][s_degree[n1-hosts]] = n2-hosts;
      adj[n2-hosts][s_degree[n2-hosts]] = n1-hosts;
      s_degree[n1-hosts]++;
      s_degree[n2-hosts]++;
    }
    else{
      if(n1 >= hosts)
        h_degree[n1-hosts]++;
      else
        h_degree[n2-hosts]++;
    }
  }

  //  for(int i=0;i<switches;i++){
  //    printf("%d : ", i+hosts);
  //    for(int j=0;j<s_degree[i];j++){
  //      printf("%3d", adj[i][j]+hosts);
  //    }
  //    printf("\n");
  //  }
  
  dfs(0, s_degree, radix, adj, hosts, check, s_output);

  int k = 0;
  for(int i=0;i<switches;i++)
    for(int j=0;j<h_degree[s_output[i]-hosts];j++)
      h_output[k++] = s_output[i];

  //  for(int i=0;i<switches;i++)
  //    printf("%d\n", s_output[i]);

  //  for(int i=0;i<hosts;i++)
  //    printf("%d\n", h_output[i]);

  k = 0;
  for(int i=0;i<lines;i++){
    if(edge[i][0] < hosts)      edge[i][1] = h_output[k++];
    else if(edge[i][1] < hosts) edge[i][0] = h_output[k++];
  }

  //  printf("%d %d %d\n", hosts, switches, radix);
  //  for(int i=0;i<lines;i++)
  //    printf("%d %d\n", edge[i][0], edge[i][1]);
  //  exit(0);

  free(check); free(s_output); free(h_output); free(s_degree); free(h_degree); free(adj);
}

int main(int argc, char *argv[])
{
  char fname[MAX_FILENAME_LENGTH];
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED, lines = NOT_DEFINED;
  
  check_args(argc, argv);
  strcpy(fname, argv[1]);
  read_property(fname, &hosts, &switches, &radix, &lines);
  fprintf(stderr, "h = %d s = %d r = %d edges = %d\n", hosts, switches, radix, lines);

  int (*edge)[2] = malloc(sizeof(int) * lines * 2); // int edge[lines][2];
  read_edge(fname, edge);

  exchange_switch(hosts, switches, radix, lines, edge);

  print_header();
  print_hosts(hosts);
  print_lines(hosts, lines, edge);
  print_routes(hosts, switches, lines, edge);
  print_footer();
  
  free(edge);
  return 0;
}
