 #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#define PERFORMANCE "500Gf"
#define BANDWIDTH   "100Gbps"
#define NIC_LATENCY "600ns"
#define SW_LATENCY  "100ns"
#define MAX_FILENAME_LENGTH (256)
#define NOT_DEFINED         (-1)
#define INF                 (1000000000)
#define DELTA               (0.000001)
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)

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

// Ref. https://blue-9.hatenadiary.com/entry/2018/01/25/190000
//      https://marycore.jp/prog/c-lang/compare-floating-point-number/
void print_link_ctn(const int src, const int dst, const int vertices, const int lines,
		    const int edge[lines][2], double dist[vertices][vertices], int count[vertices][vertices])
{
  double cost[vertices], min;
  int target, via[vertices];
  bool used[vertices];
  
  for(int i=0;i<vertices;i++){
    cost[i] = INF;
    used[i] = false;
    via[i]  = NOT_DEFINED;
  }

  cost[src] = 0;
  while(1){
    min = INF;
    for(int i=0;i<vertices;i++){
      if(!used[i] && min > cost[i]){
        min = cost[i];
        target = i;
      }
    }

    if(target == dst) break;

    for(int n=0;n<vertices;n++){
      if(cost[n] > dist[target][n] + cost[target]){
        cost[n] = dist[target][n] + cost[target];
        via[n] = target;
      }
    }
    used[target] = true;
  }

  int node = dst;
  if(node > via[node]){
    if(DELTA != 0){
      double a = dist[via[node]][node] - (int)dist[via[node]][node] - ((double)1.0 - DELTA);
      double b = dist[node][via[node]] - (int)dist[node][via[node]] - ((double)1.0 - DELTA);
      if(fmax(fabs(a), fabs(b)) <= DBL_EPSILON)
        ERROR("Make Delta small\n");
    }
    
    count[via[node]][node]++;
    dist[via[node]][node] += DELTA;
    dist[node][via[node]] += DELTA;
    printf("<link_ctn id=\"link%dto%d\"/>", via[node], node);
  }
  else{
    if(DELTA != 0){
      double a = dist[via[node]][node] - (int)dist[via[node]][node] - ((double)1.0 - DELTA);
      double b = dist[node][via[node]] - (int)dist[node][via[node]] - ((double)1.0 - DELTA);
      if(fmax(fabs(a), fabs(b)) <= DBL_EPSILON)
        ERROR("Make Delta small\n");
    }
    
    count[node][via[node]]++;
    dist[via[node]][node] += DELTA;
    dist[node][via[node]] += DELTA;
    printf("<link_ctn id=\"link%dto%d\"/>", node, via[node]);
  }
  
  while(1){
    node = via[node];
    if(node == src) break;
    if(node > via[node]){
      if(DELTA != 0){
        double a = dist[via[node]][node] - (int)dist[via[node]][node] - ((double)1.0 - DELTA);
        double b = dist[node][via[node]] - (int)dist[node][via[node]] - ((double)1.0 - DELTA);
        if(fmax(fabs(a), fabs(b)) <= DBL_EPSILON)
          ERROR("Make Delta small\n");
      }
          
      count[via[node]][node]++;
      dist[via[node]][node] += DELTA;
      dist[node][via[node]] += DELTA;
      printf("<link_ctn id=\"link%dto%d\"/>", via[node], node);
    }
    else{
      if(DELTA != 0){
        double a = dist[via[node]][node] - (int)dist[via[node]][node] - ((double)1.0 - DELTA);
        double b = dist[node][via[node]] - (int)dist[node][via[node]] - ((double)1.0 - DELTA);
        if(fmax(fabs(a), fabs(b)) <= DBL_EPSILON)
          ERROR("Make Delta small\n");
      }
          
      count[node][via[node]]++;
      dist[via[node]][node] += DELTA;
      dist[node][via[node]] += DELTA;
      printf("<link_ctn id=\"link%dto%d\"/>", node, via[node]);
    }
  }
}

void print_routes(const int hosts, const int switches, const int lines, const int edge[lines][2])
{
  int vertices = hosts + switches;
  double (*dist)[vertices] = malloc(sizeof(double) * vertices * vertices); // double dist[vertices][vertices];
  
  for(int i=0;i<vertices;i++)
    for(int j=0;j<vertices;j++)
      dist[i][j] = INF;

  for(int i=0;i<lines;i++){
    int n = edge[i][0];
    int m = edge[i][1];
    dist[n][m] = dist[m][n] = 1;
  }

  int count[vertices][vertices];
  for(int i=0;i<vertices;i++)
    for(int j=0;j<vertices;j++)
      count[i][j] = 0;
  
  for(int src=0;src<hosts;src++){
    for(int dst=src+1;dst<hosts;dst++){
      printf("    <route src=\"host%d\" dst=\"host%d\">", dst, src);
      print_link_ctn(src, dst, vertices, lines, edge, dist, count);
      printf("</route>\n");
    }
  }

  //  for(int i=0;i<vertices;i++)
  //    for(int j=0;j<vertices;j++)
  //      if(dist[i][j] >= 2 && dist[i][j] != INF){
  //	printf("Error !\n");
  //	exit(1);
  //      }


  for(int i=0;i<hosts;i++)
    for(int j=0;j<hosts;j++)
      if(count[i][j] != 0)
  	printf("AAA %d\n", count[i][j]);


  free(dist);
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

int get_random(const int max)
{
  return (int)(rand()*((double)max)/(1.0+RAND_MAX));
}

void shuffle_edge(const int hosts, const int lines, int edge[lines][2])
{
  int array[hosts];
  for(int i=0;i<hosts;i++) array[i] = i;

  for(int i=0;i<hosts*10;i++){
    int r0 = get_random(hosts);
    int	r1 = get_random(hosts);
    int tmp   = array[r0];
    array[r0] = array[r1];
    array[r1] = tmp;
  }

  for(int i=0;i<lines;i++)
    for(int j=0;j<2;j++)
      if(edge[i][j] < hosts)
        edge[i][j] = array[edge[i][j]];
}

void check_args(const int argc, char* argv[])
{
  if(argc != 3)
    ERROR("%s filename seed\n", argv[0]);
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
  
  int seed = atoi(argv[2]);
  srand(seed);
  shuffle_edge(hosts, lines, edge);
  //  for(int i=0;i<lines;i++)
  //    printf("%d %d\n", edge[i][0], edge[i][1]);
  // exit(0);

  print_header();
  print_hosts(hosts);
  print_lines(hosts, lines, edge);
  print_routes(hosts, switches, lines, edge);
  print_footer();
  
  free(edge);
  return 0;
}
