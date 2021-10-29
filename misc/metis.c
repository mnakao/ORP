#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  if(argc != 2){
    fprintf(stderr, "%s filename\n", argv[0]);
    exit(1);
  }

  FILE *fp = NULL;
  if((fp = fopen(argv[1], "r")) == NULL){
    fprintf(stderr, "File not found\n");
    exit(1);
  }

  int hosts = -1, switches = -1, radix = -1;
  fscanf(fp, "%d %d %d", &hosts, &switches, &radix);
  
  int (*adj)[radix] = malloc(sizeof(int) * switches * radix); // int adj[switches][radix];
  int *h_degree     = malloc(sizeof(int) * switches);
  int *s_degree     = malloc(sizeof(int) * switches);
  for(int i=0;i<switches;i++){
    h_degree[i] = s_degree[i] = 0;
    for(int j=0;j<radix;j++)
      adj[i][j] = 0;
  }

  int n1, n2, lines = 0;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF){
    if(n1 >= hosts && n2 >= hosts){
      adj[n1-hosts][s_degree[n1-hosts]] = n2-hosts;
      adj[n2-hosts][s_degree[n2-hosts]] = n1-hosts;
      s_degree[n1-hosts]++;
      s_degree[n2-hosts]++;
      lines++;
    }
    else{
      if(n1 >= hosts)
        h_degree[n1-hosts]++;
      else
        h_degree[n2-hosts]++;
    }
  }
  fclose(fp);

  // Output
  printf("%d %d 010 1\n", switches, lines);
  for(int i=0;i<switches;i++){
    printf("%d ", h_degree[i]);
    for(int j=0;j<s_degree[i];j++)
      printf("%d ", adj[i][j]+1);
    printf("\n");
  }

  
  free(adj);
  free(h_degree);
  return 0;
}
