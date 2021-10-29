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

  int n1, n2, lines = 0;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF)
    if(n1 >= hosts && n2 >= hosts)
      lines++;

  int (*edge)[2] = malloc(sizeof(int) * lines * 2); // int edge[lines][2];
  int *h_degree  = malloc(sizeof(int) * switches);
  for(int i=0;i<switches;i++) h_degree[i] = 0;

  rewind(fp);
  fscanf(fp, "%d %d %d", &hosts, &switches, &radix); // Skip the first line

  int i = 0;
  while(fscanf(fp, "%d %d", &n1, &n2) != EOF){
    if(n1 >= hosts && n2 >= hosts){
      edge[i][0] = n1 - hosts;
      edge[i][1] = n2 - hosts;
      i++;
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
  printf("%d %d 10\n", lines, switches);
  for(int i=0;i<lines;i++)
    printf("%d %d\n", edge[i][0]+1, edge[i][1]+1); // Vertex number in hmetis is 1-origin.

  for(int i=0;i<switches;i++)
    printf("%d\n", h_degree[i]);

  free(edge);
  free(h_degree);
  return 0;
}
