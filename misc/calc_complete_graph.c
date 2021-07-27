#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define NOT_DEFINED -1

static void print_help(char *argv)
{
  fprintf(stderr, "%s -H hosts -S switches -R radix\n", argv);
  fprintf(stderr, "  -H : Number of hosts\n");
  fprintf(stderr, "  -R : Radix\n");
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *radix)
{
  int result;
  while((result = getopt(argc,argv,"H:R:"))!=-1){
    switch(result){
    case 'H':
      *hosts = atoi(optarg);
      if(*hosts <= 3)
        ERROR("-H value > 3\n");
      break;
    case 'R':
      *radix = atoi(optarg);
      if(*radix <= 3)
        ERROR("-R value > 3\n");
      break;
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED;
  set_args(argc, argv, &hosts, &radix);
  if(hosts == NOT_DEFINED || radix == NOT_DEFINED)
    print_help(argv[0]);

  //  int lines = hosts + (switches * (switches-1))/2;
  //  if(lines > (switches * radix - hosts)/2 + hosts)
  //    printf("This is not a complete graph.\n");
  //  else
  //    printf("This is a complete graph.\n");

  // hosts + (switches * (switches-1))/2 < (switches * radix - hosts)/2 + hosts
  // -> switches * (switches-1) < switches * radix - hosts
  // -> switches^2 - (radix + 1) * switches + hosts < 0
  // switches = ((radix + 1) +- sqrt((radix)^2 - 4 * hosts)) / 2

  if(radix * radix - 4 * hosts < 0){
    printf("This graph is not a complete graph\n");
  }
  else{
    double r = sqrt(radix * radix - 4 * hosts);
    double min =  ((radix + 1) - r)/2.0;
    double max =  ((radix + 1) + r)/2.0;
    printf("%f < switches < %f\n", min, max);
    
  }
  return 0;
}
