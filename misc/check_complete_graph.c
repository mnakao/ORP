#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)
#define NOT_DEFINED -1

static void print_help(char *argv)
{
  fprintf(stderr, "%s -H hosts -S switches -R radix\n", argv);
  fprintf(stderr, "  -H : Number of hosts\n");
  fprintf(stderr, "  -S : Number of switches\n");
  fprintf(stderr, "  -R : Radix\n");
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *switches, int *radix)
{
  int result;
  while((result = getopt(argc,argv,"H:S:R:"))!=-1){
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
    default:
      print_help(argv[0]);
    }
  }
}

int main(int argc, char *argv[])
{
  int hosts = NOT_DEFINED, switches = NOT_DEFINED, radix = NOT_DEFINED;
  set_args(argc, argv, &hosts, &switches, &radix);
  if(hosts == NOT_DEFINED || radix == NOT_DEFINED || switches == NOT_DEFINED)
    print_help(argv[0]);

  int lines = hosts + (switches * (switches-1))/2;
  if(lines > (switches * radix - hosts)/2 + hosts)
    printf("This is not a complete graph.\n");
  else
    printf("This is a complete graph.\n");

  return 0;
}
