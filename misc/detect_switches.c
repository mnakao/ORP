#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#define NOT_DEFINED -1
#define ERROR(...) do{fprintf(stderr,__VA_ARGS__); exit(1);}while(0)

static void print_help(char *argv)
{
  fprintf(stderr, "%s -H hosts -R radix [-s start] [-e end]\n", argv);
  fprintf(stderr, "  -H : Number of hosts\n");
  fprintf(stderr, "  -R : Radix\n");
  fprintf(stderr, "  -s : Start number of switches\n");
  fprintf(stderr, "  -e : End number of switches (If not set, stop when the optimum number of switches is found)\n");
  exit(1);
}

static void set_args(const int argc, char **argv, int *hosts, int *radix, int *start, int *end)
{
  int result;
  while((result = getopt(argc,argv,"H:R:s:e:"))!=-1){
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
    case 's':
      *start = atoi(optarg);
      if(*start < 3)
        ERROR("-s value > 3\n");
      break;
    case 'e':
      *end = atoi(optarg);
      break;
    default:
      print_help(argv[0]);
    }
  }
}

static double moore_bound(const double nodes, const double degree)
{
  if(degree + 1 >= nodes)
    return 1;
  
  double diam = -1, n = 1, r = 1, aspl = 0.0, prev_tmp;
  while(1){
    double tmp = n + degree * pow(degree-1, r-1);
    if(tmp >= nodes || (r > 1 && prev_tmp == tmp))
      break;
    
    n = tmp;
    aspl += r * degree * pow(degree-1, r-1);
    diam = r++;
    prev_tmp = tmp;
  }

  diam++;
  aspl += diam * (nodes - n);
  aspl /= (nodes - 1);
 
  return aspl;
}

static double continuous_moore_bound(const int hosts, const int switches, const int radix)
{
  double h = hosts;
  double s = switches;
  double r = radix;
  return moore_bound(s,r-h/s)*(s*h-h)/(s*h-s)+2;
}

int main(int argc, char *argv[])
{
  bool first = true;
  int hosts = NOT_DEFINED, radix = NOT_DEFINED, start = 3, end = NOT_DEFINED, s = NOT_DEFINED;
  double prev = DBL_MAX;
  set_args(argc, argv, &hosts, &radix, &start, &end);
  if(hosts == NOT_DEFINED || radix == NOT_DEFINED)
    print_help(argv[0]);

  while(1){
    if(start*radix-2*(start-1) >= hosts){
      double tmp = continuous_moore_bound(hosts, start, radix);
      printf("%d\t%f\n", start, tmp);
      if(prev <= tmp){
        if(first){
          s = start;
          first = false;
        }
        if(end == NOT_DEFINED)
          break;
      }
      
      prev = tmp;
    }
    if(start == end) break;
    start++;
  }

  if(s != NOT_DEFINED)
    printf("Host = %d, Radix = %d, Switch = %d\n", hosts, radix, s-1);
  
  return 0;
}
