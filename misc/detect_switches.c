#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

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

int main()
{
  int hosts = 65536, radix = 64;
  
  int s = 3;
  double prev = DBL_MAX;
  while(1){
    if(s*radix-2*(s-1) >= hosts){
      double tmp = continuous_moore_bound(hosts, s, radix);
      printf(" %d %f\n", s, tmp);
      if(prev <= tmp)
        break;
      prev = tmp;
    }
    s++;
  }

  printf("Host = %d, Radix = %d, Switch = %d\n", hosts, radix, s-1);
  return 0;
}
