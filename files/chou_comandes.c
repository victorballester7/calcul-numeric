#include <stdio.h>
#include "rk78.h"

int main (int argc, char *argv[]) {
   double a, b;
   int n;
   if (
       argc!=4
       || sscanf(argv[1], "%lf", &a)!=1
       || sscanf(argv[2], "%lf", &b)!=1
       || sscanf(argv[3], "%d", &n)!=1
   ) {
      printf("./chou_comandes a b n\n");
      return -1;
   }
   printf("a %lf b %lf n %d\n", a, b, n);
   return 0;
}
