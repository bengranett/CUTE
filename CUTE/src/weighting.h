

#ifndef _WEIGHTING_H_
#define _WEIGHTING_H_


typedef struct {
  double xmin;
  double step;
  int n;
  double extrap_low;
  double extrap_high;
  double * func;
} Func;

int read_pair_weights_file(char *filename, Func * func);
double interpolate_weight(double x, Func * func);
void free_weight_func(Func * func);


#endif //_WEIGHTING_H_
