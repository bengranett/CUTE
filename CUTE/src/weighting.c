
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "define.h"
#include "common.h"
#include "weighting.h"


int read_pair_weights_file(char *filename, Func * func)
{
    /* Read the weights from an ascii file.
    Column 1: log(theta) - log of angular separation in radians.
    Column 2: pairwise weight to multiply (inverse of completeness)
    */
    FILE *fd;
    int line_count;
    int size;
    int sr;
    int j;
    double x;
    double y;
    double xprev;
    char line[1024];

    print_info("*** Reading pair weights from file %s\n\n", filename);

    fd = fopen(filename, "r");

    if(fd == NULL) error_open_file(filename);

    line_count = 0;
    size = 0;

    while((fgets(line, sizeof(line), fd))!=NULL)
    {
        line_count ++;
        sr = sscanf(line, "%lf %lf", &x, &y);
        if (sr != 2) {
            error_read_line(filename, line_count);
            continue;
        }
        size ++;
    }


    rewind(fd);

    func->n = size;
    func->func = (double *) my_malloc(size*sizeof(double));

    xprev = 0;

    j = 0;
    while((fgets(line, sizeof(line), fd))!=NULL)
    {
        sr = sscanf(line, "%lf %lf", &x, &y);
        if (sr != 2) continue;
        func->func[j] = y;

        x = log10(x*DTORAD);

        if (j == 0) func->xmin = x;
        if (j == 1) func->step = x - xprev;
        j++;
        xprev = x;
    }

    print_info("*** weighting function: xmin:%lf deg, xmax:%lf deg, size:%i\n",
        pow(10,func->xmin)/DTORAD,
        pow(10,func->xmin + func->step*(func->n-1))/DTORAD,
        func->n);

    return 0;
}



double interpolate_weight(double x, Func * func)
{
    /* Compute the function at a value by linear interpolation.
    The function is given at linearly sampled points.

    Parameters
    ----------
    x : double
    y : Func
     Function to evaluate.

    Returns
    -------
    double : interpolated value
    */
    int i;
    double f,y;

    f = (x - func->xmin) / func->step;
    i = (int) f;
    if (i < 0) return func->extrap_low;
    if (i >= func->n-1) return func->extrap_high;
    
    f -= i;  // fractional part

    y = (func->func[i+1] - func->func[i]) * f + func->func[i];

    return y;
}

void free_weight_func(Func * func){
    free(func->func);
}