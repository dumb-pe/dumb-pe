#include <math.h>

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

int
solve_cubic (double a, double b, double c,
             double *x0, double *x1, double *x2)
{
  double q = (a * a - 3 * b);
  double r = (2 * a * a * a - 9 * a * b + 27 * c);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  if (R2 < Q3)
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double ratio = sgnR * sqrt (R2 / Q3);
      double theta = acos (ratio);
      double norm = -2 * sqrt (Q);
      *x0 = norm * cos (theta / 3) - a / 3;
      *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;

      /* Sort *x0, *x1, *x2 into increasing order */

      if (*x0 > *x1)
        SWAP(*x0, *x1) ;

      if (*x1 > *x2)
        {
          SWAP(*x1, *x2) ;

          if (*x0 > *x1)
            SWAP(*x0, *x1) ;
        }

      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
      double B = Q / A ;
      *x0 = A + B - a / 3;
      return 1;
    }
}
