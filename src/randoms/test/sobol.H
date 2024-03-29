
//
//  32 bit integer routines.
//
int i4_bit_hi1 ( int n );
int i4_bit_lo0 ( int n );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4_sobol ( int dim_num, int *seed, float quasi[ ] );
int i4_uniform ( int b, int c, int *seed );
unsigned int i4_xor ( unsigned int i, unsigned int j );
//
//  64 bit integer routines.
//
int i8_bit_hi1 ( long long int n );
int i8_bit_lo0 ( long long int n );
long long int i8_max ( long long int i1, long long int i2 );
long long int i8_min ( long long int i1, long long int i2 );
void i8_sobol ( int dim_num, long long int *seed, double quasi[ ] );
long long int i8_uniform ( long long int b, long long int c, int *seed );
unsigned long long int i8_xor ( unsigned long long int i, 
  unsigned long long int j );
//
//  32 bit real rotuines.
//
float r4_abs ( float x );
int r4_nint ( float x );
float r4_uniform_01 ( int *seed );
//
//  64 bit real routines.
//
double r8_abs ( double x );
int r8_nint ( double x );
double r8_uniform_01 ( int *seed );
//
//  Utilities.
//
void timestamp ( void );
