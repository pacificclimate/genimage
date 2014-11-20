#include <math.h>
#include <assert.h>

enum OPERATOR { LT, GT };

const static int seasons[][3] = { {11, 0, 1}, {2, 3, 4}, {5, 6, 7}, {8, 9, 10} };
const static int annual[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

void precip_as_snow(const double* temp, const double* prec, double* pas, const double temp_missing, const double prec_missing, const double pas_missing, const int month, const int data_size);
void number_frost_free_days(const double* temp, double* nffd, const double temp_missing, const double nffd_missing, const int month, const int data_size);
void degree_days_pcic(const double* temp, double* dd, const double temp_missing, const double dd_missing, const int degree, const OPERATOR op, const int month, const int data_size);

