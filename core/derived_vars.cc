#include "derived_vars.h"

const static double dpm[] = { 31,      28.2425, 31, 30, 31, 30, 
			      31,      31,      30, 31, 30, 31 };

void precip_as_snow(const double* temp, const double* prec, double* pas, const double temp_missing, const double prec_missing, const double pas_missing, const int month, const int data_size) {
  const double x0[] = { -2.9901, -1.3948,  0.5473,  2.0928,   4.078,     0,     0,     0,  1.4927,  0.8099, -1.5627, -2.5909 };
  const double b[] = { -2.50353, -2.0004, -1.5719, -1.6527, -1.7428,     0,     0,     0, -2.8948, -1.6612, -2.4907, -2.2108 };
  const bool override[] = { false,   false,   false,   false,   false,  true,  true,  true,   false,   false,   false,   false };
  
  if(override[month])
    return;

  for(int j = 0; j < data_size; ++j) {
    if(temp[j] == temp_missing || prec[j] == prec_missing || pas[j] == pas_missing) {
      pas[j] = pas_missing;
    } else {
      const double v = dpm[month] * prec[j] / (1 + exp(-(temp[j] - x0[month]) / b[month]));
      pas[j] += v;
    }
  }
}

void number_frost_free_days(const double* temp, double* nffd, const double temp_missing, const double nffd_missing, const int month, const int data_size) {
  const double a = 1.15;
  const double b = 0.40;

  for(int j = 0; j < data_size; ++j) {
    if(temp[j] == temp_missing || nffd[j] == nffd_missing) {
      nffd[j] = nffd_missing;
    } else {
      const double v = dpm[month] / (1 + a * exp(-temp[j] * b));
      nffd[j] += v;
    }
  }
}

void degree_days_pcic_inner(const double* temp, double* dd, const double temp_missing, const double dd_missing, const int degree, const double exp_range_min, const double exp_range_max, const double distrib_shape, const double pow_param, const double n, const double q, const OPERATOR op, const int month, const int data_size) {
  const double mdays = dpm[month];
  if(op == GT) {
    for(int j = 0; j < data_size; ++j) {
      const double t = temp[j];
      double& d = dd[j];
      if(t == temp_missing || d == dd_missing) {
	d = dd_missing;
      } else if(t < exp_range_min) {
	d += 0;
      } else if(t < exp_range_max) {
	const double t_delta = t - q; 
	const double t_scaled = (t_delta * distrib_shape);
	d += (t_scaled * t_scaled / (t_scaled + exp(pow_param - t_delta * n))) * mdays;
      } else {
	d += (t - degree) * mdays;
      }
    }
  } else if(op == LT) {
    for(int j = 0; j < data_size; ++j) {
      const double t = temp[j];
      double& d = dd[j];
      if(t == temp_missing || d == dd_missing) {
	d = dd_missing;
      } else if(t < exp_range_min) {
	d += (degree - t) * mdays;
      } else if(t < exp_range_max) {
	const double t_delta = q - t; 
	const double t_scaled = (t_delta * distrib_shape);
	d += (t_scaled * t_scaled / (t_scaled + exp(pow_param - t_delta * n))) * mdays;
      } else {
	d += 0;
      }
    }
  }
}

void degree_days_pcic(const double* temp, double* dd, const double temp_missing, const double dd_missing, const int degree, const OPERATOR op, const int month, const int data_size) {
  if(op == LT && degree == 0) {
    degree_days_pcic_inner(temp, dd, temp_missing, dd_missing, degree, -11.7899066630181, 7.69214413729465, 0.681813434146791, 3.67164181825537, 0.161518176109162, 7.69214413729465, op, month, data_size);
  } else if(op == GT && degree == 5) {
    degree_days_pcic_inner(temp, dd, temp_missing, dd_missing, degree, -22.0649947108863, 11.2863018596233, 0.261494273957761, 12.9438413542389, 0.351606686990658, -22.0649947108863, op, month, data_size);
  } else if(op == LT && degree == 18) {
    degree_days_pcic_inner(temp, dd, temp_missing, dd_missing, degree, 11.1379738017336, 25.7802430562119, 0.59354132885639, 4.54637311062196, 0.253178333316004, 25.7802430562119, op, month, data_size);
  } else if(op == GT && degree == 18) {
    degree_days_pcic_inner(temp, dd, temp_missing, dd_missing, degree, 7.42257877561195, 23.1212731523644, 0.448357427876676, 6.17677443749065, 0.333077588715948, 7.42257877561195, op, month, data_size);
  } else {
    assert(false);
#if 0
    for(int j = 0; j < data_size; ++j) {
      if(temp[j] == temp_missing || dd[j] == dd_missing) {
	dd[j] = dd_missing;
      } else {
	const double v = (op == GT) ? (temp[j] - degree) * dpm[month] : (degree - temp[j]) * dpm[month];
	dd[j] += (v > 0) ? v : 0;
      }
    }
#endif
  }
}
