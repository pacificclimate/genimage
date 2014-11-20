#include "subset_derived_dm.h"

vector<int> get_months_for_season(int timeofyear) {
  // Build up list of months
  vector<int> months;
  if(timeofyear < DJF) {
    months.push_back(timeofyear);
  } else if(timeofyear != ANN) {
    for(int i = 0; i < 3; i++)
      months.push_back(seasons[timeofyear - DJF][i]);
  } else {
    for(int i = 0; i < 12; i++)
      months.push_back(annual[i]);
  }
  return(months);
}

DataGrid<double> DerivedSubsetDataManager::get_anomaly(const DataSpec& s) {
  DataGrid<double> dd = get_data(DataSpec(s.model, s.expt, s.timeslice, s.variable, s.timeofyear, ABSOLUTE));
  DataGrid<double> dd_base = get_data(use_alt_baseline()
				      ? DataSpec(config.baseline_model, config.baseline_expt, dd.base_period(), s.variable, s.timeofyear, ABSOLUTE)
				      : DataSpec(s.model, s.expt, dd.base_period(), s.variable, s.timeofyear, ABSOLUTE));
  double* dd_data = dd.values().get();
  double* dd_base_data = dd_base.values().get();
  double dd_missing = dd.missing();
  double base_missing = dd_base.missing();

  if(s.percent_change) {
    fprintf(stderr, "Derived percent change code activated\n");
    for(int i = 0; i < dd.grid_size(); i++)
      if(dd_data[i] == dd_missing || dd_base_data[i] == base_missing) {
	dd_data[i] = dd_missing;
      } else {
	if(dd_base_data[i] == 0)
	  dd_data[i] = 0;
	else
	  dd_data[i] = ((dd_data[i] - dd_base_data[i]) / dd_base_data[i]) * 100;
      }
  } else
    for(int i = 0; i < dd.grid_size(); i++)
      if(dd_data[i] == dd_missing || dd_base_data[i] == base_missing) {
	dd_data[i] = dd_missing;
      } else {
	dd_data[i] -= dd_base_data[i];
      }

  dd.set_data(dd.x_size(), dd.y_size(), dd.missing(), dd.values(), dd.x_grid(), dd.y_grid(), dd.proj4_string(), dd.projection(), true, dd.base_period());

  return(dd);
}

DataGrid<double> DerivedSubsetDataManager::get_data(const DataSpec& s) {
  if((s.variable == "dl00" || s.variable == "dg05" || s.variable == "dl18" || s.variable == "dg18" || s.variable == "nffd" || s.variable == "pass") && s.anom == ANOMALY)
    return(get_anomaly(s));

  if(s.variable == "dl00") {
    return(degree_days_est(s, 0, LT, 0, 1));
  } else if(s.variable == "dg05") {
    return(degree_days_est(s, 5, GT, 0, 1));
  } else if(s.variable == "dl18") {
    return(degree_days_est(s, 18, LT, 0, 1));
  } else if(s.variable == "dg18") {
    return(degree_days_est(s, 18, GT, 0, 1));
  } else if(s.variable == "nffd") {
    if(s.timeofyear == ANN) {
      return(nffd(s, 0.553488, 1.05207352));
    } else {
      return(nffd(s, 0, 1));
    }
  } else if(s.variable == "pass") {
    if(s.timeofyear == ANN) {
      return(pas(s, 19.656818, 0.8905668));
    } else {
      return(pas(s, 0, 1));
    }
  } else {
    return(DataManager::get_data(s));
  }
}

DataGrid<double> DerivedSubsetDataManager::degree_days_est(const DataSpec& s, const double degree, const OPERATOR op, const double add_factor, const double scale_factor) {
  vector<int> months = get_months_for_season(s.timeofyear);

  DataGrid<double> dd_grid(s);
  double* dd_est = 0;
  double missing = 1e20f;
  for(unsigned int i = 0; i < months.size(); i++) {
    DataGrid<double> dg = get_data(DataSpec(s.model, s.expt, s.timeslice, "temp", months[i], ABSOLUTE));
    assert(dg.has_data());
    double* data = dg.values().get();

    if(i == 0) {
      missing = dg.missing();
      dd_est = new double[dg.grid_size()];
      std::fill(dd_est, &dd_est[dg.grid_size()], 0.0);
      dd_grid.set_data(dg.x_size(), dg.y_size(), dg.missing(), boost::shared_ptr<double>(dd_est), dg.x_grid(), dg.y_grid(), dg.proj4_string(), dg.projection(), dg.anomaly(), dg.base_period());
    }

    // Add to degree days, with given operator (GT = greater than, LT = less than)
    degree_days_pcic(data, dd_est, missing, missing, degree, op, months[i], dg.grid_size());
  }

  // Apply add and scale factors
  if(add_factor != 0 && scale_factor != 1)
    for(int j = 0; j < dd_grid.grid_size(); j++)
      if(dd_est[j] != missing)
	dd_est[j] = dd_est[j] * scale_factor + add_factor;
  
  return(dd_grid);
}

DataGrid<double> DerivedSubsetDataManager::nffd(const DataSpec& s, const double add_factor, const double scale_factor) {
  vector<int> months = get_months_for_season(s.timeofyear);

  DataGrid<double> nffd_grid(s);
  double* nffd_est = 0;
  double missing = 1e20f;

  for(unsigned int i = 0; i < months.size(); i++) {
    DataGrid<double> dg = get_data(DataSpec(s.model, s.expt, s.timeslice, "temp", months[i], ABSOLUTE));
    assert(dg.has_data());
    double* data = dg.values().get();

    if(i == 0) {
      missing = dg.missing();
      nffd_est = new double[dg.grid_size()];
      std::fill(nffd_est, &nffd_est[dg.grid_size()], 0.0);
      nffd_grid.set_data(dg.x_size(), dg.y_size(), dg.missing(), boost::shared_ptr<double>(nffd_est), dg.x_grid(), dg.y_grid(), dg.proj4_string(), dg.projection(), dg.anomaly(), dg.base_period());
    }

    number_frost_free_days(data, nffd_est, missing, missing, months[i], dg.grid_size());
  }

  // Apply add and scale factors
  if(add_factor != 0 && scale_factor != 1)
    for(int j = 0; j < nffd_grid.grid_size(); j++)
      if(nffd_est[j] != missing)
	nffd_est[j] = nffd_est[j] * scale_factor + add_factor;
  
  return(nffd_grid);
}

DataGrid<double> DerivedSubsetDataManager::pas(const DataSpec& s, const double add_factor, const double scale_factor) {
  vector<int> months = get_months_for_season(s.timeofyear);

  DataGrid<double> pas_grid(s);
  double* pas_est = 0;
  double missing = 1e20f;
  
  for(unsigned int i = 0; i < months.size(); i++) {
    const int m = months[i];

    DataGrid<double> dg = get_data(DataSpec(s.model, s.expt, s.timeslice, "temp", m, ABSOLUTE));

    // FIXME: The data needs to _NOT_ be in percent. This causes fuckage left, right, and center. This is one example -- percent is specified manually here.
    DataGrid<double> dg_prec = get_data(DataSpec(s.model, s.expt, s.timeslice, "prec", m, ABSOLUTE, 1));
    assert(dg.has_data());
    assert(dg_prec.has_data());
    double* data = dg.values().get();
    double* pr_data = dg_prec.values().get();
    double temp_missing = dg.missing();
    double prec_missing = dg_prec.missing();

    if(i == 0) {
      pas_est = new double[dg.grid_size()];
      std::fill(pas_est, &pas_est[dg.grid_size()], 0.0);
      pas_grid.set_data(dg.x_size(), dg.y_size(), missing, boost::shared_ptr<double>(pas_est), dg.x_grid(), dg.y_grid(), dg.proj4_string(), dg.projection(), dg.anomaly(), dg.base_period());
    }

    precip_as_snow(data, pr_data, pas_est, temp_missing, prec_missing, missing, m, dg.grid_size());
  }

  // Apply add and scale factors
  if(add_factor != 0 && scale_factor != 1)
    for(int j = 0; j < pas_grid.grid_size(); j++)
      if(pas_est[j] != missing)
	pas_est[j] = pas_est[j] * scale_factor + add_factor;
  
  return(pas_grid);
}

