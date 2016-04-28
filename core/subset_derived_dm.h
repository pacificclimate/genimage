#ifndef __GENIMAGE_DERIVEDSUBSETDATAMANAGER_H
#define __GENIMAGE_DERIVEDSUBSETDATAMANAGER_H

#include "datamanager.h"
#include "derived_vars.h"

enum TIMESOFYEAR { JAN = 0, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC, DJF, MAM, JJA, SON, ANN };

class DerivedSubsetDataManager : public DataManager {
public:
  DerivedSubsetDataManager(Config& c) : DataManager(c) { }
  virtual DataGrid<double> get_data(const DataSpec& s);
  //DataSpec(const string model, const string expt, const string timeslice, const string variable, const int timeofyear, ANOM_TYPE anom=DEFAULT): model(model), expt(expt), timeslice(timeslice), variable(variable), timeofyear(timeofyear), anom(anom) { }

  virtual DataGrid<int> get_slmask(const DataSpec& s) { return (DataManager::get_slmask(get_slmask_spec(s))); }
  virtual DataGrid<int> get_drawmask(const DataSpec& s) { return (DataManager::get_drawmask(get_slmask_spec(s))); }
  virtual DataGrid<int> get_datamask(const DataSpec& s, const double threshold) { return (DataManager::get_datamask(get_slmask_spec(s), threshold)); }
  virtual DataGrid<double> get_areagrid(const DataSpec& s) { return (DataManager::get_areagrid(get_slmask_spec(s))); }

protected:
  DataSpec get_slmask_spec(DataSpec s) { return (DataSpec(s.model, s.expt, s.timeslice, "slmask", s.timeofyear, s.anom)); }
  DataGrid<double> get_anomaly(const DataSpec& s);
  DataGrid<double> degree_days_est(const DataSpec& s, const double degree, const OPERATOR op, const double add_factor, const double scale_factor);
  DataGrid<double> nffd(const DataSpec& s, const double add_factor, const double scale_factor);
  DataGrid<double> pas(const DataSpec& s, const double add_factor, const double scale_factor);
};

#endif
