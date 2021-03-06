// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PAYEARLYCYCLE_H_
#define _PAYEARLYCYCLE_H_

#include "coupler/PISMAtmosphere.hh"
#include "base/util/iceModelVec.hh"

namespace pism {
namespace atmosphere {

//! A class containing an incomplete implementation of an atmosphere model
//! based on a temperature parameterization using mean annual and mean July
//! (mean summer) temperatures and a cosine yearly cycle. Uses a stored
//! (constant in time) precipitation field.
class YearlyCycle : public AtmosphereModel {
public:
  YearlyCycle(IceGrid::ConstPtr g);
  virtual ~YearlyCycle();

  virtual void init();
  //! This method implements the parameterization.
  virtual void mean_precipitation(IceModelVec2S &result);
  virtual void mean_annual_temp(IceModelVec2S &result);
  virtual void begin_pointwise_access();
  virtual void end_pointwise_access();
  virtual void temp_snapshot(IceModelVec2S &result);

  virtual void init_timeseries(const std::vector<double> &ts);
  virtual void temp_time_series(int i, int j, std::vector<double> &result);
  virtual void precip_time_series(int i, int j, std::vector<double> &result);
protected:
  virtual void update_impl(double my_t, double my_dt) = 0;
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars,
                                     const PIO &nc, IO_Type nctype);
protected:
  void init_internal(const std::string &input_filename, bool regrid,
                               unsigned int start);
  double m_snow_temp_july_day;
  std::string m_reference, m_precip_filename;
  IceModelVec2S m_air_temp_mean_annual, m_air_temp_mean_july, m_precipitation;
  SpatialVariableMetadata m_air_temp_snapshot;
  std::vector<double> m_ts_times, m_cosine_cycle;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAYEARLYCYCLE_H_ */
