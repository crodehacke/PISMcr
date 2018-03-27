// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PLAPSERATES_H_
#define _PLAPSERATES_H_

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"

namespace pism {

class Geometry;

template <class Model>
class PLapseRates : public Model {
public:
  PLapseRates(IceGrid::ConstPtr g, std::shared_ptr<Model> in)
    : Model(g, in) {
    m_temp_lapse_rate = 0.0;
  }

  virtual ~PLapseRates() {
    // empty
  }

protected:

  virtual void update_impl(const Geometry &geometry, double t, double dt) {
    Model::m_input_model->update(geometry, t, dt);

    m_reference_surface.update(t, dt);

    m_reference_surface.interp(t + 0.5*dt);
  }

  virtual void init_internal() {
    IceGrid::ConstPtr grid = Model::m_grid;

    options::String file(m_option_prefix + "_file",
                         "Specifies a file with top-surface boundary conditions");

    options::Integer period(m_option_prefix + "_period",
                            "Specifies the length of the climate data period", 0.0);

    options::Real reference_year(m_option_prefix + "_reference_year",
                                 "Boundary condition reference year", 0.0);

    options::Real T_lapse_rate("-temp_lapse_rate",
                               "Elevation lapse rate for the temperature, in K per km",
                               m_temp_lapse_rate);
    m_temp_lapse_rate = T_lapse_rate;

    if (not file.is_set()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "command-line option %s_file is required.",
                                    m_option_prefix.c_str());
    }

    if (reference_year.is_set()) {
      m_bc_reference_time = units::convert(Model::m_sys, reference_year, "years", "seconds");
    } else {
      m_bc_reference_time = 0;
    }

    if (period.value() < 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid %s_period %d (period length cannot be negative)",
                                    m_option_prefix.c_str(), period.value());
    }
    m_bc_period = (unsigned int)period;

    if (not m_reference_surface.was_created()) {
      unsigned int buffer_size = (unsigned int) Model::m_config->get_double("climate_forcing.buffer_size"),
        ref_surface_n_records = 1;

      PIO nc(grid->com, "netcdf3", file, PISM_READONLY);
      ref_surface_n_records = nc.inq_nrecords("usurf", "surface_altitude",
                                              Model::m_sys);
      nc.close();

      // if -..._period is not set, make n_records the minimum of the
      // buffer size and the number of available records. Otherwise try
      // to keep all available records in memory.
      if (not period.is_set()) {
        ref_surface_n_records = std::min(ref_surface_n_records, buffer_size);
      }

      if (ref_surface_n_records == 0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "can't find reference surface elevation (usurf) in %s.\n",
                                      file->c_str());
      }

      m_reference_surface.set_n_records(ref_surface_n_records);
      m_reference_surface.create(grid, "usurf");
      m_reference_surface.set_attrs("climate_forcing",
                                  "reference surface for lapse rate corrections",
                                  "m", "surface_altitude");
      m_reference_surface.set_n_evaluations_per_year((unsigned int)Model::m_config->get_double("climate_forcing.evaluations_per_year"));
    }

    Model::m_log->message(2,
               "    reading reference surface elevation from %s ...\n",
               file->c_str());

    m_reference_surface.init(file, m_bc_period, m_bc_reference_time);
  }

  void lapse_rate_correction(const IceModelVec2S &surface,
                             const IceModelVec2S &reference_surface,
                             double lapse_rate,
                             IceModelVec2S &result) const {
    if (fabs(lapse_rate) < 1e-12) {
      return;
    }

    IceModelVec::AccessList list{&surface, &reference_surface, &result};

    for (Points p(*Model::m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      result(i, j) -= lapse_rate * (surface(i,j) - reference_surface(i, j));
    }
  }
protected:
  // "mutable" is needed here because some methods (average, interp) change the state of an
  // "IceModelVec2T"
  mutable IceModelVec2T m_reference_surface;
  unsigned int m_bc_period;
  double m_bc_reference_time,          // in seconds
    m_temp_lapse_rate;
  std::string m_option_prefix;
};


} // end of namespace pism

#endif /* _PLAPSERATES_H_ */
