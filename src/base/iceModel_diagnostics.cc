// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev
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

#include <gsl/gsl_math.h>

#include "base/basalstrength/PISMYieldStress.hh"
#include "base/energy/bedrockThermalUnit.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/rheology/flowlaws.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/SSB_Modifier.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec3Custom.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "earth/PISMBedDef.hh"
#include "enthalpyConverter.hh"
#include "iceModel_diagnostics.hh"
#include "base/util/PISMVars.hh"

namespace pism {

void IceModel::init_diagnostics() {

  // Add IceModel diagnostics:
  diagnostics["cts"]              = new IceModel_cts(this);
  diagnostics["enthalpybase"]     = new IceModel_enthalpybase(this);
  diagnostics["enthalpysurf"]     = new IceModel_enthalpysurf(this);
  diagnostics["hardav"]           = new IceModel_hardav(this);
  diagnostics["liqfrac"]          = new IceModel_liqfrac(this);
  diagnostics["proc_ice_area"]    = new IceModel_proc_ice_area(this);
  diagnostics["rank"]             = new IceModel_rank(this);
  diagnostics["temp"]             = new IceModel_temp(this);
  diagnostics["temp_pa"]          = new IceModel_temp_pa(this);
  diagnostics["tempbase"]         = new IceModel_tempbase(this);
  diagnostics["tempicethk"]       = new IceModel_tempicethk(this);
  diagnostics["tempicethk_basal"] = new IceModel_tempicethk_basal(this);
  diagnostics["temppabase"]       = new IceModel_temppabase(this);
  diagnostics["tempsurf"]         = new IceModel_tempsurf(this);
  diagnostics["dHdt"]             = new IceModel_dHdt(this);

  if (flux_divergence.was_created()) {
    diagnostics["flux_divergence"] = new IceModel_flux_divergence(this);
  }

  if (climatic_mass_balance_cumulative.was_created()) {
    diagnostics["climatic_mass_balance_cumulative"] = new IceModel_climatic_mass_balance_cumulative(this);
  }

  if (nonneg_flux_2D_cumulative.was_created()) {
    diagnostics["nonneg_flux_cumulative"] = new IceModel_nonneg_flux_2D_cumulative(this);
  }

  if (grounded_basal_flux_2D_cumulative.was_created()) {
    diagnostics["grounded_basal_flux_cumulative"] = new IceModel_grounded_basal_flux_2D_cumulative(this);
  }

  if (floating_basal_flux_2D_cumulative.was_created()) {
    diagnostics["floating_basal_flux_cumulative"] = new IceModel_floating_basal_flux_2D_cumulative(this);
  }

  if (discharge_flux_2D_cumulative.was_created()) {
    diagnostics["discharge_flux_cumulative"] = new IceModel_discharge_flux_2D_cumulative(this);
  }

  //ccr -- begin
  //ccr - land
  if (land_flux_2D.was_created()) {
    diagnostics["land_flux_2D"] = new IceModel_land_flux_2D(this);
  }

  if (land_flux_2D_cumulative.was_created()) {
    diagnostics["land_flux_2D_cumulative"] = new IceModel_land_flux_2D_cumulative(this);
  }

  //ccr - ocean
  if (ocean_flux_2D.was_created()) {
    diagnostics["ocean_flux_2D"] = new IceModel_ocean_flux_2D(this);
  }

  if (ocean_flux_2D_cumulative.was_created()) {
    diagnostics["ocean_flux_2D_cumulative"] = new IceModel_ocean_flux_2D_cumulative(this);
  }
  //ccr -- end


#if (PISM_USE_PROJ4==1)
  if (global_attributes.has_attribute("proj4")) {
    std::string proj4 = global_attributes.get_string("proj4");
    diagnostics["lat_bnds"] = new IceModel_lat_lon_bounds(this, "lat", proj4);
    diagnostics["lon_bnds"] = new IceModel_lat_lon_bounds(this, "lon", proj4);
  }
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

  ts_diagnostics["ivol"]          = new IceModel_ivol(this);
  ts_diagnostics["slvol"]         = new IceModel_slvol(this);
  ts_diagnostics["divoldt"]       = new IceModel_divoldt(this);
  ts_diagnostics["iarea"]         = new IceModel_iarea(this);
  ts_diagnostics["imass"]         = new IceModel_imass(this);
  ts_diagnostics["dimassdt"]      = new IceModel_dimassdt(this);
  ts_diagnostics["ivoltemp"]      = new IceModel_ivoltemp(this);
  ts_diagnostics["ivolcold"]      = new IceModel_ivolcold(this);
  ts_diagnostics["ivolg"]         = new IceModel_ivolg(this);
  ts_diagnostics["ivolf"]         = new IceModel_ivolf(this);
  ts_diagnostics["iareatemp"]     = new IceModel_iareatemp(this);
  ts_diagnostics["iareacold"]     = new IceModel_iareacold(this);
  ts_diagnostics["iareag"]        = new IceModel_iareag(this);
  ts_diagnostics["iareaf"]        = new IceModel_iareaf(this);
  ts_diagnostics["dt"]            = new IceModel_dt(this);
  ts_diagnostics["max_diffusivity"] = new IceModel_max_diffusivity(this);
  ts_diagnostics["ienthalpy"]     = new IceModel_ienthalpy(this);
  ts_diagnostics["max_hor_vel"]   = new IceModel_max_hor_vel(this);

  ts_diagnostics["surface_ice_flux"]   = new IceModel_surface_flux(this);
  ts_diagnostics["surface_ice_flux_cumulative"]   = new IceModel_surface_flux_cumulative(this);
  ts_diagnostics["grounded_basal_ice_flux"]     = new IceModel_grounded_basal_flux(this);
  ts_diagnostics["grounded_basal_ice_flux_cumulative"]     = new IceModel_grounded_basal_flux_cumulative(this);
  ts_diagnostics["sub_shelf_ice_flux"] = new IceModel_sub_shelf_flux(this);
  ts_diagnostics["sub_shelf_ice_flux_cumulative"] = new IceModel_sub_shelf_flux_cumulative(this);
  ts_diagnostics["nonneg_rule_flux"]   = new IceModel_nonneg_flux(this);
  ts_diagnostics["nonneg_rule_flux_cumulative"]   = new IceModel_nonneg_flux_cumulative(this);
  ts_diagnostics["discharge_flux"]    = new IceModel_discharge_flux(this);
  ts_diagnostics["discharge_flux_cumulative"]    = new IceModel_discharge_flux_cumulative(this);
  ts_diagnostics["H_to_Href_flux"] = new IceModel_H_to_Href_flux(this);
  ts_diagnostics["Href_to_H_flux"] = new IceModel_Href_to_H_flux(this);
  ts_diagnostics["sum_divQ_flux"]  = new IceModel_sum_divQ_flux(this);

  //ccr -- begin
  ts_diagnostics["land_flux"]  = new IceModel_land_flux(this);
  ts_diagnostics["land_flux_cumulative"]  = new IceModel_land_flux_cumulative(this);
  ts_diagnostics["ocean_flux"] = new IceModel_ocean_flux(this);
  ts_diagnostics["ocean_flux_cumulative"] = new IceModel_ocean_flux_cumulative(this);
  //ccr -- end

  // Get diagnostics supported by the stress balance object:
  if (stress_balance != NULL) {
    stress_balance->get_diagnostics(diagnostics, ts_diagnostics);
  }

  // Get diagnostics supported by the surface model:
  if (surface != NULL) {
    surface->get_diagnostics(diagnostics, ts_diagnostics);
  }

  // Get diagnostics supported by the ocean model:
  if (ocean != NULL) {
    ocean->get_diagnostics(diagnostics, ts_diagnostics);
  }

  // Get diagnostics supported by the bed deformation model:
  if (beddef != NULL) {
    beddef->get_diagnostics(diagnostics, ts_diagnostics);
  }

  if (basal_yield_stress_model != NULL) {
    basal_yield_stress_model->get_diagnostics(diagnostics, ts_diagnostics);
  }

  if (subglacial_hydrology != NULL) {
    subglacial_hydrology->get_diagnostics(diagnostics, ts_diagnostics);
  }
}

void IceModel::list_diagnostics() {
  PetscErrorCode ierr;

  ierr = PetscPrintf(m_grid->com, "\n");
  PISM_CHK(ierr, "PetscPrintf");

  // quantities with dedicated storage
  {
    std::set<std::string> list = m_grid->variables().keys();

    if (beddef != NULL) {
      beddef->add_vars_to_output("big", list);
    }

    if (btu != NULL) {
      btu->add_vars_to_output("big", list);
    }

    if (basal_yield_stress_model != NULL) {
      basal_yield_stress_model->add_vars_to_output("big", list);
    }

    if (subglacial_hydrology != NULL) {
      subglacial_hydrology->add_vars_to_output("big", list);
    }

    if (stress_balance != NULL) {
      stress_balance->add_vars_to_output("big", list);
    }

    if (ocean != NULL) {
      ocean->add_vars_to_output("big", list);
    }

    if (surface != NULL) {
      surface->add_vars_to_output("big", list);
    }

    for (unsigned int d = 3; d > 1; --d) {

      ierr = PetscPrintf(m_grid->com,
                         "======== Available %dD quantities with dedicated storage ========\n",
                         d);
      PISM_CHK(ierr, "PetscPrintf");

      std::set<std::string>::iterator j;
      for (j = list.begin(); j != list.end(); ++j) {
        const IceModelVec *v = NULL;

        if (m_grid->variables().is_available(*j)) {
          v = m_grid->variables().get(*j);
        }

        if (v != NULL && v->get_ndims() == d) {
          const SpatialVariableMetadata &var = v->metadata();

          std::string
            name                = var.get_name(),
            units               = var.get_string("units"),
            glaciological_units = var.get_string("glaciological_units"),
            long_name           = var.get_string("long_name");

          if (not glaciological_units.empty()) {
            units = glaciological_units;
          }

          ierr = PetscPrintf(m_grid->com,
                             "   Name: %s [%s]\n"
                             "       - %s\n\n", name.c_str(), units.c_str(), long_name.c_str());
          PISM_CHK(ierr, "PetscPrintf");
        }
      }
    }

  } // done with quantities with dedicated storage

  // 2D and 3D diagnostics
  for (unsigned int d = 3; d > 1; --d) {

    ierr = PetscPrintf(m_grid->com,
                       "======== Available %dD diagnostic quantities ========\n",
                       d);
    PISM_CHK(ierr, "PetscPrintf");

    std::map<std::string, Diagnostic*>::iterator j = diagnostics.begin();
    while (j != diagnostics.end()) {
      Diagnostic *diag = j->second;

      std::string name           = j->first,
        units               = diag->get_metadata().get_string("units"),
        glaciological_units = diag->get_metadata().get_string("glaciological_units");

      if (glaciological_units.empty() == false) {
        units = glaciological_units;
      }

      if (diag->get_metadata().get_n_spatial_dimensions() == d) {

        ierr = PetscPrintf(m_grid->com, "   Name: %s [%s]\n", name.c_str(), units.c_str());
        PISM_CHK(ierr, "PetscPrintf");

        for (int k = 0; k < diag->get_nvars(); ++k) {
          SpatialVariableMetadata var = diag->get_metadata(k);

          std::string long_name = var.get_string("long_name");

          ierr = PetscPrintf(m_grid->com, "      -  %s\n", long_name.c_str());
          PISM_CHK(ierr, "PetscPrintf");
        }

        ierr = PetscPrintf(m_grid->com, "\n");
        PISM_CHK(ierr, "PetscPrintf");
      }

      ++j;
    }
  }

  // scalar time-series
  ierr = PetscPrintf(m_grid->com, "======== Available time-series ========\n");
  PISM_CHK(ierr, "PetscPrintf");

  std::map<std::string, TSDiagnostic*>::iterator j = ts_diagnostics.begin();
  while (j != ts_diagnostics.end()) {
    TSDiagnostic *diag = j->second;

    std::string name = j->first,
      long_name = diag->get_string("long_name"),
      units = diag->get_string("units"),
      glaciological_units = diag->get_string("glaciological_units");

    if (glaciological_units.empty() == false) {
      units = glaciological_units;
    }

    ierr = PetscPrintf(m_grid->com,
                       "   Name: %s [%s]\n"
                       "      -  %s\n\n",
                       name.c_str(), units.c_str(), long_name.c_str());
    PISM_CHK(ierr, "PetscPrintf");

    ++j;
  }
}


IceModel_hardav::IceModel_hardav(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "hardav"));

  // choice to use SSA power; see #285
  const double power = 1.0 / m_grid->ctx()->config()->get_double("ssa_Glen_exponent");
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("_FillValue", m_grid->ctx()->config()->get_double("fill_value"));
}

//! \brief Computes vertically-averaged ice hardness.
IceModelVec::Ptr IceModel_hardav::compute() {
  const double fillval = m_grid->ctx()->config()->get_double("fill_value");
  double *Eij; // columns of enthalpy values

  const rheology::FlowLaw *flow_law = model->stress_balance->get_stressbalance()->flow_law();
  if (flow_law == NULL) {
    flow_law = model->stress_balance->get_ssb_modifier()->flow_law();
    if (flow_law == NULL) {
      throw RuntimeError("Can't compute vertically-averaged hardness: no flow law is used.");
    }
  }

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hardav", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(model->Enth3);
  list.add(model->ice_thickness);
  list.add(*result);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Eij = model->Enth3.get_column(i,j);
      const double H = model->ice_thickness(i,j);
      if (mask.icy(i, j)) {
        (*result)(i,j) = flow_law->averaged_hardness(H, m_grid->kBelowHeight(H),
                                                     &(m_grid->z()[0]), Eij);
      } else { // put negative value below valid range
        (*result)(i,j) = fillval;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


IceModel_rank::IceModel_rank(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "rank"));

  set_attrs("processor rank", "", "1", "", 0);
  m_vars[0].set_time_independent(true);
}

IceModelVec::Ptr IceModel_rank::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "rank", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  IceModelVec::AccessList list;
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    (*result)(p.i(),p.j()) = m_grid->rank();
  }

  return result;
}


IceModel_cts::IceModel_cts(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "cts", m_grid->z()));

  set_attrs("cts = E/E_s(p), so cold-temperate transition surface is at cts = 1", "",
            "", "", 0);
}

IceModelVec::Ptr IceModel_cts::compute() {

  // update vertical levels (in case the grid was extended
  m_vars[0].set_levels(m_grid->z());

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "cts", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->setCTSFromEnthalpy(*result);

  return result;
}

IceModel_proc_ice_area::IceModel_proc_ice_area(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "proc_ice_area"));

  set_attrs("number of cells containing ice in a processor's domain", "",
            "", "", 0);
  m_vars[0].set_time_independent(true);
}

IceModelVec::Ptr IceModel_proc_ice_area::compute() {

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2Int *ice_mask = m_grid->variables().get_2d_mask("mask");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "proc_ice_area", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  int ice_filled_cells = 0;

  MaskQuery mask(*ice_mask);

  IceModelVec::AccessList list;
  list.add(*ice_mask);
  list.add(*thickness);
  for (Points p(*m_grid); p; p.next()) {
    if (mask.icy(p.i(), p.j())) {
      ice_filled_cells += 1;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    (*result)(p.i(), p.j()) = ice_filled_cells;
  }

  return result;
}


IceModel_temp::IceModel_temp(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "temp", m_grid->z()));

  set_attrs("ice temperature", "land_ice_temperature", "K", "K", 0);
  m_vars[0].set_double("valid_min", 0);
}

IceModelVec::Ptr IceModel_temp::compute() {

  // update vertical levels (in case the grid was extended
  m_vars[0].set_levels(m_grid->z());

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "temp", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3 *enthalpy = m_grid->variables().get_3d_scalar("enthalpy");

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Tij = result->get_column(i,j);
      Enthij = enthalpy->get_column(i,j);
      for (unsigned int k=0; k <m_grid->Mz(); ++k) {
        const double depth = (*thickness)(i,j) - m_grid->z(k);
        Tij[k] = EC->temperature(Enthij[k], EC->pressure(depth));
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


IceModel_temp_pa::IceModel_temp_pa(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "temp_pa", m_grid->z()));

  set_attrs("pressure-adjusted ice temperature (degrees above pressure-melting point)", "",
            "deg_C", "deg_C", 0);
  m_vars[0].set_double("valid_max", 0);
}

IceModelVec::Ptr IceModel_temp_pa::compute() {
  bool cold_mode = m_grid->ctx()->config()->get_boolean("do_cold_ice_methods");
  double melting_point_temp = m_grid->ctx()->config()->get_double("water_melting_point_temperature");

  // update vertical levels (in case the m_grid was extended
  m_vars[0].set_levels(m_grid->z());

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "temp_pa", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3  *enthalpy  = m_grid->variables().get_3d_scalar("enthalpy");

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      Tij = result->get_column(i,j);
      Enthij = enthalpy->get_column(i,j);
      for (unsigned int k=0; k < m_grid->Mz(); ++k) {
        const double depth = (*thickness)(i,j) - m_grid->z(k),
          p = EC->pressure(depth);
        Tij[k] = EC->pressure_adjusted_temperature(Enthij[k], p);

        if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
          // is 273.15
          if (EC->is_temperate(Enthij[k],p) && ((*thickness)(i,j) > 0)) {
            Tij[k] = melting_point_temp;
          }
        }

      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->shift(-melting_point_temp);

  return result;
}

IceModel_temppabase::IceModel_temppabase(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "temppabase"));

  set_attrs("pressure-adjusted ice temperature at the base of ice", "",
            "Celsius", "Celsius", 0);
}

IceModelVec::Ptr IceModel_temppabase::compute() {

  bool cold_mode = m_grid->ctx()->config()->get_boolean("do_cold_ice_methods");
  double melting_point_temp = m_grid->ctx()->config()->get_double("water_melting_point_temperature");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "temp_pa_base", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3 *enthalpy = m_grid->variables().get_3d_scalar("enthalpy");

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  const double *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      Enthij = enthalpy->get_column(i,j);

      const double depth = (*thickness)(i,j),
        p = EC->pressure(depth);
      (*result)(i,j) = EC->pressure_adjusted_temperature(Enthij[0], p);

      if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
        // is 273.15
        if (EC->is_temperate(Enthij[0],p) && ((*thickness)(i,j) > 0)) {
          (*result)(i,j) = melting_point_temp;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->shift(-melting_point_temp);

  return result;
}

IceModel_enthalpysurf::IceModel_enthalpysurf(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "enthalpysurf"));

  set_attrs("ice enthalpy at 1m below the ice surface", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_grid->ctx()->config()->get_double("fill_value"));
}

IceModelVec::Ptr IceModel_enthalpysurf::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "enthalpysurf", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  double fill_value = m_grid->ctx()->config()->get_double("fill_value");

  // compute levels corresponding to 1 m below the ice surface:

  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) = std::max(model->ice_thickness(i,j) - 1.0, 0.0);
  }

  model->Enth3.getSurfaceValues(*result, *result);  // z=0 slice

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (model->ice_thickness(i,j) <= 1.0) {
      (*result)(i,j) = fill_value;
    }
  }

  return result;
}

IceModel_enthalpybase::IceModel_enthalpybase(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "enthalpybase"));

  set_attrs("ice enthalpy at the base of ice", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_grid->ctx()->config()->get_double("fill_value"));
}

IceModelVec::Ptr IceModel_enthalpybase::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "enthalpybase", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->Enth3.getHorSlice(*result, 0.0);  // z=0 slice

  result->mask_by(model->ice_thickness, m_grid->ctx()->config()->get_double("fill_value"));

  return result;
}


IceModel_tempbase::IceModel_tempbase(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "tempbase"));

  set_attrs("ice temperature at the base of ice", "",
            "K", "K", 0);
  m_vars[0].set_double("_FillValue", m_grid->ctx()->config()->get_double("fill_value"));
}

IceModelVec::Ptr IceModel_tempbase::compute() {

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::Ptr enth = IceModel_enthalpybase(model).compute();

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  // result contains basal enthalpy; note that it is allocated by
  // IceModel_enthalpybase::compute().

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(*result);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double depth = (*thickness)(i,j),
        pressure = EC->pressure(depth);
      if (mask.icy(i, j)) {
        (*result)(i,j) = EC->temperature((*result)(i,j), pressure);
      } else {
        (*result)(i,j) = m_grid->ctx()->config()->get_double("fill_value");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->metadata(0) = m_vars[0];
  return result;
}

IceModel_tempsurf::IceModel_tempsurf(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "tempsurf"));

  set_attrs("ice temperature at 1m below the ice surface", "",
            "K", "K", 0);
  m_vars[0].set_double("_FillValue", m_grid->ctx()->config()->get_double("fill_value"));
}

IceModelVec::Ptr IceModel_tempsurf::compute() {

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::Ptr enth = IceModel_enthalpysurf(model).compute();
  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  // result contains surface enthalpy; note that it is allocated by
  // IceModel_enthalpysurf::compute().

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*thickness);

  double depth = 1.0,
    pressure = EC->pressure(depth);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((*thickness)(i,j) > 1) {
        (*result)(i,j) = EC->temperature((*result)(i,j), pressure);
      } else {
        (*result)(i,j) = m_grid->ctx()->config()->get_double("fill_value");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->metadata(0) = m_vars[0];
  return result;
}


IceModel_liqfrac::IceModel_liqfrac(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "liqfrac", m_grid->z()));

  set_attrs("liquid water fraction in ice (between 0 and 1)", "",
            "1", "1", 0);
  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("valid_max", 1);
}

IceModelVec::Ptr IceModel_liqfrac::compute() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "liqfrac", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  bool cold_mode = m_grid->ctx()->config()->get_boolean("do_cold_ice_methods");

  if (cold_mode) {
    result->set(0.0);
  } else {
    model->compute_liquid_water_fraction(model->Enth3, *result);
  }

  return result;
}

IceModel_tempicethk::IceModel_tempicethk(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "tempicethk"));

  set_attrs("temperate ice thickness (total column content)", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_grid->ctx()->config()->get_double("fill_value"));
}

IceModelVec::Ptr IceModel_tempicethk::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "tempicethk", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  double *Enth = NULL;

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(*result);
  list.add(model->Enth3);
  list.add(model->ice_thickness);

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.icy(i, j)) {
        Enth = model->Enth3.get_column(i,j);
        double temperate_ice_thickness = 0.0;
        double ice_thickness = model->ice_thickness(i,j);
        const unsigned int ks = m_grid->kBelowHeight(ice_thickness);

        for (unsigned int k=0; k<ks; ++k) { // FIXME issue #15
          double pressure = EC->pressure(ice_thickness - m_grid->z(k));

          if (EC->is_temperate(Enth[k], pressure)) {
            temperate_ice_thickness += m_grid->z(k+1) - m_grid->z(k);
          }
        }

        double pressure = EC->pressure(ice_thickness - m_grid->z(ks));
        if (EC->is_temperate(Enth[ks], pressure)) {
          temperate_ice_thickness += ice_thickness - m_grid->z(ks);
        }

        (*result)(i,j) = temperate_ice_thickness;
      } else {
        // ice-free
        (*result)(i,j) = m_grid->ctx()->config()->get_double("fill_value");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceModel_tempicethk_basal::IceModel_tempicethk_basal(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "tempicethk_basal"));

  set_attrs("thickness of the basal layer of temperate ice", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_grid->ctx()->config()->get_double("fill_value"));
}

/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
IceModelVec::Ptr IceModel_tempicethk_basal::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "tempicethk_basal", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  double *Enth = NULL;
  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  MaskQuery mask(model->vMask);

  const double fill_value = m_grid->ctx()->config()->get_double("fill_value");

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(*result);
  list.add(model->ice_thickness);
  list.add(model->Enth3);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double thk = model->ice_thickness(i,j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (mask.ice_free(i, j)) {
        (*result)(i,j) = fill_value;
        continue;
      }

      Enth = model->Enth3.get_column(i,j);
      double pressure;
      unsigned int ks = m_grid->kBelowHeight(thk),
        k = 0;

      while (k <= ks) {         // FIXME issue #15
        pressure = EC->pressure(thk - m_grid->z(k));

        if (EC->is_temperate(Enth[k],pressure)) {
          k++;
        } else {
          break;
        }
      }
      // after this loop 'pressure' is equal to the pressure at the first level
      // that is cold

      // no temperate ice at all; go to the next grid point
      if (k == 0) {
        (*result)(i,j) = 0.0;
        continue;
      }

      // the whole column is temperate (except, possibly, some ice between
      // z(ks) and the total thickness; we ignore it)
      if (k == ks + 1) {
        (*result)(i,j) = m_grid->z(ks);
        continue;
      }

      double
        pressure_0 = EC->pressure(thk - m_grid->z(k-1)),
        dz         = m_grid->z(k) - m_grid->z(k-1),
        slope1     = (Enth[k] - Enth[k-1]) / dz,
        slope2     = (EC->enthalpy_cts(pressure) - EC->enthalpy_cts(pressure_0)) / dz;

      if (slope1 != slope2) {
        (*result)(i,j) = m_grid->z(k-1) +
          (EC->enthalpy_cts(pressure_0) - Enth[k-1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        (*result)(i,j) = std::max((*result)(i,j), m_grid->z(k-1));
        (*result)(i,j) = std::min((*result)(i,j), m_grid->z(k));
      } else {
        throw RuntimeError::formatted("Linear interpolation of the thickness of"
                                      " the basal temperate layer failed:\n"
                                      "(i=%d, j=%d, k=%d, ks=%d)\n",
                                      i, j, k, ks);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}

IceModel_flux_divergence::IceModel_flux_divergence(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "flux_divergence"));

  set_attrs("flux divergence", "", "m s-1", "m year-1", 0);
}

IceModelVec::Ptr IceModel_flux_divergence::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "flux_divergence", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->flux_divergence);

  return result;
}

IceModel_climatic_mass_balance_cumulative::IceModel_climatic_mass_balance_cumulative(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "climatic_mass_balance_cumulative"));

  set_attrs("cumulative ice-equivalent climatic mass balance", "",
            "kg m-2", "kg m-2", 0);
}

IceModelVec::Ptr IceModel_climatic_mass_balance_cumulative::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "climatic_mass_balance_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->climatic_mass_balance_cumulative);

  return result;
}

IceModel_ivol::IceModel_ivol(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ivol", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);

  m_ts->metadata().set_string("long_name", "total ice volume");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_ivol::update(double a, double b) {

  double value = model->compute_ice_volume();

  m_ts->append(value, a, b);
}

IceModel_slvol::IceModel_slvol(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "slvol", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m");
  m_ts->dimension_metadata().set_string("units", m_time_units);

  m_ts->metadata().set_string("long_name", "total sea-level relevant ice IN SEA-LEVEL EQUIVALENT");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_slvol::update(double a, double b) {

  double value = model->compute_sealevel_volume();

  m_ts->append(value, a, b);
}

IceModel_divoldt::IceModel_divoldt(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "divoldt", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3 s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->rate_of_change = true;

  m_ts->metadata().set_string("long_name", "total ice volume rate of change");
}

void IceModel_divoldt::update(double a, double b) {

  double value = model->compute_ice_volume();

  // note that "value" below *should* be the ice volume
  m_ts->append(value, a, b);
}


IceModel_iarea::IceModel_iarea(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "iarea", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total ice area");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_iarea::update(double a, double b) {

  double value = model->compute_ice_area();

  m_ts->append(value, a, b);
}

IceModel_imass::IceModel_imass(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "imass", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total ice mass");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_imass::update(double a, double b) {

  double value = model->compute_ice_volume();

  m_ts->append(value * m_grid->ctx()->config()->get_double("ice_density"), a, b);
}


IceModel_dimassdt::IceModel_dimassdt(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "dimassdt", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total ice mass rate of change");

  m_ts->rate_of_change = true;
}

void IceModel_dimassdt::update(double a, double b) {

  double value = model->compute_ice_volume();

  m_ts->append(value * m_grid->ctx()->config()->get_double("ice_density"), a, b);
}


IceModel_ivoltemp::IceModel_ivoltemp(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ivoltemp", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total volume of temperate ice");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_ivoltemp::update(double a, double b) {

  double value = model->compute_ice_volume_temperate();

  m_ts->append(value, a, b);
}


IceModel_ivolcold::IceModel_ivolcold(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ivolcold", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total volume of cold ice");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_ivolcold::update(double a, double b) {

  double value = model->compute_ice_volume_cold();

  m_ts->append(value, a, b);
}

IceModel_iareatemp::IceModel_iareatemp(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "iareatemp", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "ice-covered area where basal ice is temperate");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_iareatemp::update(double a, double b) {

  double value = model->compute_ice_area_temperate();

  m_ts->append(value, a, b);
}

IceModel_iareacold::IceModel_iareacold(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "iareacold", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "ice-covered area where basal ice is cold");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_iareacold::update(double a, double b) {

  double value = model->compute_ice_area_cold();

  m_ts->append(value, a, b);
}

IceModel_ienthalpy::IceModel_ienthalpy(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ienthalpy", m_time_dimension_name);

  m_ts->metadata().set_string("units", "J");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total ice enthalpy");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_ienthalpy::update(double a, double b) {

  double value = model->compute_ice_enthalpy();

  m_ts->append(value, a, b);
}

IceModel_iareag::IceModel_iareag(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "iareag", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total grounded ice area");
}

void IceModel_iareag::update(double a, double b) {

  double value = model->compute_ice_area_grounded();

  m_ts->append(value, a, b);
}

IceModel_iareaf::IceModel_iareaf(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "iareaf", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total floating ice area");
}

void IceModel_iareaf::update(double a, double b) {

  double value = model->compute_ice_area_floating();

  m_ts->append(value, a, b);
}

IceModel_dt::IceModel_dt(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "dt", m_time_dimension_name);

  m_ts->metadata().set_string("units", "second");
  m_ts->metadata().set_string("glaciological_units", "year");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass continuity time step");
}

void IceModel_dt::update(double a, double b) {

  m_ts->append(model->dt, a, b);
}

IceModel_max_diffusivity::IceModel_max_diffusivity(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "max_diffusivity", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2 s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "maximum diffusivity");
}

void IceModel_max_diffusivity::update(double a, double b) {
  double value = model->stress_balance->max_diffusivity();

  m_ts->append(value, a, b);
}

IceModel_surface_flux::IceModel_surface_flux(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "surface_ice_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total over ice domain of top surface ice mass flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_surface_flux::update(double a, double b) {

  double value = model->surface_ice_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_surface_flux_cumulative::IceModel_surface_flux_cumulative(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "surface_ice_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative total over ice domain of top surface ice mass flux");
}

void IceModel_surface_flux_cumulative::update(double a, double b) {

  double value = model->surface_ice_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_grounded_basal_flux::IceModel_grounded_basal_flux(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "grounded_basal_ice_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total over grounded ice domain of basal mass flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_grounded_basal_flux::update(double a, double b) {

  double value = model->grounded_basal_ice_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_grounded_basal_flux_cumulative::IceModel_grounded_basal_flux_cumulative(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "grounded_basal_ice_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative total grounded basal mass flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_grounded_basal_flux_cumulative::update(double a, double b) {

  double value = model->grounded_basal_ice_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_sub_shelf_flux::IceModel_sub_shelf_flux(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "sub_shelf_ice_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total sub-shelf ice flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_sub_shelf_flux::update(double a, double b) {

  double value = model->sub_shelf_ice_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_sub_shelf_flux_cumulative::IceModel_sub_shelf_flux_cumulative(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "sub_shelf_ice_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative total sub-shelf ice flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_sub_shelf_flux_cumulative::update(double a, double b) {

  double value = model->sub_shelf_ice_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_nonneg_flux::IceModel_nonneg_flux(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "nonneg_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_nonneg_flux::update(double a, double b) {

  double value = model->nonneg_rule_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_nonneg_flux_cumulative::IceModel_nonneg_flux_cumulative(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "nonneg_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative 'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_nonneg_flux_cumulative::update(double a, double b) {

  double value = model->nonneg_rule_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_discharge_flux::IceModel_discharge_flux(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "discharge_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "discharge (calving & icebergs) flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_discharge_flux::update(double a, double b) {
  double value = model->discharge_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_discharge_flux_cumulative::IceModel_discharge_flux_cumulative(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "discharge_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative discharge (calving etc.) flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_discharge_flux_cumulative::update(double a, double b) {
  double value = model->discharge_flux_cumulative;

  m_ts->append(value, a, b);
}

  /////////////////////////
  //ccr -- begin
IceModel_land_flux::IceModel_land_flux(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "land_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "land flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_land_flux::update(double a, double b) {
  double value = model->land_flux_cumulative;

  m_ts->append(value, a, b);
}

  /////
IceModel_land_flux_cumulative::IceModel_land_flux_cumulative(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "land_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative land flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_land_flux_cumulative::update(double a, double b) {
  double value = model->land_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_ocean_flux::IceModel_ocean_flux(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ocean_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "ocean flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_ocean_flux::update(double a, double b) {
  double value = model->ocean_flux_cumulative;

  m_ts->append(value, a, b);
}

IceModel_ocean_flux_cumulative::IceModel_ocean_flux_cumulative(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ocean_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative ocean flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_ocean_flux_cumulative::update(double a, double b) {
  double value = model->ocean_flux_cumulative;

  m_ts->append(value, a, b);
}

  //ccr -- end

IceModel_dHdt::IceModel_dHdt(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "dHdt"));

  set_attrs("ice thickness rate of change", "tendency_of_land_ice_thickness",
            "m s-1", "m year-1", 0);

  Config::ConstPtr config = m_grid->ctx()->config();
  units::System::Ptr sys = m_grid->ctx()->unit_system();
  double fill_value = units::convert(sys, config->get_double("fill_value"), "m/year", "m/s");

  m_vars[0].set_double("valid_min",  units::convert(m_sys, -1e6, "m/year", "m/s"));
  m_vars[0].set_double("valid_max",  units::convert(m_sys,  1e6, "m/year", "m/s"));
  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_string("cell_methods", "time: mean");

  last_ice_thickness.create(m_grid, "last_ice_thickness", WITHOUT_GHOSTS);
  last_ice_thickness.set_attrs("internal",
                               "ice thickness at the time of the last report of dHdt",
                               "m", "land_ice_thickness");

  last_report_time = GSL_NAN;
}

IceModelVec::Ptr IceModel_dHdt::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "dHdt", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  if (gsl_isnan(last_report_time)) {
    result->set(units::convert(m_sys, 2e6, "m/year", "m/s"));
  } else {
    IceModelVec::AccessList list;
    list.add(*result);
    list.add(last_ice_thickness);
    list.add(model->ice_thickness);

    double dt = m_grid->ctx()->time()->current() - last_report_time;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      (*result)(i, j) = (model->ice_thickness(i, j) - last_ice_thickness(i, j)) / dt;
    }
  }

  // Save the ice thickness and the corresponding time:
  this->update_cumulative();

  return result;
}

void IceModel_dHdt::update_cumulative() {
  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(last_ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    last_ice_thickness(i, j) = model->ice_thickness(i, j);
  }

  last_report_time = m_grid->ctx()->time()->current();
}

IceModel_ivolg::IceModel_ivolg(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ivolg", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total grounded ice volume");
}

void IceModel_ivolg::update(double a, double b) {
  double volume = 0.0;

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(model->vMask);
  list.add(model->cell_area);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.grounded_ice(i,j)) {
      volume += model->cell_area(i,j) * model->ice_thickness(i,j);
    }
  }

  double value = GlobalSum(m_grid->com, volume);

  m_ts->append(value, a, b);
}

IceModel_ivolf::IceModel_ivolf(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "ivolf", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total floating ice volume");
}

void IceModel_ivolf::update(double a, double b) {
  double volume = 0.0;

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(model->vMask);
  list.add(model->cell_area);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.floating_ice(i,j)) {
      volume += model->cell_area(i,j) * model->ice_thickness(i,j);
    }
  }

  double value = GlobalSum(m_grid->com, volume);

  m_ts->append(value, a, b);
}

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
/*!
 * This is the value used by the adaptive time-stepping code in the CFL condition
 * for horizontal advection (i.e. in energy and mass conservation time steps).
 *
 * This is not the maximum horizontal speed, but rather the maximum of components.
 *
 * Note that this picks up the value computed during the time-step taken at a
 * reporting time. (It is not the "average over the reporting interval computed using
 * differencing in time", as other rate-of-change diagnostics.)
 */
IceModel_max_hor_vel::IceModel_max_hor_vel(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "max_hor_vel", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m/second");
  m_ts->metadata().set_string("glaciological_units", "m/year");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "maximum abs component of horizontal ice velocity"
                                " over grid in last time step during time-series reporting interval");
}

void IceModel_max_hor_vel::update(double a, double b) {

  double gmaxu = model->gmaxu, gmaxv = model->gmaxv;

  m_ts->append(gmaxu > gmaxv ? gmaxu : gmaxv, a, b);
}

IceModel_H_to_Href_flux::IceModel_H_to_Href_flux(IceModel *m)
  : TSDiag<IceModel>(m) {
  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "H_to_Href_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass flux from thk to Href");
  m_ts->metadata().set_string("comment", "does not correspond to mass gain or loss");
  m_ts->rate_of_change = true;
}

void IceModel_H_to_Href_flux::update(double a, double b) {

  m_ts->append(model->H_to_Href_flux_cumulative, a, b);
}


IceModel_Href_to_H_flux::IceModel_Href_to_H_flux(IceModel *m)
  : TSDiag<IceModel>(m) {
  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "Href_to_H_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass flux from Href to thk");
  m_ts->metadata().set_string("comment", "does not correspond to mass gain or loss");
  m_ts->rate_of_change = true;
}

void IceModel_Href_to_H_flux::update(double a, double b) {

  m_ts->append(model->Href_to_H_flux_cumulative, a, b);
}



IceModel_sum_divQ_flux::IceModel_sum_divQ_flux(IceModel *m)
  : TSDiag<IceModel>(m) {
  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "sum_divQ_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "sum(divQ)");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_sum_divQ_flux::update(double a, double b) {

  m_ts->append(model->sum_divQ_SIA_cumulative + model->sum_divQ_SSA_cumulative,
             a, b);
}


IceModel_nonneg_flux_2D_cumulative::IceModel_nonneg_flux_2D_cumulative(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "nonneg_flux_cumulative"));

  set_attrs("cumulative non-negative rule (thk >= 0) flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_nonneg_flux_2D_cumulative::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->nonneg_flux_2D_cumulative);

  return result;
}

IceModel_grounded_basal_flux_2D_cumulative::IceModel_grounded_basal_flux_2D_cumulative(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "grounded_basal_flux_cumulative"));

  set_attrs("cumulative grounded basal flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_grounded_basal_flux_2D_cumulative::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->grounded_basal_flux_2D_cumulative);

  return result;
}

IceModel_floating_basal_flux_2D_cumulative::IceModel_floating_basal_flux_2D_cumulative(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "floating_basal_flux_cumulative"));

  set_attrs("cumulative floating basal flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_floating_basal_flux_2D_cumulative::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->floating_basal_flux_2D_cumulative);

  return result;
}


IceModel_discharge_flux_2D_cumulative::IceModel_discharge_flux_2D_cumulative(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                     "discharge_flux_cumulative"));

  set_attrs("cumulative ice discharge (calving) flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_discharge_flux_2D_cumulative::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "discharge_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->discharge_flux_2D_cumulative);

  return result;
}

  // ccr -- begin
IceModel_land_flux_2D::IceModel_land_flux_2D(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                     "land_flux_2D"));

  set_attrs("ice land flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_land_flux_2D::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "land_flux_2D", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->land_flux_2D);

  return result;
}

IceModel_land_flux_2D_cumulative::IceModel_land_flux_2D_cumulative(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                     "land_flux_2D_cumulative"));

  set_attrs("cumulative ice land flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_land_flux_2D_cumulative::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "land_flux_2D_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->land_flux_2D_cumulative);

  return result;
}
  /////////////////
IceModel_ocean_flux_2D::IceModel_ocean_flux_2D(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                     "ocean_flux_2D"));

  set_attrs("ice ocean flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_ocean_flux_2D::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "ocean_flux_2D", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->ocean_flux_2D);

  return result;
}

IceModel_ocean_flux_2D_cumulative::IceModel_ocean_flux_2D_cumulative(IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                     "ocean_flux_2D_cumulative"));

  set_attrs("cumulative ice ocean flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_ocean_flux_2D_cumulative::compute() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "ocean_flux_2D_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->ocean_flux_2D_cumulative);

  return result;
}

  // ccr -- end

#if (PISM_USE_PROJ4==1)
IceModel_lat_lon_bounds::IceModel_lat_lon_bounds(IceModel *m,
                                                 const std::string &var_name,
                                                 const std::string &proj_string)
  : Diag<IceModel>(m) {
  assert(var_name == "lat" || var_name == "lon");
  m_var_name = var_name;

  // set metadata:
  std::vector<double> levels(4);
  for (int k = 0; k < 4; ++k) {
    levels[k] = k;
  }

  m_vars.push_back(SpatialVariableMetadata(m_sys, m_var_name + "_bnds", levels));
  m_vars[0].get_z().set_name("nv4");
  m_vars[0].get_z().clear_all_strings();
  m_vars[0].get_z().clear_all_doubles();
  m_vars[0].set_time_independent(true);

  if (m_var_name == "lon") {
    set_attrs("longitude bounds", "", "degree_east", "degree_east", 0);
    m_vars[0].set_double("valid_min", -180);
    m_vars[0].set_double("valid_max", 180);
  } else {
    set_attrs("latitude bounds", "", "degree_north", "degree_north", 0);
    m_vars[0].set_double("valid_min", -90);
    m_vars[0].set_double("valid_max", 90);
  }

  lonlat = pj_init_plus("+proj=latlong +datum=WGS84 +ellps=WGS84");
  if (lonlat == NULL) {
    throw RuntimeError("projection initialization failed\n"
                       "('+proj=latlong +datum=WGS84 +ellps=WGS84').\n");
  }

  pism = pj_init_plus(proj_string.c_str());
  if (pism == NULL) {
    // if we got here, then lonlat was allocated already
    pj_free(lonlat);
    throw RuntimeError::formatted("proj.4 string '%s' is invalid.", proj_string.c_str());
  }
}

IceModel_lat_lon_bounds::~IceModel_lat_lon_bounds() {
  pj_free(pism);
  pj_free(lonlat);
}

IceModelVec::Ptr IceModel_lat_lon_bounds::compute() {

  std::map<std::string,std::string> attrs;
  std::vector<double> indices(4);

  IceModelVec3Custom::Ptr result(new IceModelVec3Custom);
  result->create(m_grid, m_var_name + "_bnds", "nv4",
                 indices, attrs);
  result->metadata(0) = m_vars[0];

  double dx2 = 0.5 * m_grid->dx(), dy2 = 0.5 * m_grid->dy();
  double x_offsets[] = {-dx2, dx2, dx2, -dx2};
  double y_offsets[] = {-dy2, -dy2, dy2, dy2};

  bool latitude = true;
  if (m_var_name == "lon") {
    latitude = false;
  }

  IceModelVec::AccessList list;
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double x0 = m_grid->x(i), y0 = m_grid->y(j);

    double *values = result->get_column(i,j);

    for (int k = 0; k < 4; ++k) {
      double
        x = x0 + x_offsets[k],
        y = y0 + y_offsets[k];

      // compute lon,lat coordinates:
      pj_transform(pism, lonlat, 1, 1, &x, &y, NULL);

      // NB! proj.4 converts x,y pairs into lon,lat pairs in *radians*.

      if (latitude) {
        values[k] = y * RAD_TO_DEG;
      } else {
        values[k] = x * RAD_TO_DEG;
      }
    }
  }

  return result;
}
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

} // end of namespace pism
