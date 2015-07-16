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

#ifndef _ICEMODEL_DIAGNOSTICS_H_
#define _ICEMODEL_DIAGNOSTICS_H_

#include "iceModel.hh"
#include "base/util/PISMDiagnostic.hh"

namespace pism {

//! \brief Computes vertically-averaged ice hardness.
class IceModel_hardav : public Diag<IceModel>
{
public:
  IceModel_hardav(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes a diagnostic field filled with processor rank values.
class IceModel_rank : public Diag<IceModel>
{
public:
  IceModel_rank(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes CTS, CTS = E/E_s(p).
class IceModel_cts : public Diag<IceModel>
{
public:
  IceModel_cts(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes the number of ice-filled cells is a processor's domain.
class IceModel_proc_ice_area : public Diag<IceModel>
{
public:
  IceModel_proc_ice_area(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes ice temperature from enthalpy.
class IceModel_temp : public Diag<IceModel>
{
public:
  IceModel_temp(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Compute the pressure-adjusted temperature in degrees C corresponding
//! to ice temperature.
class IceModel_temp_pa : public Diag<IceModel>
{
public:
  IceModel_temp_pa(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes basal values of the pressure-adjusted temperature.
class IceModel_temppabase : public Diag<IceModel>
{
public:
  IceModel_temppabase(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes surface values of ice enthalpy.
class IceModel_enthalpysurf : public Diag<IceModel>
{
public:
  IceModel_enthalpysurf(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes enthalpy at the base of the ice.
class IceModel_enthalpybase : public Diag<IceModel>
{
public:
  IceModel_enthalpybase(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes ice temperature at the base of the ice.
class IceModel_tempbase : public Diag<IceModel>
{
public:
  IceModel_tempbase(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes ice temperature at the surface of the ice.
class IceModel_tempsurf : public Diag<IceModel>
{
public:
  IceModel_tempsurf(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes the liquid water fraction.
class IceModel_liqfrac : public Diag<IceModel>
{
public:
  IceModel_liqfrac(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes the total thickness of temperate ice in a column.
class IceModel_tempicethk : public Diag<IceModel>
{
public:
  IceModel_tempicethk(IceModel *m);
  virtual IceModelVec::Ptr compute();
};
//! \brief Computes the thickness of the basal layer of temperate ice.
class IceModel_tempicethk_basal : public Diag<IceModel>
{
public:
  IceModel_tempicethk_basal(IceModel *m);
  virtual IceModelVec::Ptr compute();
};
//! \brief Computes the flux divergence.
class IceModel_flux_divergence : public Diag<IceModel>
{
public:
  IceModel_flux_divergence(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes the total ice volume.
class IceModel_ivol : public TSDiag<IceModel>
{
public:
  IceModel_ivol(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice volume, which is relevant for sea-level
class IceModel_slvol : public TSDiag<IceModel>
{
public:
  IceModel_slvol(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the rate of change of the total ice volume.
class IceModel_divoldt : public TSDiag<IceModel>
{
public:
  IceModel_divoldt(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice area.
class IceModel_iarea : public TSDiag<IceModel>
{
public:
  IceModel_iarea(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice mass.
class IceModel_imass : public TSDiag<IceModel>
{
public:
  IceModel_imass(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the rate of change of the total ice mass.
class IceModel_dimassdt : public TSDiag<IceModel>
{
public:
  IceModel_dimassdt(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total volume of the temperate ice.
class IceModel_ivoltemp : public TSDiag<IceModel>
{
public:
  IceModel_ivoltemp(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total volume of the cold ice.
class IceModel_ivolcold : public TSDiag<IceModel>
{
public:
  IceModel_ivolcold(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total area of the temperate ice.
class IceModel_iareatemp : public TSDiag<IceModel>
{
public:
  IceModel_iareatemp(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total area of the cold ice.
class IceModel_iareacold : public TSDiag<IceModel>
{
public:
  IceModel_iareacold(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total ice enthalpy.
class IceModel_ienthalpy : public TSDiag<IceModel>
{
public:
  IceModel_ienthalpy(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total grounded ice area.
class IceModel_iareag : public TSDiag<IceModel>
{
public:
  IceModel_iareag(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total floating ice area.
class IceModel_iareaf : public TSDiag<IceModel>
{
public:
  IceModel_iareaf(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total grounded ice volume.
class IceModel_ivolg : public TSDiag<IceModel>
{
public:
  IceModel_ivolg(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Computes the total floating ice volume.
class IceModel_ivolf : public TSDiag<IceModel>
{
public:
  IceModel_ivolf(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the mass continuity time step.
class IceModel_dt : public TSDiag<IceModel>
{
public:
  IceModel_dt(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports maximum diffusivity.
class IceModel_max_diffusivity : public TSDiag<IceModel>
{
public:
  IceModel_max_diffusivity(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total surface ice flux.
class IceModel_surface_flux : public TSDiag<IceModel>
{
public:
  IceModel_surface_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative total surface ice flux.
class IceModel_surface_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_surface_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total basal ice flux over the grounded region.
class IceModel_grounded_basal_flux : public TSDiag<IceModel>
{
public:
  IceModel_grounded_basal_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative total basal ice flux over the grounded region.
class IceModel_grounded_basal_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_grounded_basal_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total sub-shelf ice flux.
class IceModel_sub_shelf_flux : public TSDiag<IceModel>
{
public:
  IceModel_sub_shelf_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total sub-shelf ice flux.
class IceModel_sub_shelf_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_sub_shelf_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the 'numerical' ice flux resulting from enforcing the 'thk
//! >= 0' rule.
class IceModel_nonneg_flux : public TSDiag<IceModel>
{
public:
  IceModel_nonneg_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative 'numerical' ice flux resulting from enforcing the 'thk
//! >= 0' rule.
class IceModel_nonneg_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_nonneg_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the total discharge flux.
class IceModel_discharge_flux : public TSDiag<IceModel>
{
public:
  IceModel_discharge_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative total discharge flux.
class IceModel_discharge_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_discharge_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports cumulative surface mass balance.
class IceModel_climatic_mass_balance_cumulative : public Diag<IceModel>
{
public:
  IceModel_climatic_mass_balance_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Computes dHdt, the ice thickness rate of change.
class IceModel_dHdt : public Diag<IceModel>
{
public:
  IceModel_dHdt(IceModel *m);
  virtual IceModelVec::Ptr compute();
  virtual void update_cumulative();
protected:
  IceModelVec2S last_ice_thickness;
  double last_report_time;
};

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
class IceModel_max_hor_vel : public TSDiag<IceModel>
{
public:
  IceModel_max_hor_vel(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the mass flux from the mass tracked using ice thickness
//! (thk) to the mass tracked using Href.
class IceModel_H_to_Href_flux : public TSDiag<IceModel>
{
public:
  IceModel_H_to_Href_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the mass flux from the mass tracked using Href to the mass
//! tracked using ice thickness (thk).
class IceModel_Href_to_H_flux : public TSDiag<IceModel>
{
public:
  IceModel_Href_to_H_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the sum(div Q) flux (to diagnose issues in the mass
//! transport scheme).
class IceModel_sum_divQ_flux : public TSDiag<IceModel>
{
public:
  IceModel_sum_divQ_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the 2D cumulative (numerical) flux due to enforcing
//! non-negativity of ice thickness.
class IceModel_nonneg_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_nonneg_flux_2D_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};


//! \brief Reports the 2D cumulative grounded basal flux.
class IceModel_grounded_basal_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_grounded_basal_flux_2D_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Reports the 2D cumulative floating basal flux.
class IceModel_floating_basal_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_floating_basal_flux_2D_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Reports the 2D cumulative discharge (calving) flux.
class IceModel_discharge_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_discharge_flux_2D_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//////////////////////
//ccr --- begin
//! \brief Reports the integrated land interaction flux.
class IceModel_land_flux : public TSDiag<IceModel>
{
public:
  IceModel_land_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the integrated cumulative land interaction flux.
class IceModel_land_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_land_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the integrated ocean interaction flux.
class IceModel_ocean_flux : public TSDiag<IceModel>
{
public:
  IceModel_ocean_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the integrated cumulative ocean interaction flux.
class IceModel_ocean_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_ocean_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};
//ccr --- begin

//! \brief Reports the crevasses calving flux
class IceModel_crevasses_calv_flux : public TSDiag<IceModel>
{
public:
  IceModel_crevasses_calv_flux(IceModel *m);
  virtual void update(double a, double b);
};

//! \brief Reports the cumulative crevasses calving flux
class IceModel_crevasses_calv_flux_cumulative : public TSDiag<IceModel>
{
public:
  IceModel_crevasses_calv_flux_cumulative(IceModel *m);
  virtual void update(double a, double b);
};

//ccr --- end
//! \brief Reports the 2D crevasses calving flux
class IceModel_crevasses_calv_flux_2D : public Diag<IceModel>
{
public:
  IceModel_crevasses_calv_flux_2D(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Reports the 2D cumulative crevasses calving flux
class IceModel_crevasses_calv_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_crevasses_calv_flux_2D_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};
//ccr ---

//! \brief Reports the 2D land flux.
class IceModel_land_flux_2D : public Diag<IceModel>
{
public:
   IceModel_land_flux_2D(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Reports the 2D cumulative land flux.
class IceModel_land_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_land_flux_2D_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Reports the 2D ocean flux.
class IceModel_ocean_flux_2D : public Diag<IceModel>
{
public:
  IceModel_ocean_flux_2D(IceModel *m);
  virtual IceModelVec::Ptr compute();
};

//! \brief Reports the 2D cumulative ocean flux.
class IceModel_ocean_flux_2D_cumulative : public Diag<IceModel>
{
public:
  IceModel_ocean_flux_2D_cumulative(IceModel *m);
  virtual IceModelVec::Ptr compute();
};
//ccr --- end

} // end of namespace pism

#if (PISM_USE_PROJ4==1)

#include <proj_api.h>

namespace pism {

//! \brief Computes latitude and longitude bounds.
class IceModel_lat_lon_bounds : public Diag<IceModel>
{
public:
  IceModel_lat_lon_bounds(IceModel *m,
                          const std::string &var_name,
                          const std::string &proj_string);
  ~IceModel_lat_lon_bounds();
  virtual IceModelVec::Ptr compute();
protected:
  std::string m_var_name;
  projPJ pism, lonlat;
};

} // end of namespace pism

#elif (PISM_USE_PROJ4==0)
// do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

#endif  /* _ICEMODEL_DIAGNOSTICS_H_ */

