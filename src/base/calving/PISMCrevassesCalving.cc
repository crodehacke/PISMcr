/* Copyright (C) 2013, 2014, 2015 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PISMCrevassesCalving.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/MaxTimestep.hh"
#include "coupler/PISMOcean.hh"
#include "base/util/PISMComponent.hh"

namespace pism {
namespace calving {

CrevassesCalving::CrevassesCalving(IceGrid::ConstPtr g,
                           stressbalance::StressBalance *stress_balance)
  : Component(g), m_stencil_width(2),
    m_stress_balance(stress_balance) {


  m_strain_rates.create(m_grid, "edot", WITH_GHOSTS,
			m_stencil_width,
			2);
  m_thk_loss.create(m_grid, "temporary_storage", WITH_GHOSTS, 1);

  m_strain_rates.metadata(0).set_name("edot_1");
  m_strain_rates.set_attrs("internal",
			   "major principal component of horizontal strain-rate",
			   "1/s", "", 0);

  m_strain_rates.metadata(1).set_name("edot_2");
  m_strain_rates.set_attrs("internal",
			   "minor principal component of horizontal strain-rate",
			   "1/s", "", 1);

  m_restrict_timestep = m_config->get_boolean("cfl_crevasses_calving");

  //ccr -- here crevasses calving specific variables
  m_crevasse_dw.create(m_grid, "temporary_crevasse_dw", WITH_GHOSTS, 1);
  m_crevasse_dw.set_attrs("diagnostic", "water table in surface crevasses",
			 "m", "");

  //ccr-future: m_crevasse_dwdt.create(m_grid, "temporary_crevasse_dwdt", WITH_GHOSTS, 1);
  //ccr-future: m_crevasse_dwdt.set_attrs("diagnostic", "change water table in surface crevasses",
  //ccr-future: 			 "kg m-2", "");

  // Only
  m_surface_h.create(m_grid, "temporary_crevasses_surface_depth",
		     WITH_GHOSTS, 1);
  m_surface_h.set_attrs("diagnostic", "surface crevasses depth",
			"m", "");

  m_bottom_h.create(m_grid, "temporary_crevasses_bottom_depth",
		    WITH_GHOSTS, 1);
  m_bottom_h.set_attrs("diagnostic", "bottom crevasses depth",
		       "m", "");

  m_total_h.create(m_grid, "temporary_total_crevasses_depth",
		   WITH_GHOSTS, 1);
  m_total_h.set_attrs("diagnostic", "total crevasses depth",
		     "m", "");


  // FIXME ccr : command line switch should overwrite the pism_config values
  //bool dw_set = false, dw0_set = false;
  m_dw = m_config->get_double("crevasses_calving_dw");
  m_dw0 = m_config->get_double("crevasses_calving_dw0");
  m_ice_density = m_config->get_double("ice_density");

}

CrevassesCalving::~CrevassesCalving() {
  // empty
}

void CrevassesCalving::init() {

  m_log->message(2,
		 "* Initializing the 'crevasses_calving' mechanism... setting constants...\n"
		 "         standard crevasses water depth dw = %.2f m\n"
		 "         minimal crevasses water depth dw0 = %.2f m\n",
		 m_dw, m_dw0);


  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted("-calving crevasses_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

  m_strain_rates.set(0.0);

  m_surface_h.set(0.0);
  m_bottom_h.set(0.0);
  m_total_h.set(0.0);
  //ccr FIXME Ongoing melting shall increase the value
  //    of m_crevasse_dw while it otherwise decreases
  // NOTE: It would require to read this value during restart !!
  m_crevasse_dw.set(m_dw);
  //ccr-future: m_crevasse_dwdt.set(0.0);


}

//@article{Nicketal2010,
//author = {Nick, FM and van der Veen, CJ and Vieli, A},
//file = {:home/cr/Documents/Mendeley Desktop/2010/Nick, Veen, Vieli - 2010.pdf:pdf},
//journal = {Journal of Glaciology},
//mendeley-groups = {Glacier/Ice sheet/shelf/Ice Shelf,Glacier/Ice sheet/shelf/Calving},
//number = {199},
//pages = {781--794},
//title = {{A physically based calving model applied to marine outlet glaciers and implications for the glacier dynamics}},
//url = {http://www.ingentaconnect.com/content/igsoc/jog/2010/00000056/00000199/art00004},
//volume = {56},
//year = {2010}
//}

//! \brief Uses principal strain rates to estimate crevasses depth to calving, if applicable.
/*!
  See [\ref Nicketal2010].
*/


void CrevassesCalving::update(double dt,
			      IceModelVec2Int &pism_mask,
			      IceModelVec2S &Href,
			      IceModelVec2S &ice_thickness,
			      IceModelVec2S &surface_h,
			      IceModelVec2S &bottom_h,
			      IceModelVec2S &crevasse_dw,
			      IceModelVec2S &crevasse_flux_2D,
			      const double sea_level=0.0) {
                              //ccr-future: IceModelVec2S &crevasse_dwdt=0.0,
                              //ccr-future: const IceModelVec2S &smb=0.0) {

  const IceModelVec2S &bed = *m_grid->variables().get_2d_scalar("bedrock_altitude");
  const IceModelVec3 &T3 = *m_grid->variables().get_3d_scalar("temp");

  // Distance (grid cells) from calving front where strain rate is evaluated
  double crevasses_calving_rate_horizontal;
  const double reciprocal_dt = 1.0/dt;

  m_thk_loss.set(0.0);

  update_strain_rates();

  // Determine the actual water table in the crevasses (dw)
  // Currently we *ONLY* copy the value of 
  update_water_table_crevasses(m_crevasse_dw, crevasse_dw);//ccr-future: , crevasse_dwdt, smb);

  // Update the surface and bottom crevasses depth
  update_crevasses_depths(surface_h, bottom_h, pism_mask,
			  Href, ice_thickness, bed, T3,
			  crevasse_dw, sea_level);


  MaskQuery mask(pism_mask);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(pism_mask);
  list.add(Href);
  list.add(m_thk_loss);
  list.add(m_total_h);
  list.add(crevasse_flux_2D);
  list.add(m_crevasse_dw);
  list.add(crevasse_dw);


  // Compute where crevasses-driven calving occurs along the margin
  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    //double area = m_grid->dx() * m_grid->dy();

    crevasse_flux_2D(i, j) = 0.0; //If not calving, take zero

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.ice_free_ocean(i, j) &&
        mask.next_to_ice(i, j)) {

      // Ice margin with contact to an ocean point
      if (m_total_h(i, j) >= ice_thickness(i, j) + Href(i, j) ) {
	m_thk_loss(i, j) = ice_thickness(i, j) + Href(i, j);
      }
      // end of "if (ice_free_ocean && next_to_floating)"
    } else {
      // Icy margin with contact to an ocean point
      if (mask.ice_margin(i, j) && mask.next_to_ice_free_ocean(i, j)) {
	if (m_total_h(i, j) >= ice_thickness(i, j) + Href(i, j)) {
	  m_thk_loss(i, j) = ice_thickness(i, j) + Href(i, j);
	}
      }
    }
  }

  // FIXME : Calving also beyond the margin, when a continous line
  //         of damages starting/ending at the margin is detected
  //compute_crevasses_calving_inland(m_total_h, ice_thickness, pism_mask)

  m_thk_loss.update_ghosts();

  // Actually apply the omputed crevasses-driven calving.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j) && m_thk_loss(i, j) > 0) {
      crevasse_flux_2D(i, j) = m_thk_loss(i, j) * m_ice_density;
      
      //horizontal calving rate (m/s):dx/dt
      // if (crevasse_flux_2D(i, j) > 0.0) { //fixme Needed??
      crevasses_calving_rate_horizontal = m_grid->dx() * reciprocal_dt;

      Href(i, j) = 0.0;
      ice_thickness(i, j) = 0.0;
      pism_mask(i, j) = MASK_ICE_FREE_OCEAN;
    }
  }

  pism_mask.update_ghosts();

  remove_narrow_tongues(pism_mask, ice_thickness);

  ice_thickness.update_ghosts();
  Href.update_ghosts();
  pism_mask.update_ghosts();
}


/**
 * @brief Compute the maximum time-step length allowed by the CFL
 * condition applied to the crevasses calving rate. To active it set
 * the pism_config the cfl_crevasses_calving="yes";
 *
 * Note: this code uses the mask variable obtained from the Vars
 * dictionary. This is not the same mask that is used in the update()
 * call, since max_timestep() is called *before* the mass-continuity
 * time-step.
 *
 * @return 0 on success
 */
  MaxTimestep CrevassesCalving::max_timestep(const double sea_level=0.0) {

  if (not m_restrict_timestep) {
    return MaxTimestep();
  }

  // About 9 hours which corresponds to 10000 km/year on a 10 km grid
  const double dt_min = units::convert(m_sys, 0.001, "years", "seconds");
  const double assumed_speed_maximum = 1.0e5/(365.0*86400.0); //Unit: meter/second
  const double scaler_of_modeled_speed = 1.015;

  double my_velocity = 0.0,
    my_crevasse_calv_distance_max     = 0.0,
    my_crevasse_calv_distance_counter = 0.0;
    //my_crevasse_calv_distance_mean  = 0.0,

  const IceModelVec2Int &mask = *m_grid->variables().get_2d_mask("mask");
  const IceModelVec2S &Href = *m_grid->variables().get_2d_scalar("Href");
  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("icethickness");
  const IceModelVec2S &bed = *m_grid->variables().get_2d_scalar("bedrock_altitude");
  const IceModelVec3 &T3 = *m_grid->variables().get_3d_scalar("temp");
  const IceModelVec2S &crevasse_dw = *m_grid->variables().get_2d_scalar("crevasses_dw");
  // Take the surface speed, which is usally the highest in vertical columns of ice
  const IceModelVec2S &surf_speed = *m_grid->variables().get_2d_scalar("velsurf_mag");

  //double sea_level0 = *m_ocean->sea_level_elevation(); //FIXME How to access it??

  MaskQuery m(mask);

  update_strain_rates();

  // Consider the last water table in the crevasses (dw)
  update_crevasses_depths(m_surface_h, m_bottom_h, mask,
			  Href, ice_thickness, bed, T3,
			  crevasse_dw, sea_level);

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(m_total_h);
  list.add(Href);
  list.add(ice_thickness);
  list.add(surf_speed);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // find partially filled or empty grid boxes on the ice-free
    // ocean which have ice neighbors
    if ((m.ice_free_ocean(i, j) && m.next_to_ice(i, j)) == false) {
      continue;
    }
    double total_ice_thickness = ice_thickness(i, j) + Href(i, j);

    my_velocity = std::max(my_velocity, surf_speed(i, j));
  
    if (m_total_h(i, j) >= total_ice_thickness) {
      if (m.ice_margin(i, j)) {
	// Fully filled boxes
	double crevasse_calv_distance_horizontal = m_grid->dx();

	my_crevasse_calv_distance_counter += 1.0;
	  //my_crevasse_calv_distance_mean += crevasse_calv_distance_horizontal;
	my_crevasse_calv_distance_max = std::max(my_crevasse_calv_distance_max, crevasse_calv_distance_horizontal);

      } else if (m.next_to_ice(i, j)) {
	// Potentially partly filled boxes
	double fraction_coverage;
	if (ice_thickness(i, j) > 0.0) {
	  fraction_coverage = std::min(Href(i, j) / ice_thickness(i, j), 1.0);
	} else {
	  fraction_coverage = 0.99;
	}
	double crevasse_calv_distance_horizontal = m_grid->dx() * fraction_coverage;
  
	my_crevasse_calv_distance_counter += fraction_coverage;
	  //my_crevasse_calv_distance_mean += crevasse_calv_distance_horizontal;
	my_crevasse_calv_distance_max = std::max(my_crevasse_calv_distance_max, crevasse_calv_distance_horizontal);
      }
    }

     //    if (m.ice_margin(i, j) && 
     //	m_total_h(i, j) >= total_ice_thickness) {
     //      // Fully filled boxes
     //
     //      double crevasse_calv_distance_horizontal = m_grid->dx();
     //
     //      my_crevasse_calv_distance_counter += 1.0;
     //	//my_crevasse_calv_distance_mean += crevasse_calv_distance_horizontal;
     //      my_crevasse_calv_distance_max = std::max(my_crevasse_calv_distance_max, crevasse_calv_distance_horizontal);
     //
     //    } else if (m.next_to_ice(i, j) && ice_thickness(i, j) > 0.0 &&
     //	       m_total_h(i, j) >= total_ice_thickness) {
     //      // Potentially partly filled boxes
     //
     //      double fraction_coverage = std::min(Href(i, j) / ice_thickness(i, j), 1.0);
     //      double crevasse_calv_distance_horizontal = m_grid->dx() * fraction_coverage;
     //  
     //      my_crevasse_calv_distance_counter += fraction_coverage;
     //	//my_crevasse_calv_distance_mean += crevasse_calv_distance_horizontal;
     //      my_crevasse_calv_distance_max = std::max(my_crevasse_calv_distance_max, crevasse_calv_distance_horizontal);
     //    }
  }

  double velocity_max = 0.0,
    crevasse_calv_distance_max = 0.0, //crevasse_calv_distance_mean = 0.0,
    crevasse_calv_distance_counter = 0.0;
    //crevasse_calv_distance_mean    = GlobalSum(m_grid->com, my_crevasse_calv_distance_mean);
  crevasse_calv_distance_counter = GlobalSum(m_grid->com, my_crevasse_calv_distance_counter);
  crevasse_calv_distance_max     = GlobalMax(m_grid->com, my_crevasse_calv_distance_max);


  my_velocity = units::convert(m_sys, my_velocity, "m year-1", "m seconds-1");
  velocity_max = GlobalMax(m_grid->com, my_velocity);
  if (crevasse_calv_distance_counter > 0) {
    velocity_max = 
      std::max(velocity_max*scaler_of_modeled_speed, assumed_speed_maximum);
  }

  double dt = crevasse_calv_distance_max / velocity_max;

  m_log->message(2,
             "  crevassescalving: max velocity = %.2f m/a with distances = %.1f km gives dt=%.5f a; %d cells\n",
             units::convert(m_sys, velocity_max, "m/s", "m/year"),
             units::convert(m_sys, crevasse_calv_distance_max, "meter", "kilometer"),
             units::convert(m_sys, dt, "seconds", "years"),
             (int)crevasse_calv_distance_counter);

  return MaxTimestep(std::max(dt, dt_min));
}

void CrevassesCalving::add_vars_to_output_impl(const std::string &/*keyword*/, std::set<std::string> &/*result*/) {
  // empty
}

void CrevassesCalving::define_variables_impl(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                                  IO_Type /*nctype*/) {
  // empty
}

void CrevassesCalving::write_variables_impl(const std::set<std::string> &/*vars*/, const PIO& /*nc*/) {
  // empty
}

/**
 * Update the strain rates field.
 *
 * Note: this code uses the mask obtained from the Vars
 * dictionary, because the velocity field used to compute it need not
 * extend past the ice margin corresponding to the *beginning* of the
 * time-step.
 *
 * @return 0 on success
 */
void CrevassesCalving::update_strain_rates() {
  const IceModelVec2V &ssa_velocity = m_stress_balance->advective_velocity();
  const IceModelVec2Int &mask = *m_grid->variables().get_2d_mask("mask");
  m_stress_balance->compute_2D_principal_strain_rates(ssa_velocity, mask, m_strain_rates);
}

//ccr -------------------
/**
 * Update the water table in the crevasses
 *
 * Note: this code uses the mask obtained from the Vars
 * dictionary, because the velocity field used to compute it need not
 * extend past the ice margin corresponding to the *beginning* of the
 * time-step.
 * Note: Currently we only copy values of the field m_crevasses_dw
 *
 * @return 0 on success
 */
void CrevassesCalving::update_water_table_crevasses(IceModelVec2S &crevasse_dw_input,
						    IceModelVec2S &crevasse_dw_return) {
                                                    //ccr-future: IceModelVec2S &crevasse_dwdt,
                                                    //ccr-future: const IceModelVec2S &smb) {

  crevasse_dw_return.copy_from(crevasse_dw_input);
  
  //IceModelVec::AccessList list;
  //list.add(crevasse_dw_input);
  //list.add(crevasse_dw_return);
  ////ccr-future: list.add(crevasse_dwdt);
  //
  //double factor = exp(-dt/dt_factor);
  //const double factor = 1.0;
  //for (Points p(*m_grid); p; p.next()) {
  //  const int i = p.i(), j = p.j();
  //  crevasse_dw_return(i,j) = crevasse_dw_input(i,j) * factor; //ccr-future:  + crevasse_dwdt(i, j);
  //  crevasse_dw_return(i,j) = std::max(crevasse_dw_return(i,j),  m_dw0);
  //}
}


/**
 * Update the crevasses depth fields.
 *
 * Note: this code uses the mask obtained from the Vars
 * dictionary, because the velocity field used to compute it need not
 * extend past the ice margin corresponding to the *beginning* of the
 * time-step.
 *
 * @return 0 on success
 */
void CrevassesCalving::
  update_crevasses_depths(IceModelVec2S &surface_h,
			  IceModelVec2S &bottom_h,
			  const IceModelVec2Int &pism_mask,
			  const IceModelVec2S &Href,
			  const IceModelVec2S &ice_thickness,
			  const IceModelVec2S &bed,
			  const IceModelVec3 &T3,
			  const IceModelVec2S &crevasse_dw,
			  const double sea_level=0.0) {//fixme: future

  //int offset = m_stencil_width; //ccr fixme Needed??

  MaskQuery mask(pism_mask);

  const double A_cold = m_config->get_double("Paterson_Budd_A_cold");
  const double A_warm = m_config->get_double("Paterson_Budd_A_warm");
  const double Q_cold = m_config->get_double("Paterson_Budd_Q_cold");
  const double Q_warm = m_config->get_double("Paterson_Budd_Q_warm");
  const double crit_temp = m_config->get_double("Paterson_Budd_critical_temperature");
  const double R = m_config->get_double("ideal_gas_constant");
  const double p_air = m_config->get_double("surface_pressure");
  const double &ice_density = m_ice_density; //m_config->get_double("ice_density");
  const double fresh_water_density = m_config->get_double("fresh_water_density");
  const double sea_water_density = m_config->get_double("sea_water_density");
  const double standard_gravity = m_config->get_double("standard_gravity");
  const double Glen_n_ssa = m_config->get_double("ssa_Glen_exponent");
  const double beta = m_config->get_double("beta_CC");

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(pism_mask);
  list.add(Href);
  list.add(m_strain_rates);
  list.add(m_thk_loss);
  list.add(crevasse_dw);
  list.add(surface_h);
  list.add(bottom_h);
  list.add(m_total_h);
  list.add(T3);
  list.add(bed);

  // Compute required fields
  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    // Crevasses specific values
    double A, Re0, Re1, R_iceg;

    // Compute fields 
    m_total_h(i, j) = 0.0;
    surface_h(i, j) = 0.0;
    bottom_h(i, j) = 0.0;

    //if (mask.ice_margin(i, j)) {
    if (mask.icy(i, j)) {
      unsigned int k = 0; 
      double Hab = 0.0;
      double depth = sea_level - bed(i, j);
      double Tice = T3(i, j, k);

      if (ice_density*Href(i, j) < depth) {
	// floating ice: ice thickness below sea level is less than bedroch depth
	Hab = 0.0;
      } else {
	// ice thickness below water exceeds bedrock depth
	Hab = Href(i, j) - sea_water_density * depth;
      }
      double p = p_air + ice_density * standard_gravity * Href(i, j);
      double T_pa = Tice - beta * p;
      
      if (T_pa < crit_temp) {
	A = A_cold * exp(-Q_cold/(R*T_pa));
      } else {
	A = A_warm * exp(-Q_warm/(R*T_pa));
      }


      //?? Re0 = 2.0 * pow(( fabs(m_strain_rates(i, j, 0)) /A),(1.0/Glen_n_ssa)); // principal component 0
      //?? Re1 = 2.0 * pow(( fabs(m_strain_rates(i, j, 1)) /A),(1.0/Glen_n_ssa)); // principal component 1
      
      Re0 = 2.0 * pow(( std::max(m_strain_rates(i, j, 0),0.0) /A),(1.0/Glen_n_ssa)); // principal component 0
      Re1 = 2.0 * pow(( std::max(m_strain_rates(i, j, 1),0.0) /A),(1.0/Glen_n_ssa)); // principal component 1

      
      R_iceg = sqrt(Re0*Re0+Re1*Re1)/(ice_density/standard_gravity);
      
      surface_h(i, j) = 
	R_iceg + fresh_water_density/ice_density * crevasse_dw(i, j);
      bottom_h(i, j) = // std::max(, 0) avoids unphysical negative values
	(1.0/(sea_water_density-1.0)) * std::max(R_iceg-Hab, 0.0);
      m_total_h(i, j) = surface_h(i, j) + bottom_h(i, j);

      // Box height weighted temperature (not used)
      //   const unsigned int ks = m_grid->kBelowHeight(thk(i, j));
      //     Hab = (m_grid->z(k) - m_grid->z(k-1));
      //     Tice = T3(i, j, k) * Hab;
      //     for (k=0; k< m_grid->Mz()-1; ++k){
      //   //for (k=0; k<= ks; ++k){
      //      double dz0 = m_grid->z(k) - m_grid->z(k-1);
      //      double dz1 = m_grid->z(k+1) - m_grid->z(k);
      //      Hab = Hab + dz0; 
      //      // Check if Href is the correct variable // FIXMEE
      //      if ( Href(i, j) < Hab+dz1 ) break;  // leave the loop, when the integrated height plus the next layer exceeds the ice thickness
      //     }
    }
  }
  // Update ghosted fields
  surface_h.update_ghosts();
  bottom_h.update_ghosts();
  m_total_h.update_ghosts();

}
//ccr ------------------


/** Remove tips of one-cell-wide ice tongues ("noses").  Changes ice thickness.
 *
 * The center icy cell in ice tongues like this one (and equivalent)
 *
 * @code
   O O ?
   X X O
   O O ?
   @endcode
 *
 * where "O" is ice-free and "?" is any mask value, are removed.
 * Ice tongues like this one
 *
 * @code
   # O ?
   X X O
   # O ?
   @endcode
 * where one or two of the "#" cells are ice-filled, are not removed.
 *
 * See the code for the precise rule, which uses `ice_free_ocean()` for the "O"
 * cells if the center cell has grounded ice, and uses `ice_free()` if the
 * center cell has floating ice.
 *
 * @note We use `pism_mask` (and not ice_thickness) to make decisions.
 * This means that we can update `ice_thickness` in place without
 * introducing a dependence on the grid traversal order.
 *
 * @param[in,out] pism_mask cell type mask
 * @param[in,out] ice_thickness modeled ice thickness
 *
 * @return 0 on success
 */
void CrevassesCalving::remove_narrow_tongues(IceModelVec2Int &pism_mask,
					     IceModelVec2S &ice_thickness) {
  MaskQuery mask(pism_mask);

  IceModelVec::AccessList list;
  list.add(pism_mask);
  list.add(ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (mask.ice_free(i, j)) {
      // FIXME: it might be better to have access to bedrock elevation b(i,j)
      // and sea level SL so that the predicate can be
      //   mask.ice_free(i,j) || (mask.grounded_ice(i,j) && (b(i,j) >= SL)))
      continue;
    }

    bool ice_free_N = false,  ice_free_E = false,
      ice_free_S = false, ice_free_W = false,
      ice_free_NE = false, ice_free_NW = false,
      ice_free_SE = false, ice_free_SW = false;

    if (mask.grounded_ice(i,j)) {
      // if (i,j) is grounded ice then we will remove it if it has
      // exclusively ice-free ocean neighbors
      ice_free_N  = mask.ice_free_ocean(i, j + 1);
      ice_free_E  = mask.ice_free_ocean(i + 1, j);
      ice_free_S  = mask.ice_free_ocean(i, j - 1);
      ice_free_W  = mask.ice_free_ocean(i - 1, j);
      ice_free_NE = mask.ice_free_ocean(i + 1, j + 1);
      ice_free_NW = mask.ice_free_ocean(i - 1, j + 1);
      ice_free_SE = mask.ice_free_ocean(i + 1, j - 1);
      ice_free_SW = mask.ice_free_ocean(i - 1, j - 1);
    } else if (mask.floating_ice(i,j)) {
      // if (i,j) is floating then we will remove it if its neighbors are
      // ice-free, whether ice-free ocean or ice-free ground
      ice_free_N  = mask.ice_free(i, j + 1);
      ice_free_E  = mask.ice_free(i + 1, j);
      ice_free_S  = mask.ice_free(i, j - 1);
      ice_free_W  = mask.ice_free(i - 1, j);
      ice_free_NE = mask.ice_free(i + 1, j + 1);
      ice_free_NW = mask.ice_free(i - 1, j + 1);
      ice_free_SE = mask.ice_free(i + 1, j - 1);
      ice_free_SW = mask.ice_free(i - 1, j - 1);
    }

    if ((ice_free_W == false &&
         ice_free_NW         &&
         ice_free_SW         &&
         ice_free_N          &&
         ice_free_S          &&
         ice_free_E)         ||
        (ice_free_N == false &&
         ice_free_NW         &&
         ice_free_NE         &&
         ice_free_W          &&
         ice_free_E          &&
         ice_free_S)         ||
        (ice_free_E == false &&
         ice_free_NE         &&
         ice_free_SE         &&
         ice_free_W          &&
         ice_free_S          &&
         ice_free_N)         ||
        (ice_free_S == false &&
         ice_free_SW         &&
         ice_free_SE         &&
         ice_free_W          &&
         ice_free_E          &&
         ice_free_N)) {
      pism_mask(i, j) = MASK_ICE_FREE_OCEAN;
      ice_thickness(i, j) = 0.0;
    }
  }
}

} // end of namespace calving
} // end of namespace pism
