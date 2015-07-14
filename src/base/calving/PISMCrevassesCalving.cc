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

//#include "base/iceModel.hh" //ccr needed?? fixme
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

  m_restrict_timestep = m_config->get_boolean("cfl_eigen_calving"); //ccr FIXME

  m_K = m_config->get_double("eigen_calving_K"); //ccr FIXME

  //ccr -- here crevasses calving specific variables
  m_crevasse_dw.create(m_grid, "temporary_crevasse_dw", WITH_GHOSTS, 1);
  m_crevasse_dw.set_attrs("diagnostic", "water table in surface crevasses",
			 "m", "");

  // Only
  m_surface_h.create(m_grid, "temporary_surface_crevasses_depth",
		     WITH_GHOSTS, 1);
  m_surface_h.set_attrs("diagnostic", "surface crevasses depth",
			"m", "");

  m_bottom_h.create(m_grid, "temporary_bottom_crevasse_depth",
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


  // initialiaztion of fields
  IceModelVec::AccessList list;
  list.add(m_surface_h);
  list.add(m_bottom_h);
  list.add(m_total_h);
  list.add(m_crevasse_dw);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    m_surface_h(i, j) = 0.0;
    m_bottom_h(i, j) = 0.0;
    m_total_h(i, j) = 0.0;
    //ccr FIXME Ongoing melting shall increase the value
    //    of m_crevasse_dw while it otherwise decreases
    m_crevasse_dw(i, j) = m_dw;
  }
}

CrevassesCalving::~CrevassesCalving() {
  // empty
}

void CrevassesCalving::init() {

  m_log->message(2,
             "* Initializing the 'eigen-calving' mechanism...\n");

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted("-calving crevasses_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

  m_strain_rates.set(0.0);

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
			      const IceModelVec2S &bed,
			      const IceModelVec3 &T3,
			      const double sea_level=0.0) {

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width; //fixme Needed??
  double crevasses_calving_rate_horizontal;
  const double reciprocal_dt = 1.0/dt;

  m_thk_loss.set(0.0);

  update_strain_rates();

  // Determine the actual water table in the crevasses (dw)
  update_water_table_crevasses(m_crevasse_dw, crevasse_dw);

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

  //->update_crevasses_depths  list.add(m_strain_rates);
  //->update_crevasses_depths  list.add(T3);
  //->update_crevasses_depths  list.add(bed);
  //->update_crevasses_depths  list.add(m_surface_h);
  //->update_crevasses_depths  list.add(m_bottom_h);
  //->update_crevasses_depths  list.add(m_crevasse_dw); //ccr
  
  //
  //->update_crevasses_depths  // Compute required fields
  //->update_crevasses_depths  for (Points pt(*m_grid); pt; pt.next()) {
  //->update_crevasses_depths    const int i = pt.i(), j = pt.j();
  //->update_crevasses_depths    // Crevasses specific values
  //->update_crevasses_depths    double
  //->update_crevasses_depths      A, Re0, Re1, R_iceg;
  //->update_crevasses_depths
  //->update_crevasses_depths    // Compute fields 
  //->update_crevasses_depths    m_surface_h(i, j) = 0.0;
  //->update_crevasses_depths    m_bottom_h(i, j) = 0.0;
  //->update_crevasses_depths    m_total_h(i, j) = 0.0;
  //->update_crevasses_depths    //if (mask.ice_margin(i, j)) {
  //->update_crevasses_depths    if (mask.icy(i, j)) {
  //->update_crevasses_depths      double Tice, Hab = 0.0;
  //->update_crevasses_depths      unsigned int k=0; 
  //->update_crevasses_depths      double beta=m_config->get_double("beta_CC");
  //->update_crevasses_depths      double depth = sea_level - bed(i, j);
  //->update_crevasses_depths
  //->update_crevasses_depths      Tice = T3(i, j, k);
  //->update_crevasses_depths      // Box height weighted temperature (not used)
  //->update_crevasses_depths      //   const unsigned int ks = m_grid->kBelowHeight(thk(i, j));
  //->update_crevasses_depths      //     Hab = (m_grid->z(k) - m_grid->z(k-1));
  //->update_crevasses_depths      //     Tice = T3(i, j, k) * Hab;
  //->update_crevasses_depths      //     for (k=0; k< m_grid->Mz()-1; ++k){
  //->update_crevasses_depths      //     //for (k=0; k<= ks; ++k){
  //->update_crevasses_depths      //      //Test double dz0=10;
  //->update_crevasses_depths      ////old code; might be not correct
  //->update_crevasses_depths      ////    double dz0 = m_grid->z(k);
  //->update_crevasses_depths      ////    double dz1 = m_grid->z(k+1);
  //->update_crevasses_depths      ////fixme This should be correct or?? 
  //->update_crevasses_depths      //      double dz0 = m_grid->z(k) - m_grid->z(k-1);
  //->update_crevasses_depths      //      double dz1 = m_grid->z(k+1) - m_grid->z(k);
  //->update_crevasses_depths      //      Hab = Hab + dz0; 
  //->update_crevasses_depths      //      // Check if Href is the correct variable // FIXMEE
  //->update_crevasses_depths      //      if ( Href(i, j) < Hab+dz1 ) break;  // leave the loop, when the integrated height plus the next layer exceeds the ice thickness
  //->update_crevasses_depths      //     }
  //->update_crevasses_depths      if (m_ice_density*Href(i, j) < depth) {
  //->update_crevasses_depths	// floating ice: ice thickness below sea level is less than bedroch depth
  //->update_crevasses_depths	Hab = 0.0;
  //->update_crevasses_depths      } else {
  //->update_crevasses_depths	// ice thickness below water exceeds bedrock depth
  //->update_crevasses_depths	Hab = Href(i, j) - m_sea_water_density * depth;
  //->update_crevasses_depths      }
  //->update_crevasses_depths      double p = m_p_air + m_ice_density * m_standard_gravity * Href(i, j);
  //->update_crevasses_depths      double T_pa = Tice - beta * p;
  //->update_crevasses_depths      
  //->update_crevasses_depths      if (T_pa < m_crit_temp) {
  //->update_crevasses_depths	A = m_A_cold * exp(-m_Q_cold/(m_R*T_pa));
  //->update_crevasses_depths      } else {
  //->update_crevasses_depths	A = m_A_warm * exp(-m_Q_warm/(m_R*T_pa));
  //->update_crevasses_depths      }
  //->update_crevasses_depths      
  //->update_crevasses_depths      Re0 = 2.0 * pow(( std::max(m_strain_rates(i, j, 0),0.0) /A),(1.0/m_Glen_n_ssa)); // principal component 0, 1st eigenvalue
  //->update_crevasses_depths      Re1 = 2.0 * pow(( std::max(m_strain_rates(i, j, 1),0.0) /A),(1.0/m_Glen_n_ssa)); // principal component 1, 2nd eigenvalue
  //->update_crevasses_depths      
  //->update_crevasses_depths      R_iceg = sqrt(Re0*Re0+Re1*Re1)/(m_ice_density/m_standard_gravity);
  //->update_crevasses_depths      
  //->update_crevasses_depths      m_surface_h(i, j) = R_iceg + m_fresh_water_density/m_ice_density * m_crevasse_dw(i, j);
  //->update_crevasses_depths      m_bottom_h(i, j) = (1.0/(m_sea_water_density-1.0)) * std::max(R_iceg-Hab, 0.0); // std::max( , 0.) avoids unphysical negative values
  //->update_crevasses_depths      m_total_h(i, j) = surface_h(i, j) + bottom_h(i, j);
  //->update_crevasses_depths    }
  //->update_crevasses_depths  }
  

  // Compute where crevasses-driven calving occurs along the margin
  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    // Neighbor-averaged ice thickness:
    double H_average = 0.0;
    //double area = m_grid->dx() * m_grid->dy();
    int N_floating_neighbors = 0, M = 0;
    // M is the number of cells used to compute averages of strain
    // rate components.

    crevasse_flux_2D(i, j) = 0.0; //If not calving, take zero

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.ice_free_ocean(i, j) &&
        mask.next_to_ice(i, j)) {

      //if (mask.floating_ice(i + 1, j)) {
      //  N_floating_neighbors += 1;
      //  H_average += ice_thickness(i + 1, j);
      //}
      //
      //if (mask.floating_ice(i - 1, j)) {
      //  N_floating_neighbors += 1;
      //  H_average += ice_thickness(i - 1, j);
      //}
      //
      //if (mask.floating_ice(i, j + 1)) {
      //  N_floating_neighbors += 1;
      //  H_average += ice_thickness(i, j + 1);
      //}
      //
      //if (mask.floating_ice(i, j - 1)) {
      //  N_floating_neighbors += 1;
      //  H_average += ice_thickness(i, j - 1);
      //}
      //
      //if (N_floating_neighbors > 0) {
      //  H_average /= N_floating_neighbors;
      //}

      // Ice margin with contact to an ocean point
      if (m_total_h(i, j) >= ice_thickness(i, j) + Href(i, j) ) {
	m_thk_loss(i, j) = ice_thickness(i, j) + Href(i, j);
      }

      //// apply calving rate at partially filled or empty grid cells
      //if (crevasses_calving_rate_horizontal > 0.0) {
      //	Href(i, j) -= crevasses_calving_rate_horizontal * dt; // in m
      //	
      //	if (Href(i, j) < 0.0) {
      //	  // Partially filled grid cell became ice-free
      //	  
      //	  m_thk_loss(i, j) = -Href(i, j); // in m, corresponds to additional ice loss
      //	  Href(i, j)       = 0.0;
      //	  
      //	  // additional mass loss will be distributed among
      //	  // N_floating_neighbors:
      //	  if (N_floating_neighbors > 0) {
      //	    m_thk_loss(i, j) /= N_floating_neighbors;
      //	  }
      //	}
      //}

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

    //old double thk_loss_ij = 0.0;
    //old 
    //old if (mask.floating_ice(i, j) &&
    //old     (m_thk_loss(i + 1, j) > 0.0 || m_thk_loss(i - 1, j) > 0.0 ||
    //old      m_thk_loss(i, j + 1) > 0.0 || m_thk_loss(i, j - 1) > 0.0)) {
    //old 
    //old   thk_loss_ij = (m_thk_loss(i + 1, j) + m_thk_loss(i - 1, j) +
    //old                  m_thk_loss(i, j + 1) + m_thk_loss(i, j - 1));     // in m/(dt=s)
    //old 
    //old   // Note std::max: we do not account for further calving
    //old   // ice-inwards! Alternatively CFL criterion for time stepping
    //old   // could be adjusted to maximum of calving rate
    //old   Href(i, j) = std::max(ice_thickness(i, j) - thk_loss_ij, 0.0); // in m
    //old 
    //old   ice_thickness(i, j) = 0.0;
    //old   pism_mask(i, j) = MASK_ICE_FREE_OCEAN;
    //old }
  }


  pism_mask.update_ghosts();

  remove_narrow_tongues(pism_mask, ice_thickness);

  ice_thickness.update_ghosts();
  Href.update_ghosts();
  pism_mask.update_ghosts();
}


/**
 * @brief Compute the maximum time-step length allowed by the CFL
 * condition applied to the crevasses calving rate.
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
  double dt_min = units::convert(m_sys, 0.001, "years", "seconds");

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;
  double my_velocity = 0.0,
    my_crevasse_calv_distance_max     = 0.0,
    my_crevasse_calv_distance_counter = 0.0;
    //my_crevasse_calv_distance_mean  = 0.0,

  const double assumed_speed_maximum = 1.0e5/(365.0*86400.0); //Unit: meter/second
  const double scaler_of_modeled_speed = 1.015;

  const IceModelVec2Int &mask = *m_grid->variables().get_2d_mask("mask");
  const IceModelVec2S &Href = *m_grid->variables().get_2d_scalar("Href");
  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("icethickness");
  const IceModelVec2S &bed = *m_grid->variables().get_2d_scalar("bedrock_altitude");
  const IceModelVec3 &T3 = *m_grid->variables().get_3d_scalar("temp");
  const IceModelVec2S &crevasse_dw = *m_grid->variables().get_2d_scalar("crevasses_dw");
  // Take the surface speed, which is usally the highest in vertical columns of ice
  const IceModelVec2S &surf_speed = *m_grid->variables().get_2d_scalar("velsurf_mag");

  //double sea_level = ocean->sea_level_elevation(); //FIXME How to access it??

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
  
    if (m.ice_margin(i, j) && 
	m_total_h(i, j) >= total_ice_thickness) {
      // Fully filled boxes

      double crevasse_calv_distance_horizontal = m_grid->dx();

      my_crevasse_calv_distance_counter += 1.0;
	//my_crevasse_calv_distance_mean += crevasse_calv_distance_horizontal;
      my_crevasse_calv_distance_max = std::max(my_crevasse_calv_distance_max, crevasse_calv_distance_horizontal);

    } else if (m.next_to_ice(i, j) && ice_thickness(i, j) > 0.0 &&
	       m_total_h(i, j) >= total_ice_thickness) {
      // Potentially partly filled boxes

      double fraction_coverage = std::min(Href(i, j) / ice_thickness(i, j), 1.0);
      double crevasse_calv_distance_horizontal = m_grid->dx() * fraction_coverage;
  
      my_crevasse_calv_distance_counter += fraction_coverage;
	//my_crevasse_calv_distance_mean += crevasse_calv_distance_horizontal;
      my_crevasse_calv_distance_max = std::max(my_crevasse_calv_distance_max, crevasse_calv_distance_horizontal);
    }
  }

  double velocity_max = 0.0,
    crevasse_calv_distance_max = 0.0, //crevasse_calv_distance_mean = 0.0,
    crevasse_calv_distance_counter = 0.0;
    //crevasse_calv_distance_mean    = GlobalSum(m_grid->com, my_crevasse_calv_distance_mean);
  crevasse_calv_distance_counter = GlobalSum(m_grid->com, my_crevasse_calv_distance_counter);
  crevasse_calv_distance_max     = GlobalMax(m_grid->com, my_crevasse_calv_distance_max);


  my_velocity = units::convert(m_sys, my_velocity, "m year-1", "m seconds-1");
  velocity_max = GlobalMax(m_grid->com, my_velocity);
  velocity_max = std::max(velocity_max*scaler_of_modeled_speed, assumed_speed_maximum);

     //if (crevasse_calv_distance_counter > 0.0) {
     //  crevasse_calv_distance_mean /= crevasse_calv_distance_counter;
     //} else {
     //  crevasse_calv_distance_mean = 0.0;
     //}

  double denom = crevasse_calv_distance_max / m_grid->dx();
  const double epsilon = units::convert(m_sys, 0.001 / (m_grid->dx() + m_grid->dy()), "seconds", "years");

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
void CrevassesCalving::
  update_water_table_crevasses(IceModelVec2S &m_crevasse_dw_in,
			       IceModelVec2S &crevasse_dw) {

  IceModelVec::AccessList list;
  list.add(m_crevasse_dw_in);
  list.add(crevasse_dw);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    crevasse_dw(i,j) = m_crevasse_dw_in(i,j);
  }
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

  const double A_cold = m_config->get_double("Paterson-Budd_A_cold");
  const double A_warm = m_config->get_double("Paterson-Budd_A_warm");
  const double Q_cold = m_config->get_double("Paterson-Budd_Q_cold");
  const double Q_warm = m_config->get_double("Paterson-Budd_Q_warm");
  const double crit_temp = m_config->get_double("Paterson-Budd_critical_temperature");
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
  list.add(crevasse_dw);//list.add(m_crevasse_dw);
  list.add(surface_h);  //list.add(m_surface_h);
  list.add(bottom_h);   //list.add(m_bottom_h);
  list.add(m_total_h);
  list.add(T3);
  list.add(bed);

  // Compute required fields
  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    // Crevasses specific values
    double A, Re0, Re1, R_iceg;

    // Compute fields 
    surface_h(i, j) = 0.0; //m_surface_h(i, j) = 0.0;
    bottom_h(i, j) = 0.0;  //m_bottom_h(i, j) = 0.0;
    m_total_h(i, j) = 0.0;
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


      //?? Re0 = 2.0 * pow(( abs(m_strain_rates(i, j, 0)) /A),(1.0/Glen_n_ssa)); // principal component 0
      //?? Re1 = 2.0 * pow(( abs(m_strain_rates(i, j, 1)) /A),(1.0/Glen_n_ssa)); // principal component 1
      
      Re0 = 2.0 * pow(( std::max(m_strain_rates(i, j, 0),0.0) /A),(1.0/Glen_n_ssa)); // principal component 0
      Re1 = 2.0 * pow(( std::max(m_strain_rates(i, j, 1),0.0) /A),(1.0/Glen_n_ssa)); // principal component 1


      
      R_iceg = sqrt(Re0*Re0+Re1*Re1)/(ice_density/standard_gravity);
      
      surface_h(i, j) = R_iceg + fresh_water_density/ice_density * crevasse_dw(i, j); //m_surface_h(i, j) = 
      bottom_h(i, j) = (1.0/(sea_water_density-1.0)) * std::max(R_iceg-Hab, 0.0); //m_bottom_h(i, j) = // std::max( , 0.) avoids unphysical negative values
      m_total_h(i, j) = surface_h(i, j) + bottom_h(i, j); //m_surface_h(i, j) + m_bottom_h(i, j);

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
//ccr-m_      if (m_ice_density*Href(i, j) < depth) {
//ccr-m_	// floating ice: ice thickness below sea level is less than bedroch depth
//ccr-m_	Hab = 0.0;
//ccr-m_      } else {
//ccr-m_	// ice thickness below water exceeds bedrock depth
//ccr-m_	Hab = Href(i, j) - m_sea_water_density * depth;
//ccr-m_      }
//ccr-m_      double p = m_p_air + m_ice_density * m_standard_gravity * Href(i, j);
//ccr-m_      double T_pa = Tice - beta * p;
//ccr-m_      
//ccr-m_      if (T_pa < m_crit_temp) {
//ccr-m_	A = m_A_cold * exp(-m_Q_cold/(m_R*T_pa));
//ccr-m_      } else {
//ccr-m_	A = m_A_warm * exp(-m_Q_warm/(m_R*T_pa));
//ccr-m_      }
//ccr-m_      
//ccr-m_      Re0 = 2.0 * pow(( std::max(m_strain_rates(i, j, 0),0.0) /A),(1.0/m_Glen_n_ssa)); // principal component 0, 1st eigenvalue
//ccr-m_      Re1 = 2.0 * pow(( std::max(m_strain_rates(i, j, 1),0.0) /A),(1.0/m_Glen_n_ssa)); // principal component 1, 2nd eigenvalue
//ccr-m_      
//ccr-m_      R_iceg = sqrt(Re0*Re0+Re1*Re1)/(m_ice_density/m_standard_gravity);
//ccr-m_      
//ccr-m_      m_surface_h(i, j) = R_iceg + m_fresh_water_density/m_ice_density * crevasse_dw(i, j);
//ccr-m_      m_bottom_h(i, j) = (1.0/(m_sea_water_density-1.0)) * std::max(R_iceg-Hab, 0.0); // std::max( , 0.) avoids unphysical negative values
//ccr-m_      m_total_h(i, j) = m_surface_h(i, j) + m_bottom_h(i, j);
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
