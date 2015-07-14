// Copyright (C) 2004-2015 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __iceModel_hh
#define __iceModel_hh

//! \file iceModel.hh Definition of class IceModel.
/*! \file iceModel.hh
  IceModel is a big class which is an ice flow model.  It contains all parts that
  are not well-defined, separated components.  Such components are better places
  to put sub-models that have a clear, general interface to the rest of an ice
  sheet model.

  IceModel has pointers to well-defined components, when they exist.

  IceModel generally interprets user options, and initializes components based on
  such options.  It manages the initialization sequences (%e.g. a restart from a
  file containing a complete model state, versus bootstrapping).
*/

#include <map>
#include <set>
#include <string>
#include <vector>

// IceModel owns a bunch of fields, so we have to include this.
#include "base/util/iceModelVec.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/Context.hh"
#include "base/util/Logger.hh"
#include "base/util/PISMTime.hh"

namespace pism {

class MaxTimestep;

namespace ocean {
class OceanModel;
}

namespace surface {
class SurfaceModel;
}

namespace stressbalance {
class StressBalance;
}

namespace hydrology {
class Hydrology;
}

namespace calving {
class EigenCalving;
class CrevassesCalving;
class OceanKill;
class FloatKill;
class CalvingAtThickness;
class IcebergRemover;
}

namespace energy {
class BedThermalUnit;
}

namespace bed {
class BedDef;
}

// forward declarations
class IceGrid;
class YieldStress;
class Diagnostic;
class TSDiagnostic;

//! The base class for PISM.  Contains all essential variables, parameters, and flags for modelling an ice sheet.
class IceModel {
  // The following classes implement various diagnostic computations.
  // 2D and 3D:
  friend class IceModel_hardav;
  friend class IceModel_bwp;
  friend class IceModel_cts;
  friend class IceModel_dhdt;
  friend class IceModel_enthalpybase;
  friend class IceModel_enthalpysurf;
  friend class IceModel_tempbase;
  friend class IceModel_liqfrac;
  friend class IceModel_tempicethk;
  friend class IceModel_tempicethk_basal;
  friend class IceModel_new_mask;
  friend class IceModel_climatic_mass_balance_cumulative;
  friend class IceModel_dHdt;
  friend class IceModel_flux_divergence;
  // scalar:
  friend class IceModel_ivol;
  friend class IceModel_slvol;
  friend class IceModel_divoldt;
  friend class IceModel_iarea;
  friend class IceModel_imass;
  friend class IceModel_dimassdt;
  friend class IceModel_ivoltemp;
  friend class IceModel_ivolcold;
  friend class IceModel_ivolg;
  friend class IceModel_ivolf;
  friend class IceModel_iareatemp;
  friend class IceModel_iareacold;
  friend class IceModel_ienthalpy;
  friend class IceModel_iareag;
  friend class IceModel_iareaf;
  friend class IceModel_dt;
  friend class IceModel_max_diffusivity;
  friend class IceModel_surface_flux;
  friend class IceModel_surface_flux_cumulative;
  friend class IceModel_grounded_basal_flux;
  friend class IceModel_grounded_basal_flux_cumulative;
  friend class IceModel_sub_shelf_flux;
  friend class IceModel_sub_shelf_flux_cumulative;
  friend class IceModel_nonneg_flux;
  friend class IceModel_nonneg_flux_cumulative;
  friend class IceModel_discharge_flux;
  friend class IceModel_discharge_flux_cumulative;
  friend class IceModel_nonneg_flux_2D_cumulative;
  friend class IceModel_grounded_basal_flux_2D_cumulative;
  friend class IceModel_floating_basal_flux_2D_cumulative;
  friend class IceModel_discharge_flux_2D_cumulative;
  friend class IceModel_max_hor_vel;
  friend class IceModel_sum_divQ_flux;
  friend class IceModel_H_to_Href_flux;
  friend class IceModel_Href_to_H_flux;
  //ccr -- begin
  friend class IceModel_land_flux;                //ccr 1D, snapshot
  friend class IceModel_land_flux_cumulative;     //ccr 1D, cumulative
  friend class IceModel_land_flux_2D;             //ccr 2D, snapshot
  friend class IceModel_land_flux_2D_cumulative;  //ccr 2D, cumulative
  friend class IceModel_ocean_flux;	          //ccr 1D, snapshot
  friend class IceModel_ocean_flux_cumulative;    //ccr 1D, cumulative
  friend class IceModel_ocean_flux_2D;	          //ccr 2D, snapshot
  friend class IceModel_ocean_flux_2D_cumulative; //ccr 2D, cumulative
  //ccr -- crevasses calving
  friend class IceModel_crevasses_calv_flux;
  friend class IceModel_crevasses_calv_flux_cumulative;
  friend class IceModel_crevasses_calv_flux_2D;
  friend class IceModel_crevasses_calv_flux_2D_cumulative;
  friend class IceModel_crevasses_dw;
  //ccr -- end
public:
  // see iceModel.cc for implementation of constructor and destructor:
  IceModel(IceGrid::Ptr g, Context::Ptr context);
  virtual ~IceModel(); // must be virtual merely because some members are virtual

  // see iMinit.cc
  virtual void time_setup();

  IceGrid::Ptr grid() const;

  Context::Ptr ctx() const;

  virtual void allocate_submodels();
  virtual void allocate_stressbalance();
  virtual void allocate_bed_deformation();
  virtual void allocate_bedrock_thermal_unit();
  virtual void allocate_subglacial_hydrology();
  virtual void allocate_basal_yield_stress();
  virtual void allocate_couplers();
  virtual void allocate_iceberg_remover();

  virtual void init_couplers();
  virtual void init_step_couplers();
  virtual void model_state_setup();
  virtual void set_vars_from_options();
  virtual void allocate_internal_objects();
  virtual void misc_setup();
  virtual void init_diagnostics();
  virtual void init_calving();

  virtual void list_diagnostics();

  // see iceModel.cc
  void init();


  /** Run PISM in the "standalone" mode. */
  virtual void run();
  /** Advance the current PISM run to a specific time */
  virtual void run_to(double time);
  virtual void step(bool do_mass_continuity, bool do_energy, bool do_age, bool do_skip);
  void reset_counters();

  // see iMbootstrap.cc
  virtual void bootstrapFromFile(const std::string &fname);
  virtual void bootstrap_2d(const std::string &fname);
  virtual void putTempAtDepth();

  // see iMoptions.cc
  virtual void setFromOptions();
  virtual std::set<std::string> output_size_from_option(const std::string &option,
                                                        const std::string &description,
                                                        const std::string &default_value);
  virtual std::set<std::string> set_output_size(const std::string &keyword);
  virtual std::string get_output_size(const std::string &option);

  // see iMutil.cc
  virtual void additionalAtStartTimestep();
  virtual void additionalAtEndTimestep();
  virtual void compute_cell_areas(); // is an initialization step; should go there

  // see iMIO.cc
  virtual void initFromFile(const std::string &name);
  virtual void writeFiles(const std::string &default_filename);
  virtual void write_model_state(const PIO &nc);
  virtual void write_metadata(const PIO &nc,
                              bool write_mapping,
                              bool write_run_stats);
  virtual void write_variables(const PIO &nc, const std::set<std::string> &vars,
                               IO_Type nctype);
protected:

  //! Computational grid
  const IceGrid::Ptr m_grid;
  //! Configuration flags and parameters
  const Config::Ptr m_config;
  //! Execution context
  const Context::Ptr m_ctx;
  //! Unit system
  const units::System::Ptr m_sys;
  //! Logger
  const Logger::Ptr m_log;
  //! Time manager
  const Time::Ptr m_time;

  VariableMetadata global_attributes, //!< stores global attributes saved in a PISM output file
    mapping,                    //!< grid projection (mapping) parameters
    run_stats;                  //!< run statistics

  hydrology::Hydrology   *subglacial_hydrology;
  YieldStress *basal_yield_stress_model;

  energy::BedThermalUnit *btu;

  calving::IcebergRemover     *iceberg_remover;
  calving::OceanKill          *ocean_kill_calving;
  calving::FloatKill          *float_kill_calving;
  calving::CalvingAtThickness *thickness_threshold_calving;
  calving::EigenCalving       *eigen_calving;
  calving::CrevassesCalving   *crevasses_calving;

  surface::SurfaceModel *surface;
  ocean::OceanModel   *ocean;
  bed::BedDef       *beddef;
  bool external_surface_model, external_ocean_model;

  // state variables and some diagnostics/internals
  IceModelVec2S ice_surface_elevation,          //!< ice surface elevation; ghosted
    ice_thickness,              //!< ghosted
    basal_yield_stress,         //!< ghosted
    basal_melt_rate,           //!< rate of production of basal meltwater (ice-equivalent); no ghosts
    vLongitude, //!< Longitude; ghosted to compute cell areas
    vLatitude,  //!< Latitude; ghosted to compute cell areas
    geothermal_flux,   //!< geothermal flux; no ghosts
    vFD,    //!< fracture density
    vFG,    //!< fracture growth rate
    vFH,    //!< fracture healing rate
    vFE,    //!< fracture flow enhancement
    vFA,    //!< fracture age
    vFT,    //!< fracture toughness
    bedtoptemp,     //!< temperature seen by bedrock thermal layer, if present; no ghosts
    vHref,          //!< accumulated mass advected to a partially filled grid cell
    climatic_mass_balance,              //!< accumulation/ablation rate; no ghosts
    climatic_mass_balance_cumulative,    //!< cumulative climatic_mass_balance
    grounded_basal_flux_2D_cumulative, //!< grounded basal (melt/freeze-on) cumulative flux
    floating_basal_flux_2D_cumulative, //!< floating (sub-shelf) basal (melt/freeze-on) cumulative flux
    nonneg_flux_2D_cumulative,         //!< cumulative nonnegative-rule flux
    discharge_flux_2D_cumulative,      //!< cumulative discharge (calving) flux (2D field)
    ice_surface_temp,           //!< ice temperature at the ice surface but below firn; no ghosts
    liqfrac_surface,    //!< ice liquid water fraction at the top surface of the ice
    shelfbtemp,         //!< ice temperature at the shelf base; no ghosts
    shelfbmassflux,     //!< ice mass flux into the ocean at the shelf base; no ghosts
    cell_area,          //!< cell areas (computed using the WGS84 datum)
    flux_divergence;    //!< flux divergence

  IceModelVec2S land_flux_2D, //!< flux on land (incl. lakes)
    ocean_flux_2D,            //!< flux on ocean (excl. lakes)
    land_flux_2D_cumulative,  //!< cumulative flux on land (incl. lakes)
    ocean_flux_2D_cumulative; //!< cumulative flux on the ocean (excl. lakes)

  // some more diagnostics/internals
  //public: //??? fixme
  IceModelVec2S crevasses_calv_flux_2D,//!< crevasses calving
    crevasses_calv_flux_2D_cumulative, //!< cumulative crevasses calving
    surface_crevasses_h,  //!< surface crevasses depth
    bottom_crevasses_h,   //!< bottom crevasses depth
    crevasses_dw;         //!< water table in crevasses



public:  IceModelVec2S* get_geothermal_flux();
  void setCTSFromEnthalpy(IceModelVec3 &result);
protected:

  IceModelVec2 strain_rates; //!< major and minor principal components of horizontal strain-rate tensor
  
  IceModelVec2 deviatoric_stresses; //!< components of horizontal stress tensor along axes and shear stress

  IceModelVec2Int vMask, //!< \brief mask for flow type with values ice_free_bedrock,
  //!< grounded_ice, floating_ice, ice_free_ocean
    vBCMask; //!< mask to determine Dirichlet boundary locations
 
  IceModelVec2V vBCvel; //!< Dirichlet boundary velocities
  
  IceModelVec2S gl_mask, //!< mask to determine grounding line position
    gl_mask_x, //!< mask to determine grounding line position in x-direction
    gl_mask_y; //!< mask to determine grounding line position in y-direction

  IceModelVec3
  T3,             //!< absolute temperature of ice; K (ghosted)
    Enth3,          //!< enthalpy; J / kg (ghosted)
    age3;           //!< age of ice; s (ghosted because it is averaged onto the staggered-grid)

  // parameters
  double   dt,     //!< mass continuity time step, s
    t_TempAge,  //!< time of last update for enthalpy/temperature
    dt_TempAge,  //!< enthalpy/temperature and age time-steps
    maxdt_temporary, dt_force,
    dt_from_cfl, CFLmaxdt, CFLmaxdt2D,
    gmaxu, gmaxv, gmaxw,  // global maximums on 3D grid of abs value of vel components
    grounded_basal_ice_flux_cumulative,
    nonneg_rule_flux_cumulative,
    sub_shelf_ice_flux_cumulative,
    surface_ice_flux_cumulative,
    sum_divQ_SIA_cumulative,
    sum_divQ_SSA_cumulative,
    Href_to_H_flux_cumulative,
    H_to_Href_flux_cumulative,
    discharge_flux_cumulative;      //!< cumulative discharge (calving) flux

  double land_flux,            //!< land flux
    land_flux_cumulative,      //!< cumulative land flux
    ocean_flux,                //!< ocean flux
    ocean_flux_cumulative,     //!< cumulative ocean flux
    crevasses_calv_flux,       //!< crevasses calving
    crevasses_calv_flux_cumulative; //!< cumulative crevasses calving

  unsigned int skipCountDown,
    CFLviolcount;

  // flags
  std::string m_adaptive_timestep_reason;

  std::string stdout_flags;

protected:
  // see iceModel.cc
  virtual void createVecs();

  // see iMadaptive.cc
  virtual double max_timestep_cfl_3d();
  virtual double max_timestep_cfl_2d();
  virtual double max_timestep_diffusivity();
  virtual void max_timestep(double &dt_result, unsigned int &skip_counter);
  virtual unsigned int countCFLViolations();
  virtual unsigned int skip_counter(double input_dt, double input_dt_diffusivity);

  // see iMage.cc
  virtual void ageStep();

  // see iMenergy.cc
  virtual void energyStep();
  virtual void get_bed_top_temp(IceModelVec2S &result);
  virtual bool checkThinNeigh(const IceModelVec2S &thickness,
                              int i, int j, double threshold);

  virtual void combine_basal_melt_rate();

  // see iMenthalpy.cc
  virtual void compute_enthalpy_cold(const IceModelVec3 &temperature, IceModelVec3 &result);
  virtual void compute_enthalpy(const IceModelVec3 &temperature,
                                const IceModelVec3 &liquid_water_fraction,
                                IceModelVec3 &result);
  virtual void compute_liquid_water_fraction(const IceModelVec3 &enthalpy,
                                             IceModelVec3 &result);

  virtual void enthalpyAndDrainageStep(unsigned int *vertSacrCount,
                                       double* liquifiedVol,
                                       unsigned int *bulgeCount);

  // see iMgeometry.cc
  virtual void updateSurfaceElevationAndMask();
  virtual void update_mask(const IceModelVec2S &bed,
                           const IceModelVec2S &ice_thickness,
                           IceModelVec2Int &mask);
  virtual void update_surface_elevation(const IceModelVec2S &bed,
                                        const IceModelVec2S &ice_thickness,
                                        IceModelVec2S &result);
  virtual void cell_interface_fluxes(bool dirichlet_bc,
                                     int i, int j,
                                     StarStencil<Vector2> input_velocity,
                                     StarStencil<double> input_flux,
                                     StarStencil<double> &output_velocity,
                                     StarStencil<double> &output_flux);
  virtual void adjust_flow(StarStencil<int> mask,
                           StarStencil<double> &SSA_velocity,
                           StarStencil<double> &SIA_flux);
  virtual void massContExplicitStep();
  virtual void update_floatation_mask();
  virtual void do_calving();
  virtual void Href_cleanup();
  virtual void update_cumulative_discharge(const IceModelVec2S &thickness,
                                           const IceModelVec2S &thickness_old,
                                           const IceModelVec2S &Href,
                                           const IceModelVec2S &Href_old);


  // see iMIO.cc
  virtual void dumpToFile(const std::string &filename);
  virtual void regrid(int dimensions);
  virtual void regrid_variables(const std::string &filename,
                                const std::set<std::string> &regrid_vars,
                                unsigned int ndims);
  virtual void init_enthalpy(const std::string &filename, bool regrid, int last_record);

  // see iMfractures.cc
  virtual void calculateFractureDensity();

  // see iMpartgrid.cc
  double get_threshold_thickness(StarStencil<int> Mask,
                                 StarStencil<double> thickness,
                                 StarStencil<double> surface_elevation,
                                 double bed_elevation,
                                 bool reduce_frontal_thickness);
  virtual void residual_redistribution(IceModelVec2S &residual);
  virtual void residual_redistribution_iteration(IceModelVec2S &residual, bool &done);

  // see iMreport.cc
  virtual double compute_temperate_base_fraction(double ice_area);
  virtual double compute_original_ice_fraction(double ice_volume);
  virtual void summary(bool tempAndAge);
  virtual void summaryPrintLine(bool printPrototype, bool tempAndAge,
                                double delta_t,
                                double volume, double area,
                                double meltfrac, double max_diffusivity);

  // see iMreport.cc;  methods for computing diagnostic quantities:
  // scalar:
  virtual double compute_ice_volume();
  virtual double compute_sealevel_volume();
  virtual double compute_ice_volume_temperate();
  virtual double compute_ice_volume_cold();
  virtual double compute_ice_area();
  virtual double compute_ice_area_temperate();
  virtual double compute_ice_area_cold();
  virtual double compute_ice_area_grounded();
  virtual double compute_ice_area_floating();
  virtual double compute_ice_enthalpy();

  // see iMtemp.cc
  virtual void excessToFromBasalMeltLayer(double rho, double c, double L,
                                          double z, double dz,
                                          double *Texcess, double *bwat);
  virtual void temperatureStep(unsigned int *vertSacrCount, unsigned int *bulgeCount);

  // see iMutil.cc
  virtual int endOfTimeStepHook();
  virtual void stampHistoryCommand();
  virtual void stampHistoryEnd();
  virtual void stampHistory(const std::string &);
  virtual void update_run_stats();
  virtual void check_maximum_thickness();
  virtual void check_maximum_thickness_hook(const int old_Mz);

protected:
  // working space (a convenience)
  static const int nWork2d=2;
  IceModelVec2S vWork2d[nWork2d];

  // 3D working space
  IceModelVec3 vWork3d;

  stressbalance::StressBalance *stress_balance;

public:
  stressbalance::StressBalance* get_stress_balance();
protected:

  std::map<std::string,Diagnostic*> diagnostics;
  std::map<std::string,TSDiagnostic*> ts_diagnostics;

  // Set of variables to put in the output file:
  std::set<std::string> output_vars;

  // This is related to the snapshot saving feature
  std::string snapshots_filename;
  bool save_snapshots, snapshots_file_is_ready, split_snapshots;
  std::vector<double> snapshot_times;
  std::set<std::string> snapshot_vars;
  unsigned int current_snapshot;
  void init_snapshots();
  void write_snapshot();

  // scalar time-series
  bool save_ts;                 //! true if the user requested time-series output
  std::string ts_filename;              //! file to write time-series to
  std::vector<double> ts_times; //! times requested
  unsigned int current_ts;      //! index of the current time
  std::set<std::string> ts_vars;                //! variables requested
  void init_timeseries();
  void flush_timeseries();
  void write_timeseries();
  MaxTimestep ts_max_timestep(double my_t);

  // spatially-varying time-series
  bool save_extra, extra_file_is_ready, split_extra;
  std::string extra_filename;
  std::vector<double> extra_times;
  unsigned int next_extra;
  double last_extra;
  std::set<std::string> extra_vars;
  TimeBoundsMetadata extra_bounds;
  TimeseriesMetadata timestamp;
  void init_extras();
  void write_extras();
  MaxTimestep extras_max_timestep(double my_t);

  // automatic backups
  double backup_interval;
  std::string backup_filename;
  double last_backup_time;
  std::set<std::string> backup_vars;
  void init_backups();
  void write_backup();

  // last time at which PISM hit a multiple of X years, see the
  // timestep_hit_multiples configuration parameter
  double timestep_hit_multiples_last_time;

  // diagnostic viewers; see iMviewers.cc
  virtual void init_viewers();
  virtual void update_viewers();
  virtual void view_field(const IceModelVec *field);
  std::set<std::string> map_viewers, slice_viewers;
  int     id, jd;            // sounding indexes
  std::map<std::string,petsc::Viewer::Ptr> viewers;

private:
  double start_time;    // this is used in the wall-clock-time backup code
};

} // end of namespace pism

#endif /* __iceModel_hh */

