// Copyright (C) 2009-2015 Andreas Aschwanden and Ed Bueler and Constantine Khroulev
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

#include "enthSystem.hh"
#include <gsl/gsl_math.h>
#include "base/util/PISMConfigInterface.hh"
#include "base/util/iceModelVec.hh"
#include "base/enthalpyConverter.hh"

#include "base/util/error_handling.hh"
#include "base/util/ColumnInterpolation.hh"

namespace pism {
namespace energy {

enthSystemCtx::enthSystemCtx(const std::vector<double>& storage_grid,
                             const std::string &prefix,
                             double dx,  double dy, double dt,
                             const Config &config,
                             const IceModelVec3 &Enth3,
                             const IceModelVec3 &u3,
                             const IceModelVec3 &v3,
                             const IceModelVec3 &w3,
                             const IceModelVec3 &strain_heating3,
                             EnthalpyConverter::Ptr EC)
: columnSystemCtx(storage_grid, prefix, dx, dy, dt, u3, v3, w3),
  m_Enth3(Enth3),
  m_strain_heating3(strain_heating3),
  m_EC(EC) {

  // set some values so we can check if init was called
  m_R_cold   = -1.0;
  m_R_temp   = -1.0;
  m_lambda = -1.0;
  m_D0 = GSL_NAN;
  m_U0 = GSL_NAN;
  m_B0 = GSL_NAN;

  m_ice_density = config.get_double("ice_density");
  m_ice_c   = config.get_double("ice_specific_heat_capacity");
  m_ice_k   = config.get_double("ice_thermal_conductivity");
  m_p_air   = config.get_double("surface_pressure");

  m_ice_K  = m_ice_k / m_ice_c;
  m_ice_K0 = m_ice_K * config.get_double("enthalpy_temperate_conductivity_ratio");

  size_t Mz = m_z.size();
  m_Enth.resize(Mz);
  m_Enth_s.resize(Mz);
  m_strain_heating.resize(Mz);
  m_R.resize(Mz);

  m_E_n.resize(Mz);
  m_E_e.resize(Mz);
  m_E_s.resize(Mz);
  m_E_w.resize(Mz);

  m_nu = m_dt / m_dz;

  m_R_factor = m_dt / (PetscSqr(m_dz) * m_ice_density);
  m_R_cold = m_ice_K * m_R_factor;
  m_R_temp = m_ice_K0 * m_R_factor;

  if (config.get_boolean("use_temperature_dependent_thermal_conductivity")) {
    m_k_depends_on_T = true;
  } else {
    m_k_depends_on_T = false;
  }

  // check if c(T) is a constant function:
  if (EC->c(260) != EC->c(270)) {
    m_c_depends_on_T = true;
  } else {
    m_c_depends_on_T = false;
  }
}


enthSystemCtx::~enthSystemCtx() {
}

/*!
  In this implementation \f$k\f$ does not depend on temperature.
 */
double enthSystemCtx::k_from_T(double T) {

  if (m_k_depends_on_T) {
    return 9.828 * exp(-0.0057 * T);
  }

  return m_ice_k;
}

void enthSystemCtx::initThisColumn(int my_i, int my_j, bool my_ismarginal,
                                   double my_ice_thickness) {
  m_ice_thickness = my_ice_thickness;
  m_ismarginal    = my_ismarginal;

  init_column(my_i, my_j, m_ice_thickness);

  if (m_ks == 0) {
    return;
  }

  coarse_to_fine(m_u3, m_i, m_j, &m_u[0]);
  coarse_to_fine(m_v3, m_i, m_j, &m_v[0]);
  coarse_to_fine(m_w3, m_i, m_j, &m_w[0]);
  coarse_to_fine(m_strain_heating3, m_i, m_j, &m_strain_heating[0]);
  coarse_to_fine(m_Enth3, m_i, m_j, &m_Enth[0]);

  coarse_to_fine(m_Enth3, m_i, m_j+1, &m_E_n[0]);
  coarse_to_fine(m_Enth3, m_i+1, m_j, &m_E_e[0]);
  coarse_to_fine(m_Enth3, m_i, m_j-1, &m_E_s[0]);
  coarse_to_fine(m_Enth3, m_i-1, m_j, &m_E_w[0]);

  compute_enthalpy_CTS();

  m_lambda = compute_lambda();

  assemble_R();
}

//! Compute the CTS value of enthalpy in an ice column.
/*! m_Enth_s is set to the enthalpy value for the pressure-melting
  temperature with zero water fraction at the corresponding z level.
 */
void enthSystemCtx::compute_enthalpy_CTS() {

  for (unsigned int k = 0; k <= m_ks; k++) {
    const double
      depth = m_ice_thickness - k * m_dz,
      p = m_EC->pressure(depth); // FIXME issue #15
    m_Enth_s[k] = m_EC->enthalpy_cts(p);
  }

  const double Es_air = m_EC->enthalpy_cts(m_p_air);
  for (unsigned int k = m_ks+1; k < m_Enth_s.size(); k++) {
    m_Enth_s[k] = Es_air;
  }
}

//! Compute the lambda for BOMBPROOF.
/*!
See page \ref bombproofenth.
 */
double enthSystemCtx::compute_lambda() {
  double result = 1.0; // start with centered implicit for more accuracy
  const double epsilon = 1e-6 / 3.15569259747e7;

  for (unsigned int k = 0; k <= m_ks; k++) {
    if (m_Enth[k] > m_Enth_s[k]) { // lambda = 0 if temperate ice present in column
      result = 0.0;
    } else {
      const double denom = (fabs(m_w[k]) + epsilon) * m_ice_density * m_ice_c * m_dz;
      result = std::min(result, 2.0 * m_ice_k / denom);
    }
  }
  return result;
}


void enthSystemCtx::setDirichletSurface(double my_Enth_surface) {
#if (PISM_DEBUG==1)
  if ((m_nu < 0.0) || (m_R_cold < 0.0) || (m_R_temp < 0.0)) {
    throw RuntimeError("setDirichletSurface() should only be called after\n"
                       "initAllColumns() in enthSystemCtx");
  }
#endif
  m_Enth_ks = my_Enth_surface;
}

void enthSystemCtx::checkReadyToSolve() {
  if (m_nu < 0.0 || m_R_cold < 0.0 || m_R_temp < 0.0) {
    throw RuntimeError("not ready to solve: need initAllColumns() in enthSystemCtx");
  }
  if (m_lambda < 0.0) {
    throw RuntimeError("not ready to solve: need setSchemeParamsThisColumn() in enthSystemCtx");
  }
}


//! Set coefficients in discrete equation for \f$E = Y\f$ at base of ice.
/*!
This method should only be called if everything but the basal boundary condition
is already set.
 */
void enthSystemCtx::setDirichletBasal(double Y) {
#if (PISM_DEBUG==1)
  checkReadyToSolve();
  if (gsl_isnan(m_D0) == 0 || gsl_isnan(m_U0) == 0 || gsl_isnan(m_B0) == 0) {
    throw RuntimeError("setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  m_D0 = 1.0;
  m_U0 = 0.0;
  m_B0  = Y;
}


//! Set coefficients in discrete equation for Neumann condition at base of ice.
/*!
This method generates the Neumann boundary condition for the linear system.

The Neumann boundary condition is
   @f[ \frac{\partial E}{\partial z} = - \frac{\phi}{K} @f]
where \f$\phi\f$ is the heat flux.  Here \f$K\f$ is allowed to vary, and takes
its value from the value computed in assemble_R().

The boundary condition is combined with the partial differential equation by the
technique of introducing an imaginary point at \f$z=-\Delta z\f$ and then
eliminating it.

In other words, we combine the centered finite difference approximation
@f[ \frac{ E_{1} - E_{-1} }{2\dz}  = -\frac{\phi}{K} @f]
with

@f[ -R_{k-\frac12} E_{k-1} + \left( 1 + R_{k-\frac12} + R_{k+\frac12} \right) E_{k} - R_{k+\frac12} E_{k+1} + \text{advective terms} = \dots @f]

to get

@f{align*}{
   \frac{E_{1}-E_{-1}}{2\,\Delta z} & = -\frac{\phi}{K_{0}}, \\
   E_{1}-E_{-1} & = -\frac{2\,\Delta z\,\phi}{K_{0}}, \\
    E_{-1}\,R_{-\frac12}-R_{-\frac12}\,E_{1} & = \frac{2\,R_{-\frac12}\,\Delta z\,\phi}{K_{0}}, \\
    -R_{\frac12}\,E_{1}+E_{0}\,\left(R_{\frac12}+R_{-\frac12}+1\right)-E_{-1}\,R_{-\frac12} + \text{advective terms} & = \dots, \\
    \left(-R_{\frac12}-R_{-\frac12}\right)\,E_{1}+E_{0}\,\left(R_{\frac12}+R_{-\frac12}+1\right) + \text{advective terms} & = \frac{2\,R_{-\frac12}\,\Delta z\,\phi}{K_{0}}+\dots.
@f}

The error in the pure conductive and smooth conductivity case is @f$ O(\dz^2) @f$.

This method should only be called if everything but the basal boundary condition
is already set.

 */
void enthSystemCtx::setBasalHeatFlux(double heat_flux) {
#if (PISM_DEBUG==1)
 checkReadyToSolve();
  if (gsl_isnan(m_D0) == 0 || gsl_isnan(m_U0) == 0 || gsl_isnan(m_B0) == 0) {
    throw RuntimeError("setting basal boundary conditions twice in enthSystemCtx");
  }
#endif
  // extract K from R[0], so this code works even if K=K(T)
  // recall:   R = (ice_K / ice_density) * dt / PetscSqr(dz)
  const double
    K      = (m_ice_density * PetscSqr(m_dz) * m_R[0]) / m_dt,
    Rc     = m_R[0],
    Rr     = m_R[1],
    Rminus = Rc,
    Rplus  = 0.5 * (Rc + Rr);
  m_D0 = 1.0 + Rminus + Rplus;
  // modified upper-diagonal term:
  m_U0 = - Rminus - Rplus;
  // m_Enth[0] (below) is there due to the fully-implicit discretization
  // in time, the second term is the modification of the right-hand
  // side implementing the Neumann B.C. (see the doxygen comment)
  m_B0 = m_Enth[0] + 2.0 * Rminus * m_dz * heat_flux / K;
  // treat vertical velocity using first-order upwinding:
  if (not m_ismarginal) {
    const double UpEnthu = (m_u[0] < 0 ?
                            m_u[0] * (m_E_e[0] -  m_Enth[0]) / m_dx :
                            m_u[0] * (m_Enth[0]  - m_E_w[0]) / m_dx);
    const double UpEnthv = (m_v[0] < 0 ?
                            m_v[0] * (m_E_n[0] -  m_Enth[0]) / m_dy :
                            m_v[0] * (m_Enth[0]  - m_E_s[0]) / m_dy);

    m_B0 += m_dt * ((m_strain_heating[0] / m_ice_density) - UpEnthu - UpEnthv);  // = rhs[0]
  }
}


//! \brief Assemble the R array.  The R value switches at the CTS.
/*!  In a simple abstract diffusion
  \f[ \frac{\partial u}{\partial t} = D \frac{\partial^2 u}{\partial z^2}, \f]
with time steps \f$\Delta t\f$ and spatial steps \f$\Delta z\f$ we define
  \f[ R = \frac{D \Delta t}{\Delta z^2}. \f]
This is used in an implicit method to write each line in the linear system, for
example [\ref MortonMayers]:
  \f[ -R U_{j-1}^{n+1} + (1+2R) U_j^{n+1} - R U_{j+1}^{n+1} = U_j^n. \f]
  
In the case of conservation of energy [\ref AschwandenBuelerKhroulevBlatter],
  \f[ u=E \qquad \text{ and } \qquad D = \frac{K}{\rho} \qquad \text{ and } \qquad K = \frac{k}{c}. \f]
Thus
  \f[ R = \frac{k \Delta t}{\rho c \Delta z^2}. \f]
 */
void enthSystemCtx::assemble_R() {
  if (m_k_depends_on_T == false && m_c_depends_on_T == false) {

    for (unsigned int k = 1; k <= m_ks; k++) {
      m_R[k] = (m_Enth[k] < m_Enth_s[k]) ? m_R_cold : m_R_temp;
    }
    //still the cold ice value, if no temperate layer above
    m_R[0] = (m_Enth[1] < m_Enth_s[1]) ? m_R_cold : m_R_temp;

  } else {

    for (unsigned int k = 1; k <= m_ks; k++) {
      if (m_Enth[k] < m_Enth_s[k]) {
        // cold case
        const double depth = m_ice_thickness - k * m_dz;
        double T = m_EC->temperature(m_Enth[k],
                                   m_EC->pressure(depth)); // FIXME: issue #15

        m_R[k] = ((m_k_depends_on_T ? k_from_T(T) : m_ice_k) / m_EC->c(T)) * m_R_factor;
      } else {
        // temperate case
        m_R[k] = m_R_temp;
      }
    }
    // still the cold ice value, if no temperate layer above
    if (m_Enth[1] < m_Enth_s[1]) {
      double T = m_EC->temperature(m_Enth[0],
                                   m_EC->pressure(m_ice_thickness)); // FIXME: issue #15
      m_R[0] = ((m_k_depends_on_T ? k_from_T(T) : m_ice_k) / m_EC->c(T)) * m_R_factor;
    } else {
      // temperate layer case
      m_R[0] = m_R_temp;
    }

  }

  // R[k] for k > m_ks are never used
#if (PISM_DEBUG==1)
  for (unsigned int k = m_ks + 1; k < m_R.size(); ++k) {
    m_R[k] = GSL_NAN;
  }
#endif
}


/*! \brief Solve the tridiagonal system, in a single column, which
 *  determines the new values of the ice enthalpy.
 *
 * We are solving a convection-diffusion equation, treating the @f$ z @f$ direction implicitly and
 * @f$ x, y @f$ directions explicitly. See @ref bombproofenth for the documentation of the treatment
 * of convection terms. The notes below document the way we treat diffusion using a simplified PDE.
 *
 * To discretize
 * @f[ \diff{}{z} \left( K(E) \diff{E}{z}\right) = \diff{E}{t} @f]
 *
 * at a location @f$ k @f$ of the vertical grid, we use centered finite differences in space,
 * backward differences in time, and evaluate @f$ K(E) @f$ at staggered-grid locations:
 *
 * @f[ \frac{K_{k+\frac12}\frac{E_{k+1} - E_{k}}{\dz} - K_{k-\frac12}\frac{E_{k} - E_{k-1}}{\dz}}{\dz} = \frac{E_{k} - E^{n-1}_{k}}{\dt}, @f]
 *
 * where @f$ E = E^{n} @f$ for clarity and @f$ K_{k\pm \frac12} = K(E^{n-1}_{k\pm \frac12}) @f$,
 * %i.e. we linearize this PDE by evaluating @f$ K(E) @f$ at the _previous_ time-step.
 *
 * We define
 *
 * @f[ R_i = \frac{\dt\, K_i}{\dz^2}, @f]
 *
 * and the discretization takes form
 *
 * @f[ -R_{k-\frac12} E_{k-1} + \left( 1 + R_{k-\frac12} + R_{k+\frac12} \right) E_{k} - R_{k+\frac12} E_{k+1} = E^{n-1}_{k}. @f]
 *
 * In the assembly of the tridiagonal system this corresponds to
 *
 * @f{align*}{
 * L_i &= - \frac12 (R_{i} + R_{i-1}),\\
 * D_i &= 1 + \frac12 (R_{i} + R_{i-1}) + \frac12 (R_{i} + R_{i+1}),\\
 * U_i &= - \frac12 (R_{i} + R_{i+1}),\\
 * b_i &= E^{n-1}_{i},
 * @f}
 *
 * where @f$ L_i, D_i, U_i @f$ are lower-diagonal, diagonal, and upper-diagonal entries
 * corresponding to an equation @f$ i @f$ and @f$ b_i @f$ is the corresponding right-hand side.
 * (Staggered-grid values are approximated by interpolating from the regular grid).
 *
 * This method is _unconditionally stable_ and has a maximum principle (see [@ref MortonMayers,
 * section 2.11]).
 */
void enthSystemCtx::solveThisColumn(std::vector<double> &x) {

  TridiagonalSystem &S = *m_solver;

#if (PISM_DEBUG==1)
  checkReadyToSolve();
  if (gsl_isnan(m_D0) || gsl_isnan(m_U0) || gsl_isnan(m_B0)) {
    throw RuntimeError("solveThisColumn() should only be called after\n"
                       "  setting basal boundary condition in enthSystemCtx");
  }
#endif

  // k=0 equation is already established
  // L[0] = 0.0;  // not used
  S.D(0)   = m_D0;
  S.U(0)   = m_U0;
  S.RHS(0) = m_B0;

  // generic ice segment in k location (if any; only runs if m_ks >= 2)
  for (unsigned int k = 1; k < m_ks; k++) {
    const double
        Rminus = 0.5 * (m_R[k-1] + m_R[k]),
        Rplus  = 0.5 * (m_R[k]   + m_R[k+1]);
    S.L(k) = - Rminus;
    S.D(k) = 1.0 + Rminus + Rplus;
    S.U(k) = - Rplus;
    const double AA = m_nu * m_w[k];
    if (m_w[k] >= 0.0) {  // velocity upward
      S.L(k) -= AA * (1.0 - m_lambda/2.0);
      S.D(k) += AA * (1.0 - m_lambda);
      S.U(k) += AA * (m_lambda/2.0);
    } else {            // velocity downward
      S.L(k) -= AA * (m_lambda/2.0);
      S.D(k) -= AA * (1.0 - m_lambda);
      S.U(k) += AA * (1.0 - m_lambda/2.0);
    }
    S.RHS(k) = m_Enth[k];
    if (not m_ismarginal) {
      const double UpEnthu = (m_u[k] < 0 ?
                              m_u[k] * (m_E_e[k] -  m_Enth[k]) / m_dx :
                              m_u[k] * (m_Enth[k]  - m_E_w[k]) / m_dx);
      const double UpEnthv = (m_v[k] < 0 ?
                              m_v[k] * (m_E_n[k] -  m_Enth[k]) / m_dy :
                              m_v[k] * (m_Enth[k]  - m_E_s[k]) / m_dy);

      S.RHS(k) += m_dt * ((m_strain_heating[k] / m_ice_density) - UpEnthu - UpEnthv);
    }
  }

  // set Dirichlet boundary condition at top
  if (m_ks > 0) {
    S.L(m_ks) = 0.0;
  }
  S.D(m_ks) = 1.0;
  if (m_ks < m_z.size() - 1) {
    S.U(m_ks) = 0.0;
  }
  S.RHS(m_ks) = m_Enth_ks;

  // Solve it; note drainage is not addressed yet and post-processing may occur
  try {
    S.solve(m_ks + 1, x);
  }
  catch (RuntimeError &e) {
    e.add_context("solving the tri-diagonal system (enthSystemCtx) at (%d,%d)\n"
                  "saving system to m-file... ", m_i, m_j);
    reportColumnZeroPivotErrorMFile(m_ks + 1);
    throw;
  }

  // air above
  for (unsigned int k = m_ks+1; k < x.size(); k++) {
    x[k] = m_Enth_ks;
  }

#if (PISM_DEBUG==1)
  // if success, mark column as done by making scheme params and b.c. coeffs invalid
  m_lambda = -1.0;
  m_D0     = GSL_NAN;
  m_U0     = GSL_NAN;
  m_B0     = GSL_NAN;
#endif
}

void enthSystemCtx::save_system(std::ostream &output, unsigned int system_size) const {
  m_solver->save_system(output, system_size);
  m_solver->save_vector(output, m_R, system_size, m_solver->prefix() + "_R");
}

} // end of namespace energy
} // end of namespace pism
