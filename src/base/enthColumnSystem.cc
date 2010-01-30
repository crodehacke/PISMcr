// Copyright (C) 2009-2010 Andreas Aschwanden and Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "pism_const.hh"        // e.g. for MASK_FLOATING
#include "enthalpyConverter.hh"
#include "enthColumnSystem.hh"


enthSystemCtx::enthSystemCtx(int my_Mz, int my_Mbz)
      : columnSystemCtx(my_Mz + my_Mbz - 1) {
  Mz = my_Mz;
  Mbz = my_Mbz;
  k0 = Mbz - 1; // max size nmax of system is Mz + k0 = Mz + Mbz - 1
  // set flags to indicate nothing yet set
  initAllDone = false;
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;
  // set values so we can check if init was called on all
  dx = -1;
  dy = -1;
  dtTemp = -1;
  dzEQ = -1;
  dzbEQ = -1;
  ice_rho = -1;
  ice_c   = -1;
  ice_k   = -1;
  ice_nu  = -1;
  bed_rho = -1;
  bed_c   = -1;
  bed_k   = -1;
  Enth = NULL;
  Enth_s = NULL;
  Tb = NULL;
  u = NULL;
  v = NULL;
  w = NULL;
  Sigma = NULL;
  Enth3 = NULL;
}


PetscErrorCode enthSystemCtx::initAllColumns() {
  // check whether each parameter & pointer got set
  if (dx <= 0.0) { SETERRQ(2,"un-initialized dx in enthSystemCtx"); }
  if (dy <= 0.0) { SETERRQ(3,"un-initialized dy in enthSystemCtx"); }
  if (dtTemp <= 0.0) { SETERRQ(4,"un-initialized dtTemp in enthSystemCtx"); }
  if (dzEQ <= 0.0) { SETERRQ(5,"un-initialized dzEQ in enthSystemCtx"); }
  if (dzbEQ <= 0.0) { SETERRQ(6,"un-initialized dzbEQ in enthSystemCtx"); }
  if (ice_rho <= 0.0) { SETERRQ(7,"un-initialized ice_rho in enthSystemCtx"); }
  if (ice_c <= 0.0) { SETERRQ(8,"un-initialized ice_c in enthSystemCtx"); }
  if (ice_k <= 0.0) { SETERRQ(9,"un-initialized ice_k in enthSystemCtx"); }
  if (ice_nu < 0.0) { SETERRQ(21,"un-initialized ice_nu in enthSystemCtx"); } // 0.0 is allowed
  if (bed_rho <= 0.0) { SETERRQ(10,"un-initialized bed_rho in enthSystemCtx"); }
  if (bed_c <= 0.0) { SETERRQ(11,"un-initialized bed_c in enthSystemCtx"); }
  if (bed_k <= 0.0) { SETERRQ(12,"un-initialized bed_k in enthSystemCtx"); }
  if (Enth == NULL) { SETERRQ(13,"un-initialized pointer Enth in enthSystemCtx"); }
  if (Enth_s == NULL) { SETERRQ(20,"un-initialized pointer Enth_s in enthSystemCtx"); }
  if (Tb == NULL) { SETERRQ(14,"un-initialized pointer Tb in enthSystemCtx"); }
  if (u == NULL) { SETERRQ(15,"un-initialized pointer u in enthSystemCtx"); }
  if (v == NULL) { SETERRQ(16,"un-initialized pointer v in enthSystemCtx"); }
  if (w == NULL) { SETERRQ(17,"un-initialized pointer w in enthSystemCtx"); }
  if (Sigma == NULL) { SETERRQ(18,"un-initialized pointer Sigma in enthSystemCtx"); }
  if (Enth3 == NULL) { SETERRQ(19,"un-initialized pointer Enth3 in enthSystemCtx"); }
  // set derived constants
  nuEQ = dtTemp / dzEQ;
  iceK = ice_k / (ice_rho * ice_c);
  iceRcold = iceK * dtTemp / PetscSqr(dzEQ);
  iceRtemp = ice_nu * dtTemp / PetscSqr(dzEQ);
  bedK = bed_k / (bed_rho * bed_c);
  bedR = bedK * dtTemp / PetscSqr(dzbEQ);
  dzav = 0.5 * (dzEQ + dzbEQ);
  // done
  initAllDone = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setSchemeParamsThisColumn(
                     PismMask my_mask, bool my_isMarginal, PetscScalar my_lambda) {
  if (!initAllDone) {  SETERRQ(2,
     "setSchemeParamsThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (schemeParamsValid) {  SETERRQ(3,
     "setSchemeParamsThisColumn() called twice (?) in enthSystemCtx"); }
  mask = my_mask;
  isMarginal = my_isMarginal;
  lambda = my_lambda;
  schemeParamsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setSurfaceBoundaryValuesThisColumn(PetscScalar my_Enth_surface) {
  if (!initAllDone) {  SETERRQ(2,
     "setSurfaceBoundaryValuesThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (surfBCsValid) {  SETERRQ(3,
     "setSurfaceBoundaryValuesThisColumn() called twice (?) in enthSystemCtx"); }
  Enth_ks = my_Enth_surface;
  surfBCsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Enth_shelfbase, PetscScalar my_Rb) {
  if (!initAllDone) {  SETERRQ(2,
     "setBasalBoundaryValuesThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (basalBCsValid) {  SETERRQ(3,
     "setBasalBoundaryValuesThisColumn() called twice (?) in enthSystemCtx"); }
  Ghf = my_Ghf;
  Enth_shelfbase = my_Enth_shelfbase;
  Rb = my_Rb;
  basalBCsValid = true;
  return 0;
}


PetscErrorCode enthSystemCtx::viewConstants(PetscViewer viewer) {
  PetscErrorCode ierr;

  if (!viewer) {
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_SELF,&viewer); CHKERRQ(ierr);
  }

  PetscTruth iascii;
  ierr = PetscTypeCompare((PetscObject)viewer,PETSC_VIEWER_ASCII,&iascii); CHKERRQ(ierr);
  if (!iascii) { SETERRQ(1,"Only ASCII viewer for enthSystemCtx::viewConstants()\n"); }
  
  ierr = PetscViewerASCIIPrintf(viewer,"\n<<VIEWING enthSystemCtx:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "initAllDone = %d\n"
                     "for ALL columns:\n",
                     initAllDone); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  dx,dy,dtTemp,dzEQ,dzbEQ = %8.2f,%8.2f,%10.3e,%8.2f,%8.2f\n",
                     dx,dy,dtTemp,dzEQ,dzbEQ); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  ice_rho,ice_c,ice_k = %10.3e,%10.3e,%10.3e\n",
                     ice_rho,ice_c,ice_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  bed_rho,bed_c,bed_k = %10.3e,%10.3e,%10.3e\n",
                     bed_rho,bed_c,bed_k); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  nuEQ,dzav = %10.3e,%10.3e\n",
                     nuEQ,dzav); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  iceK,iceRcold,iceRtemp = %10.3e,%10.3e,%10.3e,\n",
		     iceK,iceRcold,iceRtemp); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  bedK,bedR = %10.3e,%10.3e\n",
		     bedK,bedR); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "schemeParamsValid,surfBCsValid,basalBCsValid = (%d,%d,%d)\n"
                     "for THIS column:\n",
                     schemeParamsValid,surfBCsValid,basalBCsValid); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  i,j,ks = %d,%d,%d\n",
                     i,j,ks); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
                     "  Enth_ks,Ghf,Enth_shelfbase,Rb = %10.3e,%10.3e,%10.3e,%10.3e\n",
                     Enth_ks,Ghf,Enth_shelfbase,Rb); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,">>\n\n"); CHKERRQ(ierr);
  return 0;
}


//! Solve the tridiagonal system which determines the new values, in a single
//! \brief column, of the ice enthalpy and the bedrock temperature.
/*!
See the page documenting BOMBPROOF.  We implement equations FIXME: 
(\ref bombtwo), (\ref bedrockeqn),
(\ref geothermalbedeqn), (\ref icebedfinalcold), (\ref icebedfinaltemperate),
(\ref neartopofbedrock), and (\ref icebasenobedrock).
 */
PetscErrorCode enthSystemCtx::solveThisColumn(PetscScalar **x) {
  PetscErrorCode ierr;
  if (!initAllDone) {  SETERRQ(2,
     "solveThisColumn() should only be called after initAllColumns() in enthSystemCtx"); }
  if (!schemeParamsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setSchemeParamsThisColumn() in enthSystemCtx"); }
  if (!surfBCsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setSurfaceBoundaryValuesThisColumn() in enthSystemCtx"); }
  if (!basalBCsValid) {  SETERRQ(3,
     "solveThisColumn() should only be called after setBasalBoundaryValuesThisColumn() in enthSystemCtx"); }

/* PRINCIPLES ABOUT THESE MODIFICATIONS OF tempSystemCtx::solveThisColumn(): 
1)  coefficients in system are unitless and therefore most D,L,U expressions are not altered
2)  old temperature equation had units of Kelvin on each side; new equation has units
    of enthalpy, namely J kg-1, on each side
3)  item 2) means all rhs[] expressions must be modified, and expressions not proportional
    to a temperature or enthalpy generally are multiplied by c (which has units J kg-1 K-1)
*/

  if (Mbz > 1) { // bedrock present: build k=0:Mbz-2 eqns

    // case where Mbz==2 is different; it is stopped by IceEnthalpy model

    // should gives O(\Delta t,\Delta z^2) convergence
    // note L[0] not an allocated location:
    D[0] = (1.0 + 2.0 * bedR);
    U[0] = - 2.0 * bedR;  
    rhs[0] = Tb[0] + 2.0 * dtTemp * Ghf / (bed_rho * bed_c * dzbEQ);
    //rhs[0] = Enth_b[0] + 2.0 * dtTemp * Ghf / (bed_thermal_rho * dzbEQ);

    // bedrock only; pure vertical conduction problem
    // because Mbz >=3, k0 = Mbz-1 >= 2, and the loop goes at least once
    for (PetscInt k = 1; k < k0; k++) {
      L[k] = -bedR;
      D[k] = 1.0 + 2.0 * bedR;
      if (k == k0-1) {
        U[k] = -bedR * (bed_c / ice_c);
      } else {
        U[k] = -bedR;
      }
      rhs[k] = Tb[k];
      //rhs[k] = Enth_b[k];
    }
  }

  // bottom part of ice (and top of bedrock in some cases): k=k0=Mbz-1 eqn
  if (ks == 0) { // no ice; set Enth[0] to surface enthalpy if grounded
    if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
    D[k0] = 1.0;
    U[k0] = 0.0;
    // if floating and no ice then worry only about bedrock temps
    if (mask >= MASK_FLOATING) {
      // essentially no ice but floating ... ask PISMOceanCoupler
      rhs[k0] = Enth_shelfbase;
    } else { // top of bedrock sees atmosphere
      rhs[k0] = Enth_ks; 
    }
  } else { // ks > 0; there is ice
    if (mask >= MASK_FLOATING) {
      // just apply Dirichlet condition to base of column of ice in an ice shelf
      if (k0 > 0) { L[k0] = 0.0; } // note L[0] not allocated 
      D[k0] = 1.0;
      U[k0] = 0.0;
      rhs[k0] = Enth_shelfbase;
    } else if (Mbz == 1) {
      // grounded but no bedrock layer; apply geothermal flux here
      // WARNING: subtle consequences of finite volume argument for basal segment;
      // see BOMBPROOF docs
      const PetscScalar R = (Enth[k0] > Enth_s[k0]) ? iceRtemp : iceRcold;
      // L[k0] = 0.0;  (note this is not an allocated location!) 
      D[k0] = 1.0 + 2.0 * R;
      U[k0] = - 2.0 * R;
      if (w[0] < 0.0) { // velocity downward: upwind vertical
        const PetscScalar AA = nuEQ * w[0];
        D[k0] -= AA;
        U[k0] += AA;
      }
      rhs[k0] = Enth[k0];
      rhs[k0] += (2.0 * nuEQ / ice_rho) * (Ghf + 0.5 * Rb); // geothermal and half of frictional heat
      if (!isMarginal) {
        planeStar ss;
        ierr = Enth3->getPlaneStarZ(i,j,0.0,&ss);
        const PetscScalar UpEnthu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                                 u[0] * (ss.ij  - ss.im1) / dx;
        const PetscScalar UpEnthv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                                 v[0] * (ss.ij  - ss.jm1) / dy;
        rhs[k0] -= dtTemp * (UpEnthu + UpEnthv);     // hor. advection
        rhs[k0] += dtTemp * Sigma[0] / ice_rho;      // strain heat
      }
    } else { 
      // there is *grounded* ice AND there is bedrock at interface
      // WARNING: subtle consequences of finite volume argument across interface;
      // in this segment even temperate ice is conductive (same value as cold ice)
      // see BOMBPROOF docs
      // note sure how to proceed here, keep cold value or replace with temperate?
      const PetscScalar rho_ratio  = ice_rho / bed_rho,
                        c_ratioINV = bed_c / ice_c;
      L[k0] = - 2.0 * bedR;
      if (Enth[k0] > Enth_s[k0]) {
        D[k0] = rho_ratio * (1.0 + 2.0 * iceRcold);
      } else {
        D[k0] = (1.0 + 2.0 * bedR) * c_ratioINV + rho_ratio * (1.0 + 2.0 * iceRcold);
      }
      U[k0] = - rho_ratio * 2.0 * iceRcold;
      if (w[0] < 0.0) { // velocity downward: upwind vertical
        const PetscScalar AA = rho_ratio * dtTemp * w[0] / dzav;
        D[k0] -= AA;
        U[k0] += AA;
      }
      if (Enth[k0] > Enth_s[k0]) { // decide cold vs temperate using prev enthalpy
        rhs[k0] = rho_ratio * Enth[0] - 2 * bedR * c_ratioINV * Enth_s[k0];
      } else {
        rhs[k0] = (rho_ratio + c_ratioINV) * Enth[0];
      }
      rhs[k0] += 2.0 * dtTemp * Rb / (bed_rho * dzav); // frictional heat
      if (!isMarginal) {
        planeStar ss;
        ierr = Enth3->getPlaneStarZ(i,j,0.0,&ss);
        const PetscScalar UpEnthu = (u[0] < 0) ? u[0] * (ss.ip1 -  ss.ij) / dx :
                                                 u[0] * (ss.ij  - ss.im1) / dx;
        const PetscScalar UpEnthv = (v[0] < 0) ? v[0] * (ss.jp1 -  ss.ij) / dy :
                                                 v[0] * (ss.ij  - ss.jm1) / dy;
        rhs[k0] -= dtTemp * rho_ratio * (UpEnthu + UpEnthv); // hor. advection
        rhs[k0] += dtTemp * Sigma[0] / bed_rho;      // strain heat
      }
    }
  }

  // generic ice segment: build k0+1:k0+ks-1 eqns
  for (PetscInt k = 1; k < ks; k++) {
    const PetscScalar AA = nuEQ * w[k],
                      // zero conduction if omega > 0  <==>  E > E_s(p) :
                      R = (Enth[k] > Enth_s[k]) ? iceRtemp : iceRcold;
    if (w[k] >= 0.0) {  // velocity upward
      L[k0+k] = - R - AA * (1.0 - lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * R + AA * (1.0 - lambda);
      U[k0+k] = - R + AA * (lambda/2.0);
    } else {  // velocity downward
      L[k0+k] = - R - AA * (lambda/2.0);
      D[k0+k] = 1.0 + 2.0 * R - AA * (1.0 - lambda);
      U[k0+k] = - R + AA * (1.0 - lambda/2.0);
    }
    rhs[k0+k] = Enth[k];
    if (!isMarginal) {
      planeStar ss;
      ierr = Enth3->getPlaneStarZ(i,j,k * dzEQ,&ss);
      const PetscScalar UpEnthu = (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx :
                                               u[k] * (ss.ij  - ss.im1) / dx;
      const PetscScalar UpEnthv = (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy :
                                               v[k] * (ss.ij  - ss.jm1) / dy;
      rhs[k0+k] += dtTemp * ((Sigma[k] / ice_rho) - UpEnthu - UpEnthv);
    }
  }
      
  // surface b.c.
  if (ks>0) {
    L[k0+ks] = 0.0;
    D[k0+ks] = 1.0;
    // ignore U[k0+ks]
    rhs[k0+ks] = Enth_ks;
  }

  // mark column as done
  schemeParamsValid = false;
  surfBCsValid = false;
  basalBCsValid = false;

  // solve it; note melting not addressed yet
  return solveTridiagonalSystem(k0+ks+1,x);
}

