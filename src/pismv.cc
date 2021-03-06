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

static char help[] =
"Ice sheet driver for PISM (SIA and SSA) verification.  Uses exact solutions\n"
"  to various coupled subsystems.  Computes difference between exact solution\n"
"  and numerical solution.  Can also just compute exact solution (-eo).\n"
"  Currently implements tests A, B, C, D, E, F, G, H, K, L.\n\n";

#include <string>

#include "base/util/IceGrid.hh"
#include "base/util/PISMConfig.hh"
#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_options.hh"
#include "verif/iceCompModel.hh"
#include "base/util/Context.hh"
#include "base/util/Logger.hh"
#include "base/util/PISMTime.hh"
#include "base/enthalpyConverter.hh"

using namespace pism;

//! Allocate the PISMV (verification) context. Uses ColdEnthalpyConverter.
Context::Ptr pismv_context(MPI_Comm com, const std::string &prefix) {
  // unit system
  units::System::Ptr sys(new units::System);

  // logger
  Logger::Ptr logger = logger_from_options(com);

  // configuration parameters
  Config::Ptr config = config_from_options(com, *logger, sys);

  config->set_string("calendar", "none");

  set_config_from_options(*config);

  print_config(*logger, 3, *config);

  Time::Ptr time = time_from_options(com, config, sys);

  EnthalpyConverter::Ptr EC = EnthalpyConverter::Ptr(new ColdEnthalpyConverter(*config));

  return Context::Ptr(new Context(com, sys, config, EC, time, logger, prefix));
}

GridParameters pismv_grid_defaults(Config::Ptr config,
                                   char testname) {
  // This sets the defaults for each test; command-line options can override this.

  GridParameters P;

  // use the non-periodic grid:
  P.periodicity = NOT_PERIODIC;
  // equal spacing is the default for all the tests except K
  P.Lx = config->get_double("grid_Lx");
  P.Ly = config->get_double("grid_Ly");

  P.Mx = config->get_double("grid_Mx");
  P.My = config->get_double("grid_My");

  SpacingType spacing = EQUAL;
  double Lz = config->get_double("grid_Lz");
  unsigned int Mz = config->get_double("grid_Mz");

  switch (testname) {
  case 'A':
  case 'E':
    // use 1600km by 1600km by 4000m rectangular domain
    P.Lx = 800e3;
    P.Ly = P.Lx;
    Lz = 4000;
    break;
  case 'B':
  case 'H':
    // use 2400km by 2400km by 4000m rectangular domain
    P.Lx = 1200e3;
    P.Ly = P.Lx;
    Lz = 4000;
    break;
  case 'C':
  case 'D':
    // use 2000km by 2000km by 4000m rectangular domain
    P.Lx = 1000e3;
    P.Ly = P.Lx;
    Lz = 4000;
    break;
  case 'F':
  case 'G':
  case 'L':
    // use 1800km by 1800km by 4000m rectangular domain
    P.Lx = 900e3;
    P.Ly = P.Lx;
    Lz = 4000;
    break;
  case 'K':
  case 'O':
    // use 2000km by 2000km by 4000m rectangular domain, but make truely periodic
    config->set_double("grid_Mbz", 2);
    config->set_double("grid_Lbz", 1000);
    P.Lx = 1000e3;
    P.Ly = P.Lx;
    Lz = 4000;
    P.periodicity = XY_PERIODIC;
    spacing = QUADRATIC;
    break;
  case 'V':
    P.My = 3;             // it's a flow-line setup
    P.Lx = 500e3;            // 500 km long
    P.periodicity = Y_PERIODIC;
    break;
  default:
    throw RuntimeError("desired test not implemented\n");
  }

  P.z = IceGrid::compute_vertical_levels(Lz, Mz, spacing,
                                         config->get_double("grid_lambda"));
  return P;
}

IceGrid::Ptr pismv_grid(Context::Ptr ctx, char testname) {
  options::String input_file("-i", "Specifies a PISM input file");
  options::forbidden("-bootstrap");

  if (input_file.is_set()) {
    Periodicity p = string_to_periodicity(ctx->config()->get_string("grid_periodicity"));

    // get grid from a PISM input file
    std::vector<std::string> names;
    names.push_back("enthalpy");
    names.push_back("temp");

    return IceGrid::FromFile(ctx, input_file, names, p);
  } else {
    // use defaults set by pismv_grid_defaults()
    GridParameters P = pismv_grid_defaults(ctx->config(), testname);
    P.horizontal_size_from_options();
    P.horizontal_extent_from_options();
    P.vertical_grid_from_options(ctx->config());
    P.ownership_ranges_from_options(ctx->size());

    return IceGrid::Ptr(new IceGrid(ctx, P));
  }
}

int main(int argc, char *argv[]) {
  MPI_Comm com = MPI_COMM_WORLD;

  petsc::Initializer petsc(argc, argv, help);
  com = PETSC_COMM_WORLD;
      
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();

    verbPrintf(2, com, "PISMV %s (verification mode)\n",
               PISM_Revision);

    if (options::Bool("-version", "stop after printing print PISM version")) {
      return 0;
    }

    std::string usage =
      "  pismv -test x [-no_report] [-eo] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -test x     SIA-type verification test (x = A|B|C|D|E|F|G|H|K|L)\n"
      "  -no_report  do not give error report at end of run\n"
      "  -eo         do not do numerical run; exact solution only\n"
      "(see User's Manual for tests I and J).\n";

    std::vector<std::string> required;
    required.push_back("-test");

    bool done = show_usage_check_req_opts(com, "pismv", required, usage);
    if (done) {
      return 0;
    }

    Context::Ptr ctx = pismv_context(com, "pismv");
    Config::Ptr config = ctx->config();

    // determine test (and whether to report error)
    std::string testname = options::Keyword("-test", "Specifies PISM verification test",
                                            "A,B,C,D,E,F,G,H,K,L", "A");

    IceGrid::Ptr g = pismv_grid(ctx, testname[0]);

    // actually construct and run one of the derived classes of IceModel
    // run derived class for compensatory source SIA solutions
    // (i.e. compensatory accumulation or compensatory heating)
    IceCompModel m(g, ctx, testname[0]);

    m.init();

    m.run();
    verbPrintf(2,com, "done with run\n");

    m.reportErrors();

    // provide a default output file name if no -o option is given.
    m.writeFiles("unnamed.nc");

    print_unused_parameters(*ctx->log(), 3, *config);
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}

