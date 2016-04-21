PISMdev
=======

Parallel Ice Sheet Model (PISM) work of crodehacke

**Note**
* Original code could be found at [github.com/pism/pism](https://github.com/pism/pism)
* Related full documentations at [www.pism-docs.org](http://www.pism-docs.org)
* More details are given below.

The development presented here is supported by the
* [European Research Council](http://erc.europa.eu) under the [European Community's Seventh Framework Programme (FP7/2007-2013)](http://ec.europa.eu/research/fp7/index_en.cfm?pg=documents)/ERC grant agreement 610055 as part of the [Ice2Ice project](http://www.ice2ice.eu)
* Nordic Center of Excellence [eSTICC](https://esticc.net) (_eScience Tool for Investigating Climate Change in northern high latitudes_) funded by [Nordforsk](http://www.nordforsk.org) (grant 57001)

Thanks to various collaborators 

## Below the original information from the PISM source


PISM, a Parallel Ice Sheet Model
================================

The Parallel Ice Sheet Model is an open source, parallel, high-resolution ice sheet model:

* hierarchy of available stress balances
* marine ice sheet physics, dynamic calving fronts
* polythermal, enthalpy-based conservation of energy scheme
* extensible coupling to atmospheric and ocean models
* verification and validation tools
* complete [documentation](http://www.pism-docs.org/) for users and developers
* uses [MPI](http://www-unix.mcs.anl.gov/mpi/) and [PETSc](http://www-unix.mcs.anl.gov/petsc/petsc-as/) for parallel simulations
* reads and writes [CF-compliant](http://cf-pcmdi.llnl.gov/) [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) files

PISM is jointly developed at the [University of Alaska, Fairbanks (UAF)](http://www.uaf.edu/) and the [Potsdam Institute for Climate Impact Research (PIK)](http://www.pik-potsdam.de/).  UAF developers are based in the [Glaciers Group](http://www.gi.alaska.edu/snowice/glaciers/) at the [Geophysical Institute](http://www.gi.alaska.edu).

PISM development is supported by the [NASA Modeling, Analysis, and Prediction program](http://map.nasa.gov/) (grant #NNX13AM16G) and the [NASA Cryospheric Sciences program](http://ice.nasa.gov/) (grant #NNX13AK27G).


Homepage
--------

[www.pism-docs.org](http://www.pism-docs.org/)


Download and Install
--------------------

See [instructions for getting the latest release](http://www.pism-docs.org/wiki/doku.php?id=stable_version).


Generating Documentation
------------------------

See the [INSTALL.md](INSTALL.md) file in this directory.

Contributing
------------

Want to contribute? Great! See [Committing to PISM](http://www.pism-docs.org/wiki/doku.php?id=committing).
