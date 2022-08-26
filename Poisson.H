#ifndef POISSON_H_
#define POISSON_H_

#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_Geometry.H>

void InitData (amrex::Geometry const& a_geom, amrex::MultiFab& State);

#endif
