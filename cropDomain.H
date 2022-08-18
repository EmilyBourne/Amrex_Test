#pragma once
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>


amrex::Vector<amrex::Dim3>
getCorners(const amrex::Box& box);


bool
inDomain(const amrex::Box& test_box, const amrex::Geometry& simulation_geom,
        amrex::Real rmin, amrex::Real rmax);


amrex::BoxArray
cropDomain(const amrex::Box& simulation_domain, const amrex::RealBox& rb,
        const amrex::Array<int,AMREX_SPACEDIM>& is_periodic,
        int nrefine, int ncoarsen, amrex::Real rmin, amrex::Real rmax,
        int refRatio = 2);
