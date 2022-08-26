
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>

#include "Poisson.H"
#include "cropDomain.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        int verbose = 1;
        int n_cell = 128;
        int max_grid_size = 32;
        int n_levels = 2;
        int ref_ratio = 2;
        int n_coarsen = 2;
        int solve = 1;
        double rmin = -1;
        double rmax = -1;

        // read parameters
        {
            ParmParse pp;
            pp.query("verbose", verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("ref_ratio", ref_ratio);
            pp.query("n_levels", n_levels);
            pp.query("n_coarsen", n_coarsen);
            pp.query("solve", solve);
            pp.query("rmin", rmin);
            pp.query("rmax", rmax);
        }
        AMREX_ASSERT(n_levels >= 1);
        AMREX_ASSERT(n_coarsen >= 1);

        Vector<Geometry> geom(n_levels);
        Vector<BoxArray> grids(n_levels);
        Vector<DistributionMapping> dmap(n_levels);
        RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};
        Geometry::Setup(&rb, 0, is_periodic.data());
        Box domain_coarse(IntVect{AMREX_D_DECL(0,0,0)},
                   IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
        Box domain = domain_coarse;
        geom[0].define(domain);
        domain.refine(ref_ratio);
        geom[1].define(domain);

        domain = domain_coarse;
        grids[0].define(domain);
        grids[0].maxSize(max_grid_size);
        domain.growLo(Direction::x, -5*n_cell/8);
        domain.refine(ref_ratio);
        grids[1].define(domain);
        grids[1].maxSize(max_grid_size);

        dmap[0].define(grids[0]);
        dmap[1].define(grids[1]);
        //{
        //    Geometry geom_coarse;
        //    geom_coarse.define(domain_coarse, rb, CoordSys::cartesian, is_periodic);
        //    geom.push_back(geom_coarse);

        //    BoxArray grid_coarse(domain_coarse); // define the BoxArray to be a single grid
        //    grid_coarse.maxSize(max_grid_size); // chop domain up into boxes with length max_Grid_size
        //    grids.push_back(grid_coarse);

        //    dmap.push_back(DistributionMapping(grids[0])); // create a processor distribution mapping given the BoxARray
        //}
        //for (int refine_level(1); refine_level < n_levels; ++refine_level)
        //{
        //    Geometry geom_fine;
        //    DistributionMapping dmap_fine;
        //    domain.refine(ref_ratio);
        //    geom_fine.define(domain, rb, CoordSys::cartesian, is_periodic);
        //    BoxArray grids_fine(cropDomain(domain_coarse, rb, is_periodic,
        //                refine_level, n_coarsen-refine_level, rmin, rmax));
        //    grids_fine.maxSize(max_grid_size);
        //    dmap_fine.define(grids_fine);
        //    geom.push_back(geom_fine);
        //    grids.push_back(grids_fine);
        //    dmap.push_back(dmap_fine);
        //}

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible
        // build a simple geometry using the "eb2." parameters in the inputs file
        for (int lev(0); lev < n_levels; ++lev)
        {
            EB2::Build(geom[lev], required_coarsening_level, max_coarsening_level);
        }

        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();

        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;

        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};

        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        Vector<const EBFArrayBoxFactory*> factory;
        // charge density and electric potential
        Vector<MultiFab> q;
        Vector<MultiFab> phi;
        for (int i(0); i < n_levels; ++i)
        {
            const EB2::Level& eb_level = eb_is.getLevel(geom[i]);
            factory.push_back(new EBFArrayBoxFactory(eb_level, geom[i], grids[i], dmap[i], ng_ebs, ebs));
            q  .push_back(MultiFab(grids[i], dmap[i], 1, 0, MFInfo(), *factory[i]));
            phi.push_back(MultiFab(grids[i], dmap[i], 1, 0, MFInfo(), *factory[i]));
            InitData(q[i]);
        }

        LPInfo info;

        MLEBABecLap mlebabec(geom,grids,dmap,info,factory);

        // define array of LinOpBCType for domain boundary conditions
        std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
        std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bc_lo[idim] = LinOpBCType::Periodic;
            bc_hi[idim] = LinOpBCType::Periodic;
        }

        // Boundary of the whole domain. This functions must be called,
        // and must be called before other bc functions.
        mlebabec.setDomainBC(bc_lo,bc_hi);

        // operator looks like (ACoef - div BCoef grad) phi = rhs

        Vector<MultiFab> acoef;
        Vector<Array<MultiFab,AMREX_SPACEDIM>> bcoef;
        acoef.resize(n_levels);
        bcoef.resize(n_levels);

        for (int i(0); i < n_levels; ++i)
        {
            // see AMReX_MLLinOp.H for an explanation
            mlebabec.setLevelBC(i, nullptr);

            // set ACoef to zero
            acoef[i].define(grids[i], dmap[i], 1, 0, MFInfo(), *factory[i]);
            acoef[i].setVal(0.);
            mlebabec.setACoeffs(i, acoef[i]);
            // set BCoef to 1.0 (and array of face-centered coefficients)
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                bcoef[i][idim].define(amrex::convert(grids[i],IntVect::TheDimensionVector(idim)), dmap[i], 1, 0, MFInfo(), *factory[i]);
                bcoef[i][idim].setVal(1.0);
            }
            mlebabec.setBCoeffs(i, amrex::GetArrOfConstPtrs(bcoef[i]));

            // think of this beta as the "BCoef" associated with an EB face
            MultiFab beta(grids[i], dmap[i], 1, 0, MFInfo(), *factory[i]);
            beta.setVal(1.);

            // set homogeneous Dirichlet BC for EB
            mlebabec.setEBHomogDirichlet(i,beta);
        }


        // scaling factors; these multiply ACoef and BCoef
        Real ascalar = 0.0;
        Real bscalar = 1.0;
        mlebabec.setScalars(ascalar, bscalar);

        MLMG mlmg(mlebabec);

        // relative and absolute tolerances for linear solve
        const Real tol_rel = 1.e-10;
        const Real tol_abs = 0.0;

        mlmg.setVerbose(verbose);

        // Solve linear system
        for (int i(0); i < n_levels; ++i)
        {
            phi[i].setVal(0.0); // initial guess for phi
        }
        if (solve)
            mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(q), tol_rel, tol_abs);

        // store plotfile variables; q and phi
        Vector<MultiFab> plotfile_mf(n_levels);

        for (int lev(0); lev < n_levels; ++lev)
        {
            plotfile_mf[lev].define(grids[lev], dmap[lev], 2, 0, MFInfo(), *factory[lev]);
            MultiFab::Copy(plotfile_mf[lev],  q[lev],0,0,1,0);
            MultiFab::Copy(plotfile_mf[lev],phi[lev],0,1,1,0);
        }

        EB_WriteMultiLevelPlotfile("plt", n_levels, GetVecOfConstPtrs(plotfile_mf),
                {"q", "phi"}, geom, 0.0, Vector<int>(n_levels, 0),
                            Vector<IntVect>(n_levels, IntVect{ref_ratio}));


        for (int i(0); i < n_levels; ++i)
        {
            delete factory[i];
        }
    }

    amrex::Finalize();
}
