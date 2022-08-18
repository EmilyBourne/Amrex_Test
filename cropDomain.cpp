#include "cropDomain.H"

/*
 * Return a list of the coordinates of the corners of the box
 * provided
 * */
amrex::Vector<amrex::Dim3>
getCorners(const amrex::Box& box)
{
    amrex::Dim3 lbound = amrex::lbound(box);
    amrex::Dim3 ubound = amrex::ubound(box);
    amrex::Vector<amrex::Dim3> corners;
    corners.push_back({lbound.x    , lbound.y    , lbound.z    });
    corners.push_back({ubound.x+1.0, lbound.y    , lbound.z    });
    corners.push_back({lbound.x    , ubound.y+1.0, lbound.z    });
    corners.push_back({ubound.x+1.0, ubound.y+1.0, lbound.z    });
#if (AMREX_SPACEDIM > 2)
    corners.push_back({lbound.x    , lbound.y    , ubound.z+1.0});
    corners.push_back({ubound.x+1.0, lbound.y    , ubound.z+1.0});
    corners.push_back({lbound.x    , ubound.y+1.0, ubound.z+1.0});
    corners.push_back({ubound.x+1.0, ubound.y+1.0, ubound.z+1.0});
#endif
    return corners;
}

/* Determine whether any part of test_box lies within the radial bounds [rmin, rmax]
 */
bool
inDomain(const amrex::Box& test_box, const amrex::Geometry& simulation_geom,
        amrex::Real rmin, amrex::Real rmax)
{
    const auto prob_lo = simulation_geom.ProbLoArray();
    const auto dx      = simulation_geom.CellSizeArray();
    amrex::Vector<amrex::Dim3> corners (getCorners(test_box));
    int position = 0;

    for (auto c : corners)
    {
        // Find the radius at the corner
        amrex::Real x = prob_lo[0] + dx[0]*c.x;
        amrex::Real y = prob_lo[1] + dx[1]*c.y;
#if (AMREX_SPACEDIM > 2)
        amrex::Real z = prob_lo[2] + dx[2]*c.z;
#endif
        amrex::Real xi[AMREX_SPACEDIM];
        amrex::Real corner_radius_squared = x*x + y*y;
#if (AMREX_SPACEDIM > 2)
        corner_radius_squared += z*z;
#endif
        amrex::Real corner_radius = sqrt(corner_radius_squared);

        // Update the position data
        if (corner_radius < rmin)
            position -= 1;
        if (corner_radius > rmax)
            position += 1;
    }

    if (rmin < 0)
        // position is 4 if all corners are > rmax
        // (less than rmin is not possible)
        return position != 4;
    else if (rmax < 0)
        // position is 0 if all corners are < rmin
        // (everything is greater than rmax)
        return position != 0;
    else
        // position is 4 if all corners are > rmax
        // position is -4 if all corners < rmin
        return position != 4 and position != -4;
}

/*
 * Starting from the box described by simulation_domain, refine
 * then crop any resulting boxes which are ourside the requested
 * rmin/rmax
 *
 * Parameters
 * ----------
 *  simulation_domain : The original domain
 *  rb                : The real coordinates of the original domain
 *  is_periodic       : Indicates along which dimensions the domain is periodic
 *  nrefine           : The number of times the final domain should be refined
 *                       beyond the original refinement of the simulation_domain
 *  ncoarsen          : The number of times the simulation_domain should be coarsened
 *                       before deciding which boxes should be kept in the domain
 *  rmin              : The lower radius of the domain of interest
 *  rmax              : The upper radius of the domain of interest
 *
 * */
amrex::BoxArray
cropDomain(const amrex::Box& simulation_domain, const amrex::RealBox& rb,
        const amrex::Array<int,AMREX_SPACEDIM>& is_periodic,
        int nrefine, int ncoarsen, amrex::Real rmin, amrex::Real rmax,
        int ref_ratio)
{
    // If no refinement then exit
    if (rmin < 0 and rmax < 0)
        return amrex::BoxArray({simulation_domain});

    // Find the number of coarsening levels, restricted by the number
    // of cells in the box
    amrex::IntVect sizes(simulation_domain.size());
    int ncells = sizes[0] < sizes[1] ? sizes[0] : sizes[1];
    int max_nlevel = 0;
    while (ncells >>= 1)
        ++max_nlevel;
    ncoarsen = ncoarsen < max_nlevel ? ncoarsen : max_nlevel;

    // Coarsen the domain the requested number of times
    amrex::BoxList final_domain;
    amrex::Box domain = simulation_domain;
    if (ncoarsen > 0)
    {
        for (int i(0); i<ncoarsen; ++i)
            domain.coarsen(ref_ratio);
    }
    else
    {
        for (int i(0); i<-ncoarsen; ++i)
            domain.refine(ref_ratio);
    }

    // Split the coarsened domain into boxes of size 1 which are
    // kept or removed depending on how it is cropped
    amrex::BoxArray possible_boxes(domain);
    possible_boxes.maxSize(1);
    amrex::Geometry simulation_geom(domain, rb, amrex::CoordSys::cartesian, is_periodic);

    // Loop over the boxes and decide if they are in the domain or not
    AMREX_ASSERT(rmin > 0 or rmax > 0);
    for (amrex::Long i(0); i < possible_boxes.size(); ++i)
    {
        const amrex::Box box = possible_boxes[i];
        if (inDomain(box, simulation_geom, rmin, rmax))
            final_domain.push_back(box);
    }

    // For each of the chosen boxes remove the coarsening used to select
    // boxes and refine it the requested amount
    for (auto& f_domain : final_domain)
    {
        for (int j(0); j < (ncoarsen+nrefine); ++j)
        {
            f_domain.refine(ref_ratio);
        }
    }

    return amrex::BoxArray(final_domain);
}
