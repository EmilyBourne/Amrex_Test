#include "Poisson.H"

using namespace amrex;

void InitData (Geometry const& a_geom, MultiFab& State)
{
    const auto prob_lo = a_geom.ProbLoArray();
    const auto dx = a_geom.CellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        const Array4<Real>& q = State.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = prob_lo[0] + dx[0]*(i + 0.5);
            Real y = prob_lo[1] + dx[1]*(j + 0.5);
#if (AMREX_SPACEDIM > 2)
            Real z = prob_lo[2] + dx[2]*(k + 0.5);
#endif

            if (x < 0.5)
                q(i,j,k) = 1.0;
            else
                q(i,j,k) = 0.0;
            //if (i==70 && j==70 && k==0) {
            //    q(i,j,k) = 1.0;
            //} else {
            //    q(i,j,k) = 0.0;
            //}
        });
    }
}
