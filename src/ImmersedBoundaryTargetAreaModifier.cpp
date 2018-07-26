/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "ImmersedBoundaryTargetAreaModifier.hpp"

#include <algorithm>
#include <numeric>

#include "Exception.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryElement.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
void ImmersedBoundaryTargetAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateTargetAreas(rCellPopulation);
}

template<unsigned DIM>
void ImmersedBoundaryTargetAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    if (mMinTargetArea == DOUBLE_UNSET || mMaxTargetArea == DOUBLE_UNSET)
    {
        EXCEPTION("You must set the max and min target areas when using this modifier.");
    }

    if (dynamic_cast<ImmersedBoundaryCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("This modifier is only for use with Immersed Boundary cell populations.");
    }

    // Set initial target area with the current cell volume
    for (const auto& p_cell : rCellPopulation.rGetCells())
    {
        p_cell->GetCellData()->SetItem("target area", rCellPopulation.GetVolumeOfCell(p_cell));
    }

    UpdateTargetAreas(rCellPopulation);
}
#include "Debug.hpp"
template<unsigned DIM>
void ImmersedBoundaryTargetAreaModifier<DIM>::UpdateTargetAreas(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    /*
     * Method outline:
     *
     * Determine how "crowded" each IB element is: this is the ratio of the element's perimeter to its voronoi
     * perimeter.  The voronoi perimeter is the perimeter of the union of the voronoi cells of every node in the
     * element.  The closer this ratio is to 1.0, the more "crowded" this element is, and we assume it should then
     * reduce in volume.
     */

    auto p_cell_population = dynamic_cast<ImmersedBoundaryCellPopulation<DIM>*>(&rCellPopulation);
    ImmersedBoundaryMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();

    std::vector<double> crowding;
    crowding.reserve(r_mesh.GetNumElements());

    for (unsigned elem_idx = 0; elem_idx < r_mesh.GetNumElements(); ++elem_idx)
    {
        crowding.emplace_back(r_mesh.GetVoronoiSurfaceAreaOfElement(elem_idx) / r_mesh.GetSurfaceAreaOfElement(elem_idx));
    }

    const double crowding_mean = std::accumulate(crowding.begin(), crowding.end(), 0.0) / crowding.size();
    const double crowding_std = std::sqrt(
            std::inner_product(crowding.begin(), crowding.end(), crowding.begin(), 0.0) / crowding.size() - crowding_mean * crowding_mean);

    std::vector<double> deviations(crowding.size());
    std::transform(crowding.begin(), crowding.end(), deviations.begin(), [&](const double c){
        return (c - crowding_mean) / crowding_std;
    });

    const double old_total_target_area = GetTotalTargetArea(rCellPopulation);

    // The response rate ought to be independent of timestep
    const double response = mResponseSpeed * SimulationTime::Instance()->GetTimeStep();
    assert(response < 1.0);

    for (unsigned elem_idx = 0; elem_idx < r_mesh.GetNumElements(); ++elem_idx)
    {
        auto p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_idx);

        const double t_area = p_cell->GetCellData()->GetItem("target area") * (1.0 + response * deviations[elem_idx]);
        const double actual_area = t_area < mMinTargetArea ? mMinTargetArea : t_area > mMaxTargetArea ? mMaxTargetArea : t_area;
        p_cell->GetCellData()->SetItem("target area", actual_area);
    }

    const double tweak = (old_total_target_area - GetTotalTargetArea(rCellPopulation)) / p_cell_population->GetNumAllCells();

    for (const auto& p_cell : p_cell_population->rGetCells())
    {
        p_cell->GetCellData()->SetItem("target area", tweak + p_cell->GetCellData()->GetItem("target area"));
    }
}

template<unsigned DIM>
double ImmersedBoundaryTargetAreaModifier<DIM>::GetTotalTargetArea(AbstractCellPopulation<DIM,DIM>& rCellPopulation) const noexcept
{
    return std::accumulate(
        rCellPopulation.rGetCells().begin(), rCellPopulation.rGetCells().end(), 0.0, [](double val, CellPtr p_cell) {
            return val + p_cell->GetCellData()->GetItem("target area");
        });
}

template<unsigned DIM>
double ImmersedBoundaryTargetAreaModifier<DIM>::GetMinTargetArea() const noexcept
{
    return mMinTargetArea;
}

template<unsigned DIM>
double ImmersedBoundaryTargetAreaModifier<DIM>::GetMaxTargetArea() const noexcept
{
    return mMaxTargetArea;
}

template<unsigned DIM>
void ImmersedBoundaryTargetAreaModifier<DIM>::SetMinTargetArea(double minTargetArea) noexcept
{
    assert(minTargetArea >= 0.0);
    mMinTargetArea = minTargetArea;
}

template<unsigned DIM>
void ImmersedBoundaryTargetAreaModifier<DIM>::SetMaxTargetArea(double maxTargetArea) noexcept
{
    assert(maxTargetArea >= 0.0);
    mMaxTargetArea = maxTargetArea;
}

template<unsigned DIM>
void ImmersedBoundaryTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinTargetArea>" << mMinTargetArea << "</MinTargetArea>\n";
    *rParamsFile << "\t\t\t<MaxTargetArea>" << mMaxTargetArea << "</MaxTargetArea>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ImmersedBoundaryTargetAreaModifier<1>;
template class ImmersedBoundaryTargetAreaModifier<2>;
template class ImmersedBoundaryTargetAreaModifier<3>;
