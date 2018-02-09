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

#ifndef TESTHONEYCOMBSIMULATION_HPP_
#define TESTHONEYCOMBSIMULATION_HPP_

// Needed for the test environment
#include "AbstractCellBasedTestSuite.hpp"

// From Chaste
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FluidSource.hpp"
#include "ImmersedBoundaryBoundaryCellWriter.hpp"
#include "ImmersedBoundaryMorseMembraneForce.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"

// Needed due to design on numerical methods (which aren't used with IB simulations)
#include <boost/make_shared.hpp>
#include "ForwardEulerNumericalMethod.hpp"

// From this user project
#include "VoronoiImmersedBoundaryMeshGenerator.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestThreeRegionSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestVoronoiImmersedBoundarySimulation()
    {
        /*
         * @param numElementsX  The number of elements requested across the mesh
         * @param numElementsY  The number of elements requested up the mesh
         * @param numRelaxationSteps  The number of Lloyd's Relaxation steps in the Voronoi iteration
         * @param numFluidGridPoints  The number of fluid mesh points, which determines the node spacing (with targetNodeSpacingRatio)
         * @param maxWidthOrHeightOfMesh The maximum width or height the mesh may be (default 0.9)
         * @param proportionalGapBetweenElements The gap between elements as a proportion of cell size (default 5%)
         * @param targetNodeSpacingRatio The target ratio of node spacing to fluid mesh spacing (default 0.5)
         */
        VoronoiImmersedBoundaryMeshGenerator gen(10, 10, 3, 128, 0.9, 0.08, 0.5);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

        auto vertex_dist = gen.GetVertexMeshPolygonDistribution();
        auto ib_____dist = p_mesh->GetPolygonDistribution();

        std::cout << p_mesh->GetSpacingRatio() << std::endl;

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(false);
        cell_population.SetInteractionDistance(2.0 * cell_population.GetInteractionDistance());

        cell_population.SetReMeshFrequency(UINT_MAX);
        cell_population.AddCellWriter<ImmersedBoundaryBoundaryCellWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR(ImmersedBoundaryMorseMembraneForce<2>, p_boundary_force);
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementWellDepth(1e7);
        p_boundary_force->SetElementRestLength(0.25);

        // Set simulation properties
        double dt = 0.01;
        simulator.SetOutputDirectory("VertexIbComp/TestVoroniMesh");
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(200.0 * dt);

        simulator.Solve();
    }
};

#endif /*TESTHONEYCOMBSIMULATION_HPP_*/