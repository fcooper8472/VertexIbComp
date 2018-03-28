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

#ifndef TESTVORONOIVSHONEYCOMBMESHGEN_HPP_
#define TESTVORONOIVSHONEYCOMBMESHGEN_HPP_

#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "FarhadifarForce.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "PolygonDistributionWriter.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#include "FakePetscSetup.hpp"

class TestVoronoiVsHoneycombMeshGen : public AbstractCellBasedTestSuite
{
public:



    void TestHoneycombSims() throw(Exception)
    {
        for (unsigned sim = 0; sim < mNumSims; ++sim)
        {
            SimulationTime::Instance()->SetStartTime(0.0);
            CellId::ResetMaxCellId();
            RandomNumberGenerator::Instance()->Reseed(time(NULL));

            RunHoneycombSimulation(sim);

            SimulationTime::Destroy();
            RandomNumberGenerator::Destroy();
            CellPropertyRegistry::Instance()->Clear();
        }
    }

    void TestVoronoiSims() throw(Exception)
    {
        for (unsigned sim = 0; sim < mNumSims; ++sim)
        {
            SimulationTime::Instance()->SetStartTime(0.0);
            CellId::ResetMaxCellId();
            RandomNumberGenerator::Instance()->Reseed(time(NULL));

            RunVoronoiSimulation(sim);

            SimulationTime::Destroy();
            RandomNumberGenerator::Destroy();
            CellPropertyRegistry::Instance()->Clear();
        }
    }

private:

    const unsigned mNumSims = 50u;

    void RunHoneycombSimulation(unsigned simIdx)
    {
        HoneycombVertexMeshGenerator generator(6, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        auto p_proliferative_type = boost::make_shared<TransitCellProliferativeType>();
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_proliferative_type);


        for (auto& p_cell : cells)
        {
            // Reset the cell cycle model for each cell
            auto p_ccm = static_cast<UniformG1GenerationalCellCycleModel*>(p_cell->GetCellCycleModel());
            p_ccm->SetMaxTransitGenerations(0u);

            p_cell->SetCellCycleModel(p_ccm);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Add the polygon distribution writer
        auto p_polygon_distribution_writer = boost::make_shared<PolygonDistributionWriter<2, 2>>();
        cell_population.AddPopulationWriter(p_polygon_distribution_writer);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VoronoiVsHoneycomb/honeycomb" + std::to_string(simIdx));
        simulator.SetEndTime(15.0);

        simulator.SetSamplingTimestepMultiple(10);

        auto p_force = boost::make_shared<FarhadifarForce<2>>();
        simulator.AddForce(p_force);

        auto p_growth_modifier = boost::make_shared<SimpleTargetAreaModifier<2>>();
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();
    }

    void RunVoronoiSimulation(unsigned simIdx)
    {
        VoronoiVertexMeshGenerator generator(6, 6, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        auto p_proliferative_type = boost::make_shared<TransitCellProliferativeType>();
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_proliferative_type);


        for (auto& p_cell : cells)
        {
            // Reset the cell cycle model for each cell
            auto p_ccm = static_cast<UniformG1GenerationalCellCycleModel*>(p_cell->GetCellCycleModel());
            p_ccm->SetMaxTransitGenerations(0u);

            p_cell->SetCellCycleModel(p_ccm);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Add the polygon distribution writer
        auto p_polygon_distribution_writer = boost::make_shared<PolygonDistributionWriter<2, 2>>();
        cell_population.AddPopulationWriter(p_polygon_distribution_writer);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VoronoiVsHoneycomb/voronoi" + std::to_string(simIdx));
        simulator.SetEndTime(15.0);

        simulator.SetSamplingTimestepMultiple(10);

        auto p_force = boost::make_shared<FarhadifarForce<2>>();
        simulator.AddForce(p_force);

        auto p_growth_modifier = boost::make_shared<SimpleTargetAreaModifier<2>>();
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();
    }
};

#endif /* TESTVORONOIVSHONEYCOMBMESHGEN_HPP_ */
