/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTPOLYGONDISTRIBUTIONROBUSTNESS_HPP_
#define TESTPOLYGONDISTRIBUTIONROBUSTNESS_HPP_

// Needed for the test environment
#include "AbstractCellBasedTestSuite.hpp"

// From Chaste
#include "ImmersedBoundaryMesh.hpp"
#include "RandomNumberGenerator.hpp"

// From this user project
#include "VoronoiImmersedBoundaryMeshGenerator.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestPolygonDistributionRobustness : public AbstractCellBasedTestSuite
{
public:

    void TestRobustnessWithElementGap() throw(Exception)
    {
        // Parameters for investigation
        const unsigned num_runs_per_gap = 10u;
        const unsigned num_gaps = 10u;
        const std::array<double> min_max_gap = {{0.01, 0.1}};  // 1% to 10%

        // Parameters for mesh
        const unsigned num_elements_x = 20u;
        const unsigned num_elements_y = 20u;
        const unsigned num_lloyd_steps = 3u;
        const unsigned num_fluid_mesh_pts = 128u;
        const double max_mesh_size = 0.9;

        // Vector of seeds
        std::vector<unsigned> seeds(num_runs_per_gap);
        for (unsigned i = 0; i < seeds.size(); ++i)
        {
            seeds[i] = i;
        }

        // Vector of gaps
        std::vector<double> gaps(num_gaps);
        for (unsigned i = 0; i < gaps.size(); ++i)
        {
            gaps[i] = min_max_gap[0] + i * (min_max_gap[1] - min_max_gap[0]) / gaps.size();
        }

        PRINT_VECTOR(seeds);
        PRINT_VECTOR(gaps);


        for (const unsigned& seed : seeds)
        {
            for (const double& gap : gaps)
            {
                RandomNumberGenerator::Instance()->Reseed(seed);

                VoronoiImmersedBoundaryMeshGenerator gen(num_elements_x,
                                                         num_elements_y,
                                                         num_lloyd_steps,
                                                         num_fluid_mesh_pts,
                                                         max_mesh_size,
                                                         gap);

                ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

                auto vertex_dist = gen.GetVertexMeshPolygonDistribution();
                auto ib_dist = p_mesh->GetPolygonDistribution();

                PRINT_VECTOR(vertex_dist);
                PRINT_VECTOR(ib_dist);
            }
        }
    }
};

#endif /*TESTPOLYGONDISTRIBUTIONROBUSTNESS_HPP_*/