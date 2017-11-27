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

#include "Timer.hpp"

class TestPolygonDistributionRobustness : public AbstractCellBasedTestSuite
{
public:

    void TestRobustnessWithElementGap() throw(Exception)
    {
        // Parameters for investigation
        const unsigned num_runs_per_gap = 20u;
        const unsigned num_gaps = 11u;
        const std::array<double, 2> min_max_gap = {{0.1, 0.2}};  // 3% to 20%

        // Parameters for mesh
        const unsigned num_elements_x = 15u;
        const unsigned num_elements_y = 15u;
        const unsigned num_lloyd_steps = 3u;
        const unsigned num_fluid_mesh_pts = 256u;
        const double max_mesh_size = 0.9;

        // Comparison lambda for difference between two polygon distributions
        auto abs_difference = [](const std::array<unsigned, 13>& truth, const std::array<unsigned, 13>& compare) -> double
        {
            const unsigned total_elems = std::accumulate(truth.begin(), truth.end(), 0u);
            double cumulative_difference = 0.0;
            for (unsigned i = 0; i < truth.size(); ++i)
            {
                // Being very (probably more than necessary) careful about subtracting unsigned values
                auto t = static_cast<double>(truth[i]);
                auto c = static_cast<double>(compare[i]);
                cumulative_difference += std::fabs(t - c);
            }
            return cumulative_difference / total_elems;
        };

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
            gaps[i] = min_max_gap[0] + i * (min_max_gap[1] - min_max_gap[0]) / (gaps.size() - 1);
        }

        PRINT_VECTOR(seeds);
        PRINT_VECTOR(gaps);


        Timer timer;
        timer.Reset();

        double gen_time = 0.0;
        double dist_time = 0.0;
        double total_time = timer.GetWallTime();

        for (const double& gap : gaps)
        {
            double average_difference = 0.0;
            for (const unsigned& seed : seeds)
            {
                RandomNumberGenerator::Instance()->Reseed(seed);

                timer.Reset();
                VoronoiImmersedBoundaryMeshGenerator gen(num_elements_x,
                                                         num_elements_y,
                                                         num_lloyd_steps,
                                                         num_fluid_mesh_pts,
                                                         max_mesh_size,
                                                         gap);
                gen_time += timer.GetElapsedTime();

                ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
                p_mesh->SetNeighbourDist(0.1);

                auto vertex_dist = gen.GetVertexMeshPolygonDistribution();

                timer.Reset();
                auto ib_dist = p_mesh->GetPolygonDistribution();
                dist_time += timer.GetElapsedTime();
                timer.Reset();

                average_difference += abs_difference(vertex_dist, ib_dist);
            }
            average_difference /= seeds.size();
            PRINT_2_VARIABLES(gap, average_difference);
        }

        total_time = timer.GetWallTime() - total_time;

        PRINT_VARIABLE(gen_time);
        PRINT_VARIABLE(dist_time);
        PRINT_VARIABLE(total_time);
        PRINT_VARIABLE(100.0 * gen_time / total_time);
        PRINT_VARIABLE(100.0 * dist_time / total_time);
    }
};

#endif /*TESTPOLYGONDISTRIBUTIONROBUSTNESS_HPP_*/