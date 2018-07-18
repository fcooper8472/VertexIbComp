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

#ifndef TESTGENERATERANDOMFIELDFIGS_HPP_
#define TESTGENERATERANDOMFIELDFIGS_HPP_

#include <vector>

// Needed for the test environment
#include <cxxtest/TestSuite.h>

// From Chaste
#include "ChasteMakeUnique.hpp"
#include "Node.hpp"
#include "OffLatticeRandomFieldGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "SuperellipseGenerator.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"


class TestIbNodeVoronoiFig : public CxxTest::TestSuite
{
public:

    void TestGenerateVoronoiDiagram()
    {
        unsigned num_points = 37;
        double exponent = 8.0;
        double width = 1.0;
        double height = 1.0;
        double bot_left_x = -0.5;
        double bot_left_y = -0.5;

        SuperellipseGenerator gen(num_points, exponent, width, height, bot_left_x, bot_left_y);

        // Get vectors
        std::vector<c_vector<double, 2> > vector_points = gen.GetPointsAsVectors();

        std::vector<c_vector<double, 2>> modified_points;

        for (auto vec : vector_points)
        {
            double dist = 0.05;

            if (vec[0] > 0.0 && vec[1] > 0.0)
            {
                vec[0] += dist;
                vec[1] += dist;
            }
            else if (vec[0] > 0.0 && vec[1] < 0.0)
            {
                vec[0] += dist;
                vec[1] -= dist;
            }
            else if (vec[0] < 0.0 && vec[1] > 0.0)
            {
                vec[0] -= dist;
                vec[1] += dist;
            }
            else
            {
                vec[0] -= dist;
                vec[1] -= dist;
            }

            vec[0] += RandomNumberGenerator::Instance()->NormalRandomDeviate(0.0, 0.01);
            vec[1] += RandomNumberGenerator::Instance()->NormalRandomDeviate(0.0, 0.01);

            modified_points.emplace_back(Create_c_vector(vec[0], vec[1]));
        }


        using boost::polygon::voronoi_builder;
        using boost::polygon::voronoi_diagram;
        using boost_point = boost::polygon::point_data<int>;

        // Convert points to integers for generating the voronoi diagram
        std::vector<boost_point> int_points(modified_points.size());
        std::transform(modified_points.begin(), modified_points.end(), int_points.begin(),
        [](const c_vector<double, 2>& p)->boost_point
        {
            auto x = static_cast<int>(p[0] * INT_MAX);
            auto y = static_cast<int>(p[1] * INT_MAX);
            return boost_point(x, y);
        });

        // Construct the Voronoi tessellation of these 9 x mTotalNumElements points
        voronoi_diagram<double> vd;
        construct_voronoi(int_points.begin(), int_points.end(), &vd);

        // Output the original node locations
        {
            OutputFileHandler results_handler("VertexIbComp/IbNodeVoronoiFigure", false);
            out_stream results_file = results_handler.OpenOutputFile("node_locations.csv");

            for (const auto& node_location : modified_points)
            {
                (*results_file) << node_location[0] << "," << node_location[1] << "\n";
            }

            results_file->close();
        }

        // Output the voronoi vertex locations
        for (const auto& cell : vd.cells())
        {
            const std::string f_name = "vertex_locations_" + std::to_string(cell.source_index()) + ".csv";

            OutputFileHandler results_handler("VertexIbComp/IbNodeVoronoiFigure", false);
            out_stream results_file = results_handler.OpenOutputFile(f_name);

            const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();

            do
            {
                if (edge->is_primary() && edge->is_finite())
                {
                    const double x = edge->vertex0()->x() / INT_MAX;
                    const double y = edge->vertex0()->y() / INT_MAX;
                    std::cout << x << "," << y << "\n";
                    (*results_file) << x << "," << y << "\n";
                }

                edge = edge->next();

            } while (edge != cell.incident_edge());

            results_file->close();
        }
    }
};

#endif /*TESTGENERATERANDOMFIELDFIGS_HPP_*/