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
#include "Toroidal2dVertexMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"


class TestPolygonDistributionRobustness : public CxxTest::TestSuite
{
private:

    void GenerateAndSampleFromRandomField(const std::vector<Node<2>*>& rNodes,
                                          const double domainSize,
                                          const double lengthScale,
                                          const double proportionOfTrace)
    {
        // Global params for random field generator
        const std::array<double, 2> lower_corner {{0.0, 0.0}};
        const std::array<double, 2> upper_corner {{domainSize, domainSize}};
        const std::array<bool, 2> periodicity {{true, true}};
        const unsigned num_eigenvals = std::lround(rNodes.size() / 2.0);

        OffLatticeRandomFieldGenerator<2> gen(
                lower_corner,
                upper_corner,
                periodicity,
                num_eigenvals,
                lengthScale * domainSize,
                domainSize
        );

        gen.TuneNumEigenvals(rNodes, proportionOfTrace);
        gen.Update(rNodes);

        // Output to file
        const std::string dir_name = "VertexIbComp/RandomFields";
        const std::string file_name = std::to_string(lengthScale) + "-" + std::to_string(proportionOfTrace) + ".rfs";

        OutputFileHandler results_handler(dir_name, false);
        out_stream results_file = results_handler.OpenOutputFile(file_name);

        const unsigned num_fields_to_sample = 500u;
        for (unsigned field_num = 0; field_num < num_fields_to_sample; ++field_num)
        {
            const std::vector<double> field = gen.SampleRandomField();

            for (const auto& sample : field)
            {
                (*results_file) << sample << '\n';
            }
        }

        results_file->close();
    }

    void GenerateAndVisualiseRandomField(const std::vector<Node<2>*>& rNodes,
                                         const double domainSize,
                                         const double lengthScale)
    {
        // Global params for random field generator
        const std::array<double, 2> lower_corner {{0.0, 0.0}};
        const std::array<double, 2> upper_corner {{domainSize, domainSize}};
        const std::array<bool, 2> periodicity {{true, true}};
        const unsigned num_eigenvals = std::lround(rNodes.size() / 2.0);

        OffLatticeRandomFieldGenerator<2> gen(
                lower_corner,
                upper_corner,
                periodicity,
                num_eigenvals,
                lengthScale * domainSize,
                domainSize
        );

        gen.TuneNumEigenvals(rNodes, 0.95);
        gen.Update(rNodes);

        // Output to file
        const std::string dir_name = "VertexIbComp/RandomFields";

        for (const auto& field_num : {0, 1, 2})
        {
            const std::string file_name = std::to_string(lengthScale) + "-" + std::to_string(field_num) + ".rfv";

            OutputFileHandler results_handler(dir_name, false);
            out_stream results_file = results_handler.OpenOutputFile(file_name);

            const std::vector<double> field = gen.SampleRandomField();

            for (unsigned field_val = 0; field_val < field.size(); ++field_val)
            {
                const double x_loc = rNodes[field_val]->rGetLocation()[0];
                const double y_loc = rNodes[field_val]->rGetLocation()[1];
                const double value = field[field_val];

                (*results_file) << x_loc << ',' << y_loc << ',' << value << '\n';
            }

            results_file->close();
        }
    }

    void CalculateNumEigenvalsAfterTuning(const std::vector<std::vector<Node<2>*>>& rNodes,
                                          const double domainSize,
                                          const double lengthScale,
                                          const std::vector<double>& rProportionsOfTrace)
    {
        // Global params for random field generator
        const std::array<double, 2> lower_corner {{0.0, 0.0}};
        const std::array<double, 2> upper_corner {{domainSize, domainSize}};
        const std::array<bool, 2> periodicity {{true, true}};
        const unsigned num_eigenvals = std::lround(rNodes[0].size() / 2.0);

        OffLatticeRandomFieldGenerator<2> gen(
                lower_corner,
                upper_corner,
                periodicity,
                num_eigenvals,
                lengthScale * domainSize,
                domainSize
        );

        // Output to file
        const std::string dir_name = "VertexIbComp/RandomFields";
        const std::string file_name = std::to_string(lengthScale) + ".rfe";

        OutputFileHandler results_handler(dir_name, false);
        out_stream results_file = results_handler.OpenOutputFile(file_name);


        for (const auto& proportion_of_trace : rProportionsOfTrace)
        {
            std::vector<double> eigenval_proportions;

            for (const auto& node_vec : rNodes)
            {
                const unsigned num_evals = gen.TuneNumEigenvals(node_vec, proportion_of_trace);
                eigenval_proportions.push_back(static_cast<double>(num_evals) / node_vec.size());
            }

            const double mean = std::accumulate(eigenval_proportions.begin(),
                                                eigenval_proportions.end(),
                                                0.0) / eigenval_proportions.size();

            const double var = std::inner_product(eigenval_proportions.begin(),
                                                  eigenval_proportions.end(),
                                                  eigenval_proportions.begin(),
                                                  0.0) / eigenval_proportions.size() - mean * mean;

            (*results_file) << std::to_string(proportion_of_trace) << ","
                            << std::to_string(mean) << ","
                            << std::to_string(std::sqrt(std::abs(var))) << '\n';
        }

        results_file->close();
    }

public:

    void TestGenerateAndOuputDataToFiles() throw(Exception)
    {
        // Global params
        const double domain_size = 20.0;
        // Generating multiple domains turned out to have no effect on the number of eigenvalues that need to be calculated
        const double num_meshes_to_create = 1u;

        // Generate some nodes from a realistic vertex mesh geometry
        std::vector<std::vector<Node<2>*>> nodes_vec_vec(num_meshes_to_create);
        {
            auto p_mesh_gen = our::make_unique<VoronoiVertexMeshGenerator>(20, 20, 3, 1.0);

            for(auto& nodes_vec : nodes_vec_vec)
            {
                p_mesh_gen->RefreshSeedsAndRegenerateMesh();
                auto p_mesh = p_mesh_gen->GetToroidalMesh();
                for (unsigned node_idx = 0; node_idx < p_mesh->GetNumNodes(); ++node_idx)
                {
                    nodes_vec.emplace_back(new Node<2>(node_idx, p_mesh->GetNode(node_idx)->rGetLocation()));
                }
            }
        }

        std::cout << "There were " << nodes_vec_vec[0].size() << " nodes generated.\n";

        /*
         * Generate the samples from the random fields
         */
        for (const auto& length_scale : {0.5, 0.1, 0.01})
        {
            for (const auto& trace_prop : {0.5, 0.8, 0.95})
            {
                GenerateAndSampleFromRandomField(nodes[0], domain_size, length_scale, trace_prop);
            }
        }

        /*
         * Generate the visualisations of the random fields
         */
        for (const auto& length_scale : {0.5, 0.1, 0.01})
        {
            GenerateAndVisualiseRandomField(nodes[0], domain_size, length_scale);
        }

        /*
         * Calculate the dependence on trace proportion of the number of eigenvalues
         */
        std::vector<double> proportions_of_trace;
        for (unsigned trace_proportion = 0; trace_proportion < 20u; ++trace_proportion)
        {
            proportions_of_trace.emplace_back(0.98 - 0.02 * trace_proportion);
        }

        std::vector<double> length_scales;
        for (unsigned length_scale = 0; length_scale < 9u; ++length_scale)
        {
            length_scales.emplace_back(0.01 + 0.01 * length_scale);
        }

        for (const auto& length_scale : length_scales)
        {
            CalculateNumEigenvalsAfterTuning(nodes_vec_vec, domain_size, length_scale, proportions_of_trace);
        }

        for (const auto& nodes_vec : nodes_vec_vec)
        {
            for (const auto& node : nodes_vec)
            {
                delete (node);
            }
        }
    }
};

#endif /*TESTGENERATERANDOMFIELDFIGS_HPP_*/