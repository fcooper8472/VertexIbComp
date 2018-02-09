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
#include "Node.hpp"
#include "OffLatticeRandomFieldGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"


class TestPolygonDistributionRobustness : public CxxTest::TestSuite
{
private:

    void GenerateRandomField(const std::vector<Node<2>*>& rNodes,
                             const double lengthScale,
                             const double proportionOfTrace)
    {
        const std::string dir_name = "VertexIbComp/RandomFields/";
        const std::string file_name = std::to_string(lengthScale) + "-" + std::to_string(proportionOfTrace) + ".rf";

        OutputFileHandler results_handler("", false);
        out_stream results_file = results_handler.OpenOutputFile(dir_name + file_name);

//        (*results_file) << "hi\n";

        results_file->close();
    }

public:

    void TestGenerateAndOuputDataToFile()
    {
        // Global params
        RandomNumberGenerator* p_rand = RandomNumberGenerator::Instance();
        const unsigned num_nodes = 500u;

        // Generate some nodes
        std::vector<Node<2>*> nodes(num_nodes);
        for (unsigned node_idx = 0; node_idx < num_nodes; ++node_idx)
        {
            nodes[node_idx] = new Node<2>(node_idx, Create_c_vector(10.0 * p_rand->ranf(), 10.0 * p_rand->ranf()));
        }


        for (const auto& length : {5.0, 1.0, 0.1})
        {
            for (const auto& trace : {0.5, 0.8, 0.95})
            {
                GenerateRandomField(nodes, length, trace);
            }
        }





        for (const auto& node : nodes)
        {
            delete(node);
        }





//        // Change this vaiable to your cached random field (relative to $CHASTE_TEST_OUTPUT)
//        const std::string file_path = "CachedRandomFields/xy_0.000_1.000_2.000_5.000_10_12_0_0_100_0.500.rfg";
//
//        // Change this variable to the number of instances of the random field you want to sample
//        const unsigned num_fields_to_sample = 100u;
//
//        // Generate random field from cache and sample from it
//        UniformGridRandomFieldGenerator<2> gen(file_path);
//
//        OutputFileHandler results_handler("", false);
//        out_stream results_file = results_handler.OpenOutputFile(file_path + ".check");
//
//        // Get a number of random field instances and dump the values to file
//        for (unsigned i = 0; i < num_fields_to_sample; ++i)
//        {
//            const std::vector<double> grf = gen.SampleRandomField();
//            for (const double& var : grf)
//            {
//                (*results_file) << var << '\n';
//            }
//        }
//
//        // Tidy up
//        results_file->close();
    }
};

#endif /*TESTGENERATERANDOMFIELDFIGS_HPP_*/