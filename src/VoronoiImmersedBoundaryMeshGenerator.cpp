/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "VoronoiImmersedBoundaryMeshGenerator.hpp"

#include "VoronoiVertexMeshGenerator.hpp"
#include "VertexElement.hpp"
#include "Node.hpp"
#include "MutableVertexMesh.hpp"
#include "UblasCustomFunctions.hpp"

VoronoiImmersedBoundaryMeshGenerator::VoronoiImmersedBoundaryMeshGenerator(unsigned numElementsX,
                                                                           unsigned numElementsY,
                                                                           unsigned numRelaxationSteps,
                                                                           unsigned numFluidGridPoints,
                                                                           double maxWidthOrHeightOfMesh,
                                                                           double proportionalGapBetweenElements,
                                                                           double targetNodeSpacingRatio)
        : mpIbMesh(nullptr),
          mpVertexMesh(nullptr),
          mNumElementsX(numElementsX),
          mNumElementsY(numElementsY),
          mNumRelaxationSteps(numRelaxationSteps),
          mNumFluidGridPoints(numFluidGridPoints),
          mMaxWidthOrHeightOfMesh(maxWidthOrHeightOfMesh),
          mProportionalGapBetweenElements(proportionalGapBetweenElements),
          mTargetNodeSpacingRatio(targetNodeSpacingRatio)
{
    assert(maxWidthOrHeightOfMesh > 0.0);
    assert(maxWidthOrHeightOfMesh < 1.0);

    assert(proportionalGapBetweenElements > 0.0);
    assert(proportionalGapBetweenElements < 1.0);

    VoronoiVertexMeshGenerator vertex_mesh_gen(mNumElementsX, mNumElementsY, mNumRelaxationSteps);
    DeepCopyVertexMesh(vertex_mesh_gen.GetMeshAfterReMesh());

    GenerateImmersedBoundaryMesh();
}

VoronoiImmersedBoundaryMeshGenerator::~VoronoiImmersedBoundaryMeshGenerator()
{
    delete mpIbMesh;
    delete mpVertexMesh;
}

void VoronoiImmersedBoundaryMeshGenerator::DeepCopyVertexMesh(MutableVertexMesh<2, 2>* pMeshToCopy)
{
    assert(pMeshToCopy != nullptr);

    // We need to construct new nodes and elements so we don't have mpVertexMesh sharing data with pMeshToCopy
    std::vector<Node<2>*> new_nodes;
    std::vector<VertexElement<2,2>*> new_elems;

    // Copy nodes
    for (unsigned node_idx = 0 ; node_idx < pMeshToCopy->GetNumNodes() ; node_idx++)
    {
        const Node<2>* p_node_to_copy = pMeshToCopy->GetNode(node_idx);

        // Get all the information about the node we are copying
        const unsigned copy_index = p_node_to_copy->GetIndex();
        const bool copy_is_boundary = p_node_to_copy->IsBoundaryNode();
        const c_vector<double, 2>& copy_location = p_node_to_copy->rGetLocation();

        // We assume the nodes are sequentially ordered. This should be the case as the mesh is from a generator.
        assert(copy_index == node_idx);

        // Create a new node and emplace it at the back of new_nodes
        new_nodes.emplace_back(new Node<2>(copy_index, copy_location, copy_is_boundary));
    }

    // Copy elements
    for (unsigned elem_idx = 0; elem_idx < pMeshToCopy->GetNumElements(); elem_idx++)
    {
        const VertexElement<2,2>* p_elem_to_copy = pMeshToCopy->GetElement(elem_idx);

        // Get the information relating to the element we are copying
        const unsigned copy_index     = p_elem_to_copy->GetIndex();
        const unsigned copy_num_nodes = p_elem_to_copy->GetNumNodes();

        // The vertex element is created from a vector of nodes
        std::vector<Node<2>*> nodes_this_elem;

        // Loop through the nodes in p_elem_to_copy and add the corresponding nodes that we have already copied
        for (unsigned node_local_idx = 0 ; node_local_idx < copy_num_nodes ; node_local_idx++)
        {
            const Node<2>* p_local_node = p_elem_to_copy->GetNode(node_local_idx);
            const unsigned local_node_global_idx = p_local_node->GetIndex();

            nodes_this_elem.emplace_back(new_nodes[local_node_global_idx]);
        }

        // Create a new element and emplace it at the back of new_nodes
        new_elems.emplace_back(new VertexElement<2,2>(copy_index, nodes_this_elem));
    }

    mpVertexMesh = new MutableVertexMesh<2,2>(new_nodes, new_elems);
}

void VoronoiImmersedBoundaryMeshGenerator::GenerateImmersedBoundaryMesh()
{
    assert(mpVertexMesh != nullptr);

    ChasteCuboid<2> bounding_box = mpVertexMesh->CalculateBoundingBox();

    const double vertex_mesh_width = bounding_box.GetWidth(0u);
    const double vertex_mesh_height = bounding_box.GetWidth(1u);
    const double max_size_of_vertex_mesh = std::max(vertex_mesh_width, vertex_mesh_height);
    const double scale_factor = mMaxWidthOrHeightOfMesh / max_size_of_vertex_mesh;

    // Determine the gap needed in order to centre the scaled IB mesh within the bounding box [0,1]x[0,1]
    const double left_gap = 0.5 * (1.0 - (vertex_mesh_width * scale_factor));
    const double bottom_gap = 0.5 * (1.0 - (vertex_mesh_height * scale_factor));

    // The displacement which, when added to the scaled node locations, will centre the new IB mesh in [0,1]x[0,1]
    const auto displacement = Create_c_vector(left_gap + bounding_box.rGetLowerCorner()[0] * scale_factor,
                                              bottom_gap + bounding_box.rGetLowerCorner()[1] * scale_factor);

    // Containers in which to store the new nodes and elements
    std::vector<Node<2>*> new_nodes;
    std::vector<ImmersedBoundaryElement<2,2>*> new_elems;

    for (unsigned elem_idx = 0; elem_idx < mpVertexMesh->GetNumElements(); ++elem_idx)
    {
        VertexElement<2, 2>* const p_elem = mpVertexMesh->GetElement(elem_idx);
        const unsigned num_nodes_elem = p_elem->GetNumNodes();

        // Calculate the scaled locations (not taking into account shrinking to provide a gap between elements)
        std::vector<c_vector<double, 2>> scaled_locations(num_nodes_elem);
        for (unsigned node_local_idx = 0; node_local_idx < num_nodes_elem; ++node_local_idx)
        {
            scaled_locations[node_local_idx] = (scale_factor * p_elem->GetNode(node_local_idx)->rGetLocation()) + displacement;
        }

        const c_vector<double, 2> centre_of_mass = std::accumulate(scaled_locations.begin(),
                                                                   scaled_locations.end(),
                                                                   zero_vector<double>(2)) / scaled_locations.size();

        // Rescale the locations, taking the gap into account
        for (auto& location : scaled_locations)
        {
            location = centre_of_mass + (1.0 - mProportionalGapBetweenElements) * (location - centre_of_mass);
        }

        const double perimeter_estimate = mpVertexMesh->GetSurfaceAreaOfElement(elem_idx) * scale_factor * (1.0 - mProportionalGapBetweenElements);
    }
}

ImmersedBoundaryMesh<2,2>* VoronoiImmersedBoundaryMeshGenerator::GetMesh()
{
    return mpIbMesh;
}

MutableVertexMesh<2,2>* VoronoiImmersedBoundaryMeshGenerator::GetMutableVertexMesh()
{
    return mpVertexMesh;
}

std::vector<double> VoronoiImmersedBoundaryMeshGenerator::GetPolygonDistribution()
{
//    assert(mpMesh != nullptr);
//
//    // Number of elements in the mesh
//    unsigned num_elems = mpMesh->GetNumElements();
//
//    // Store the number of each class of polygons
//    std::vector<unsigned> num_polygons(mMaxExpectedNumSidesPerPolygon - 2, 0);
//
//    // Container to return the polygon distribution
//    std::vector<double> polygon_dist;
//
//    // Loop over elements in the mesh to get the number of each class of polygon
//    for (unsigned elem_idx = 0; elem_idx < num_elems; elem_idx++)
//    {
//        unsigned num_nodes_this_elem = mpMesh->GetElement(elem_idx)->GetNumNodes();
//
//        // All polygons are assumed to have 3, 4, 5, ..., mMaxExpectedNumSidesPerPolygon sides
//        assert(num_nodes_this_elem > 2);
//        assert(num_nodes_this_elem <= mMaxExpectedNumSidesPerPolygon);
//
//        // Increment correct place in counter - triangles in place 0, squares in 1, etc
//        num_polygons[num_nodes_this_elem - 3]++;
//    }
//
//    // Loop over the vector of polygon numbers and calculate the distribution vector to return
//    unsigned elems_accounted_for = 0;
//    for (unsigned polygon = 0 ; polygon < num_polygons.size() ; polygon++)
//    {
//        elems_accounted_for += num_polygons[polygon];
//
//        polygon_dist.push_back(static_cast<double>(num_polygons[polygon]) / static_cast<double>(num_elems));
//
//        // Only fill the vector of polygon distributions to the point where there are none of higher class in the mesh
//        if (elems_accounted_for == num_elems)
//        {
//            break;
//        }
//    }
//
//    return polygon_dist;
    return 0;
}

double VoronoiImmersedBoundaryMeshGenerator::GetAreaCoefficientOfVariation()
{
//    assert(mpMesh != nullptr);
//
//    // Number of elements in the mesh, and check there are at least two
//    unsigned num_elems = mpMesh->GetNumElements();
//    assert(num_elems > 1);
//
//    double var = 0.0;
//
//    // Loop over elements in the mesh to get the contributions to the variance
//    for (unsigned elem_idx = 0 ; elem_idx < num_elems ; elem_idx++)
//    {
//        double deviation = mpMesh->GetVolumeOfElement(elem_idx) - mElementTargetArea;
//        var += deviation * deviation;
//    }
//
//    var /= static_cast<double>(num_elems - 1);
//
//    return sqrt(var) / mElementTargetArea;
    return 0;
}
