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

#include "VoronoiImmersedBoundaryMeshGenerator.hpp"

#include "ChasteMakeUnique.hpp"
#include "ImmersedBoundaryElement.hpp"
#include "MeshUtilityFunctions.hpp"
#include "MutableVertexMesh.hpp"
#include "Node.hpp"
#include "UblasCustomFunctions.hpp"
#include "VertexElement.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/point.hpp>

VoronoiImmersedBoundaryMeshGenerator::VoronoiImmersedBoundaryMeshGenerator(unsigned numElementsX,
                                                                           unsigned numElementsY,
                                                                           unsigned numRelaxationSteps,
                                                                           unsigned numFluidGridPoints,
                                                                           double maxWidthOrHeightOfMesh,
                                                                           double absoluteGapBetweenElements,
                                                                           double targetNodeSpacingRatio)
        : mpIbMesh(nullptr),
          mpVertexMesh(nullptr),
          mNumElementsX(numElementsX),
          mNumElementsY(numElementsY),
          mNumRelaxationSteps(numRelaxationSteps),
          mNumFluidGridPoints(numFluidGridPoints),
          mMaxWidthOrHeightOfMesh(maxWidthOrHeightOfMesh),
          mAbsoluteGapBetweenElements(absoluteGapBetweenElements),
          mTargetNodeSpacingRatio(targetNodeSpacingRatio)
{
    assert(maxWidthOrHeightOfMesh > 0.0);
    assert(maxWidthOrHeightOfMesh < 1.0);

    assert(absoluteGapBetweenElements > 0.0);
    assert(absoluteGapBetweenElements < 1.0);

    // Scaling necessary to correctly pre-size the vertex mesh to lie in a subset of [0, maxWidthOrHeightOfMesh]^2
    const auto max_x_y = std::max(numElementsX, numElementsY);
    const double target_area = (mMaxWidthOrHeightOfMesh * mMaxWidthOrHeightOfMesh) / (max_x_y * max_x_y);

    // Get a Mutable Vertex Mesh from the existing voronoi generator, and perform a deep copy into mpVertexMesh
    {
        VoronoiVertexMeshGenerator vertex_mesh_gen(mNumElementsX, mNumElementsY, mNumRelaxationSteps, target_area);
        DeepCopyAndRepositionVertexMesh(vertex_mesh_gen.GetMesh());
    }

    // Use mpVertexMesh to create an IB mesh in mpIbMesh with the same geometry
    GenerateImmersedBoundaryMesh();
}

void VoronoiImmersedBoundaryMeshGenerator::DeepCopyAndRepositionVertexMesh(MutableVertexMesh<2, 2>* pMeshToCopy)
{
    assert(pMeshToCopy != nullptr);

    // We need to construct new nodes and elements so we don't have mpVertexMesh sharing data with pMeshToCopy
    std::vector<Node<2>*> new_nodes;
    std::vector<VertexElement<2,2>*> new_elems;

    // Determine the vector necessary to reposition the mesh to the centre of [0,1]x[0,1]
    const bool longer_in_x = mNumElementsY < mNumElementsX;
    const double length_ratio = longer_in_x ? (double) mNumElementsY / mNumElementsX : (double) mNumElementsX / mNumElementsY;

    const double long_margin = 0.5 * (1.0 - mMaxWidthOrHeightOfMesh);
    const double short_margin = 0.5 * (1.0 - (mMaxWidthOrHeightOfMesh * length_ratio));

    const auto correction = longer_in_x ? Create_c_vector(short_margin, long_margin) : Create_c_vector(long_margin, short_margin);

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

        // Create a new node and emplace it at the back of new_nodes, with corrected location
        new_nodes.emplace_back(new Node<2>(copy_index, copy_location + correction, copy_is_boundary));
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

    mpVertexMesh = our::make_unique<MutableVertexMesh<2,2>>(new_nodes, new_elems);
}

void VoronoiImmersedBoundaryMeshGenerator::GenerateImmersedBoundaryMesh()
{
    assert(mpVertexMesh != nullptr);

    // Containers in which to store the new nodes and elements
    std::vector<Node<2>*> new_nodes;
    std::vector<ImmersedBoundaryElement<2,2>*> new_elems;

    // Reposition a c_vector to lie within the unit square
    auto RepositionToUnitSquare = [](const c_vector<double, 2>& a)->c_vector<double, 2>
    {
        double x = a[0];
        double y = a[1];

        while(x < 0.0){x += 1.0;}
        while(x >= 1.0){x -= 1.0;}
        while(y < 0.0){y += 1.0;}
        while(y >= 1.0){y -= 1.0;}

        return Create_c_vector(x, y);
    };

    for (unsigned elem_idx = 0; elem_idx < mpVertexMesh->GetNumElements(); ++elem_idx)
    {
        VertexElement<2, 2>* const p_vertex_elem = mpVertexMesh->GetElement(elem_idx);

        // Get the original locations of nodes in the current element
        std::vector<c_vector<double, 2>> vertex_locations;
        for (unsigned local_idx = 0; local_idx < p_vertex_elem->GetNumNodes(); ++local_idx)
        {
            vertex_locations.emplace_back(p_vertex_elem->GetNode(local_idx)->rGetLocation());
        }

        // Determine the shortest length edge to get a sensible number of steps for performing the reduction
        double shortest_edge_length = DBL_MAX;
        for (unsigned local_idx = 0; local_idx < vertex_locations.size(); ++local_idx)
        {
            const unsigned next_idx = AdvanceMod(local_idx, 1, vertex_locations.size());
            const double edge_length = norm_2(vertex_locations[next_idx] - vertex_locations[local_idx]);

            if (edge_length < shortest_edge_length)
            {
                shortest_edge_length = edge_length;
            }
        }

        // Proceed with repositioning in small increments so as not to run in to problems of multiple intersections.
        // We want to move no more than half the shortest edge length in any given step.
        const auto num_steps = std::lround(std::max(1.0, 2.0 * mAbsoluteGapBetweenElements / shortest_edge_length));
        const double step_dist = mAbsoluteGapBetweenElements / num_steps;

        for (unsigned step = 0; step < num_steps; ++step)
        {
            for (unsigned node_local_idx = 0; node_local_idx < vertex_locations.size(); ++node_local_idx)
            {
                const unsigned next_idx = AdvanceMod(node_local_idx, 1, vertex_locations.size());
                const unsigned prev_idx = AdvanceMod(node_local_idx, -1, vertex_locations.size());

                const c_vector<double, 2>& this_pos = vertex_locations[node_local_idx];
                const c_vector<double, 2>& next_pos = vertex_locations[next_idx];
                const c_vector<double, 2>& prev_pos = vertex_locations[prev_idx];

                // We need to bisect the angle that this node makes with the previous and the next, and walk the node in
                // along the line of bisection
                const c_vector<double, 2> unit_vec_to_next = (next_pos - this_pos) / norm_2(next_pos - this_pos);
                const c_vector<double, 2> unit_vec_to_prev = (prev_pos - this_pos) / norm_2(prev_pos - this_pos);
                const c_vector<double, 2> unit_bisecting_vec = 0.5 * (unit_vec_to_next + unit_vec_to_prev);

                const double angle_of_bisection = std::acos(inner_prod(unit_bisecting_vec, unit_vec_to_prev));
                const double length = step_dist / std::sin(angle_of_bisection);

                vertex_locations[node_local_idx] = this_pos + length * unit_bisecting_vec;
            }

            // Now check for intersections caused by small edges becoming inverted due to movement
            using geom_point = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
            using geom_segment = boost::geometry::model::segment<geom_point>;

            std::vector<geom_segment> line_segments;
            line_segments.reserve(vertex_locations.size());

            for (unsigned idx = 0; idx < vertex_locations.size(); ++idx)
            {
                const c_vector<double, 2>& this_pos = vertex_locations[idx];
                const c_vector<double, 2>& next_pos = vertex_locations[AdvanceMod(idx, 1, vertex_locations.size())];

                line_segments.emplace_back(geom_point(this_pos[0], this_pos[1]),
                                           geom_point(next_pos[0], next_pos[1]));
            }

            bool intersections_this_elem = false;
            for (unsigned segment = 0; segment < line_segments.size(); ++segment)
            {
                unsigned next_next = AdvanceMod(segment, 2, line_segments.size());

                if (boost::geometry::intersects(line_segments[segment], line_segments[next_next]))
                {
                    intersections_this_elem = true;

                    std::vector<geom_point> intersection;
                    boost::geometry::intersection(line_segments[segment], line_segments[next_next], intersection);
                    assert(intersection.size() == 1);

                    // Intersection between segment i and i+2, so segment i+1 needs merging.
                    // Segment i+1 has scaled_location[i+1] & scaled_location[i+2]
                    const unsigned start_idx = AdvanceMod(segment, 1, line_segments.size());
                    const unsigned end_idx = AdvanceMod(segment, 2, line_segments.size());


                    const c_vector<double, 2> merged_location = Create_c_vector(intersection[0].get<0>(),
                                                                                intersection[0].get<1>());
                    vertex_locations[start_idx] = merged_location;
                    vertex_locations[end_idx] = merged_location;
                }
            }

            // Clear out the vertex_locations vector of duplicate values
            if (intersections_this_elem)
            {
                // Erase all but the first consecutive identical element
                auto last = std::unique(vertex_locations.begin(), vertex_locations.end(),
                                        [](const c_vector<double, 2>& a, const c_vector<double, 2>& b)
                                        {
                                            return a[0] == b[0] && a[1] == b[1];
                                        });
                vertex_locations.erase(last, vertex_locations.end());
            }
        }

        // Calculate the IB node locations by evenly spacing along the path defined by the vertex locations
        const bool closed_path = true;
        const bool permute_path = true;  // permute the path to remove any bias from consistent starting location
        const double ideal_node_spacing = mTargetNodeSpacingRatio / mNumFluidGridPoints;
        std::vector<c_vector<double, 2>> ib_node_locations = EvenlySpaceAlongPath(vertex_locations,
                                                                                  closed_path,
                                                                                  permute_path,
                                                                                  0u,
                                                                                  ideal_node_spacing);


        // Evaluate the splines at equally-spaced points to create new nodes for the IB mesh at required locations
        std::vector<Node<2>*> nodes_this_elem;
        nodes_this_elem.reserve(ib_node_locations.size());
        for (const auto& location : ib_node_locations)
        {
            new_nodes.emplace_back(new Node<2>(new_nodes.size(), RepositionToUnitSquare(location), true));
            nodes_this_elem.emplace_back(new_nodes.back());
        }

        // Create the element and set whether it is on the boundary or not
        new_elems.emplace_back(new ImmersedBoundaryElement<2, 2>(new_elems.size(), nodes_this_elem));
        new_elems.back()->SetIsBoundaryElement(p_vertex_elem->IsElementOnBoundary());
    }

    // No use case yet for laminas in this type of simulation
    std::vector<ImmersedBoundaryElement<1,2>*> empty_laminas_vec{};

    mpIbMesh = our::make_unique<ImmersedBoundaryMesh<2, 2>>(new_nodes, new_elems, empty_laminas_vec, mNumFluidGridPoints, mNumFluidGridPoints);

    // Replace the default balancing fluid sources with a source at each vertex
    auto& r_balancing_sources = mpIbMesh->rGetBalancingFluidSources();

    for (auto& source : r_balancing_sources)
    {
        delete(source);
    }
    r_balancing_sources.clear();

    r_balancing_sources.reserve(mpVertexMesh->GetNumNodes());
    for (unsigned node_idx = 0; node_idx < mpVertexMesh->GetNumNodes(); ++node_idx)
    {
        const auto idx = r_balancing_sources.size();
        const auto location = mpVertexMesh->GetNode(node_idx)->rGetLocation();

        r_balancing_sources.emplace_back(new FluidSource<2>(idx, RepositionToUnitSquare(location)));
        r_balancing_sources.back()->SetStrength(0.0);
    }
}

ImmersedBoundaryMesh<2,2>* VoronoiImmersedBoundaryMeshGenerator::GetMesh()
{
    return mpIbMesh.get();
}

MutableVertexMesh<2,2>* VoronoiImmersedBoundaryMeshGenerator::GetMutableVertexMesh()
{
    return mpVertexMesh.get();
}

std::array<unsigned, 13> VoronoiImmersedBoundaryMeshGenerator::GetVertexMeshPolygonDistribution()
{
    std::array<unsigned, 13> polygon_dist = {{0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u}};

    for (auto elem_it = mpVertexMesh->GetElementIteratorBegin(); elem_it != mpVertexMesh->GetElementIteratorEnd(); ++elem_it)
    {
        if (!elem_it->IsElementOnBoundary())
        {
            // Accumulate all 12+ sided shapes
            unsigned num_neighbours = std::min<unsigned>(12u, mpVertexMesh->GetNeighbouringElementIndices(elem_it->GetIndex()).size());
            polygon_dist[num_neighbours]++;
        }
    }

    return polygon_dist;
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
