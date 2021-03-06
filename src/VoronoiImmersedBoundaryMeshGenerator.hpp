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

#ifndef VORONOIIMMERSEDBOUNDARYMESHGENERATOR_HPP_
#define VORONOIIMMERSEDBOUNDARYMESHGENERATOR_HPP_

#include <cmath>
#include <memory>
#include <vector>

#include "ImmersedBoundaryMesh.hpp"
#include "MutableVertexMesh.hpp"

#include <boost/polygon/voronoi.hpp>

/**
 * Mesh generator that creates a 2D Voronoi tessellation using a number of Lloyd's relaxation steps
 * (http://en.wikipedia.org/wiki/Lloyd%27s_algorithm).
 *
 * This generator calls the corresponding vertex generator, then converts the mesh into an Immersed Boundary mesh.
 */
class VoronoiImmersedBoundaryMeshGenerator
{
    friend class TestVoronoiImmersedBoundaryMeshGenerator;

protected:

    /** A pointer to the mesh generated by this class */
    std::unique_ptr<ImmersedBoundaryMesh<2, 2>> mpIbMesh;

    /** A pointer to the mesh underlying vertex mesh from which the immersed boundary mesh is derived */
    std::unique_ptr<MutableVertexMesh<2, 2>> mpVertexMesh;

    /** The number of elements requested across the mesh. */
    unsigned mNumElementsX;

    /** The number of elements requested up the mesh. */
    unsigned mNumElementsY;

    /** The number of Lloyd's relaxation steps requested in the Voronoi iteration. */
    unsigned mNumRelaxationSteps;

    /** The number of fluid grid points in the X and Y directions */
    unsigned mNumFluidGridPoints;

    /** The requested average target area of elements in the mesh. */
    double mMaxWidthOrHeightOfMesh;

    /** The absolute gap distance between elements */
    double mAbsoluteGapBetweenElements;

    /** The target ratio between fluid-mesh spacing and node spacing */
    double mTargetNodeSpacingRatio;

    /**
     * Helper method for the constructor.
     *
     * Perform a deep copy of a mesh generated by a VoronoiVertexMeshGenerator, into mpVertexMesh, in order for the
     * mesh to persist once the generator goes out of scope and deletes its mesh.
     * The mesh is also repositioned to be central within [0,1]x[0,1].
     *
     * @param pMeshToCopy the MutableVertexMesh generated by a VoronoiVertexMeshGenerator.
     */
    void DeepCopyAndRepositionVertexMesh(MutableVertexMesh<2, 2>* pMeshToCopy);

    /**
     * Helper method for the constructor.
     *
     * Use the copied mpVertexMesh to generate a corresponding Immersed Boundary Mesh.
     */
    void GenerateImmersedBoundaryMesh();

public:

    /**
     * Constructor.
     *
     * @param numElementsX  The number of elements requested across the mesh
     * @param numElementsY  The number of elements requested up the mesh
     * @param numRelaxationSteps  The number of Lloyd's Relaxation steps in the Voronoi iteration
     * @param numFluidGridPoints  The number of fluid mesh points, which determines the node spacing (with targetNodeSpacingRatio)
     * @param maxWidthOrHeightOfMesh The maximum width or height the mesh may be (default 0.9)
     * @param absoluteGapBetweenElements The gap between elements (default 0.01)
     * @param targetNodeSpacingRatio The target ratio of node spacing to fluid mesh spacing (default 0.5)
     */
    VoronoiImmersedBoundaryMeshGenerator(unsigned numElementsX,
                                         unsigned numElementsY,
                                         unsigned numRelaxationSteps,
                                         unsigned numFluidGridPoints,
                                         double maxWidthOrHeightOfMesh=0.9,
                                         double absoluteGapBetweenElements=0.01,
                                         double targetNodeSpacingRatio=0.5);

    /**
     * Null constructor for derived classes to call.
     */
    VoronoiImmersedBoundaryMeshGenerator() = default;

    /**
     * Default destructor
     */
    virtual ~VoronoiImmersedBoundaryMeshGenerator() = default;

    /**
     * Return a pointer to the Mutable Vertex Mesh on which the geometry for the Immersed Boundary Mesh is based.
     *
     * @return mpVertexMesh
     */
    virtual MutableVertexMesh<2,2>* GetMutableVertexMesh();

    /**
     * Return a pointer to the Immersed Boundary Mesh generated by this class.
     *
     * @return mpIbMesh
     */
    virtual ImmersedBoundaryMesh<2,2>* GetMesh();

    /**
     * Calculate the polygon distribution for the underlying vertex mesh: number of {0, 1, 2, 3, 4, 5,..., 12+}-gons.
     * Note that the vector will always begin {0, 0, 0, ...} as there can be no 0, 1, or 2-gons, but this choice means
     * that accessing the nth element of the vector gives you the number of n-gons which seems to be most natural.
     * All 12-sided and higher order polygons are accumulated in the array[12] position.
     *
     * @return an array of length 13 representing the polygon distribution.
     */
    std::array<unsigned, 13> GetVertexMeshPolygonDistribution();

    /**
     * Computes the coefficient of variation of the areas of elements in the mesh, defined to be the sample standard
     * deviation in area divided by the mean area.
     *
     * @return The coefficient of variation of the area of elements in the mesh
     */
    double GetAreaCoefficientOfVariation();

};

#endif /*VORONOIIMMERSEDBOUNDARYMESHGENERATOR_HPP_*/
