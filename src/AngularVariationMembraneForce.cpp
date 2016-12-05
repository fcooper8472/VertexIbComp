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

#include "AngularVariationMembraneForce.hpp"

template <unsigned DIM>
AngularVariationMembraneForce<DIM>::AngularVariationMembraneForce()
        : AbstractImmersedBoundaryForce<DIM>(),
          mpMesh(NULL),
          mSpringConstant(1e6),
          mRestLengthMultiplier(0.5)
{
}

template <unsigned DIM>
AngularVariationMembraneForce<DIM>::~AngularVariationMembraneForce()
{
}

template <unsigned DIM>
void AngularVariationMembraneForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                                                              ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // This is only run once, on the first pass, and performs a set-up of all necessary class members
    if (mpMesh == NULL)
    {
        mpMesh = &(rCellPopulation.rGetMesh());
    }

    /*
     * The following is called each time step, and calculates the forces acting on each element and lamina
     */

    // Used in the calculation of the spring constant
    mIntrinsicSpacingSquared = rCellPopulation.GetIntrinsicSpacing() * rCellPopulation.GetIntrinsicSpacing();

    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = mpMesh->GetElementIteratorBegin();
         elem_it != mpMesh->GetElementIteratorEnd();
         ++elem_it)
    {
        CalculateForcesOnElement(*elem_it);
    }
}

template <unsigned DIM>
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AngularVariationMembraneForce<DIM>::CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>& rElement)
{
    // Get index and number of nodes of current element
    unsigned elem_idx = rElement.GetIndex();
    unsigned num_nodes = rElement.GetNumNodes();

    // Helper variables
    double normed_dist;
    c_vector<double, DIM> aggregate_force;

    // Make a vector to store the force on node i+1 from node i
    std::vector<c_vector<double, DIM> > elastic_force_to_next_node(num_nodes);

    /*
     * Get the node spacing ratio for this element.  The rest length and spring constant are derived from this
     * characteristic length.
     *
     * The spring constant is derived with reference to the intrinsic spacing, so that with different node spacings
     * the user-defined parameters do not have to be updated.
     *
     * The correct factor to increase the spring constant by is (intrinsic spacing / node_spacing)^2.  One factor
     * takes into account the energy considerations of the elastic springs, and the other takes account of the
     * factor of node_spacing used in discretising the force relation.
     */

    double node_spacing;
    double spring_constant;
    double rest_length;

    node_spacing = mpMesh->GetAverageNodeSpacingOfElement(elem_idx, false);

    spring_constant = mSpringConstant * mIntrinsicSpacingSquared / (node_spacing * node_spacing);
    rest_length = mRestLengthMultiplier * node_spacing;

    // Loop over nodes and calculate the force exerted on node i+1 by node i
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Index of the next node, calculated modulo number of nodes in this element
        unsigned next_idx = (node_idx + 1) % num_nodes;

        double modified_spring_constant = spring_constant;
        double modified_rest_length = rest_length;

        // Hooke's law linear spring force
        elastic_force_to_next_node[node_idx] = mpMesh->GetVectorFromAtoB(rElement.GetNodeLocation(node_idx), rElement.GetNodeLocation(next_idx));
        normed_dist = norm_2(elastic_force_to_next_node[node_idx]);
        elastic_force_to_next_node[node_idx] *= modified_spring_constant * (normed_dist - modified_rest_length) / normed_dist;
    }

    // Add the contributions of springs adjacent to each node
    for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // Get index of previous node
        unsigned prev_idx = (node_idx + num_nodes - 1) % num_nodes;

        aggregate_force = elastic_force_to_next_node[node_idx] - elastic_force_to_next_node[prev_idx];

        // Add the aggregate force contribution to the node
        rElement.GetNode(node_idx)->AddAppliedForceContribution(aggregate_force);
    }
}

template <unsigned DIM>
void AngularVariationMembraneForce<DIM>::SetSpringConstant(double springConstant)
{
    mSpringConstant = springConstant;
}

template <unsigned DIM>
double AngularVariationMembraneForce<DIM>::GetSpringConstant()
{
    return mSpringConstant;
}

template <unsigned DIM>
void AngularVariationMembraneForce<DIM>::SetRestLengthMultiplier(double restLengthMultiplier)
{
    mRestLengthMultiplier = restLengthMultiplier;
}

template <unsigned DIM>
double AngularVariationMembraneForce<DIM>::GetRestLengthMultiplier()
{
    return mRestLengthMultiplier;
}

template <unsigned DIM>
void AngularVariationMembraneForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringConstant>" << mSpringConstant << "</SpringConstant>\n";
    *rParamsFile << "\t\t\t<RestLengthMultiplier>" << mRestLengthMultiplier << "</RestLengthMultiplier>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

// Explicit instantiation
template class AngularVariationMembraneForce<1>;
template class AngularVariationMembraneForce<2>;
template class AngularVariationMembraneForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AngularVariationMembraneForce)
