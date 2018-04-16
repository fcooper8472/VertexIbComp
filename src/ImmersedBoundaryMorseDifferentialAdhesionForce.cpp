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

#include "ImmersedBoundaryMorseDifferentialAdhesionForce.hpp"

#include "CellLabel.hpp"

template <unsigned DIM>
ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::ImmersedBoundaryMorseDifferentialAdhesionForce()
        : AbstractImmersedBoundaryForce<DIM>(),
          mRepulsionWellDepth(1e3),
          mAdhesionAtoAWellDepth(1e3),
          mAdhesionAtoBWellDepth(1e3),
          mAdhesionBtoBWellDepth(1e3),
          mRestLength(0.25),
          mWellWidth(0.25)
{
}

template <unsigned DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::AddImmersedBoundaryForceContribution(
        std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // We assume this force is not used with laminas; modifications necessary if laminas are present in the mesh
    assert(rCellPopulation.rGetMesh().GetNumLaminas() == 0u);

    for (const auto& node_pair : rNodePairs)
    {
        // Interactions only exist between pairs of nodes that are not in the same boundary / lamina
        if (rCellPopulation.rGetMesh().NodesInDifferentElementOrLamina(node_pair.first, node_pair.second))
        {
            Node<DIM>* const p_node_a = node_pair.first;
            Node<DIM>* const p_node_b = node_pair.second;

            c_vector<double, DIM> vec_a2b = rCellPopulation.rGetMesh().GetVectorFromAtoB(p_node_a->rGetLocation(),
                                                                                         p_node_b->rGetLocation());
            const double normed_dist = norm_2(vec_a2b);

            // Force non-zero only within interaction distance, by definition
            if (normed_dist < rCellPopulation.GetInteractionDistance())
            {
                const unsigned a_elem_idx = *(p_node_a->ContainingElementsBegin());
                const unsigned b_elem_idx = *(p_node_b->ContainingElementsBegin());

                const bool a_elem_labelled = rCellPopulation.GetCellUsingLocationIndex(a_elem_idx)->template HasCellProperty<CellLabel>();
                const bool b_elem_labelled = rCellPopulation.GetCellUsingLocationIndex(b_elem_idx)->template HasCellProperty<CellLabel>();


                const double node_a_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(a_elem_idx, false);
                const double node_b_elem_spacing = rCellPopulation.rGetMesh().GetAverageNodeSpacingOfElement(b_elem_idx, false);

                const double elem_spacing = 0.5 * (node_a_elem_spacing + node_b_elem_spacing);

                const double eff_rest_length = mRestLength * rCellPopulation.GetInteractionDistance();
                const double eff_well_width = mWellWidth * rCellPopulation.GetInteractionDistance();


                double eff_well_depth = elem_spacing / rCellPopulation.GetIntrinsicSpacing();

                if (normed_dist < eff_rest_length)
                {
                    eff_well_depth *= mRepulsionWellDepth;
                }
                else if (a_elem_labelled && b_elem_labelled)
                {
                    eff_well_depth *= mAdhesionBtoBWellDepth;
                }
                else if (a_elem_labelled || b_elem_labelled)
                {
                    eff_well_depth *= mAdhesionAtoBWellDepth;
                }
                else
                {
                    eff_well_depth *= mAdhesionAtoAWellDepth;
                }

                const double morse_exp = std::exp((eff_rest_length - normed_dist) / eff_well_width);

                /*
                 * We must scale each applied force by a factor of elem_spacing / local spacing, so that forces
                 * balance when spread to the grid later (where the multiplicative factor is the local spacing)
                 */
                vec_a2b *= 2.0 * eff_well_width * eff_well_depth * morse_exp * (1.0 - morse_exp) / normed_dist;

                c_vector<double, DIM> force_a2b = vec_a2b * (elem_spacing / node_a_elem_spacing);
                p_node_a->AddAppliedForceContribution(force_a2b);

                c_vector<double, DIM> force_b2a = vec_a2b * (-1.0 * elem_spacing / node_b_elem_spacing);
                p_node_b->AddAppliedForceContribution(force_b2a);
            }
        }
    }
}

template <unsigned DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RepulsionWellDepth>" << mRepulsionWellDepth << "</RepulsionWellDepth>\n";
    *rParamsFile << "\t\t\t<AdhesionAtoAWellDepth>" << mAdhesionAtoAWellDepth << "</AdhesionAtoAWellDepth>\n";
    *rParamsFile << "\t\t\t<AdhesionAtoBWellDepth>" << mAdhesionAtoBWellDepth << "</AdhesionAtoBWellDepth>\n";
    *rParamsFile << "\t\t\t<AdhesionBtoBWellDepth>" << mAdhesionBtoBWellDepth << "</AdhesionBtoBWellDepth>\n";
    *rParamsFile << "\t\t\t<RestLength>" << mRestLength << "</RestLength>\n";
    *rParamsFile << "\t\t\t<WellWidth>" << mWellWidth << "</WellWidth>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template<unsigned int DIM>
double ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::GetRepulsionWellDepth() const
{
    return mRepulsionWellDepth;
}

template<unsigned int DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::SetRepulsionWellDepth(double repulsionWellDepth)
{
    mRepulsionWellDepth = repulsionWellDepth;
}

template<unsigned int DIM>
double ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::GetAdhesionAtoAWellDepth() const
{
    return mAdhesionAtoAWellDepth;
}

template<unsigned int DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::SetAdhesionAtoAWellDepth(double adhesionAtoAWellDepth)
{
    mAdhesionAtoAWellDepth = adhesionAtoAWellDepth;
}

template<unsigned int DIM>
double ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::GetAdhesionAtoBWellDepth() const
{
    return mAdhesionAtoBWellDepth;
}

template<unsigned int DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::SetAdhesionAtoBWellDepth(double adhesionAtoBWellDepth)
{
    mAdhesionAtoBWellDepth = adhesionAtoBWellDepth;
}

template<unsigned int DIM>
double ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::GetAdhesionBtoBWellDepth() const
{
    return mAdhesionBtoBWellDepth;
}

template<unsigned int DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::SetAdhesionBtoBWellDepth(double adhesionBtoBWellDepth)
{
    mAdhesionBtoBWellDepth = adhesionBtoBWellDepth;
}

template <unsigned DIM>
double ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::GetRestLength() const
{
    return mRestLength;
}

template <unsigned DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::SetRestLength(double restLength)
{
    mRestLength = restLength;
}

template <unsigned DIM>
double ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::GetWellWidth() const
{
    return mWellWidth;
}

template <unsigned DIM>
void ImmersedBoundaryMorseDifferentialAdhesionForce<DIM>::SetWellWidth(double wellWidth)
{
    mWellWidth = wellWidth;
}

// Explicit instantiation
template class ImmersedBoundaryMorseDifferentialAdhesionForce<1>;
template class ImmersedBoundaryMorseDifferentialAdhesionForce<2>;
template class ImmersedBoundaryMorseDifferentialAdhesionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMorseDifferentialAdhesionForce)
