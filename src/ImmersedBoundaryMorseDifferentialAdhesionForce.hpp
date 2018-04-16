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

#ifndef IMMERSEDBOUNDARYMORSEDIFFERENTIALADHESIONFORCE_HPP_
#define IMMERSEDBOUNDARYMORSEDIFFERENTIALADHESIONFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractImmersedBoundaryForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"

#include <iostream>

/**
 * A force class for use in immersed boundary simulations. This force implements Morse-potential-like links between
 * nodes in adjacent immersed boundaries. https://en.wikipedia.org/wiki/Morse_potential
 * The well depth is a constant interaction strength, the rest length is an equilibrium bond distance, and the well
 * width is a parameter governing the profile of the curve.
 *
 * Further, the class behaves differently in the presence of labelled cells, where different adhesion parameters govern
 * the interactions between labelled/non-labelled cells.
 */
template <unsigned DIM>
class ImmersedBoundaryMorseDifferentialAdhesionForce : public AbstractImmersedBoundaryForce<DIM>
{
private:
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractImmersedBoundaryForce<DIM> >(*this);
        archive& mRepulsionWellDepth;
        archive& mAdhesionAtoAWellDepth;
        archive& mAdhesionAtoBWellDepth;
        archive& mAdhesionBtoBWellDepth;
        archive& mRestLength;
        archive& mWellWidth;
    }

    /** The basic interaction strength for interactions closer than the rest length */
    double mRepulsionWellDepth;

    /** The basic interaction strength between two non-labelled cells */
    double mAdhesionAtoAWellDepth;

    /** The basic interaction strength between one labelled and one non-labelled cell */
    double mAdhesionAtoBWellDepth;

    /** The basic interaction strength between two labelled cells */
    double mAdhesionBtoBWellDepth;

    /** The basic rest length associated with interactions, as a fraction of cell population's interaction distance */
    double mRestLength;

    /** The well width as a fraction of the cell population's interaction distance */
    double mWellWidth;

public:
    /**
     * Constructor.
     */
    ImmersedBoundaryMorseDifferentialAdhesionForce();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryMorseDifferentialAdhesionForce() = default;

    /**
     * Overridden AddImmersedBoundaryForceContribution() method.
     *
     * Calculates the force on each node in the immersed boundary cell population as a result of cell-cell interactions.
     *
     * @param rNodePairs reference to a vector set of node pairs between which to contribute the force
     * @param rCellPopulation reference to the cell population
     */
    void AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
                                              ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputImmersedBoundaryForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile);

    /** @return mRepulsionWellDepth */
    double GetRepulsionWellDepth() const;

    /** @param repulsionWellDepth the new value of mRepulsionWellDepth */
    void SetRepulsionWellDepth(double repulsionWellDepth);

    /** @return mAdhesionAtoAWellDepth */
    double GetAdhesionAtoAWellDepth() const;

    /** @param adhesionAtoAWellDepth the new value of mAdhesionAtoAWellDepth */
    void SetAdhesionAtoAWellDepth(double adhesionAtoAWellDepth);

    /** @return mAdhesionAtoBWellDepth */
    double GetAdhesionAtoBWellDepth() const;

    /** @param adhesionAtoBWellDepth the new value of mAdhesionAtoBWellDepth */
    void SetAdhesionAtoBWellDepth(double adhesionAtoBWellDepth);

    /** @return mAdhesionBtoBWellDepth */
    double GetAdhesionBtoBWellDepth() const;

    /** @param adhesionBtoBWellDepth the new value of mAdhesionBtoBWellDepth */
    void SetAdhesionBtoBWellDepth(double adhesionBtoBWellDepth);

    /** @return mRestLength */
    double GetRestLength() const;

    /** @param restLength the new value of mRestLength */
    void SetRestLength(double restLength);

    /** @return mWellWidth */
    double GetWellWidth() const;

    /** @param wellWidth the new value of mWellWidth */
    void SetWellWidth(double wellWidth);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ImmersedBoundaryMorseDifferentialAdhesionForce)

#endif /*IMMERSEDBOUNDARYMORSEDIFFERENTIALADHESIONFORCE_HPP_*/
