#include "OffLatticeRandomFieldForce.hpp"

template<unsigned DIM>
OffLatticeRandomFieldForce<DIM>::OffLatticeRandomFieldForce()
    : AbstractForce<DIM>(),
      mMovementParameter(0.01)
{
}

template <unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::SetUpRandomFieldGenerator(const std::vector<Node<DIM>*>& rNodes,
                                                                const double lengthScale,
                                                                const double traceProportion,
                                                                const double domainGrowthProportion)
{
    std::array<double, DIM> lower_corner;
    std::array<double, DIM> upper_corner;

    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        auto compare_dim = [&](Node<DIM>* a, Node<DIM>* b){return a->rGetLocation()[dim] < b->rGetLocation()[dim];};
        auto min_max = std::minmax_element(rNodes.begin(), rNodes.end(), compare_dim);

        const double lower = min_max.first->rGetLocation[dim];
        const double upper = min_max.second->rGetLocation[dim];
        const double extra_each_end = 0.5 * (domainGrowthProportion - 1.0) * (upper - lower);

        lower_corner[dim] = lower - extra_each_end;
        upper_corner[dim] = upper + extra_each_end;
    }

    // We assume there is no periodicity, for now
    std::array<bool, DIM> periodicity;
    periodicity.fill(false);

    mpRandomFieldGenerator = our::make_unique<OffLatticeRandomFieldGenerator<DIM>(
            lower_corner, upper_corner, periodicity, 1u, lengthScale
    );

    mpRandomFieldGenerator->TuneNumEigenvals(rNodes, traceProportion);
}


template<unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (mpRandomFieldGenerator == nullptr)
    {
        EXCEPTION("You must set up the random field generator: call SetUpRandomFieldGenerator()");
    }

    double dt = SimulationTime::Instance()->GetTimeStep();

    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
                c_vector<double, DIM> force_contribution;
        for (unsigned i=0; i<DIM; i++)
        {
            /*
             * The force on this cell is scaled with the timestep such that when it is
             * used in the discretised equation of motion for the cell, we obtain the
             * correct formula
             *
             * x_new = x_old + sqrt(2*D*dt)*W
             *
             * where W is a standard normal random variable.
             */
            double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

            force_contribution[i] = (sqrt(2.0*mMovementParameter*dt)/dt)*xi;
        }
        node_iter->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OffLatticeRandomFieldForce<1>;
template class OffLatticeRandomFieldForce<2>;
template class OffLatticeRandomFieldForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeRandomFieldForce)
