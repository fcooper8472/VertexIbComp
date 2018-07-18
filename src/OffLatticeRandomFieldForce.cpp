#include "OffLatticeRandomFieldForce.hpp"

#include <algorithm>

#include "ChasteMakeUnique.hpp"

#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"


template<unsigned DIM>
OffLatticeRandomFieldForce<DIM>::OffLatticeRandomFieldForce()
    : AbstractForce<DIM>(),
      mDiffusionStrength(0.01)
{
}

template <unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::SetUpRandomFieldGenerator(const std::string cachedFieldName)
{
    if (cachedFieldName.empty())
    {
        mpRandomFieldGenerator = nullptr;
    }
    else
    {
        mpRandomFieldGenerator = our::make_unique<UniformGridRandomFieldGenerator<DIM>>(cachedFieldName);
    }
}

template<unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // If the field is null, the noise lengthscale is zero and we add uncorrelated random noise
    if (mpRandomFieldGenerator == nullptr)
    {
        AddUncorrelatedForce(rCellPopulation);
    }
    else
    {
        AddCorrelatedForce(rCellPopulation);
    }
}

template<unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::AddUncorrelatedForce(AbstractCellPopulation<DIM>& rCellPopulation) const noexcept
{
    // The multiplicative pre-factor applied to each element from the random field
    const double force_scale_factor = std::sqrt(2.0 * mDiffusionStrength / SimulationTime::Instance()->GetTimeStep());

    for (auto& p_node : rCellPopulation.rGetMesh().rGetNodes())
    {
        c_vector<double, DIM> force;
        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            force[dim] = force_scale_factor * RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
        }

        p_node->AddAppliedForceContribution(force);
    }
}

template<unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::AddCorrelatedForce(AbstractCellPopulation<DIM>& rCellPopulation) const noexcept
{
    // We need to sample a random field for each dimension
    std::vector<std::vector<double>> random_fields(DIM);
    std::generate(random_fields.begin(), random_fields.end(), [&](){return mpRandomFieldGenerator->SampleRandomField();});

    // The multiplicative pre-factor applied to each element from the random field
    const double force_scale_factor = std::sqrt(2.0 * mDiffusionStrength / SimulationTime::Instance()->GetTimeStep());

    for (auto& p_node : rCellPopulation.rGetMesh().rGetNodes())
    {
        c_vector<double, DIM> force;
        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            force[dim] = force_scale_factor * mpRandomFieldGenerator->Interpolate(random_fields[dim], p_node->rGetLocation());
        }

        p_node->AddAppliedForceContribution(force);
    }
}

template<unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DiffusionStrength>" << mDiffusionStrength << "</DiffusionStrength> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned int DIM>
double OffLatticeRandomFieldForce<DIM>::GetDiffusionStrength() const
{
    return mDiffusionStrength;
}

template<unsigned int DIM>
void OffLatticeRandomFieldForce<DIM>::SetDiffusionStrength(double diffusionStrength)
{
    mDiffusionStrength = diffusionStrength;
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
