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
    mpRandomFieldGenerator = our::make_unique<UniformGridRandomFieldGenerator<DIM>>(cachedFieldName);
}
#include "Debug.hpp"
template<unsigned DIM>
void OffLatticeRandomFieldForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (mpRandomFieldGenerator == nullptr)
    {
        EXCEPTION("You must set up the random field generator: call SetUpRandomFieldGenerator()");
    }

    const auto& r_nodes_vec = rCellPopulation.rGetMesh().rGetNodes();

    Timer t;
    t.Reset();

    // We need to sample a random field for each dimension
    std::vector<std::vector<double>> random_fields(DIM);
    std::generate(random_fields.begin(), random_fields.end(), [&](){return mpRandomFieldGenerator->SampleRandomField();});

    const double time_generate = t.GetElapsedTime();
    t.Reset();

    // The multiplicative pre-factor applied to each element from the random field
    const double force_scale_factor = std::sqrt(2.0 * mDiffusionStrength / SimulationTime::Instance()->GetTimeStep());

    for (auto& p_node : r_nodes_vec)
    {
        c_vector<double, DIM> force;
        for (unsigned dim = 0; dim < DIM; ++dim)
        {
            force[dim] = mpRandomFieldGenerator->Interpolate(random_fields[dim], p_node->rGetLocation());
        }

        p_node->AddAppliedForceContribution(force_scale_factor * force);
    }

    const double time_interpolate = t.GetElapsedTime();

    const double pct_generate = 100.0 * time_generate / (time_generate + time_interpolate);
    const double pct_interpolate = 100.0 * time_interpolate / (time_generate + time_interpolate);

    PRINT_2_VARIABLES(pct_generate, pct_interpolate);
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
