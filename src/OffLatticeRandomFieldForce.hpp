#ifndef OFFLATTICERANDOMFIELDFORCE_HPP_
#define OFFLATTICERANDOMFIELDFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "UniformGridRandomFieldGenerator.hpp"

/**
 * A force class to model random cell movement, making use of Gaussian Random Fields.
 */
template<unsigned DIM>
class OffLatticeRandomFieldForce : public AbstractForce<DIM>
{
private:

    /** The strength of the force on each node such that F = sqrt(2 * mDiffusionStrength / dt) * x, for x~N(0,1) */
    double mDiffusionStrength;

    /** An owning pointer to the random field generator that creates appropriate correlation between nodes */
    std::unique_ptr<UniformGridRandomFieldGenerator<DIM>> mpRandomFieldGenerator;

    /**
     * Helper method for AddForceContribution.  Add uncorrelated noise in the case of no random field.
     * @param rCellPopulation the cell population to add noise to
     */
    void AddUncorrelatedForce(AbstractCellPopulation<DIM>& rCellPopulation) const noexcept;

    /**
     * Helper method for AddForceContribution.  Add correlated noise sampled from the random field.
     * @param rCellPopulation the cell population to add noise to
     */
    void AddCorrelatedForce(AbstractCellPopulation<DIM>& rCellPopulation) const noexcept;

    /** Archiving */
    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mDiffusionStrength;
    }

public:

    /** Constructor */
    OffLatticeRandomFieldForce();

    /** Destructor */
    ~OffLatticeRandomFieldForce() = default;

    /**
     * Set up the random field generator, from a cached field.
     *
     * @param cachedFieldName the filename of a cached random field, relative to $CHASTE_TEST_OUTPUT
     */
    void SetUpRandomFieldGenerator(const std::string cachedFieldName);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /** @return mDiffusionStrength */
    double GetDiffusionStrength() const;

    /** @paramm diffusionStrength the new value of mDiffusionStrength */
    void SetDiffusionStrength(double diffusionStrength);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeRandomFieldForce)

#endif /*OFFLATTICERANDOMFIELDFORCE_HPP_*/
