#ifndef OFFLATTICERANDOMFIELDFORCE_HPP_
#define OFFLATTICERANDOMFIELDFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "OffLatticeRandomFieldGenerator.hpp"

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
    std::unique_ptr<OffLatticeRandomFieldGenerator<SPACE_DIM>> mpRandomFieldGenerator;

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
     * Set up the random field generator, using the current nodes as a template.
     *
     * @param rNodes the nodes in the mesh
     * @param lengthScale the length scale for spatial correlation in the random field
     * @param domainGrowthProportion best guess of what proportion larger the mesh will eventually reach
     * @param traceProportion how well the random numbers will be distributed. Defaults to 0.95.
     */
    void SetUpRandomFieldGenerator(const std::vector<Node<DIM>*>& rNodes,
                                   const double lengthScale,
                                   const double domainGrowthProportion,
                                   const double traceProportion=0.95);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    

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
