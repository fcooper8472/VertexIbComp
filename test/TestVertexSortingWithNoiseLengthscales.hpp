#ifndef TESTVERTEXSORTINGWITHNOISELENGTHSCALES_HPP
#define TESTVERTEXSORTINGWITHNOISELENGTHSCALES_HPP


#include "AbstractCellBasedTestSuite.hpp"

#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeRandomFieldForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "ProgressReporter.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UblasCustomFunctions.hpp"
#include "UniformGridRandomFieldGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#include "FakePetscSetup.hpp"


#include "Debug.hpp"

// Global idx to ensure random seed is different for every simulation
static unsigned global_sim_idx = 0u;

class TestCellSortingLiteratePaper : public CxxTest::TestSuite
{
private:

    static constexpr double m_time_to_steady_state = 1.0; //10
    static constexpr double m_time_for_simulation = 100.0; //100
    static constexpr unsigned m_num_cells_across = 10u; //20 // this ^2 cells


    /**
     * Helper method to reset singletons.  Needed as we run multiple simulations in loops.
     */
    void SetupSingletons()
    {
        // Set up what the test suite would do
        SimulationTime::Instance()->SetStartTime(0.0);

        // Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
        RandomNumberGenerator::Instance()->Reseed(++global_sim_idx);
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();
    }

    /**
     * Helper method destroy singletons.  Needed as we run multiple simulations in loops.
     */
    void DestroySingletons()
    {
        // This is from the tearDown method of the test suite
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
        CellPropertyRegistry::Instance()->Clear();
    }

    /**
     * Helper method to randomly label cells in the cell population.
     *
     * @param rCells the list of cells to randomly label
     * @param pLabel the label to give to selected cells
     * @param labelledRatio the target fraction of cells to label
     */
    void RandomlyLabelCells(std::list<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (const auto& p_cell : rCells)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
               p_cell->AddCellProperty(pLabel);
            }
        }
    }

    /**
     * Helper method to generate a suitable random field.  The field is then cached, so for a given set of parameters
     * this method will only take a significant amount of time once.
     *
     * @param lowerCorner the lower corner of the mesh
     * @param upperCorner the upper corner of the mesh
     * @param length_scale the correlation length for the random field
     * @return the path to the cached random field
     */
    const std::string GenerateSuitableRandomField(
        const std::array<double, 2> lowerCorner,
        const std::array<double, 2> upperCorner,
        const double length_scale)
    {
        if (length_scale == 0.0)
        {
            return "";
        }

        const double x_delta = upperCorner[0] - lowerCorner[0];
        const double y_delta = upperCorner[1] - lowerCorner[1];

        TS_ASSERT_DELTA(x_delta, y_delta, DBL_EPSILON);

        const std::array<unsigned, 2> num_grid_pts = {{64u, 64u}};
        const std::array<bool, 2> periodicity = {{true, true}};
        const double trace_proportion = 0.8;

        const double grid_spacing = x_delta / num_grid_pts[0];
        PRINT_2_VARIABLES(grid_spacing, length_scale);

        // Generate and cache the random field
        UniformGridRandomFieldGenerator<2> gen(
                lowerCorner,
                upperCorner,
                num_grid_pts,
                periodicity,
                trace_proportion,
                length_scale
        );

        return gen.SaveToCache();
    }

    /**
     * The actual simulation.  This resets and destroys singletons so that this method can be called multiple times
     * in a loop.
     *
     * @param outputDir the output directory for this simulation
     * @param lengthscale the lengthscale for the random field
     * @param diffusionStrength the strength of the random force added to the simulation
     * @param rearrangementThreshold the cell rearrangement threshold (default 0.01)
     */
    void RunSimulation(const std::string outputDir, const double lengthscale, const double diffusionStrength,
                       const double rearrangementThreshold=0.01)
    {
        SetupSingletons();

        // Create a simple periodic 2D MutableVertexMesh
        VoronoiVertexMeshGenerator generator(m_num_cells_across, m_num_cells_across, 5u, 0.5 * std::sqrt(3.0));

        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        const std::array<double, 2> lower_corner = {0.0, 0.0};
        const std::array<double, 2> upper_corner = {p_mesh->GetWidth(0), p_mesh->GetWidth(1)};

        p_mesh->SetCellRearrangementThreshold(rearrangementThreshold);

        // Slows things down but can use a larger timestep and diffusion forces
        p_mesh->SetCheckForInternalIntersections(false);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        for (const auto& p_cell : cells)
        {
            // Set a target area rather than setting a growth modifier.
            // The modifiers don't work correctly as making very long G1 phases)
            p_cell->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(outputDir);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(UINT_MAX);
        simulator.SetEndTime(m_time_to_steady_state);

        // Set up force law and pass it to the simulation
        auto p_force = boost::make_shared<NagaiHondaDifferentialAdhesionForce<2>>();
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(2.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(20.0);
        simulator.AddForce(p_force);

        // Add some noise to avoid local minimum
        auto p_random_force = boost::make_shared<OffLatticeRandomFieldForce<2>>();
        p_random_force->SetDiffusionStrength(diffusionStrength);

        const std::string cached_field = GenerateSuitableRandomField(lower_corner, upper_corner, lengthscale);
        p_random_force->SetUpRandomFieldGenerator(cached_field);

        simulator.AddForce(p_random_force);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust parameters
        p_random_force->SetDiffusionStrength(diffusionStrength);

        // Run simulation
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(m_time_to_steady_state + m_time_for_simulation);

        // Set the progress reporter
        ProgressReporter& r_progress = simulator.rSetUpAndGetProgressReporter();
        r_progress.SetOutputToConsole(true);

        simulator.Solve();


        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), m_num_cells_across*m_num_cells_across);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        DestroySingletons();
    }

    std::vector<double> Range(const double low, const double high, const double inc) const noexcept
    {
        TS_ASSERT(low < high && inc > 0.0);
        std::vector<double> ret = {low};
        while(ret.back() < high)
        {
            ret.push_back(low + ret.size() * inc);
        }
        return ret;
    }

public:

    /**
     * == Single snapshot image of sorting ==
     */
    void xTestSingleSnapshot()
    {
        const double rearrangement_threshold = 0.1;

        std::stringstream sim_name;
        sim_name << "VertexIbComp/CellSorting/VertexSingleSnapshot/";

        // Note it may be necessary to lengthen the simulation runtime to get better total sorting
        RunSimulation(sim_name.str(), 0.0, 0.5, rearrangement_threshold);
    }

    /**
     * == Zero lengthscale, variable diffusion strength, default rearrangement threshold ==
     */
    void xTestDiffusionStrengthDefaultThreshold()
    {
        const unsigned num_reruns = 20u;
        const double rearrangement_threshold = 0.01;

        for (const auto& diff_str : Range(0.1, 0.9, 0.2))
        {
            for (unsigned i = 0; i < num_reruns; ++i)
            {
                std::stringstream sim_name;
                sim_name << std::setprecision(1) << std::fixed;
                sim_name << "VertexIbComp/CellSorting/DiffStrDefault/" << diff_str << "/" << i;

                RunSimulation(sim_name.str(), 0.0, diff_str, rearrangement_threshold);
            }
        }
    }

    /**
     * == Zero lengthscale, variable diffusion strength, paper rearrangement threshold ==
     */
    void xTestDiffusionStrengthPaperThreshold()
    {
        const unsigned num_reruns = 20u;
        const double rearrangement_threshold = 0.1;

        for (const auto& diff_str : Range(0.1, 0.9, 0.2))
        {
            for (unsigned i = 0; i < num_reruns; ++i)
            {
                std::stringstream sim_name;
                sim_name << std::setprecision(1) << std::fixed;
                sim_name << "VertexIbComp/CellSorting/DiffStrPaper/" << diff_str << "/" << i;

                RunSimulation(sim_name.str(), 0.0, diff_str, rearrangement_threshold);
            }
        }
    }

    /**
     * == Varying lengthscale, fixed diffusion strength, paper rearrangement threshold ==
     */
    void xTestDiffusionLengthscale()
    {
        const unsigned num_reruns = 5u;
        const double rearrangement_threshold = 0.05;
        const double diff_str = 1.0;

        for (const auto& lengthscale : Range(0.0, 2.0, 0.4))
        {
            for (unsigned i = 0; i < num_reruns; ++i)
            {
                std::stringstream sim_name;
                sim_name << std::setprecision(1) << std::fixed;
                sim_name << "VertexIbComp/CellSorting/DiffLengthscaleVM/" << lengthscale << "/" << i;

                RunSimulation(sim_name.str(), lengthscale, diff_str, rearrangement_threshold);
            }
        }
    }
};

#endif /* TESTVERTEXSORTINGWITHNOISELENGTHSCALES_HPP */
