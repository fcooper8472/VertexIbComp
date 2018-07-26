#ifndef TESTVERTEXSORTINGWITHNOISELENGTHSCALES_HPP
#define TESTVERTEXSORTINGWITHNOISELENGTHSCALES_HPP


#include <cxxtest/TestSuite.h>

#include "CellId.hpp"
#include "CellLabel.hpp"
#include "CellsGenerator.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMorseDifferentialAdhesionForce.hpp"
#include "ImmersedBoundaryMorseMembraneForce.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryTargetAreaModifier.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "ProgressReporter.hpp"
#include "SmartPointers.hpp"
#include "VoronoiImmersedBoundaryMeshGenerator.hpp"

#include "FakePetscSetup.hpp"

#include "Debug.hpp"

// Global idx to ensure random seed is different for every simulation
static unsigned global_sim_idx = 0u;

class TestCellSortingLiteratePaper : public CxxTest::TestSuite
{
private:

    static constexpr double m_time_to_steady_state = 10.0; //10
    static constexpr double m_time_for_simulation = 250.0; //100
    static constexpr unsigned m_num_cells_across = 6u; //20 // this ^2 cells


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
        const unsigned num_to_label = static_cast<unsigned>(std::round(rCells.size() * labelledRatio));

        std::vector<unsigned> cell_ids;
        for (const auto& p_cell : rCells)
        {
            cell_ids.push_back(p_cell->GetCellId());
        }

        std::random_shuffle(cell_ids.begin(), cell_ids.end());
        std::vector<unsigned> selected_ids {cell_ids.begin(), cell_ids.begin() + num_to_label};
        TS_ASSERT(selected_ids.size() == num_to_label);

        for (const auto& p_cell : rCells)
        {
            if (std::find(selected_ids.begin(), selected_ids.end(), p_cell->GetCellId()) != selected_ids.end())
            {
                p_cell->AddCellProperty(pLabel);
            }
        }
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
    void RunSimulation(const std::string outputDir, const double lengthscale, const double diffusionStrength)
    {
        SetupSingletons();

        // Create a simple 2D Immersed Boundary mesh
        const double dist_between_cells = 0.03;
        const double interaction_dist_multiple = 2.0;

        VoronoiImmersedBoundaryMeshGenerator generator(m_num_cells_across, m_num_cells_across, 20u, 128u, 1.0, dist_between_cells, 0.5);

        ImmersedBoundaryMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each element
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        // Create cell population
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(true);
        cell_population.SetReMeshFrequency(50u);
        cell_population.SetInteractionDistance(interaction_dist_multiple * dist_between_cells);
        cell_population.SetIfPopulationHasActiveSources(true);

        auto GetTotalArea = [&]() {
            double total_area = 0.0;
            for (unsigned i = 0; i < p_mesh->GetNumElements(); ++i)
            {
                total_area += p_mesh->GetVolumeOfElement(i);
            }
            return total_area;
        };

        const double surface_at_start = GetTotalArea();
        PRINT_VARIABLE(surface_at_start);

        // Set the neighbour distance to the same as the cell population interaction distance
        p_mesh->SetNeighbourDist(cell_population.GetInteractionDistance());

        // Set population to output all data to results files
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        std::vector<double> vols;
        for (const auto& p_cell : cell_population.rGetCells())
        {
            vols.push_back(cell_population.GetVolumeOfCell(p_cell));
        }

        auto p_area_modifier = boost::make_shared<ImmersedBoundaryTargetAreaModifier<2>>();
        const double vol_mean = std::accumulate(vols.begin(), vols.end(), 0.0) / vols.size();
//        const double vol_var = std::inner_product(vols.begin(), vols.end(), vols.begin(), 0.0) / vols.size() - vol_mean * vol_mean;
        p_area_modifier->SetMinTargetArea(0.5 * vol_mean);
        p_area_modifier->SetMaxTargetArea(1.5 * vol_mean);
        simulator.AddSimulationModifier(p_area_modifier);

        // Add main immersed boundary simulation modifier and random noise
        auto p_main_modifier = boost::make_shared<ImmersedBoundarySimulationModifier<2>>();
        p_main_modifier->SetNoiseLengthScale(lengthscale);
        p_main_modifier->SetNoiseSkip(2u);
        p_main_modifier->SetNoiseStrength(diffusionStrength);
        p_main_modifier->SetAdditiveNormalNoise(true);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        auto p_boundary_force = boost::make_shared<ImmersedBoundaryMorseMembraneForce<2>>();
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementWellDepth(1.5 * 1e7);

        const double basic_strength = 1.2 * 1e5;
        auto p_cell_cell_force = boost::make_shared<ImmersedBoundaryMorseDifferentialAdhesionForce<2>>();
        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetRepulsionWellDepth(100.0 * basic_strength);
        p_cell_cell_force->SetAdhesionAtoAWellDepth(basic_strength);
        p_cell_cell_force->SetAdhesionAtoBWellDepth(0.25 * basic_strength);
        p_cell_cell_force->SetAdhesionBtoBWellDepth(basic_strength);
        p_cell_cell_force->SetRestLength(0.5 * 1.0 / interaction_dist_multiple);

        simulator.SetOutputDirectory(outputDir);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/30.0);
        simulator.SetSamplingTimestepMultiple(UINT_MAX);
        simulator.SetEndTime(m_time_to_steady_state);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Run simulation
        simulator.SetSamplingTimestepMultiple(60u);
        simulator.SetEndTime(m_time_to_steady_state + m_time_for_simulation);

        // Set the progress reporter
        ProgressReporter& r_progress = simulator.rSetUpAndGetProgressReporter();
        r_progress.SetOutputToConsole(true);

        try
        {
            simulator.Solve();
        }
        catch(const Exception& e)
        {
            std::cout << e.GetMessage() << std::endl;
        }

        const double surface_at_end = GetTotalArea();

        PRINT_VARIABLE(surface_at_end / surface_at_start);

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), m_num_cells_across*m_num_cells_across);

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
     * == Varying lengthscale, fixed diffusion strength, paper rearrangement threshold ==
     */
    void TestDiffusionLengthscale() throw(Exception)
    {
        const unsigned num_reruns = 5u;
        const double diff_str = 5.0 * 1e8;

        for (const auto& lengthscale : Range(0.03, 0.07, 0.02))//{0.01, 0.03, 0.05, 0.07, 0.09})
        {
            for (unsigned rerun = 0; rerun < num_reruns; ++rerun)
            {
                PRINT_3_VARIABLES(lengthscale, rerun, diff_str);
                std::stringstream sim_name;
                sim_name << std::setprecision(4) << std::fixed;
                sim_name << "VertexIbComp/CellSorting/DiffLengthscaleIB/" << lengthscale << "/" << rerun;

                RunSimulation(sim_name.str(), lengthscale, diff_str);
            }
        }
    }
};

#endif /* TESTVERTEXSORTINGWITHNOISELENGTHSCALES_HPP */
