#ifndef TESTCELLSORTINGSIMULATIONS_HPP_
#define TESTCELLSORTINGSIMULATIONS_HPP_


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
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryBoundaryCellWriter.hpp"
#include "ImmersedBoundaryNeighbourNumberWriter.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryMorseMembraneForce.hpp"
#include "ImmersedBoundaryMorseDifferentialAdhesionForce.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeRandomFieldForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "ProgressReporter.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformGridRandomFieldGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "VoronoiImmersedBoundaryMeshGenerator.hpp"

#include "FarhadifarForce.hpp"

#include "FakePetscSetup.hpp"


#include "Debug.hpp"
/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_TIME_TO_STEADY_STATE = 2; //10
static const double M_TIME_FOR_SIMULATION = 20; //100
static const double M_NUM_CELLS_ACROSS = 10; //20 // this ^2 cells
static const double M_CELL_FLUCTUATION = 1.0;

class TestCellSortingLiteratePaper : public AbstractCellBasedTestSuite
{
private:

    /*
     * This is a helper method to randomly label cells add is used in all simulations.
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

    const std::string GenerateSuitableRandomField(const ChasteCuboid<2>& bounding_box, const double length_scale)
    {
        // The box is roughly square, so we make it perfectly so
        const double lower = std::min(bounding_box.rGetLowerCorner()[0], bounding_box.rGetLowerCorner()[1]);
        const double upper = std::min(bounding_box.rGetUpperCorner()[0], bounding_box.rGetUpperCorner()[1]);

        // Create a 20% margin, and round the values to the nearest integer for neatness
        const double extra = 0.2 * (upper - lower);
        const std::array<double, 2> lower_corner = {{std::round(lower - extra), std::round(lower - extra)}};
        const std::array<double, 2> upper_corner = {{std::round(upper + extra), std::round(upper + extra)}};

        const std::array<unsigned, 2> num_grid_pts = {{64u, 64u}};
        const std::array<bool, 2> periodicity = {{false, false}};
        const double trace_proportion = 0.8;

        const double grid_spacing = ((upper - lower) + 2.0 * extra) / num_grid_pts[0];
        PRINT_2_VARIABLES(grid_spacing, length_scale);

        // Generate and cache the random field
        UniformGridRandomFieldGenerator<2> gen(lower_corner, upper_corner, num_grid_pts, periodicity, trace_proportion, length_scale);
        return gen.SaveToCache();
    }

public:

    /*
     * == VM ==
     *
     * Simulate a population of cells exhibiting cell sorting using the
     * Cell Vertex model.
     */
    void xTestVertexMonolayerCellSorting()
    {
        // Create a simple 2D MutableVertexMesh
        VoronoiVertexMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, 3u, 0.5 * std::sqrt(3.0));
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Slows things down but can use a larger timestep and diffusion forces
        //p_mesh->SetCheckForInternalIntersections(true);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
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
        simulator.SetOutputDirectory("VertexIbComp/CellSorting/Vertex");

        // Set time step and end time for simulation
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

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
        p_random_force->SetDiffusionStrength(0.1);
MARK;
        const std::string cached_field = GenerateSuitableRandomField(p_mesh->CalculateBoundingBox(), 0.5);
        p_random_force->SetUpRandomFieldGenerator(cached_field);
MARK;
        simulator.AddForce(p_random_force);
MARK;
        ProgressReporter& r_progress = simulator.rSetUpAndGetProgressReporter();
        r_progress.SetOutputToConsole(true);

        // Run simulation
        simulator.Solve();
MARK;
        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust parameters
        p_random_force->SetDiffusionStrength(0.1*M_CELL_FLUCTUATION);

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);

        // Set the progress reporter
        r_progress = simulator.rSetUpAndGetProgressReporter();
        r_progress.SetOutputToConsole(true);

        simulator.Solve();


        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }


    /*
     * == IB ==
     *
     * Simulate a population of cells exhibiting cell sorting using the
     * Immersed Boundary method.
     */
    void TestIbMonolayerCellSorting()
    {
        // Create a simple 2D Immersed Boundary mesh
        const double dist_between_cells = 0.02;
        const double interaction_dist_multiple = 2.0;

//        const double t_to_steady_state = 10.0;
//        const double t_sim = 100.0;

        VoronoiImmersedBoundaryMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, 20u, 128u, 1.0, dist_between_cells, 0.5);

        ImmersedBoundaryMesh<2,2>* p_mesh = generator.GetMesh();

        // Get average volume of element
        auto GetAverageVolumeOfElement = [&]() -> double
        {
            double a = 0.0;
            for (unsigned i=0; i < p_mesh->GetNumElements(); i++)
            {
                a += p_mesh->GetVolumeOfElement(i);
            }
            return a / p_mesh->GetNumElements();
        };

        // Get total perimeter of elements
        auto GetTotalPerimeterOfElements = [&]() -> double
        {
            double a = 0.0;
            for (unsigned i=0; i < p_mesh->GetNumElements(); i++)
            {
                a += p_mesh->GetSurfaceAreaOfElement(i);
            }
            return a;
        };

        // Set up cells, one for each element
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0 * p_mesh->GetVolumeOfElement(i));
        }

        // Create cell population
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(true);
        cell_population.SetReMeshFrequency(50u);
        cell_population.SetInteractionDistance(interaction_dist_multiple * dist_between_cells);
        cell_population.SetIfPopulationHasActiveSources(true);

        // Set the neighbour distance to the same as the cell population interaction distance
        p_mesh->SetNeighbourDist(cell_population.GetInteractionDistance());

        PRINT_VARIABLE(cell_population.GetInteractionDistance());

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<ImmersedBoundaryBoundaryCellWriter>();
        cell_population.AddCellWriter<ImmersedBoundaryNeighbourNumberWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2, 2>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        // Add main immersed boundary simulation modifier and random noise
        auto p_main_modifier = boost::make_shared<ImmersedBoundarySimulationModifier<2>>();
        p_main_modifier->SetNoiseLengthScale(0.01);
        p_main_modifier->SetNoiseSkip(2u);
        p_main_modifier->SetNoiseStrength(1.0 * 1e10);
        p_main_modifier->SetAdditiveNormalNoise(true);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        auto p_boundary_force = boost::make_shared<ImmersedBoundaryMorseMembraneForce<2>>();
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementWellDepth(1.0 * 1e7);

        const double basic_strength = 1e5;
        auto p_cell_cell_force = boost::make_shared<ImmersedBoundaryMorseDifferentialAdhesionForce<2>>();
        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetRepulsionWellDepth(basic_strength);
        p_cell_cell_force->SetAdhesionAtoAWellDepth(basic_strength);
        p_cell_cell_force->SetAdhesionAtoBWellDepth(0.0);
        p_cell_cell_force->SetAdhesionBtoBWellDepth(basic_strength);
        p_cell_cell_force->SetRestLength(0.5 * 1.0 / interaction_dist_multiple);

        simulator.SetOutputDirectory("VertexIbComp/CellSorting/Ib");

        // Set time step and end time for simulation
        const double dt = 0.01;

        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(1u);
        simulator.SetEndTime(100.0 * dt);


        ProgressReporter& r_progress = simulator.rSetUpAndGetProgressReporter();
        r_progress.SetOutputToConsole(true);

        double vol_at_start = GetAverageVolumeOfElement();
        auto pd = p_mesh->GetPolygonDistribution();

        PRINT_VARIABLE(GetAverageVolumeOfElement() / vol_at_start);
        PRINT_VARIABLE(GetTotalPerimeterOfElements());
        PRINT_VARIABLE(std::accumulate(pd.begin(), pd.end(), 0u));

        // Run simulation
        simulator.Solve();

//        // Now label some cells
//        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
//        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);
//
//        pd = p_mesh->GetPolygonDistribution();
//        PRINT_VARIABLE(GetAverageVolumeOfElement() / vol_at_start);
//        PRINT_VARIABLE(GetTotalPerimeterOfElements());
//        PRINT_VARIABLE(std::accumulate(pd.begin(), pd.end(), 0u));
//
//        auto set_0 = p_mesh->GetNeighbouringElementIndices(0);
//        auto set_5 = p_mesh->GetNeighbouringElementIndices(5);
//
//        PRINT_VARIABLE(set_0.size());
//        PRINT_VARIABLE(set_5.size());
//
//        MARK;
//        for (auto x : set_0)
//        {
//            PRINT_VARIABLE(x);
//        }
//        MARK;
//
////
////        // Adjust parameters
////        p_random_force->SetDiffusionStrength(0.1*M_CELL_FLUCTUATION);
////
//        // Run simulation
//        simulator.SetEndTime(t_to_steady_state + t_sim);
//
//        // Set the progress reporter
//        r_progress = simulator.rSetUpAndGetProgressReporter();
//        r_progress.SetOutputToConsole(true);
//
//        simulator.Solve();
//
//        PRINT_VARIABLE(GetAverageVolumeOfElement() / vol_at_start);
//
//        // Check that the same number of cells
//        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), m_num_cells_across*m_num_cells_across);
//
//        // Test no births or deaths
//        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
//        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }
};

#endif /* TESTCELLSORTINGLITERATEPAPER_HPP_ */
