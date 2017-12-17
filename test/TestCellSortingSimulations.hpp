#ifndef TESTCELLSORTINGSIMULATIONS_HPP_
#define TESTCELLSORTINGSIMULATIONS_HPP_


#include "AbstractCellBasedTestSuite.hpp"

#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeRandomFieldForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "ProgressReporter.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VoronoiVertexMeshGenerator.hpp"

#include "FakePetscSetup.hpp"

/*
 *  This is where you can set parameters to be used in all the simulations.
 */

static const double M_TIME_TO_STEADY_STATE = 8; //10
static const double M_TIME_FOR_SIMULATION = 8; //100
static const double M_NUM_CELLS_ACROSS = 20; //20 // this ^2 cells
static const double M_CELL_FLUCTUATION = 1.0;

class TestCellSortingLiteratePaper : public AbstractCellBasedTestSuite
{
private:

    /*
     * This is a helper method to randomly label cells add is used in all simulations.
     */ 

    void RandomlyLabelCells(std::list<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (std::list<CellPtr>::iterator cell_iter = rCells.begin();
             cell_iter != rCells.end();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
               (*cell_iter)->AddCellProperty(pLabel);
            }
        }
    }

public:

    /*
     * == VM ==
     *
     * Simulate a population of cells exhibiting cell sorting using the
     * Cell Vertex model.
     */
    void TestVertexMonolayerCellSorting() throw (Exception)
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
        simulator.SetOutputDirectory("CellSorting/Vertex");

        // Set time step and end time for simulation
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(20);
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
        p_random_force->SetMovementParameter(0.1);
        simulator.AddForce(p_random_force);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust parameters
        p_random_force->SetMovementParameter(0.1*M_CELL_FLUCTUATION);

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);

        // Set the progress reporter
        ProgressReporter& r_progress = simulator.rSetUpAndGetProgressReporter();
        r_progress.SetOutputToConsole(true);

        simulator.Solve();


        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }
};

#endif /* TESTCELLSORTINGLITERATEPAPER_HPP_ */
