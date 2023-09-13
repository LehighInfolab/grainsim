// (c) 2023 Benjamin Zalatan Productions
// Keep circulating the tapes

#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

#include "vtk.h"
#include "lattice.h"
#include "octree3.h"
#include "debug_timer.h"
#include "config.h"
#include "analysis.h"

// cd C:\Stuff\School\summer 2023\grainsim
// g++ -O3 CPPGrainSim/main.cpp -o grainsim.out -static

int main(int argc, char *argv[])
{
	// Load the config file.
	config_t cfg;
	cfg.load_config();

	// Create the lattice from file.
	lattice_t *cube;

	// If the scale multiplier is not 1, scale the lattice.
	if (cfg.scale_multiplier != 1)
	{
		lattice_t *temp = vtk::from_file(cfg.initial_state_path.c_str(), false);
		cube = vtk::scale_lattice(temp, cfg.scale_multiplier, false);
		delete temp;
	}
	else
	{
		cube = vtk::from_file(cfg.initial_state_path.c_str(), false);
	}

	lattice_analyzer_t analyze;

	cube->default_mobility = cfg.default_mobility;
	cube->transitioned_mobility = cfg.transitioned_mobility;
	cube->grain_count = cfg.const_grain_count;
	cube->init();

	// Generate the checkpoint list.
	std::vector<double> checkpoints;
	if (!cfg.checkpoints.empty())
	{
		cfg.checkpoints_to_vector(&checkpoints);
	}
	unsigned char curr_checkpoint = 0;

	// Start global timer.
	debug_timer_t timer;
	timer.start();

	double timestep = 0, curr_step, log_duration = 0, transition_duration = 0, next_checkpoint = cfg.checkpoint_interval;
	int vtkcount = 0;

	if (cfg.log_transitions) cube->begin_logging_transitions(cfg.output_folder);

	// Main simulation loop.
	while (true)
	{
		// Flip a voxel and store elapsed timesteps.
		curr_step = cube->step();

		// Update current timestep.
		timestep += curr_step;
		log_duration += curr_step;
		transition_duration += curr_step;

		// Debug logging.
		if (log_duration >= 20000)
		{
			std::cout << "T = " << timestep << ", dT = " << curr_step << ", A = " << cube->system_activity() << ", Flips = " << cube->total_flips << ", tFlips = " << cube->transformed_flips << ", dTime = " << timer.lap() << " sec, tTime = " << timer.total() << " sec" << std::endl;
			log_duration = 0;
		}

		// Transition some boundaries if applicable.
		if (transition_duration >= cfg.transition_interval && cfg.transition_count > 0)
		{
			if (cfg.log_transitions) cube->set_log_timestep(timestep);

			cube->transition_boundaries(cfg.transition_count, cfg.propagation_chance, cfg.propagation_ratio, cfg.use_potential_energy);

			transition_duration = 0;
		}

		// Check if VTK should be generated.
		if (checkpoints.size() > 0 && curr_checkpoint < checkpoints.size() && timestep >= checkpoints[curr_checkpoint]) // The current timestep is an explicit checkpoint.
		{
			std::stringstream ss;
			ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << ".vtk";
			vtk::to_vtk(ss.str().c_str(), cube);
			if (cfg.log_transitions) cube->flush_log_file();

			if (cfg.generate_analysis_files)
			{
				std::cout << "Beginning analysis..." << std::endl;
				analyze.load_lattice(cube);
				ss.str(std::string());
				ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << "_analysis.txt";
				analyze.save_analysis_to_file(ss.str().c_str());
			}

			++vtkcount;
			++curr_checkpoint;

			if (cfg.max_timestep <= 0 && curr_checkpoint >= checkpoints.size()) break;
		}
		else if (cfg.checkpoint_interval > 0 && timestep >= next_checkpoint) // The current timestep surpasses the interval threshold.
		{
			std::stringstream ss;
			ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << ".vtk";
			vtk::to_vtk(ss.str().c_str(), cube);
			if (cfg.log_transitions) cube->flush_log_file();

			if (cfg.generate_analysis_files)
			{
				std::cout << "Beginning analysis..." << std::endl;
				analyze.load_lattice(cube);
				ss.str(std::string());
				ss << cfg.output_folder << cfg.identifier << "_" << std::setw(4) << std::setfill('0') << std::to_string(vtkcount + 1) << '_' << std::to_string((size_t)timestep) << "_analysis.txt";
				analyze.save_analysis_to_file(ss.str().c_str());
			}

			++vtkcount;

			next_checkpoint += cfg.checkpoint_interval;
		}

		// Break if the max timestep is reached.
		if (cfg.max_timestep > 0 && timestep >= cfg.max_timestep) break;
	}

	if (cfg.log_transitions) cube->stop_logging_transitions();










	// Commented out code for octree testing

	/*lattice_t cube(1);
	size_t size = 100;
	unsigned char height = 6;
	voxel_t *vlist = new voxel_t[size * size * size];
	octree3_t tree(128, log2(128));
	for (size_t z = 0; z < size; ++z)
	{
		for (size_t y = 0; y < size; ++y)
			for (size_t x = 0; x < size; ++x) {
				tree.delta(x, y, z, 1);
				vlist[x + y * size + z * size * size].activity = 1;
			}
	}

	for (size_t i = 0; i < 10000; ++i)
	{
		coord_t x, y, z;
		x = cube.rng(0, size);
		y = cube.rng(0, size);
		z = cube.rng(0, size);

		std::cout << "=================" << std::endl;
		std::cout << "real: " << x << ", " << y << ", " << z << std::endl;

		tree.delta(x, y, z, 1);
		vlist[x + y * size + z * size * size].activity = 1;
		std::cout << "setting vindex " << x + y * size + z * size * size << " to 1 " << std::endl;

		coord_t vx, vy, vz;
		tree.get_voxel_from_sum_activity(&vx, &vy, &vz, 1, vlist, 100);

		std::cout << "pred: " << vx << ", " << vy << ", " << vz << std::endl;

		if (x != vx || y != vy || z != vz)
		{
			std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!! MISMATCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}

		tree.delta(x, y, z, -1);
		vlist[x + y * size + z * size * size].activity = 0;
	}*/

	//tree.dump();
	//tree.dump_level(height - 1);

	/*coord_t vx, vy, vz;
	tree.get_voxel_from_sum_activity(&vx, &vy, &vz, size * size * size / 2, vlist);
	std::cout << vx << ", " << vy << ", " << vz << std::endl;*/
}