#pragma once

#include <vector>
#include <unistd.h>
#include <sstream>
#include <fstream>

#include "types.h"

struct config_t
{
	std::string initial_state_path;
	std::string identifier;
	std::string output_folder;
	std::string checkpoints;
	double checkpoint_interval = -1;
	double max_timestep = -1;
	activ_t default_mobility = 0.002, transitioned_mobility = 0.04;
	double transition_interval = 0;
	size_t transition_count = 0;
	double scale_multiplier = 1;
	double propagation_chance = 0.95;
	bool use_potential_energy = false;
	int const_grain_count = 0;
	bool log_transitions = false;
	double propagation_ratio = 0;
	bool generate_analysis_files = false;

	void checkpoints_to_vector(std::vector<double> *checkpoint_vector)
	{
		std::istringstream ss(checkpoints);
		std::string word;
		while (ss >> word)
		{
			checkpoint_vector->push_back(std::stoul(word));
		}
	}

	void load_config()
	{
		std::ifstream cfgfile("grainsim_config.txt");
		std::string line;
		while (std::getline(cfgfile, line))
		{
			std::istringstream ss(line);
			std::string key, equals, value, tempv;

			if(!(ss >> key)) continue;
			if (key[0] == '#') continue;

			if (!(ss >> equals)) continue;
			if (!(ss >> value)) continue;
			while (ss >> tempv)
			{
				value.append(" ");
				value.append(tempv);
			}

			if (key == "INITIAL_STATE_FILE")
			{
				initial_state_path = value;
			}
			else if (key == "OUTPUT_FOLDER")
			{
				output_folder = value;
			}
			else if (key == "IDENTIFIER")
			{
				identifier = value;
			}
			else if (key == "CHECKPOINTS")
			{
				checkpoints = value;
			}
			else if (key == "PERIODIC_CHECKPOINT_INTERVAL")
			{
				checkpoint_interval = std::stod(value);
			}
			else if (key == "MAX_TIMESTEP")
			{
				max_timestep = std::stod(value);
			}
			else if (key == "DEFAULT_MOBILITY")
			{
				default_mobility = std::stod(value);
			}
			else if (key == "TRANSITIONED_MOBILITY")
			{
				transitioned_mobility = std::stod(value);
			}
			else if (key == "TRANSITION_INTERVAL")
			{
				transition_interval = std::stod(value);
			}
			else if (key == "TRANSITION_COUNT")
			{
				transition_count = std::stoi(value);
			}
			else if (key == "PROPAGATION_CHANCE")
			{
				propagation_chance = std::stod(value);
			}
			else if (key == "USE_POTENTIAL_ENERGY")
			{
				use_potential_energy = value == "true";
			}
			else if (key == "SCALE_MULTIPLIER")
			{
				scale_multiplier = std::stod(value);
			}
			else if (key == "LOG_BOUNDARY_TRANSITIONS")
			{
				log_transitions = value == "true";
			}
			else if (key == "CONST_GRAIN_COUNT")
			{
				const_grain_count = std::stoi(value);
			}
			else if (key == "PROPAGATION_RATIO")
			{
				propagation_ratio = std::stod(value);
			}
			else if (key == "GENERATE_ANALYSIS_FILES")
			{
				generate_analysis_files = value == "true";
			}
			else
			{
				std::cout << "Warning: Unknown config key \"" << key << "\"." << std::endl;
			}
		}
	}
};