
#include <iostream>
#include <fstream>

#include "types.h"
#include "lattice.h"

class lattice_analyzer_t
{
private:

	lattice_t *curr_cube;

	spin_t max_grains;
	spin_t calculate_max_grains()
	{
		spin_t max_so_far = 0;
		for (coord_t z = 0; z < curr_cube->side_length; ++z)
			for (coord_t y = 0; y < curr_cube->side_length; ++y)
				for (coord_t x = 0; x < curr_cube->side_length; ++x)
				{
					spin_t curr_spin = curr_cube->voxel_at(x, y, z)->spin;
					if (curr_spin > max_so_far) max_so_far = curr_spin;
				}
		return max_so_far;
	}

	void memfill(int *arr, size_t size, int val)
	{
		for (size_t i = 0; i < size; ++i)
		{
			arr[i] = val;
		}
	}

	void check_edge(
		coord_t rx, coord_t ry, coord_t rz,
		coord_t x1, coord_t y1, coord_t z1,
		coord_t x2, coord_t y2, coord_t z2,
		coord_t x3, coord_t y3, coord_t z3,
		coord_t x4, coord_t y4, coord_t z4)
	{
		spin_t id1, id2, id3, id4;

		id1 = curr_cube->voxel_at(rx + x1, ry + y1, rz + z1)->spin;
		id2 = curr_cube->voxel_at(rx + x2, ry + y2, rz + z2)->spin;
		id3 = curr_cube->voxel_at(rx + x3, ry + y3, rz + z3)->spin;
		id4 = curr_cube->voxel_at(rx + x4, ry + y4, rz + z4)->spin;

		if (id1 != id2 && id2 == id3 && id2 == id4)
		{
			++outie_matrix[id1 * matrix_dim + id2];
		}
		else if (id2 != id1 && id1 == id3 && id1 == id4)
		{
			++outie_matrix[id2 * matrix_dim + id1];
		}
		else if (id3 != id1 && id1 == id2 && id1 == id4)
		{
			++outie_matrix[id3 * matrix_dim + id1];
		}
		else if (id4 != id1 && id1 == id2 && id1 == id3)
		{
			++outie_matrix[id4 * matrix_dim + id1];
		}
	}

	int *outie_matrix = nullptr, *vol_vector = nullptr;
	size_t matrix_dim = 0;
	void generate_matrices(bool realloc=true)
	{
		if (realloc)
		{
			if (outie_matrix != nullptr && vol_vector != nullptr)
			{
				delete outie_matrix;
				delete vol_vector;
			}
			outie_matrix = new int[matrix_dim * matrix_dim];
			vol_vector = new int[matrix_dim];
		}

		memfill(outie_matrix, matrix_dim * matrix_dim, 0);
		memfill(vol_vector, matrix_dim, 0);

		for(coord_t z = 0; z < curr_cube->side_length; ++z)
			for (coord_t y = 0; y < curr_cube->side_length; ++y)
				for (coord_t x = 0; x < curr_cube->side_length; ++x)
				{
					// back bottom
					check_edge(
						x, y, z,
						0, 0, -1, 
						0, 0, 0,
						0, -1, 0,
						0, -1, -1);
					// back left
					check_edge(
						x, y, z,
						-1, 0, 0,
						0, 0, 0, 
						0, 0, -1, 
						-1, 0, -1);
					// top left
					check_edge(
						x, y, z,
						-1, 1, 0, 
						0, 1, 0,
						0, 0, 0,
						-1, 0, 0);

					++vol_vector[curr_cube->voxel_at(x, y, z)->spin];
				}
	}

public:

	void load_lattice(lattice_t *cube)
	{
		curr_cube = cube;
		max_grains = calculate_max_grains();

		bool realloc = matrix_dim != max_grains + 1;
		matrix_dim = max_grains + 1;

		generate_matrices(realloc);
	}

	double get_curvature(spin_t a, spin_t b) // DOES NOT CHECK IF BOUNDARY EXISTS!!!!!
	{
		return (3.141592653589793 / 4.0) * (outie_matrix[a * matrix_dim + b] - outie_matrix[b * matrix_dim + a]);
	}

	void save_analysis_to_file(const char *fname)
	{
		std::ofstream afile(fname);

		std::cout << "Creating analysis file " << fname << std::endl;

		// Volumes
		afile << "VOLUMES\n";
		for (size_t i = 0; i < matrix_dim; ++i)
		{
			if (vol_vector[i] == 0) continue;
			afile << i << ' ' << vol_vector[i] << '\n';
		}

		// Curvatures
		afile << "CURVATURES\n";
		for (auto sm_iter = curr_cube->boundary_tracker.boundary_map.begin(); sm_iter != curr_cube->boundary_tracker.boundary_map.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				boundary_t *boundary = lg_iter->second;

				if (boundary->area() == 0) continue;

				afile << boundary->a_spin << ' ' << boundary->b_spin << ' ' << get_curvature(boundary->a_spin, boundary->b_spin) << '\n';
				afile << boundary->b_spin << ' ' << boundary->a_spin << ' ' << get_curvature(boundary->b_spin, boundary->a_spin) << '\n';
			}
		}

		// Curvatures
		afile << "SURFACE_AREAS\n";
		for (auto sm_iter = curr_cube->boundary_tracker.boundary_map.begin(); sm_iter != curr_cube->boundary_tracker.boundary_map.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				boundary_t *boundary = lg_iter->second;

				if (boundary->area() == 0) continue;

				afile << boundary->a_spin << ' ' << boundary->b_spin << ' ' << boundary->area() << '\n';
				afile << boundary->b_spin << ' ' << boundary->a_spin << ' ' << boundary->area() << '\n';
			}
		}

		afile.close();
	}
};