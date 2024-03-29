
#include <iostream>
#include <fstream>
#include <unordered_map>

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

#pragma pack(push, 1)
	struct boundary_info_t
	{
		int sm_to_lg_outies = 0,
			sm_to_lg_delta = 0,
			lg_to_sm_outies = 0,
			lg_to_sm_delta = 0,
			surface_area = 0;
	};
#pragma pack(pop)
	std::unordered_map<spin_t, std::unordered_map<spin_t, boundary_info_t> > sparse_info_matrix;
	std::unordered_map<spin_t, size_t> vol_map;

	void incr_sparse_outies(spin_t a, spin_t b)
	{
		if (a > b)
		{
			++sparse_info_matrix[b][a].lg_to_sm_outies;
		}
		else if (a < b)
		{
			++sparse_info_matrix[a][b].sm_to_lg_outies;
		}
	}
	void incr_sparse_delta(spin_t a, spin_t b)
	{
		if (a > b)
		{
			++sparse_info_matrix[b][a].lg_to_sm_delta;
		}
		else if (a < b)
		{
			++sparse_info_matrix[a][b].sm_to_lg_delta;
		}
	}
	void incr_sparse_sa(spin_t a, spin_t b)
	{
		++sparse_info_matrix[a > b ? b : a][a > b ? a : b].surface_area;
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
			incr_sparse_outies(id1, id2);
		}
		else if (id2 != id1 && id1 == id3 && id1 == id4)
		{
			incr_sparse_outies(id2, id1);
		}
		else if (id3 != id1 && id1 == id2 && id1 == id4)
		{
			incr_sparse_outies(id3, id1);
		}
		else if (id4 != id1 && id1 == id2 && id1 == id3)
		{
			incr_sparse_outies(id4, id1);
		}
	}

	size_t matrix_dim = 0;
	void generate_matrices()
	{
		sparse_info_matrix.clear();
		vol_map.clear();

		for(coord_t z = 0; z < curr_cube->side_length; ++z)
			for (coord_t y = 0; y < curr_cube->side_length; ++y)
				for (coord_t x = 0; x < curr_cube->side_length; ++x)
				{
					spin_t
						curr_id = curr_cube->voxel_at(x, y, z)->spin,
						fwd_id = curr_cube->voxel_at(x, y, z + 1)->spin,
						right_id = curr_cube->voxel_at(x + 1, y, z)->spin,
						up_id = curr_cube->voxel_at(x, y + 1, z)->spin;

					if (curr_id != right_id)
					{
						incr_sparse_sa(curr_id, right_id);
					}
					if (curr_id != fwd_id)
					{
						incr_sparse_sa(curr_id, fwd_id);
					}
					if (curr_id != up_id)
					{
						incr_sparse_sa(curr_id, up_id);
					}

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

					++vol_map[curr_id];
				}

	}

public:

	void load_lattice(lattice_t *cube)
	{
		curr_cube = cube;
		max_grains = calculate_max_grains();
		matrix_dim = max_grains + 1;

		generate_matrices();
	}

	double get_curvature(spin_t a, spin_t b) // DOES NOT VERIFY THAT BOUNDARY EXISTS!!!!!
	{
		if (a > b)
		{
			return (3.141592653589793 / 4.0) * (sparse_info_matrix[b][a].lg_to_sm_outies - sparse_info_matrix[b][a].sm_to_lg_outies);
		}
		else if (a < b)
		{
			return (3.141592653589793 / 4.0) * (sparse_info_matrix[a][b].sm_to_lg_outies - sparse_info_matrix[a][b].lg_to_sm_outies);
		}
		return 0;
	}

	void save_analysis_to_file(const char *fname)
	{
		std::ofstream afile(fname);

		std::cout << "Creating analysis file " << fname << std::endl;

		// Volumes
		afile << "VOLUMES\n";
		for (auto vol_iter = vol_map.begin(); vol_iter != vol_map.end(); ++vol_iter)
		{
			afile << vol_iter->first << ' ' << vol_iter->second << '\n';
		}

		// Curvatures
		afile << "CURVATURES\n";
		for (auto sm_iter = sparse_info_matrix.begin(); sm_iter != sparse_info_matrix.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				if (lg_iter->second.surface_area == 0) continue;

				afile << sm_iter->first << ' ' << lg_iter->first << ' ' << get_curvature(sm_iter->first, lg_iter->first) << '\n';
				afile << lg_iter->first << ' ' << sm_iter->first << ' ' << get_curvature(lg_iter->first, sm_iter->first) << '\n';
			}
		}

		// Curvatures
		afile << "SURFACE_AREAS\n";
		for (auto sm_iter = sparse_info_matrix.begin(); sm_iter != sparse_info_matrix.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				if (lg_iter->second.surface_area == 0) continue;

				afile << sm_iter->first << ' ' << lg_iter->first << ' ' << lg_iter->second.surface_area << '\n';
				afile << lg_iter->first << ' ' << sm_iter->first << ' ' << lg_iter->second.surface_area << '\n';
			}
		}

		// Velocities
		afile << "VELOCITIES\n";
		for (auto sm_iter = curr_cube->boundary_tracker.velocity_tracker.begin(); sm_iter != curr_cube->boundary_tracker.velocity_tracker.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				std::pair<int, int> delta = lg_iter->second;
				
				afile << sm_iter->first << ' ' << lg_iter->first << ' ' << (delta.first - delta.second) << '\n';
				afile << lg_iter->first << ' ' << sm_iter->first << ' ' << (delta.second - delta.first) << '\n';
			}
		}

		afile << "ADJACENT_BOUNDARIES\n";
		for (auto sm_iter = curr_cube->boundary_tracker.boundary_map.begin(); sm_iter != curr_cube->boundary_tracker.boundary_map.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				boundary_t *boundary = lg_iter->second;

				if (boundary->area() == 0) continue;

				afile << boundary->a_spin << '/' << boundary->b_spin;
				for (auto junc_iter = boundary->junctions.begin(); junc_iter != boundary->junctions.end(); ++junc_iter)
				{
					afile << ' ' << junc_iter->first->a_spin << '/' << junc_iter->first->b_spin;
				}

				afile << '\n';
			}
		}

		curr_cube->boundary_tracker.reset_flip_tracker();

		afile.close();
	}
};