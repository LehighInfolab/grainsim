#pragma once

#include <unordered_map>
#include <unordered_set>
#include <list>

#include "types.h"
#include "config.h"

#pragma pack(push, 1)
struct boundary_t
{
	spin_t a_spin, b_spin; // storing these is actually not necessary, but is used for logging

	bool marked_for_deletion = false;

	bool transformed = false;
	std::unordered_set<size_t> boundary_voxel_indices;

	size_t previous_surface_area = 0;
	int potential_energy = 0;

	// Stores the adjacent boundaries and the number of voxels that they share.
	std::unordered_map<boundary_t *, long> junctions;

	void delta_junction(boundary_t *b, char dArea)
	{
		junctions[b] += dArea;
	}
	void delete_junction(boundary_t *b)
	{
		if(junctions.find(b) != junctions.end())
			junctions.erase(b);
	}

	size_t area()
	{
		return boundary_voxel_indices.size();
	}
};
#pragma pack(pop)

struct boundary_tracker_t
{
	// The boundary map is actually a map of maps. When trying to find the boundary object for
	// the boundary between two grains, the root map takes the smaller spin of the two as the key,
	// and the child map takes the larger spin of the two as the key. (e.g. if trying to find the boundary
	// between grains 5 and 10, we would do it via "boundary_map.at(5).at(10)" ). This ordering prevents
	// redundancies.

	std::unordered_map<spin_t, std::unordered_map<spin_t, boundary_t *> > boundary_map;
	size_t transformed_boundary_count = 0, total_boundary_count = 0;

	// Find the boundary between two grains, or create it if it does not yet exist.
	boundary_t *find_or_create_boundary(spin_t a, spin_t b)
	{
		boundary_t *output = boundary_map[a < b ? a : b][a < b ? b : a];
		if (!output)
		{
			output = new boundary_t();
			output->a_spin = a;
			output->b_spin = b;
			boundary_map[a < b ? a : b][a < b ? b : a] = output;
			++total_boundary_count;
		}
		else if (output->marked_for_deletion)
		{
			unmark_for_deletion(output);
		}
		return output;
	}

	// Forcefully delete the boundary between two grains.
	void delete_boundary(spin_t a, spin_t b)
	{
		std::unordered_map<spin_t, boundary_t *> *sm_bucket = &boundary_map.at(a < b ? a : b);

		boundary_t *boundary = sm_bucket->at(a < b ? b : a);
		sm_bucket->erase(a < b ? b : a);

		if (boundary->transformed)
		{
			--transformed_boundary_count;
		}
		if (sm_bucket->size() == 0)
		{
			boundary_map.erase(a < b ? a : b);
		}

		--total_boundary_count;

		// just give potential energy to a random boundary...
		if (boundary->junctions.size() > 0)
		{
			boundary_t *transfer_boundary = nullptr;
			for (auto junc_iter = boundary->junctions.begin(); junc_iter != boundary->junctions.end(); ++junc_iter)
			{
				if (junc_iter->first->transformed)
				{
					if (junc_iter->first->potential_energy > 0)
					{
						transfer_boundary = junc_iter->first;
						break;
					}
					else if (transfer_boundary == nullptr)
					{
						transfer_boundary = junc_iter->first;
					}
				}
			}
			if (transfer_boundary == nullptr) transfer_boundary = boundary->junctions.begin()->first;
			transfer_boundary->potential_energy += boundary->potential_energy;
		}

		delete boundary;
	}

	// Mark a boundary for deletion.
	void mark_for_deletion(boundary_t *boundary)
	{
		boundary->marked_for_deletion = true;
	}
	// Remove the mark for deletion from a boundary.
	void unmark_for_deletion(boundary_t *boundary)
	{
		boundary->marked_for_deletion = false;
	}

	// Check if the boundary between two grains is transformed.
	bool is_transformed(spin_t a, spin_t b)
	{
		return find_or_create_boundary(a, b)->transformed;
	}

	// Add a voxel to a boundary and update that boundary's junctions.
	void add_to_boundary(spin_t a, spin_t b, size_t index, spin_t *voxel_neighbor_spins = nullptr)
	{
		boundary_t *boundary = find_or_create_boundary(a, b);
		boundary->boundary_voxel_indices.insert(index);
		if (voxel_neighbor_spins != nullptr)
		{
			for (char i = 0; i < NEIGH_COUNT; ++i)
			{
				if (voxel_neighbor_spins[i] != 0 && voxel_neighbor_spins[i] != a && voxel_neighbor_spins[i] != b)
				{
					boundary->delta_junction(find_or_create_boundary(a, voxel_neighbor_spins[i]), 1); // assume that spin "a" is the root voxel
				}
			}
		}
	}
	// Remove a voxel from a boundary and update that boundary's junctions.
	void remove_from_boundary(spin_t a, spin_t b, size_t index, spin_t *voxel_neighbor_spins = nullptr)
	{
		boundary_t *boundary = find_or_create_boundary(a, b);
		boundary->boundary_voxel_indices.erase(index);

		if (boundary->area() == 0)
		{
			mark_for_deletion(boundary);
		}
		else if (voxel_neighbor_spins != nullptr)
		{
			for (char i = 0; i < NEIGH_COUNT; ++i)
			{
				if (voxel_neighbor_spins[i] != 0 && voxel_neighbor_spins[i] != a && voxel_neighbor_spins[i] != b)
				{
					boundary->delta_junction(find_or_create_boundary(a, voxel_neighbor_spins[i]), -1); // assume that spin "a" is the root grain
				}
			}
		}
	}

	// Mark the boundary between two grains as transformed.
	void mark_transformed(spin_t a, spin_t b)
	{
		boundary_t *boundary = find_or_create_boundary(a, b);
		if (boundary->transformed) return;

		boundary->transformed = true;
		++transformed_boundary_count;
	}
	// Mark a boundary that you already have the reference for as transformed.
	void mark_transformed(boundary_t *boundary)
	{
		if (boundary->transformed) return;

		boundary->transformed = true;
		++transformed_boundary_count;
	}

	// Delete all marked boundaries from the boundary map and remove all invalid junctions.
	void remove_marked_boundaries()
	{
		std::list<boundary_t *> delete_list;

		for (auto sm_iter = boundary_map.begin(); sm_iter != boundary_map.end(); ++sm_iter)
		{
			for (auto lg_iter = sm_iter->second.begin(); lg_iter != sm_iter->second.end(); ++lg_iter)
			{
				boundary_t *boundary = lg_iter->second;
				if (boundary->marked_for_deletion || boundary->area() == 0)
				{
					delete_list.push_back(boundary);
					continue;
				}

				std::list<boundary_t *> remove_from_junctions_list;
				for (auto junc_iter = boundary->junctions.begin(); junc_iter != boundary->junctions.end(); ++junc_iter)
				{
					boundary_t *jbound = junc_iter->first;
					if (jbound->marked_for_deletion || junc_iter->second <= 0)
					{
						remove_from_junctions_list.push_back(jbound);
					}
				}

				for (auto rm_junc_iter = remove_from_junctions_list.begin(); rm_junc_iter != remove_from_junctions_list.end(); ++rm_junc_iter)
				{
					boundary->junctions.erase(*rm_junc_iter);
				}
			}
		}

		for (auto delete_iter = delete_list.begin(); delete_iter != delete_list.end(); ++delete_iter)
		{
			delete_boundary((*delete_iter)->a_spin, (*delete_iter)->b_spin);
		}
	}
};