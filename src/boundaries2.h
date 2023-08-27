#pragma once

#include <unordered_map>
#include <unordered_set>

#include "types.h"
#include "config.h"

#pragma pack(push, 1)
struct boundary_t
{
	spin_t a_spin, b_spin; // storing these is actually not necessary, but is used for logging

	bool transformed = false;
	std::unordered_set<size_t> boundary_voxel_indices;

	size_t previous_surface_area = 0;
	int potential_energy = 0;

	// Stores the adjacent boundaries and the number of junctions that they share.
	std::unordered_map<boundary_t *, size_t> junctions;

	// I may be wrong here but the junction map may cause a memory leak...
	// Is it possible that a boundary may be deleted before its removal from the junction list? Not sure...
	// Even if this is the case, I can't imagine that it causes that much of a problem, as these boundaries are generally short-lived.
	// However, on a huge cube with an abnormally long-lasting boundary I'm not so sure...

	void delta_junction(boundary_t *b, char dArea)
	{
		junctions[b] += dArea;
		if (junctions[b] == 0)
		{
			delete_junction(b);
		}
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

	std::unordered_map<spin_t, std::unordered_map<spin_t, boundary_t *>> boundary_map;
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

		/*for (auto junc_iter = boundary->junctions.begin(); junc_iter != boundary->junctions.end(); ++junc_iter)
		{
			junc_iter->first->delete_junction(boundary);
		}*/

		delete boundary;
	}

	// Check if the boundary between two grains is transformed.
	bool is_transformed(spin_t a, spin_t b)
	{
		return find_or_create_boundary(a, b)->transformed;
	}

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
	void remove_from_boundary(spin_t a, spin_t b, size_t index, spin_t *voxel_neighbor_spins = nullptr)
	{
		boundary_t *boundary = find_or_create_boundary(a, b);
		boundary->boundary_voxel_indices.erase(index);
		if (boundary->boundary_voxel_indices.size() == 0)
		{
			delete_boundary(a, b);
		}
		else if (voxel_neighbor_spins != nullptr)
		{
			for (char i = 0; i < NEIGH_COUNT; ++i)
			{
				if (voxel_neighbor_spins[i] != 0 && voxel_neighbor_spins[i] != a && voxel_neighbor_spins[i] != b)
				{
					boundary->delta_junction(find_or_create_boundary(a, voxel_neighbor_spins[i]), -1); // assume that spin "a" is the root voxel
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
};