#pragma once

#include "types.h"
#include "boundaries2.h"

// An object that represents a single voxel within the lattice.
#pragma pack(push, 1)
struct voxel_t
{
private:
	// A const that signifies that no neighbor is present within that slot.
	const spin_t NO_NEIGHBOR = 0;
	// A list of the neighboring grains that this voxel is touching (note that this list is UNIQUE, so only one entry for a grain can exist at once, and the voxel's current grain is not present).
	spin_t neighbor_spins[NEIGH_COUNT];
	// A list of the probabilities that this voxel has for flipping to each neighboring grain.
	activ_t neighbor_probs[NEIGH_COUNT];

public:
	// The spin (grain ID) of this voxel.
	spin_t spin;
	// The activity of this voxel.
	activ_t activity;
	// The index of this voxel within the lattice.
	size_t index;

	voxel_t()
	{
		for (char i = 0; i < NEIGH_COUNT; ++i)
		{
			neighbor_spins[i] = neighbor_probs[i] = 0;
		}
		spin = activity = 0;
	}

	// Set the probability that this voxel will flip to a certain grain (returns the resulting change in voxel activity).
	activ_t set_neighbor(spin_t nspin, activ_t prob, boundary_tracker_t *blist)
	{
		if (prob == 0)
		{
			return remove_neighbor(nspin, blist);
		}

		char nindex = -1;
		bool new_neighbor = true;
		for (char i = 0; i < NEIGH_COUNT; ++i)
		{
			if (neighbor_spins[i] == nspin)
			{
				nindex = i;
				new_neighbor = false;
				break;
			}
			else if (nindex < 0 && neighbor_spins[i] == NO_NEIGHBOR)
			{
				nindex = i;
			}
		}

		if (nindex < 0)
		{
			std::cout << "Error: Voxel-wise adjacent grain list overflow." << std::endl;
			exit(0);
		}

		if (new_neighbor)
		{
			neighbor_spins[nindex] = nspin;
			neighbor_probs[nindex] = prob;
			activity += prob;
			if(blist != nullptr) blist->add_to_boundary(spin, nspin, index, neighbor_spins, spin);
			return prob;
		}
		else
		{
			prob -= neighbor_probs[nindex];
			neighbor_probs[nindex] += prob;
			activity += prob;
			return prob;
		}
	}

	bool has_neighbor(spin_t nspin)
	{
		for (char i = 0; i < NEIGH_COUNT; ++i)
		{
			if (neighbor_spins[i] == nspin)
			{
				return true;
			}
		}
		return false;
	}

	// Remove a certain grain from the neighbor list (returns the resulting change in voxel activity).
	activ_t remove_neighbor(spin_t nspin, boundary_tracker_t *blist)
	{
		for (char i = 0; i < NEIGH_COUNT; ++i)
		{
			if (neighbor_spins[i] == nspin)
			{
				neighbor_spins[i] = NO_NEIGHBOR;
				activity -= neighbor_probs[i];
				if (blist != nullptr) blist->remove_from_boundary(spin, nspin, index, neighbor_spins, spin);
				return -neighbor_probs[i];
			}
		}
		return 0;
	}

	// Remove all neighbors from the list (returns the resulting change in voxel activity).
	activ_t reset(boundary_tracker_t *blist)
	{
		activ_t delta = 0;
		for (char i = 0; i < NEIGH_COUNT; ++i)
		{
			if (neighbor_spins[i] != NO_NEIGHBOR)
			{
				delta -= neighbor_probs[i];
				if (blist != nullptr) blist->remove_from_boundary(spin, neighbor_spins[i], index, neighbor_spins, spin);
			}
			neighbor_spins[i] = NO_NEIGHBOR;
		}
		activity = 0;
		return delta;
	}

	// Choose a neighbor based on a random desired activity value (0..voxel_activity).
	spin_t choose_neighbor(activ_t desired_activ)
	{
		for (char i = 0; i < NEIGH_COUNT; ++i)
		{
			if (neighbor_spins[i] == NO_NEIGHBOR) continue;

			desired_activ -= neighbor_probs[i];
			if (desired_activ <= 0) return neighbor_spins[i];
		}
		return NO_NEIGHBOR;
	}
};
#pragma pack(pop)