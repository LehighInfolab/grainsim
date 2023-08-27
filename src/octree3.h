#pragma once

#include <cmath>
#include <iostream>

#include "types.h"
#include "voxel.h"

struct octree3_t
{
private:
	// The length of one side of the area that the octree represents (generally a power of 2).
	coord_t root_size;
	// The index of the lowest level in the tree.
	unsigned char max_level;
	// The total number of activities stored in the octree.
	size_t activity_count;
	// The activity array (stored in level-order).
	activ_t *activities;
	// A table that stores powers of 8 for efficient access.
	size_t *pow_table;

public:
	octree3_t(coord_t side_length, unsigned char height)
	{
		root_size = side_length;
		max_level = height - 1;

		reset_pos();

		pow_table = new size_t[height];
		
		activity_count = 0;
		for (unsigned char i = 0; i < height; ++i)
		{
			pow_table[i] = (size_t)pow(8, i);
			activity_count += pow_table[i];
		}
		activities = new activ_t[activity_count];
		for (size_t i = 0; i < activity_count; ++i) activities[i] = 0;
	}

	// Shift the activity of a certain voxel by the specified amount.
	void delta(coord_t x, coord_t y, coord_t z, activ_t dA)
	{
		if (dA == 0) return;

		reset_pos();

		while (true) // not big on this while(true)... maybe change function returns slightly
		{
			jump_to_positional_sibling(x, y, z);
			delta_current_node(dA);
			if (!first_child())
			{
				break;
			}
		}
	}

	// Returns xyz position of the voxel where the sum of all previous voxel activities (when walking the lattice) is equal to rand_activ.
	// TODO: Try to get rid of voxel_list dependancy.
	void get_voxel_from_sum_activity(coord_t *x, coord_t *y, coord_t *z, activ_t rand_activ, voxel_t *voxel_list, coord_t true_side_length)
	{
		// Reset positional pointer to root node.
		reset_pos();

		// Walk the tree.
		while (true)
		{
			while (current_node_activity() < rand_activ)
			{
				rand_activ -= current_node_activity();
				next_on_level();
			}
			// If first_child() returns false then the current node was a leaf; walk is done.
			if (!first_child())
			{
				break;
			}
		}

		// Iterate over all voxels contained within leaf node.
		// In some cases we may want leaf nodes to contain more than one voxel, which would make these loops necessary.
		size_t vindex = 0;
		for(*z = parent_z + offset_z; *z < std::min(parent_z + offset_z + node_size, true_side_length); ++*z)
			for (*y = parent_y + offset_y; *y < std::min(parent_y + offset_y + node_size, true_side_length); ++*y)
				for (*x = parent_x + offset_x; *x < std::min(parent_x + offset_x + node_size, true_side_length); ++*x)
				{
					vindex = *x + (*y * true_side_length) + (*z * true_side_length * true_side_length);
					if (voxel_list[vindex].activity >= rand_activ)
					{
						return;
					}
					rand_activ -= voxel_list[vindex].activity;
				}
	}

	// Print out all of the activities stored on a level.
	void dump_level(unsigned char level)
	{
		size_t rindex = 0;
		for (unsigned char l = 0; l < level; ++l)
		{
			rindex += pow_table[l];
		}
		activ_t total = 0;
		for (size_t i = 0; i < pow_table[level]; ++i)
		{
			std::cout << (rindex + i) << ": " << activities[rindex + i] << std::endl;
			total += activities[rindex + i];
		}
		std::cout << "TOTAL: " << total << std::endl;
	}

	// Get the overall system activity.
	activ_t system_activity()
	{
		return activities[0];
	}

private:

	// These helper functions are used to navigate through the tree.

	size_t curr_index;
	unsigned char curr_level, curr_sibling;
	coord_t parent_x, parent_y, parent_z;
	coord_t offset_x, offset_y, offset_z;
	coord_t node_size;

	// Reset all navigation parameters and return to the root node.
	void reset_pos()
	{
		parent_x = parent_y = parent_z = 0;
		offset_x = offset_y = offset_z = 0;
		curr_level = curr_sibling = curr_index = 0;
		node_size = root_size;
	}
	// Change the activity for the current node.
	void delta_current_node(activ_t dA)
	{
		activities[curr_index] += dA;
	}
	// Get the activity for the current node.
	activ_t current_node_activity()
	{
		return activities[curr_index];
	}
	// Move to the first child of the current node (returns false if the current node is a leaf).
	bool first_child()
	{
		if (curr_level == max_level) return false;
		if (curr_level == 0)
		{
			++curr_level;
			++curr_index;
			node_size /= 2;
			return true;
		}

		parent_x += offset_x;
		parent_y += offset_y;
		parent_z += offset_z;

		offset_x = offset_y = offset_z = 0;
		node_size /= 2;

		size_t rindex = 0;
		for (unsigned char level = 0; level < curr_level; ++level)
		{
			rindex += pow_table[level];
		}
		curr_index = (rindex + pow_table[curr_level]) + ((curr_index - rindex) * 8);

		curr_sibling = 0;
		++curr_level;

		return true;
	}
	// Moves to a sibling node under the same parent (sibling can be a value from 0 to 7).
	void jump_to_sibling(unsigned char sibling)
	{
		curr_index += (sibling - curr_sibling);
		curr_sibling = sibling;

		if (sibling >= 4)
		{
			offset_z = node_size;
			sibling -= 4;
		}
		else offset_z = 0;

		if (sibling >= 2)
		{
			offset_y = node_size;
			sibling -= 2;
		}
		else offset_y = 0;

		if (sibling >= 1)
		{
			offset_x = node_size;
			sibling -= 1;
		}
		else offset_x = 0;
	}
	// Moves to the next sibling with the same parent (returns false if current node is the 8th sibling).
	bool next_on_level()
	{
		if (curr_sibling >= 7) return false;

		jump_to_sibling(curr_sibling + 1);

		return true;
	}
	// Moves to the sibling that contains the specified point (this function assumes that the point lies within the parent node's region!)
	void jump_to_positional_sibling(coord_t x, coord_t y, coord_t z)
	{
		x -= parent_x;
		y -= parent_y;
		z -= parent_z;

		unsigned char sibling = 0;

		if (z >= node_size) sibling += 4;
		if (y >= node_size) sibling += 2;
		if (x >= node_size) sibling += 1;

		jump_to_sibling(sibling);
	}
};