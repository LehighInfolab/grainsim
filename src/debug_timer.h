#pragma once

#include <chrono>

class debug_timer_t
{
private:

	std::chrono::time_point<std::chrono::system_clock> tstart, tlap, temp;
	std::chrono::duration<double> elapsed;

public:

	void start()
	{
		tstart = tlap = std::chrono::system_clock::now();
	}

	double lap()
	{
		temp = std::chrono::system_clock::now();
		elapsed = temp - tlap;
		tlap = temp;
		return elapsed.count();
	}

	double total()
	{
		temp = std::chrono::system_clock::now();
		elapsed = temp - tstart;
		return elapsed.count();
	}
};