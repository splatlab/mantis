#include "util.h"

std::string last_part(std::string str, char c) {
	uint64_t found = str.find_last_of(c);
	return str.substr(found + 1);
}

std::string first_part(std::string str, char c) {
	uint64_t found = str.find_first_of(c);
	return str.substr(0, found);
}

/* Print elapsed time using the start and end timeval */
void print_time_elapsed(std::string desc, struct timeval* start, struct
												timeval* end)
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	std::cout << desc << "Total Time Elapsed: " << std::to_string(time_elapsed) << 
		"seconds" << std::endl;
}
