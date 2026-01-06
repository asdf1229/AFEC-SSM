#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <set>
#include <map>
#include <unordered_set>
#include <sys/time.h>
#include "configuration/types.h"
#include "configuration/config.h"

#define pb push_back
#define mp make_pair
#define mmax(a,b) ((a)>(b)?(a):(b))
#define mmin(a,b) ((a)<(b)?(a):(b))

const int INF = 10000000;
const double EPS = 10e-6;
const ui block_size = 1024;

#include <cassert>

class Utility {
public:
	static FILE *open_file(const char *file_name, const char *mode) {
		FILE *f = fopen(file_name, mode);
		if(f == nullptr) {
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number) {
		std::vector<ui> sequence;
		if(number == 0) sequence.pb(0);
		while(number > 0) {
			sequence.pb(number%1000);
			number /= 1000;
		}

		char buf[5];
		std::string res;
		for(ui i = sequence.size();i > 0;i --) {
			if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
			else sprintf(buf, ",%03u", sequence[i-1]);
			res += std::string(buf);
		}
		return res;
	}
};

struct PairHash {
    size_t operator()(const std::pair<ui, ui>& p) const {
        return (size_t(p.first) << 32) ^ p.second;
    }
};

class Timer {
public:
	Timer() { m_start = timestamp(); }
	void restart() { m_start = timestamp(); }
	long long elapsed() { return timestamp() - m_start; }

private:
	long long m_start;

	// Returns a timestamp ('now') in microseconds
	long long timestamp() {
		struct timeval tp;
		gettimeofday(&tp, nullptr);
		return ((long long)(tp.tv_sec))*1000000 + tp.tv_usec;
	}
};

#endif
