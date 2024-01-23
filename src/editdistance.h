#ifndef EDITDISTANCE_H
#define EDITDISTANCE_H

#include <stdint.h>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <bitset>
#include <set>
using namespace std;

unsigned int edit_distance(const char *a, const unsigned int asize, const char *b, const unsigned int bsize);
// void create_patternmap(struct PatternMap *pm, const int64_t *a, const unsigned int size);
// unsigned int edit_distance_by_patternmap(struct PatternMap *mp, const int64_t *b, const unsigned int size);

unsigned int edit_distance(string a, string b);

std::pair<bool, std::set<int>> edit_distance2(string a, string b);

unsigned int hamming_distance(const char *a, const unsigned int asize, const char *b, const unsigned int bsize);

bool editdistance_test();

#endif /* EDITDISTANCE_H */

