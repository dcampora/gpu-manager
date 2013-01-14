//pragma once
#ifndef SEARCHKERNEL_H_
#define SEARCHKERNEL_H_
#include <vector>

void launchKernel(std::vector<char> &inputVector, int eventsNumber, int *events, unsigned int allHits, char **results);
inline int getBestDevice();
void findCudaDevice();

#endif /* SEARCHKERNEL_H_ */
