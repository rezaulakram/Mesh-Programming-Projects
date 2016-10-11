#ifndef __CORRESPONDANCE_HEADER__
#define __CORRESPONDANCE_HEADER__

#include "cholmod.h"
#include "umfpack.h"
#include "Matrix.hh"
#define LOADMODE 1
#define SAVEMODE 2
#define NORMALMODE 3

extern cholmod_common global;
extern int phase2mode;
extern int state;


typedef struct __CorrespondanceEntry
{
	unsigned long src_triangle_id;
	unsigned long target_triangle_id;
	double distance;
} CorrespondanceEntry;

void Initialize();
void Deinitialize();

void ComputeCorrespondance();

#endif