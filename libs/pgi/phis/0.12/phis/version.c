#include "config.h"
#include "phis_guk.h"

void 
phis_version(int *major, int *minor, int *patch)
{
	*major = MAJOR;
	*minor = MINOR;
	*patch = PATCH;

	return;
}
