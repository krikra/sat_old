#include <stdlib.h>
#include "sat_sys.h"
#include "sat_mem.h"

int main(int argc, char **argv)
{
	SAT *sat;

	sat = malloc(sizeof(SAT));

	sat_setdim(sat, 4);

	sat->nd_whole[0] = 10;
	sat->nd_whole[1] = 10;
	sat->nd_whole[2] = 10;
	sat->nd_whole[3] = 10;
	sat->dd_whole = 10 * 10 * 10 * 10;
	sat->dd_init = 16;

	sat_setvec(sat, 100);

	//sat_destroy(sat);
	return(0);
}
