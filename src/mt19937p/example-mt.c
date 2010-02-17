/**
 * Parallel MT19937 example.
 */

#include <time.h>
#include <stdio.h>
#include "mt19937p.h"

int main(int argc, char** argv) {
	int i;
	struct mt19937p mt;
	sgenrand(time(NULL), &mt);
	for(i=0;i<100;i++) {
		printf("%f\n",genrand(&mt));
	}
	return 0;
}
