#ifndef __H_OMPTHREADER
#define __H_OMPTHREADER

#include <omp.h>

/*
 * A header with functions for setting up OpenMP's number of threads
 * appropriately for the current system (i.e., at runtime). Call
 * setupOpenMP() elsewhere to change the settings and print a message
 * about it.
 */

int getOurThreads(int requestedThreads);
void setupOpenMP(int numThreads);

// run instructions to setup number of OpenMP threads appropriately
void setupOpenMP(int numThreads = -1) {
#ifdef crism_omp_on
    // a check on OpenMP status borrowed from Soysal's code
    int nProcessors=omp_get_num_procs();
    int nOurThreads = getOurThreads(numThreads);
    omp_set_num_threads(nOurThreads);
    printf("Using %d out of %d threads\n",nOurThreads,nProcessors);
#else
    printf("Running single-threaded.\n");
#endif
}

// chooses a number of threads for us to use. Attempt to provide requestedThreads threads,
// but use default method if its value is -1.
int getOurThreads(int requestedThreads) {
    int nProcessors=omp_get_num_procs();
    
	if (requestedThreads > 0) {
		return requestedThreads > nProcessors ? nProcessors : requestedThreads; // if more than we have requested, use them all
	}
	else {
		if (nProcessors <= 3) {
			return nProcessors; // if there are few enough, use all of them
		}
		else {
			return nProcessors - 2; // if there are many, leave 2 in background
		}
	}
}

#endif /* __H_OMPTHREADER */