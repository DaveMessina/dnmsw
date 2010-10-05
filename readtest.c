#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "mpi.h"

#define MAXLINELEN 100002
#define THREADSPERNODE 2

/* remove newline from string */
void chomp (char * string) {
    char * ptr;
    if( (ptr = strchr(string, '\n')) != NULL)
        *ptr = '\0';
}

int main (int argc, char * argv[]) {

    // mpi variables
    int rank, omp_rank, mpisupport;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpisupport);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_procs   = 4;
    int num_threads = 2;
    int tot_threads = num_procs * num_threads;
    omp_set_num_threads(num_threads);

    // open file and initialize vars
    FILE *file;
    char lineAs[tot_threads][MAXLINELEN];
    char lineBs[tot_threads][MAXLINELEN];
    file = fopen(argv[1], "r");
    if (file == NULL) {
        perror ("Error reading file");
        exit(1);
    }

    int c;
    int num_lines;
    num_lines = tot_threads * 2;

    for (c = 0; c < tot_threads; c++) {
	// grab two lines and store them in our line arrays A & B
	char lineA [MAXLINELEN];
	char lineB [MAXLINELEN];
	fgets(lineA, MAXLINELEN, file);
	fgets(lineB, MAXLINELEN, file);
	chomp(lineA);
	chomp(lineB);
	strcpy(lineAs[c], lineA);
	strcpy(lineBs[c], lineB);                
    }
	
    #pragma omp parallel private(omp_rank)
    {
	omp_rank = omp_get_thread_num();
	char a [MAXLINELEN];
	char b [MAXLINELEN];
	int line_num = (rank * THREADSPERNODE)+omp_rank;
	strcpy(a, lineAs[line_num]);
	strcpy(b, lineBs[line_num]);
	printf("process %d, thread %d, seq %s and %s\n", rank, omp_rank, a, b);
    }

    MPI_Finalize();
}
