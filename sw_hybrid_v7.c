#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "mpi.h"

#define MAXLINELEN 100002
#define ALIGNMENTS 32

// define scoring scheme
#define MATCH     1
#define MISMATCH -1
#define GAP      -1

double gettime(void) {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

void show_args(int argc, const char * argv[]) {
    int count;
    printf( "\nCommand-line arguments:\n" );
    for( count = 0; count < argc; count++ ) {
        printf( "  argv[%d]   %s\n", count, argv[count] );
    }
}

void fill_matrix(int * matrix, int x, int y, int value) {
    int i, j;
    
    for (i = 0; i <= x; i++) {
        for (j = 0; j <= y; j++) {
            //printf("before %d*%d + %d = matrix[%d]\t%d\n", x, i, j, x*i+j, matrix[x*i+j]);
            matrix[x*i + j] = value;
            //printf(" after %d*%d + %d = matrix[%d]\t%d\n", x, i, j, x*i+j, matrix[x*i+j]);
        }
    }
}

void print_matrix(int * matrix, int x, int y) {
    int i, j;
    
    for (i = 0; i <= x; i++) {
        for (j = 0; j <= y; j++) {
            int value = matrix[x*i + j];
            printf("\t(%d,%d)=%d", i, j, value);
        }
        printf("\n");
    }
    printf("\n");
}

void walk_matrix(int * score_matrix, int * ptr_matrix, char * seq1, char * seq2, int * max_i, int * max_j, int * max_score) {
    int i, j;
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);
    
    // fill every cell in the matrix
    for (i = 1; i <= seq1len; i++) {
        for (j = 1; j <= seq2len; j++) {
            int diagonal_score, left_score, up_score;
            
            // calculate match score
            char a = seq1[i-1];
            char b = seq2[j-1];
            if (a == b) {
                diagonal_score = score_matrix[seq1len*(i-1) + (j-1)] + MATCH;
            }
            else {
                diagonal_score = score_matrix[seq1len*(i-1) + (j-1)] + MISMATCH;
            }
            
            // calculate gap scores
            up_score   = score_matrix[seq1len*(i-1) +  j   ] + GAP;
            left_score = score_matrix[seq1len*i     + (j-1)] + GAP;
            
            if (diagonal_score <= 0 && up_score <= 0 && left_score <= 0) {
                score_matrix[seq1len*i + j] = 0;
                ptr_matrix[seq1len*i + j]   = 0;
                //printf("(%c %d, %c %d) = %f **CONT\n", a, i, b, j, score_matrix[i*seq1len +j]);
                continue;
            }
            
            // choose best score
            // ptr_matrix values:
            // 0 : no pointer
            // 1 : diagonal
            // 2 : left
            // 3 : up
            // (these could be replaced with a 2 bit value instead)
            if (diagonal_score >= up_score) {
                if (diagonal_score >= left_score) {
                    score_matrix[seq1len*i + j] = diagonal_score;
                    ptr_matrix[seq1len*i + j]   = 1;
                }
                else {
                    score_matrix[seq1len*i + j] = left_score;
                    ptr_matrix[seq1len*i + j]   = 2;
                }
            }
            else {
                if (up_score >= left_score) {
                    score_matrix[seq1len*i + j] = up_score;
                    ptr_matrix[seq1len*i + j]   = 3;
                    
                }
                else {
                    score_matrix[seq1len*i + j] = left_score;
                    ptr_matrix[seq1len*i + j]   = 2;
                    
                }
                
            }
            
            // set maximum score
            if (score_matrix[seq1len*i + j] > *max_score) {
                *max_i     = i;
                *max_j     = j;
                *max_score = score_matrix[seq1len*i + j];
                // printf("** now max_score is %d\n", max_score);
            }
            
            //printf("(%c %d, %c %d) = %f\n", a, i, b, j, score_matrix[i*seq1len +j]);
        }
    }
    
}

/* append a char to a string */
void append(char * string, char c)
{
    int len = strlen(string);
    string[len] = c;
    string[len+1] = '\0';
}

/* remove newline from string */
void chomp (char * string) {
    char * ptr;
    if( (ptr = strchr(string, '\n')) != NULL)
        *ptr = '\0';
}

/* reverse: reverse string s in place */
void reverse(char s[]) {
    int c, i, j;
    for (i = 0, j = strlen(s)-1; i < j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

void traceback(int * score_matrix, int * ptr_matrix, char * seq1, char * seq2, int max_i, int max_j) {
    int seq1len = strlen(seq1);
    int seq2len = strlen(seq2);
    int max_align_len = seq1len + seq2len;
    int flag = 1;   // used to exit the while loop
    int ii = max_i; // ii and jj are our position in the matrix as we traceback
    int jj = max_j; // starting from the i,j where the greatest score was obtained

    // initialize output strings
    char align1[max_align_len+1];
    strcpy(align1, "");
    char align2[max_align_len+1];
    strcpy(align2, "");
    
    while (flag) {
        int tb = ptr_matrix[ii*seq1len + jj];
        char a, b;
        
        switch (tb) {
            // if we reach a 0 traceback ptr, we're done
            case 0:
                flag = 0;
                break;
            // diagonal
            case 1:
                a = seq1[ii-1];
                b = seq2[jj-1];
                append(align1, a);
                append(align2, b);
                ii--;
                jj--;
                break;
            // left
            case 2:
                a = seq1[ii-1];
                b = '-';
                append(align1, a);
                append(align2, b);
                jj--;
                break;
            // up
            case 3:
                a = '-';
                b = seq2[jj-1];
                append(align1, a);
                append(align2, b);
                ii--;
                break;
            default:
                perror("invalid value in ptr_matrix!\n");
                break;
            break;
        }
    }

    // reverse the strings (traceback starts at the end)
    reverse(align1);
    reverse(align2);
    printf("%s\n%s\n", align1, align2);
}

void do_alignment(char * line1, char * line2) {

    // get start time
    double align_t = gettime();

    // prep seqs
    int seq1len = strlen(line1);
    int seq2len = strlen(line2);
    char seq1[seq1len];
    strcpy(seq1, line1);
    char seq2[seq2len];
    strcpy(seq2, line2);

    // initialize score and pointer matrices
    int * score_matrix = (int *) malloc((10000+1) * (10000+1) * sizeof(int));
    int * ptr_matrix   = (int *) malloc((10000+1) * (10000+1) * sizeof(int));

    double fill1_t = gettime();
    fill_matrix(score_matrix, seq1len, seq2len, 0);
    fill1_t = gettime() - fill1_t;

    double fill2_t = gettime();
    fill_matrix(  ptr_matrix, seq1len, seq2len, 0); 
    fill2_t = gettime() - fill2_t;

    //print_matrix(ptr_matrix, seq1len, seq2len);
    //print_matrix(score_matrix, seq1len, seq2len);

    // do the score calculation for each cell in the matrix
    int max_i = 0, max_j = 0, max_score = 0;
    walk_matrix(score_matrix, ptr_matrix, seq1, seq2, &max_i, &max_j, &max_score);

    // find the alignment by following the pointers back through the matrix
    traceback(score_matrix, ptr_matrix, seq1, seq2, max_i, max_j);

    // release matrix memory
    free(score_matrix);
    free(ptr_matrix);

    // print score
    printf("max score: %d\n", max_score);

    // print num of cells, time in secs, and flops
    align_t = gettime() - align_t;
    fprintf(stderr, "%d\t%f\t%E f(%f, %f)\n", seq1len*seq2len,
        align_t, (seq1len*seq2len)/align_t), fill1_t, fill2_t;
}


int main (int argc, char * argv[]) {
    
    // mpi variables
    int rank, omp_rank, mpisupport, num_procs;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpisupport);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // get start time
    double t = gettime();
    
    // Display each command-line argument
    // show_args(argc, argv);
    if (argc <= 1) {
        fprintf(stderr, "sw <seqfile> [threads]\n");
        exit(1);
    }

    // get number of threads per processor (default is 1)
    // and set that number for OpenMP
    int num_threads = 1;
    if (argc == 3) {
        num_threads = atoi(argv[2]);
    }
    omp_set_num_threads(num_threads);
    int tot_threads = num_procs * num_threads;
    printf("running with %d procs and %d threads, %d total threads\n", num_procs, omp_get_max_threads(), tot_threads );

    
    // open file and initialize vars
    FILE *file;
    char lineAs[ALIGNMENTS][MAXLINELEN];
    char lineBs[ALIGNMENTS][MAXLINELEN];
    file = fopen(argv[1], "r");
    if (file == NULL) {
        perror ("Error reading file");
        exit(1);
    }

    int c;
    int num_lines;
    num_lines = ALIGNMENTS * 2;

    // read file into arrays
    for (c = 0; c < ALIGNMENTS; c++) {

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
    fclose(file);

    // set up chunking
    int chunk_count, chunk;
    chunk_count = ALIGNMENTS/tot_threads;

    // each thread does its own alignments
    #pragma omp parallel private(omp_rank, chunk)
    {
	for (chunk = 0; chunk < chunk_count; chunk++) {
	    int chunk_start = chunk * tot_threads;

	    omp_rank = omp_get_thread_num();
	    char line1 [MAXLINELEN];
	    char line2 [MAXLINELEN];
	    int line_num = (rank * num_threads) + omp_rank + chunk_start;
	    strcpy(line1, lineAs[line_num]);
	    strcpy(line2, lineBs[line_num]);                
	    
	    //fprintf(stderr, "process %d, thread %d - chunk %d/%d - seq %s\n", rank, omp_rank, chunk, chunk_count-1, line1);
	    do_alignment(line1, line2);
	}
    }


    // get and print elapsed time
    t = gettime() - t;  
    printf("process %d elapsed time: %f secs\n", rank, t);
    
    MPI_Finalize();
    return 0;
}
