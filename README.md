README for dnmsw
----------------

Usage
-----

To compile serial versions:
`% gcc -O3 -lm -o sw_serial_v2 sw_serial_v2.c`

To run serial versions:
`% sw_serial_v2 64seqs.txt`


To compile parallel versions:
`% module add gcc openmpi/1.3-gcc
% mpicc -O3 -fopenmp -o sw_hybrid_v7 sw_hybrid_v7.c`

To run parallel versions:
`% spattach -i -n 4 time ./run_job sw_hybrid_v7 4 8 2> v7_4_8.128_1.out`
 OR
`% spattach -j 100415520913 time ./run_job sw_hybrid_v7 4 8 2> v7_4_8.128_1.out`
