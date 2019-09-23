Based on graph500 generator. It will genrate **text** graph files. If you need to generate binary graph file, either convert the text file to binary file, or comment out the text generation (line no. 108-129) and enable binary generation (line no 130-155) in the file https://github.com/the-data-lab/gstore/blob/master/graph500-generator/generator_test_mpi.c

Alternately, use the generator https://github.com/pradeep-k/gConv for binary file generation. We plan to migrate to this generator anyway as it is scalable.

Compilation require presence of MPI softwares. 
But generation can be done in single machine.

Steps:
=====
- Run `make`
- Run `./generator_test_mpi 25 16 1 1` for single threaded single machine generation of kronecker graph of 2^25 vertices and (2^25)*16 edges. It will generate a single text graph file.
- or `mpiexec -np 16 ./generator_test_mpi 25 16 1 1` for multi-process version in a single machine. It will generate 16 files for the same graph. 
- Content will look like this:

```
865705  1136538
566992  1463521
281809  1825156
1947818 1597918
598735  1213461
1634626 1213392
1075126 82992
548811  1561989
1198311 1122171
1305322 1327260
1515989 512979
1008789 1345790
1479024 2043758
1371049 136926
977545  1825573
1199337 947657
294976  29456
1747046 1850680
1008153 1643820
1177859 1636115
356390  1199737
919832  560391
.... many more
```
