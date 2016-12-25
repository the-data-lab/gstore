# gstore
*Semi-external graph engine for trillion-edge graphs.*

## Help
We will be updating this file as and when required and will serve as help file.

### Building:
  `make gstored` and `make gstoreu` or `make all`
  
`gstored` is for directed graphs. `gstoreu` is for undirected graphs. 
`libaio` is required before you build the software.

### How to run
`gstore` has proposed a new storage format called *tile based represenation* which takes advantage of *symmetry* and *smallest number of bits (SNB)* format. So, a graph need to be converted in that format before you can run.  

* `Graph generation`: We have modified Graph500 generator little bit and have added the source code to generate a kronecker graph. Go inside `graph500-generator` directory and run `make`. You need `mpi` to be installed. Thereafter run `./generator_test_mpi 25 16 1 1` (single-threaded) or  `mpirun ./generator_test_mpi 16 16 1 1` (multi-process) to generate an kronecker-25-16 graph. This will have 2^25 vertices and 2x16x2^25 edges (assuming undirected). If you run multi-process one, than you need to concatenate all the generated files in one file. At the end you will get a binary edge-list file. Lets call the above generated file as *kron_25_16b.dat*

* `Graph conversion`: Run following command to convert the *kron_25_16b.dat* to *tile based representation*:

`./gstoreu -s 25 -i kron_25_16b.dat -c 1 -o tile_25_16b.dat`

It will generate some files named tile_25_16b.dat.start, tile_25_16b.dat.grid etc.

* `Running page rank`: To run 5 iteration of pagerank job on scale 25 kronecker graph run following command.

    `./gstoreu -s 25 -i tile_25_16b.dat -j 1 -a5` 

Please see gstore.cpp main() function for meaning of different parameter such as i, j, o, c, a etc.


### Things to tune:
#### gstore.h
*  `NUM_THDS` : Total number of threads to spawn using openMP. We have set its value as 56. Change it prior to compilation.
  Also, pinning the openmp threads further helps in achieving good performance. This generally requires defining and exporting some openMP environment variable before running gstored/gstoreu. Consult the openMP documentation.
  
#### gstore.cpp  
*  `total_size`: total memory size that you want to set for streaming and caching. Default value is set as 8GB. There is a run time input for this variable.
* `memory`: Size of segment. There is two segment used for IO and processing alternatly. Make sure that you set memory and total_size in such a way that total_size/memory gives you an integer.
*  `IO_THDS`: Number of IO threads to deploy for doing IO. Ideally, you should have just one thread for one SSD but more threads for multiple SSD raid. Feel free to experiment with it.
* `AIO_MAXIO` 16384 
* `AIO_BATCHIO` 256

The above two variables decide the IO depth and IO batching crietria. AIO_BATCHIO count IO will be submitted in each system call of libaio. And the IO will behave as asynchronous till the count of submitted IO (through multiple submission call) reaches beyond AIO_MAXIO.  After that, we start reaping the IO and sending thus maintaining the IO depth. Set it to smaller number if you have fewer SSDs. Such as 16 or 32 batching criteria per-SSD.
    
#### Other tunables:
*  `Huge Page support`: We have enabled huge page (2MB huge page is only supported, 1GB pages can be supported easily with slight change) support for many algorithms. To take advantage of huge page support, we use anonymous mmap calls, hence, user has to enable to 2MB page support through boot time parameter(That's all). 
* `Trasnparent huge pages (THP)` : Transparent huge pages reduces IO throuhgput that we can get out of a single thread. Please disable it for better performance.
  
There are plently of other settings for IO optimizations such as number of request containers in old block-layer (The new block-layer is called blk-mq). Also the thread completion affinity for IO completions and many more.
  
