/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>
#include <stdio.h>
#include <assert.h>
#include <stdio.h>
#include <mpi.h>
#include "make_graph.h"
#include <assert.h>

off_t fsize(const char *filename) {
    struct stat st;
    if (stat(filename, &st) == 0)
        return st.st_size;
//    fprintf(stderr, "Cannot determine size of %s: %s\n",
//            filename, strerror(errno));
    return -1;
}

int main(int argc, char** argv) 
{
  int64_t log_numverts, degree;
  int size, rank;
  unsigned long my_edges;
  unsigned long global_edges;
  double start, stop;
  size_t i;
	
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(rank==0)fprintf(stdout,"Input: /path/to/exe log_numverts degree row-partitions column-partitions\n");
	
	//fprintf(stdout,"%d\n",argc);
	if(argc!= 5){ 
		printf("Wrong input\n");
		exit(-1);
	}
	
  log_numverts = (int64_t) atoi(argv[1]);
  degree = (int64_t) atoi(argv[2]);
	int row_par=atoi(argv[3]);
	int col_par=atoi(argv[4]);
	unsigned long sz_row_par=(INT64_C(1)<<log_numverts)/row_par;
	unsigned long sz_col_par=(INT64_C(1)<<log_numverts)/col_par;
	//if(rank == 0) printf("row-par-sz: %lu, col-par-sz: %lu\n",sz_row_par,sz_col_par);
	//fprintf(stdout,"%d reach here\n",rank);
	if(row_par>256 || col_par>256){fprintf(stderr,"Too many partions\n");}

	//exit(-1);
  if (rank == 0) fprintf(stderr, "Graph size is %" PRId64 " vertices and %" PRId64 " edges\n", INT64_C(1) << log_numverts, degree << log_numverts);

  /* Start of graph generation timing */
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  int64_t nedges;
  packed_edge* result;
  make_graph(log_numverts, degree << log_numverts, 1, 2, &nedges, &result);
  MPI_Barrier(MPI_COMM_WORLD);
  stop = MPI_Wtime();
  /* End of graph generation timing */

  my_edges = nedges;

  for (i = 0; i < my_edges; ++i) {
    assert ((get_v0_from_edge(&result[i]) >> log_numverts) == 0);
    assert ((get_v1_from_edge(&result[i]) >> log_numverts) == 0);
  }
  

  MPI_Reduce(&my_edges, &global_edges, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    fprintf(stderr, "%lu edge%s generated in %fs (%f Medges/s on %d processor%s)\n", global_edges, \
						(global_edges == 1 ? "" : "s"), (stop - start), global_edges / (stop - start) * 1.e-6, size, (size == 1 ? "" : "s"));
  }

	//support 256 partions in row and col, each file name 256 chars at most
	FILE ***fid = (FILE***)malloc(sizeof(FILE **)*256);
	char ***filename=(char ***)malloc(sizeof(char **)*256);

	for(int i=0;i<256;i++)
	{
		fid[i] = (FILE **)malloc(sizeof(FILE *)*256);
		filename[i]=(char **)malloc(sizeof(char *)*256);
		for(int j=0;j<256;j++)
			filename[i][j]=(char *)malloc(sizeof(char)*256);
	}

//TEXT	
	for(int i=0;i<row_par;i++)
		for(int j=0;j<col_par;j++)
		{
			sprintf(filename[i][j],"scale-%lu-rank-%d-of-%d-par-rowcol-%d-%d.dat",log_numverts,rank,size,i,j);
			fid[i][j]=fopen(filename[i][j],"w");
			if(fid[i][j]==NULL){fprintf(stderr, "Wrong fopen\n");}
		}
	
	for(i=0;i<my_edges;++i)
	{
		vertex_t src=result[i].v0;
		vertex_t des=result[i].v1;
		
		//store current edge
		fprintf(fid[src/sz_row_par][des/sz_col_par],"%lu\t%lu\n",src,des);

		//XXX store reverse edge
		//fprintf(fid[des/sz_col_par][src/sz_row_par],"%lu\t%lu\n",des,src);
	}

//--- end of text
/*
//BINARY
	for(int i=0;i<row_par;i++)
		for(int j=0;j<col_par;j++)
		{
			sprintf(filename[i][j],"scale-%lu-degree-%lu-rank-%d-of-%d-par-rowcol-%d-%d.bin",log_numverts,degree,rank,size,i,j);
			fid[i][j]=fopen(filename[i][j],"wb+");
			if(fid[i][j]==NULL){fprintf(stderr, "Wrong fopen\n");}
		}
	
	//binary format
	for(i=0;i<my_edges;++i)
	{
		vertex_t src=result[i].v0;
		vertex_t des=result[i].v1;
	
		//store current edge
		fwrite(result+i,sizeof(packed_edge),1,fid[src/sz_row_par][des/sz_col_par]);
	
		//store reverse edge
        //XXX: uncommet it.
		//result[i].v0=des;
		//result[i].v1=src;
		//fwrite(result+i,sizeof(packed_edge),1,fid[des/sz_row_par][src/sz_col_par]);
	}
*/
//End of binary
	for(int i=0;i<row_par;i++)
		for(int j=0;j<col_par;j++)
			fclose(fid[i][j]);
//VERIFY
//	for(int i=0;i<row_par;i++)
//		for(int j=0;j<col_par;j++)
//		{
//			fid[i][j]=fopen(filename[i][j],"rb");
//			if(fid[i][j]==NULL){fprintf(stderr, "Wrong fopen\n");}
//			
//			off_t size_offset = fsize((const char *)filename[i][j]);
//			packed_edge *check_res = (packed_edge*) malloc(size_offset);
//			fread(check_res, sizeof(packed_edge), size_offset/sizeof(packed_edge), fid[i][j]);
//			
//			char tempfilename[256];
//			sprintf(tempfilename,"scale-%lu-rank-%d-of-%d-par-rowcol-%d-%d.txt",log_numverts,rank,size,i,j);
//			FILE *tempfile=fopen(tempfilename,"w");
//			for(int k=0;k<size_offset/sizeof(packed_edge);k++)
//				fprintf(tempfile,"%lu,%lu\n",check_res[k].v0,check_res[k].v1);
//			
//			fclose(tempfile);
//			free(check_res);
//		}

	
	/*FILE *fid=fopen(filename,"wb");
		fwrite(result,sizeof(packed_edge),nedges,fid);
	fclose(fid);
	*/
	
	/*int64_t src_min, src_max;
	int64_t dest_min,dest_max;
		
	src_min=result[0].v0;
	src_max=result[0].v0;
 	dest_min=result[0].v1;
	dest_max=result[0].v1;

	for(int64_t i=1;i<my_edges;i++)
	{
		if(result[i].v0<src_min) src_min=result[i].v0;
		if(result[i].v0>src_max) src_max=result[i].v0;
		if(result[i].v1<dest_min) dest_min=result[i].v1;
		if(result[i].v1>dest_max) dest_max=result[i].v1;
	}
	
	if(rank==0)
		printf("tid0-Source: %lu ~ %lu. Dest: %lu ~ %lu\n",src_min,src_max,dest_min,dest_max);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==1)
		printf("tid0-Source: %lu ~ %lu. Dest: %lu ~ %lu\n",src_min,src_max,dest_min,dest_max);
	MPI_Barrier(MPI_COMM_WORLD);
	*/


	/*packed_edge *check_res;
	fid=fopen(filename,"rb");
	off_t size_offset = fsize((const char *)filename);
	check_res = (packed_edge*) malloc(size_offset);
	fread(check_res, sizeof(packed_edge), nedges, fid);
	assert(0 == memcmp(result, check_res, nedges));
	*///verify whether the write out is correct.


  free(result);
	MPI_Finalize();
  return 0;
}
