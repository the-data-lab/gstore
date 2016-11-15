

/*
 * Copyright 2016 The George Washington University
 * Written by Pradeep Kumar 
 * Directed by Prof. Howie Huang
 *
 * https://www.seas.gwu.edu/~howie/
 * Contact: iheartgraph@gmail.com
 *
 * 
 * Please cite the following paper:
 * 
 * Pradeep Kumar and H. Howie Huang. 2016. G-Store: High-Performance Graph Store for Trillion-Edge Processing. In Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC '16).
 
 *
 * This file is part of G-Store.
 *
 * G-Store is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * G-Store is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with G-Store.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <omp.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <algorithm>
#include <errno.h>
#include <cmath>
#include <fstream>
#include <math.h>

#include "wtime.h"
#include "gstore.h"

inline off_t fsize(const char *filename) {
    struct stat st; 
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1; 
}
void text_to_bin(string textfile)
{
	int fd;
	char* ss_head;
	char* ss;

	size_t file_size = fsize(textfile.c_str());
	fd = open( textfile.c_str(), O_CREAT|O_RDWR, 00777);

	ss_head = (char*)mmap(NULL,file_size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);

	size_t head_offset=0;
	while(ss_head[head_offset]=='%'){
		while(ss_head[head_offset]!='\n'){
			head_offset++;
		}
		head_offset++;
	}
	ss = &ss_head[head_offset];
	file_size -= head_offset;

	size_t curr=0;
	size_t next=0;

	//step 1. vert_count,edge_count,
	size_t edge_count=0;
	size_t vert_count;
	uint32_t v_max = 0;
	uint32_t v_min = 999999;//as infinity
	vertex_t a;
	while(next<file_size){
		char* sss=ss+curr;
		a = atoi(sss);

		if(v_max<a){
			v_max = a;
		}
		if(v_min>a){
			v_min = a;
		}

		while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
			next++;
		}
		while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
			next++;
		}
		curr = next;
        edge_count++; 
	}
	edge_count /=2;
	vert_count = v_max - v_min + 1;

	cerr<<"max vertex id: "<<v_max<<endl;
	cerr<<"min vertex id: "<<v_min<<endl;

	cerr<<"edge count: "<<edge_count<<endl;
	cerr<<"vert count: "<<vert_count<<endl;

	//step 4: write adjacent list 
	uint32_t v0;
    uint32_t v1;
	size_t offset =0;
	next = 0;
	curr = 0;

    string edgefile = textfile + ".edge";    
	//int fd4 = open( edgefile.c_str(), O_CREAT|O_RDWR,00777 );
	//ftruncate(fd4, edge_count*sizeof(gedge_t));
	//gedge_t* adj = (gedge_t*)mmap(NULL,edge_count*sizeof(gedge_t),
    //                                PROT_READ|PROT_WRITE,MAP_SHARED,fd4,0);
	gedge_t* adj = (gedge_t*)malloc(edge_count*sizeof(gedge_t));
	
	while(next < file_size) {
	    char* sss = ss+curr;
	    v0 = atoi(sss)-v_min;//first end of pair
	    while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
		    next++;
	    }
	    while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
		    next++;
	    }
	    curr = next;

	    char* sss1=ss+curr;
	    v1 = atoi(sss1)-v_min;//second end of pair
         
	    adj[offset].set_v0(v0);
	    adj[offset].set_v1(v1);
		

	    while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
		    next++;
	    }
	    while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
		    next++;
	    }
	    curr = next;

	    offset++;
	}
	
	munmap( ss,sizeof(char)*file_size );
	//munmap( adj,sizeof(vertex_t)*edge_count );
	close(fd);
	FILE* fd4 = fopen(edgefile.c_str(), "wb");
    assert(fd4 != 0);
    fwrite(adj, sizeof(gedge_t), edge_count, fd4);
    
    fclose(fd4);
}

void remove_dup(string edgefile)
{
    //read the binary edge file
    int fid_edge = open(edgefile.c_str(), O_RDONLY);
    struct stat st_edge;
    fstat(fid_edge, &st_edge);

    assert(st_edge.st_size != 0);

    index_t nedges = st_edge.st_size/sizeof(gedge_t);
    gedge_t* edges = (gedge_t*)mmap(0, st_edge.st_size, PROT_READ,
                                  MAP_PRIVATE, fid_edge, 0);
        
    if (MAP_FAILED == edges) {
        handle_error("failed to open file");
    }

    //index_t* edge_count = (index_t*)calloc(sizeof(index_t), p*p);
    index_t* edge_count[NUM_THDS] = {0};
    index_t*  edge_start = (index_t*)calloc(sizeof(index_t), p*p);
    index_t*  edge_cnt = (index_t*)calloc(sizeof(index_t), p*p);
    index_t prefix_sum = 0;
    index_t total_part = p*p;
    
    #pragma omp parallel num_threads(NUM_THDS) 
    {
        edge_count[omp_get_thread_num()] = (index_t*)calloc(sizeof(index_t), p*p);
    
        //---classification: dry run
        uint32_t p1, p2 ;
        for(uint64_t k = 0; k < nedges; ++k) {
            //if (edges[k].is_self_loop()) continue;
            p1 = (edges[k].get_v0() >> bit_shift2);
            p2 = (edges[k].get_v1() >> bit_shift2);
            ++edge_count[omp_get_thread_num()][p1*p + p2];
        }   

        #pragma omp for
        for (index_t ipart = 0; ipart < total_part; ++ipart) {
            for (index_t ithd = 0; ithd < NUM_THDS; ++ithd) {
                edge_cnt[ipart] += edge_count[ithd][ipart];
            }
        }
        free(edge_count[omp_get_thread_num()]);
    }

    for (index_t ipart = 0; ipart < total_part; ++ipart) {
        edge_start[ipart] = prefix_sum;
        prefix_sum += edge_cnt[ipart];
        edge_cnt[ipart] = 0; 
    }  
    assert(prefix_sum == nedges);

    //---classification-- actual run
    gedge_t* new_edge = (gedge_t*)malloc(sizeof(gedge_t)*prefix_sum);

    #pragma omp parallel for num_threads(NUM_THDS)
    for(index_t k = 0; k < nedges; ++k) {
        //if (edges[k].is_self_loop()) continue;
        part_t p1 = (edges[k].get_v0() >> bit_shift2);
        part_t p2 = (edges[k].get_v1() >> bit_shift2);
        index_t n = p1*p + p2;
        index_t m = __sync_fetch_and_add(edge_cnt + n, 1);
        new_edge[edge_start[n]+ m] = edges[k];
    }
    free(edge_cnt);
    //---remove duplicate edges
    //Kron graph has lots of duplicate edges.
    #pragma omp parallel for collapse(2) //schedule (dynamic, 1)
    for (uint32_t i = 0; i < p; ++i) {
        for (uint32_t j = 0; j < p; ++j) {
            index_t n = i*p + j;
            uint64_t cedge = edge_start[n + 1] - edge_start[n]; 
            gedge_t* part_edge = new_edge + edge_start[n];

            for (uint64_t k = 0 ; k < cedge; ++k) { //by this time all genuine self loop has been removed.
                if ((part_edge[k].get_v0() == part_edge[k].get_v1())) continue;
                for (uint64_t l = k+1; l < cedge; ++l) {
                    if (part_edge[k] == part_edge[l]) {
                        //convert it to self loop so that it
                        //doesn't interfere with the further 
                        //comparison.
                        part_edge[l].set_v0(part_edge[l].get_v1());
                    } 
                }
            }
        }
    }

    //---Write the processed file
    string file = edgefile + ".nodup";
    FILE* f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(new_edge, sizeof(gedge_t), nedges, f);
    
}

void grid::analyze_grid_size(string edgefile)
{
	read_start_in_mem(edgefile);
	
    index_t total_s_part = calc_total_part(p_s);
	index_t total_b_part = calc_total_part(p);

	index_t edge_count = _s_start_edge[total_s_part];

	index_t* b_start_edge = (index_t*)malloc((total_b_part + 1)* sizeof(index_t));
	
	index_t offset;

	for(part_t i = 0; i < p; ++i) {
		index_t n = 0;
		#ifdef HALF_GRID
		n = i;
		#endif

		for(index_t j = n; j < p; ++j) {
			#ifdef HALF_GRID
			if (i == j) {
				offset = beg_edge_offset1(i);
			} else {
				offset = beg_edge_offset2(i, j);
			}
			#else 
				offset = beg_edge_offset2(i, j);
			#endif
			index_t index = calc_index(i, j, p);
			b_start_edge[index] = _s_start_edge[offset];
		}	
	}
	b_start_edge[total_b_part] = edge_count;

	//convert to size 
	for (index_t i = 0; i < total_s_part; ++i ) 
	{
		_s_start_edge[i] = _s_start_edge[i+1] - _s_start_edge[i];
	}

	for (index_t i = 0; i < total_b_part; ++i) {
		 b_start_edge[i] = b_start_edge[i+ 1] -  b_start_edge[i]; 
	}

	sort(b_start_edge, b_start_edge + total_b_part);
	sort(_s_start_edge, _s_start_edge + total_s_part);

	//Write in a txt file
	string file = edgefile + ".ssize";
	FILE* f = fopen(file.c_str(), "w");
	index_t k = 0;
	index_t value = 0;
	index_t sum0 = 0;
	//index_t sum1 = 0;
	//index_t sum2 = 0;

	for (index_t i = 0; i < total_s_part; ++i) {
		//if (value != _s_start_edge[i]) {
			value = _s_start_edge[i];
			fprintf(f, "%lu\t%lu\n", k, value);
			++k;
		//}
		if (value > 100000) sum0 += value;
	}
	cout << "sum0 = " << sum0 << endl;
	fclose(f);

	file = edgefile + ".bsize";
	f = fopen(file.c_str(), "w");
	k = 0;
	value = 0;
	for (index_t i = 0 ; i < total_b_part; ++i) {
		if (value != b_start_edge[i]) {
			value = b_start_edge[i];
			fprintf(f, "%lu\t%lu\n", k, value);
			++k;
		}
	}
	fclose(f);
}

void conv_to_text(string edgefile)
{
    struct stat st_edge;
    stat(edgefile.c_str(), &st_edge);
    assert(st_edge.st_size != 0);
    
    FILE* f = fopen(edgefile.c_str(), "rb");
    assert(f != 0);
    gedge_t* edges = (gedge_t*) malloc(st_edge.st_size);
    assert(edges);
	index_t cedge = st_edge.st_size/sizeof(gedge_t);
    fread(edges, sizeof(gedge_t), cedge, f);
    fclose(f);
	
	string file = edgefile + ".txt";
    f = fopen(file.c_str(), "w");
	vertex_t v0, v1;	
	for(index_t i = 0; i < cedge; ++i) {
		v0 = edges[i].get_v0();
		v1 = edges[i].get_v1();
		fprintf(f, "%lu\t%lu\n", v0, v1);	
	}
	fclose(f);
	return;
}
