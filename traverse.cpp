
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
#include <asm/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <algorithm>
#include <errno.h>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "wtime.h"
#include "traverse.h"


extern cid_t invalid_cid;
extern grid* g;

void traverse_t::init(vertex_t a_vert_count, vertex_t pivot, cid_t a_color )
{
    vert_count = a_vert_count;
    memset(front_count, 0, sizeof(vertex_t)* NUM_THDS);
    
    vert_cid = (cid_t*)mmap(NULL, sizeof(cid_t)*vert_count, 
                           PROT_READ|PROT_WRITE,
                           MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, 0 , 0);
    
    if (MAP_FAILED == vert_cid) {
        vert_cid = (cid_t*)calloc(sizeof(cid_t), vert_count);
        memset(vert_cid, invalid_cid, sizeof(cid_t)*vert_count);
    }
    
    vert_cid[pivot] = a_color;
    color = a_color;
    iteration = 0;
}

void traverse_t::algo_mem_part(segment* seg)
{
    index_t t_front_count = 0; 

	index_t b_i, b_j;
    part_t big_i, big_j;
    spart_t i, j, i_end, j_end;

    #ifdef HALF_GRID
	matrix<spart_t, index_t> start_edge_half;
	start_edge_half.part_count = p_p;
    #endif

	matrix_f<spart_t, index_t> start_edge_full;
	start_edge_full.part_count = p_p;
    matrix<spart_t, index_t>* start_edge = 0;
	    
    part_meta_t* meta = seg->meta;
	index_t ctx_count = seg->ctx_count;
    
    //#pragma omp for schedule (dynamic, 1) nowait
    for (index_t l = 0; l < ctx_count; ++l) {
        get_ij(meta[l].start, big_i, big_j, i, j);
        get_s_ij(meta[l].end, i_end, j_end);

		#ifdef HALF_GRID	
        if (big_i == big_j) {
            start_edge = &start_edge_half;
            start_edge->val = g->_s_start_edge + 
                beg_edge_offset1(big_i);
        } else {
            start_edge = &start_edge_full;
            start_edge->val = g->_s_start_edge + 
                beg_edge_offset2(big_i, big_j);
        }
		#else
		start_edge = &start_edge_full;
        start_edge->val = g->_s_start_edge + 
        beg_edge_offset2(big_i, big_j);
	    #endif
        
        b_i = (big_i << bit_shift3);
	    b_j = (big_j << bit_shift3);

	    char* buf = seg->buf;
	    
	    // Align new offset. Add the offset from start edge of i,j 
	    char* new_offset = buf + meta[l].offset
                + ((start_edge->get(i,j) << bytes_in_edge_shift) & 0x1FF);

	    char* edges = new_offset 
                      - (start_edge->get(i, j) << bytes_in_edge_shift);

	    part_t j2 = j;
	    part_t j_end1 = p_p - 1;
        
        #pragma omp for schedule (dynamic, 1) nowait
	    for (part_t i1 = i; i1 <= i_end; ++i1) {
            j2 = 0;
		    if (i1 == i_end) j_end1 = j_end; 
            if (i1 == i) j2 = j; 
			#ifdef HALF_GRID
            else if(b_i == b_j)  j2 = i1; 
			#endif
		
            for (part_t j1 = j2; j1 <= j_end1; ++j1) {
                edge_t* part_edge = (edge_t*)(edges + 
                    (start_edge->get(i1, j1) << bytes_in_edge_shift));
                index_t cedge = start_edge->get_count(i1,j1);
                t_front_count += traverse_onepart(part_edge, cedge,
                               i1 + b_i, j1 + b_j);
                   
            }
	    }
    }
    
    front_count[omp_get_thread_num()] += t_front_count;
    return;
}

index_t traverse_t::traverse_onepart(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    index_t t_front_count = 0;
    
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    cid_t* vert_cid0 = vert_cid + offset0;
    cid_t* vert_cid1 = vert_cid + offset1;
    #else
    cid_t* vert_cid0 = vert_cid;
    cid_t* vert_cid1 = vert_cid;
    #endif

    vertex_t v0,v1;
    cid_t c0, c1;    

    for (uint64_t k = 0 ; k < cedge; ++k) {
        v0 = part_edge[k].v0;
        v1 = part_edge[k].v1;
  
		c0 = vert_cid0[v0];		
		c1 = vert_cid1[v1];

        if (c0 == color && c1 == invalid_cid) {
            vert_cid1[v1] = c0;
            ++t_front_count;
            if (iteration > 50) {
                cout << ", "<< offset0 + v0 << ":"  << offset1 + v1 << "(" 
                     << g->vert_degree[offset1 + v1] << ") ";
            }
        }
    }
    return t_front_count;
}

int traverse_t::iteration_finalize() 
{
    index_t t_front_count = 0;
    
    for(int ithd = 0; ithd < NUM_THDS; ++ithd) {
        t_front_count += front_count[ithd];
        front_count[ithd] = 0;
    }

    cout << "level: " << iteration  
         << " frontier count: " << t_front_count <<endl
         << endl;
    
    if (t_front_count == 0) {
        return 0;
    }
    
    ++iteration;
    return t_front_count; 
}
    
