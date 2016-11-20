
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
#include "pr.h"

extern grid* g;

using namespace std;

// Credits to :
// http://www.memoryhole.net/kyle/2012/06/a_use_for_volatile_in_multithr.html
rank_t qthread_dincr(rank_t *operand, rank_t incr)
{
    //*operand = *operand + incr;
    //return incr;
    
    union {
       rank_t   d;
       uint32_t i;
    } oldval, newval, retval;
    do {
         oldval.d = *(volatile rank_t *)operand;
         newval.d = oldval.d + incr;
         //__asm__ __volatile__ ("lock; cmpxchgq %1, (%2)"
         __asm__ __volatile__ ("lock; cmpxchg %1, (%2)"
                                : "=a" (retval.i)
                                : "r" (newval.i), "r" (operand),
                                 "0" (oldval.i)
                                : "memory");
    } while (retval.i != oldval.i);
    return oldval.d;
}

void 
pr_t::init(vertex_t a_vert_count, int a_iteration_count) 
{
    vert_count = a_vert_count;
    iteration_count = a_iteration_count;
    
    vert_rank = (rank_t*)mmap(NULL, sizeof(rank_t)*vert_count, 
                           PROT_READ|PROT_WRITE,
                           MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, 0 , 0);
    if (MAP_FAILED == vert_rank) {
        vert_rank = (rank_t*)calloc(sizeof(rank_t), vert_count);
        vert_rank_prior = (rank_t*)calloc(sizeof(rank_t), vert_count);
    } else {
    
        vert_rank_prior = (rank_t*)mmap(NULL, sizeof(rank_t)*vert_count, 
                           PROT_READ|PROT_WRITE,
                           MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, 0 , 0);
    
        if (MAP_FAILED == vert_rank_prior) {
            vert_rank_prior = (rank_t*)calloc(sizeof(rank_t), vert_count);
        }
    }
    
    double inv_vert_count = 1.0/vert_count;

    #pragma omp parallel for num_threads (NUM_THDS) schedule(dynamic, 1024)
    for(vertex_t i = 0 ;i < vert_count; i++)
    {
		vert_rank_prior[i] = inv_vert_count;
		//Initialize the rank
		if (g->vert_degree[i] != 0) { 
            g->vert_degree[i] = 1.0/g->vert_degree[i];
        }
    }
    memset(tmax, 0, sizeof(rank_t)*NUM_THDS);
}

void 
pr_t::pagerank_onepart(edge_t* part_edge, index_t cedge,  part_t i, part_t j)
{
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    rank_t* pri_rank0 = vert_rank_prior + offset0;
    rank_t*  rank1    = vert_rank + offset1; 
    #ifdef HALF_GRID
    rank_t* pri_rank1 = vert_rank_prior + offset1;
    rank_t*  rank0    = vert_rank + offset0; 
    #endif
    #else
    rank_t* pri_rank0 = vert_rank_prior;
    rank_t*  rank1    = vert_rank; 
    #ifdef HALF_GRID
    rank_t* pri_rank1 = vert_rank_prior;
    rank_t*  rank0    = vert_rank; 
    #endif
    #endif
                    
    //#pragma omp for nowait //2. bad lock, wb
    for (uint64_t k = 0 ; k < cedge; ++k) {
        vertex_t v0 = part_edge[k].v0;
        vertex_t v1 = part_edge[k].v1;
        //XXX rank_t d0 = vert_rank_prior[v0]/g->vert_degree[v0];
        rank_t d0 = pri_rank0[v0];
        qthread_dincr(rank1 + v1, d0); 
        #ifdef HALF_GRID
        //XXX rank_t d1 = vert_rank_prior[v1]/g->vert_degree[v1];
        rank_t d1 = pri_rank1[v1];
        qthread_dincr(rank0 + v0, d1); 
        #endif
    }

}

void pr_t::algo_mem_part(segment* seg)
{
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
            //#pragma omp for schedule (dynamic, 8) nowait 
            for (part_t j1 = j2; j1 <= j_end1; ++j1) {
                edge_t* part_edge = (edge_t*)(edges + (start_edge->get(i1, j1) 
                                                           << bytes_in_edge_shift));
                index_t cedge = start_edge->get_count(i1,j1);
                pagerank_onepart(part_edge, cedge, i1 + b_i, j1 + b_j);
            }
	    }
	}
}

void
pr_t::iteration_finalize(int current_iteration)
{
    /*
    #pragma omp for 
    for(vertex_t i = 0; i < vert_count; i++) {
        vert_rank[i] = 0.15 + 0.85*vert_rank[i];
		
        sdegree_t m = g->svert_degree[i];
		assert(m);
		if (m > 0) {
			tmax[omp_get_thread_num()] = max(tmax[omp_get_thread_num()],
                     abs(vert_rank[i] - m*vert_rank_prior[i]));
		
		} else {
			tmax[omp_get_thread_num()] = max(tmax[omp_get_thread_num()],
                     abs(vert_rank[i] - g->bvert_degree[-m]*vert_rank_prior[i]));
		}

        vert_rank_prior[i] = 0;
    }
    
    if (0 == omp_get_thread_num()) {
		rank_t max_diff = tmax[0];
		for(int i = 1 ; i < NUM_THDS; i++) {
			max_diff = max(max_diff, tmax[i]);
			tmax[i] = 0;
		}
		tmax[0] = 0;
		std::cout << " diff: "
			  << max_diff
			  << endl;
    }
        
    if (current_iteration < iteration_count) {
        #pragma omp for 
        for(vertex_t i = 0; i < vert_count; i++) {
            if (g->svert_degree[i] > 0) {
                vert_rank[i] = vert_rank[i]/g->svert_degree[i];
            } else {
                vert_rank[i] = vert_rank[i]/g->bvert_degree[-g->svert_degree[i]];
                assert(g->bvert_degree[-g->svert_degree[i]]);
            }
        }
    
        if (omp_get_thread_num() == 0) {
            //memset(vert_rank_prior, 0, vert_count*sizeof(rank_t));
            swap(vert_rank_prior, vert_rank);
        }
    }
    */
    
    if ((omp_get_thread_num() == 0) 
        && (current_iteration < iteration_count)) {
            memset(vert_rank_prior, 0, vert_count*sizeof(rank_t));
    }

    if (current_iteration < iteration_count) {
        #pragma omp for schedule (dynamic, 4096) 
        for(vertex_t i = 0; i < vert_count; i++) {
            vert_rank[i] = (0.15+0.85*vert_rank[i])*g->vert_degree[i];
        }
    } else {
        #pragma omp for schedule (dynamic, 4096)
        for(vertex_t i = 0; i < vert_count; i++) {
            vert_rank[i] = 0.15 + 0.85*vert_rank[i];
        }
    }
    
    if ((omp_get_thread_num() == 0) 
        && (current_iteration < iteration_count)) {
            //memset(vert_rank_prior, 0, vert_count*sizeof(rank_t));
            swap(vert_rank_prior, vert_rank);
    }
}
