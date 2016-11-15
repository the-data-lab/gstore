
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
#include <algorithm>

#include "wtime.h"
#include "kcore.h"

extern grid* g;

void kcore_t::init(vertex_t a_vert_count, int k, int32_t* pdegree, bitmap_t* a_read_part)
{
    read_part = a_read_part;
    vert_count = a_vert_count;
    kc = k;
	active_vert = new bitmap_t(vert_count);
}

void kcore_t::kcore_onepart(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    #endif
    for (uint64_t k = 0 ; k < cedge; ++k) {
        #ifdef COMPACT_GRID
        vertex_t v0 = offset0 + part_edge[k].v0;
        vertex_t v1 = offset1 + part_edge[k].v1;
        #else 
        vertex_t v0 = part_edge[k].v0;
        vertex_t v1 = part_edge[k].v1;
        #endif
        if(active_vert->get_bit(v0)) { 
			__sync_fetch_and_sub(g->vert_degree + v1, 1);
		/*	
            sdegree_t m = g->svert_degree[v1];
			if (m > 0) {
				__sync_fetch_and_sub(g->svert_degree + v1, 1);
			} 
			else if (m < 0) {
				__sync_fetch_and_sub(g->bvert_degree - m, 1);
			}
        */
		}
        #ifdef HALF_GRID
        if (active_vert->get_bit(v1)) {
			__sync_fetch_and_sub(g->vert_degree + v0, 1);
		/*
			sdegree_t m = g->svert_degree[v0];
			if (m > 0) {
				__sync_fetch_and_sub(g->svert_degree + v0, 1);
			} 
			else if (m < 0) {
				__sync_fetch_and_sub(g->bvert_degree - m, 1);
			}
         */
        }
        #endif
    }
    #ifdef HALF_GRID
    g->read_part_next->set_bit_atomic(i);
    #endif
    g->read_part_next->set_bit_atomic(j);
}

void kcore_t::kcore_onepart_col(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    #ifdef HALF_GRID
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    #endif
    for (uint64_t k = 0 ; k < cedge; ++k) {
        #ifdef COMPACT_GRID
        vertex_t v0 = offset0 + part_edge[k].v0;
        vertex_t v1 = offset1 + part_edge[k].v1;
        #else 
        vertex_t v0 = part_edge[k].v0;
        vertex_t v1 = part_edge[k].v1;
        #endif
        /*
        vert_degree[v0] -= ((prior_vert_degree[v1] < kc) 
                           && (prior_vert_degree[v1] > 0) 
                           && (prior_vert_degree[v0] > 0));
                           */
        if (active_vert->get_bit(v1)) { 
			__sync_fetch_and_sub(g->vert_degree + v0, 1);
            /*
			sdegree_t m = g->svert_degree[v0];
			if (m > 0) {
				__sync_fetch_and_sub(g->svert_degree + v0, 1);
			} 
			else if (m < 0){
				__sync_fetch_and_sub(g->bvert_degree - m, 1);
			}
             */
		}
        //__sync_fetch_and_sub(vert_degree + v0, 1);
    }
    g->read_part_next->set_bit_atomic(i);
    #endif
}

void kcore_t::kcore_onepart_row(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    #ifdef HALF_GRID
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    #endif
    for (uint64_t k = 0 ; k < cedge; ++k) {
        #ifdef COMPACT_GRID
        vertex_t v0 = offset0 + part_edge[k].v0;
        vertex_t v1 = offset1 + part_edge[k].v1;
        #else 
        vertex_t v0 = part_edge[k].v0;
        vertex_t v1 = part_edge[k].v1;
        #endif
        if(active_vert->get_bit(v0)) { 
			__sync_fetch_and_sub(g->vert_degree + v1, 1);
            /*
			if (m > 0) {
				__sync_fetch_and_sub(g->svert_degree + v1, 1);
			} 
			else if (m < 0) {
				__sync_fetch_and_sub(g->bvert_degree - m, 1);
			}
             */
		}
    }
    g->read_part_next->set_bit_atomic(j);
    #endif
}

int kcore_t::iteration_finalize() 
{
    //Find out if we need to load this partition.
    active_vert->reset(); 
    vertex_t t_front_count = 0;
    #pragma omp parallel num_threads(NUM_THDS)
    {
        #pragma omp for reduction(+:t_front_count)
        for (vertex_t ivert = 0; ivert < vert_count; ++ivert) {
			//sdegree_t m = g->svert_degree[ivert];
			sdegree_t m = g->vert_degree[ivert];
            if(m < kc && m > 0) {
                //g->svert_degree[ivert] = GONE;
                g->vert_degree[ivert] = GONE;
                ++t_front_count;
				active_vert->set_bit(ivert);
            } 
            /*
			else if (m < 0) {
				if(g->bvert_degree[-m] <= 32767) {
					g->svert_degree[ivert] = g->bvert_degree[-m];
					m = g->svert_degree[ivert];
					if(m < kc && m > 0) {
						g->svert_degree[ivert] = GONE;
						++t_front_count;
						active_vert->set_bit(ivert);
					}
				}
			}
             */
        }
    }

    cout <<" frontier count: "<<t_front_count << endl;

    if (t_front_count == 0) return 0;
    calc_part_needed();
    return t_front_count;
}
    
void kcore_t::algo_mem_part(segment* seg)
{
    //double start = mywtime();
	index_t b_i, b_j;
    part_t big_i, big_j;
    spart_t i, j, i_end, j_end;

	matrix<spart_t, index_t> start_edge_half;
	start_edge_half.part_count = p_p;
	matrix_f<spart_t, index_t> start_edge_full;
	start_edge_full.part_count = p_p;
    matrix<spart_t, index_t>* start_edge = 0;
	    
    part_meta_t* meta = seg->meta;
	index_t ctx_count = seg->ctx_count;
    
	#pragma omp for schedule (dynamic, 1) nowait 
    for (index_t l = 0; l < ctx_count; ++l) {
        get_ij(meta[l].start, big_i, big_j, i, j);
        get_s_ij(meta[l].end, i_end, j_end);

        if (big_i == big_j) {
            start_edge = &start_edge_half;
            start_edge->val = g->_s_start_edge + 
                beg_edge_offset1(big_i);
        } else {
            start_edge = &start_edge_full;
            start_edge->val = g->_s_start_edge + 
                beg_edge_offset2(big_i, big_j);
        }
	    
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

	    for (part_t i1 = i; i1 <= i_end; ++i1) {
			if (i1 == i_end) j_end1 = j_end; 
			
			if (read_part->get_bit(i1 + b_i)) {
				//#pragma omp for schedule (dynamic, 8) nowait 
				for (part_t j1 = j2; j1 <= j_end1; ++j1) {
				edge_t* part_edge = (edge_t*)(edges + 
					(start_edge->get(i1, j1) << bytes_in_edge_shift));
				index_t cedge = start_edge->get_count(i1,j1);
				kcore_onepart(part_edge, cedge,
							   i1 + b_i, j1 + b_j);
				}   
			}
			#ifdef  HALF_GRID
			else {
				//#pragma omp for schedule (dynamic, 8) nowait 
				for (part_t j1 = j2; j1 <= j_end1; ++j1) {
				if (read_part->get_bit(j1 + b_j)) { 
					edge_t* part_edge = (edge_t*)(edges + 
							(start_edge->get(i1, j1) << bytes_in_edge_shift));
					index_t cedge = start_edge->get_count(i1,j1);
					kcore_onepart_col(part_edge,
							cedge, i1 + b_i, j1 + b_j );
				}
				}   
			}
			#endif
			if(b_i == b_j) {
				j2 = i1 + 1;
			} else {
				j2 = 0;
			}
	    }
	}
    return;
}
    
