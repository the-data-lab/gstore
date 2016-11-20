
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
#include "bfs2.h"

extern grid* g;
extern const depth_t INFTY;


void bfs2_t::init(vertex_t vert_count, vertex_t root, bitmap_t* a_read_part)
{
    memset(front_count, 0, sizeof(vertex_t)* NUM_THDS);
    
    depth = (depth_t*)mmap(NULL, sizeof(depth_t)*vert_count, 
                           PROT_READ|PROT_WRITE,
                           MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, 0 , 0);
    
    if (MAP_FAILED == depth) {
        if(posix_memalign((void**)&depth, 2097152 , sizeof(depth_t)*vert_count)) {
			perror("posix_memalign");
            assert(0);
			return ;
        }
        memset(depth, INFTY, sizeof(depth_t)*vert_count);
    }
    
    depth[root] = 1;
    level = 1;
    cout << "Root = " << root << endl; 
    //calc_part_needed();
	g->read_part->set();
	ioskip = 1;
}

index_t bfs2_t::bfs_onepart(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    index_t t_front_count = 0;
    index_t t_f_count = 0;
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    depth_t* d0 = depth + offset0;
    depth_t* d1 = depth + offset1;
    #else
    depth_t* d0 = depth;
    depth_t* d1 = depth;
    #endif

	vertex_t v0, v1;
	depth_t l0, l1;
        
    //traverse edges.
    //#pragma omp for nowait 
    for (uint64_t k = 0 ; k < cedge; ++k) {
        v0 = part_edge[k].v0;
        v1 = part_edge[k].v1;
		
		
		l0 = d0[v0];
		l1 = d1[v1];
	   	
        if ((l0 != INFTY) &&  ((l1 == INFTY) || (l1 - l0 > 1) )) {
            d1[v1] = l0 + 1;
            ++t_f_count;
		}
        #ifdef HALF_GRID
        else if ((l1 != INFTY) &&  ((l0 == INFTY) || (l0 - l1 > 1) )) {
            d0[v0] = l1 + 1;
            ++t_front_count;
		}
        #endif
    }

	if (t_f_count != 0) {
		g->read_part_next->set_bit_atomic(j);
		g->read_part->set_bit_atomic(j);
    }
    #ifdef HALF_GRID
    if (t_front_count != 0) {
		g->read_part_next->set_bit_atomic(i);
	}
    #endif
    return t_front_count + t_f_count;
}

int bfs2_t::iteration_finalize()
{
	index_t total_count = 0;
    for(int ithd = 0; ithd < NUM_THDS; ++ithd) {
        total_count += front_count[ithd];
        front_count[ithd] = 0;
    }
    cout << "level: " << (int)level;
    cout << " frontier count: " << total_count <<endl;
        
    if (total_count == 0) {
		/*	
        sort(depth, depth+g->vert_count);
		vertex_t count = 0;
		int lev = 0;
		for (vertex_t i = 0; i < g->vert_count; ++i) {
			if (lev == depth[i]) {
				++count;
			} else {
				cout << "level = " << lev << " count = " << count << endl;
				lev = depth[i];
				count = 1;
			}
		}
		cout << "level = " << lev << " count = " << count << endl;
		*/	
        return 0;
	}
    ++level;
    //one level done
	swap(g->read_part_next, g->read_part);
	g->read_part_next->reset();

	#ifdef HALF_GRID
	if (total_count <= 16*p) {
		ioskip = 0;
		//cout << "io_skip = " << (int)io_skip << endl;
	}
	#else
	if (total_count <= 4*p_s) {
		ioskip = 0;
		//cout << "io_skip = " << (int)io_skip << endl;
	}
	#endif

    return total_count;
}

//The BFS for this memory partition.
void bfs2_t::algo_mem_part(segment* seg)
{
    //double start = mywtime();
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
			
			if (g->read_part->get_bit(i1+b_i) || ioskip) {	
				for (part_t j1 = j2; j1 <= j_end1; ++j1) {
					edge_t* part_edge = (edge_t*)(edges + 
						(start_edge->get(i1, j1) << bytes_in_edge_shift));
					index_t cedge = start_edge->get_count(i1,j1);
					t_front_count += bfs_onepart(part_edge, cedge,
								   i1 + b_i, j1 + b_j);
				}   
			}
			#ifdef HALF_GRID	
			else {
				for (part_t j1 = j2; j1 <= j_end1; ++j1) {
					if (g->read_part->get_bit(j1+b_j)|| ioskip) {	
						edge_t* part_edge = (edge_t*)(edges + 
							(start_edge->get(i1, j1) << bytes_in_edge_shift));
						index_t cedge = start_edge->get_count(i1,j1);
						t_front_count += bfs_onepart(part_edge, cedge,
									   i1 + b_i, j1 + b_j);
					}
				}
			}
			#endif
	    }
	}
    
    front_count[omp_get_thread_num()] += t_front_count;
    //double end = mywtime();
    //cout << "mem_iteration time = " << end - start << endl;
    return;
}

void bfs2_t::algo_mem_part2(segment* seg)
{
    //double start = mywtime();
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
    
    #pragma omp for schedule (dynamic, 1) nowait 
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

		for (part_t i1 = i; i1 <= i_end; ++i1) {
            j2 = 0;
		    if (i1 == i_end) j_end1 = j_end; 
            if (i1 == i) j2 = j;
			#ifdef HALF_GRID	
            else if(b_i == b_j)  j2 = i1; 
			#endif
			
			if (g->read_part_next->get_bit(i1 + b_i)) {	
				for (part_t j1 = j2; j1 <= j_end1; ++j1) {
					edge_t* part_edge = (edge_t*)(edges + 
						(start_edge->get(i1, j1) << bytes_in_edge_shift));
					index_t cedge = start_edge->get_count(i1,j1);
					t_front_count += bfs_onepart(part_edge, cedge,
								   i1 + b_i, j1 + b_j);
				}   
			}
			#ifdef HALF_GRID	
			else {
				for (part_t j1 = j2; j1 <= j_end1; ++j1) {
					if (g->read_part_next->get_bit(j1 + b_j)) {	
						edge_t* part_edge = (edge_t*)(edges + 
							(start_edge->get(i1, j1) << bytes_in_edge_shift));
						index_t cedge = start_edge->get_count(i1,j1);
						t_front_count += bfs_onepart(part_edge, cedge,
									   i1 + b_i, j1 + b_j);
					}
				}   
			}
			#endif
	    }
	}
    
    front_count[omp_get_thread_num()] += t_front_count;
    //double end = mywtime();
    //cout << "mem_iteration time = " << end - start << endl;
    return;
}
