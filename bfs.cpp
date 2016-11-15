
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
#include "bfs.h"

extern grid* g;

const depth_t INFTY = 0;

void bfs_t::init(vertex_t vert_count, vertex_t a_root, bitmap_t* a_read_part)
{
    read_part = a_read_part;
    root = a_root;
    depth = (depth_t*)calloc(sizeof(depth_t), vert_count);
    memset(depth, INFTY, sizeof(depth_t)*vert_count);

    depth[root] = 1;
    level = 1;
    cout << "Root = " << root << endl; 
    //#pragma omp parallel 
    calc_part_needed();
}

index_t bfs_t::bfs_onepart_row(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    index_t t_front_count = 0; 
    #ifdef HALF_GRID
    bool next_j = false; 
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    depth_t* d0 = depth + offset0;
    depth_t* d1 = depth + offset1;
    #else
    depth_t* d0 = depth;
    depth_t* d1 = depth;
    #endif
        
    //traverse edges.
    for (uint64_t k = 0 ; k < cedge; ++k) {
        vertex_t v0 = part_edge[k].v0;
        vertex_t v1 = part_edge[k].v1;

        if ((d0[v0] == level) && (d1[v1] == INFTY))  {
            d1[v1] = level + 1;
            ++t_front_count;
            next_j = true;
            //cout << "front " << v0 << endl;
        }
    }
    if (next_j) g->read_part_next->set_bit_atomic(j);
    #endif
    return t_front_count;
}

index_t bfs_t::bfs_onepart_col(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    index_t t_front_count = 0; 
    bool next_i = false; 
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    depth_t* d0 = depth + offset0;
    depth_t* d1 = depth + offset1;
    #else
    depth_t* d0 = depth;
    depth_t* d1 = depth;
    #endif
    vertex_t v0;
    vertex_t v1;    
    //traverse edges.
    for (uint64_t k = 0 ; k < cedge; ++k) {
        v0 = part_edge[k].v0;
        v1 = part_edge[k].v1;
        if ((d1[v1] == level) && (d0[v0] == INFTY)) { //|| (255 == depth[v0]))) {
            d0[v0] = level + 1;
            ++t_front_count;
            next_i = true;
        }
    }
    if (next_i) g->read_part_next->set_bit_atomic(i);
    return t_front_count;
}

index_t bfs_t::bfs_onepart(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    index_t t_front_count = 0;
    bool next_j = false;  
    
	#ifdef HALF_GRID
	bool next_i = false;
	#endif
    #ifdef COMPACT_GRID
    vertex_t offset0 = ((index_t)i << bit_shift2);
    vertex_t offset1 = ((index_t)j << bit_shift2);
    depth_t* d0 = depth + offset0;
    depth_t* d1 = depth + offset1;
    #else
    depth_t* d0 = depth;
    depth_t* d1 = depth;
    #endif
        
    //traverse edges.
    //#pragma omp for nowait 
    for (uint64_t k = 0 ; k < cedge; ++k) {
        vertex_t v0 = part_edge[k].v0;
        vertex_t v1 = part_edge[k].v1;
    
        if ((d0[v0] == level) && (d1[v1] == INFTY)) {// || (255 == depth[v1]))) {
            d1[v1] = level + 1;
            ++t_front_count;
            next_j = true;
            //cout << "front " << v1 << endl;
        } 
        #ifdef HALF_GRID
        else if ((d1[v1] == level) && (d0[v0] == INFTY)) {
                // || (255 == depth[v0]))) {
            d0[v0] = level + 1;
            ++t_front_count;
            next_i = true;
            //cout << "front " << v0 << endl;
        }
        #endif
    }
    
    #ifdef HALF_GRID
    if(next_i) g->read_part_next->set_bit_atomic(i);
    #endif
    if(next_j) g->read_part_next->set_bit_atomic(j);
    
    return t_front_count;
}

int bfs_t::iteration_finalize()
{
    for(int ithd = 1; ithd < NUM_THDS; ++ithd) {
        front_count[0] += front_count[ithd];
        front_count[ithd] = 0;
    }
    cout << "level: " << (int)level;
    cout << " frontier count: " << front_count[0] <<endl;
        
    if (front_count[0] == 0) {
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
	
		
		return 1;
	}
    ++level;
    front_count[0] = 0;
    //one level done
    swap(g->read_part, g->read_part_next);
    read_part = g->read_part;
    g->read_part_next->reset();
    //calc_part_needed();
    return 0;
}

//The BFS for this memory partition.
void bfs_t::algo_mem_part(segment* seg)
{
    //double start = mywtime();
    index_t t_front_count = 0; 
    
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
			
			if (read_part->get_bit(i1 + b_i)) {
				//#pragma omp for schedule (dynamic, 8) nowait 
				for (part_t j1 = j2; j1 <= j_end1; ++j1) {
				edge_t* part_edge = (edge_t*)(edges + 
					(start_edge->get(i1, j1) << bytes_in_edge_shift));
				index_t cedge = start_edge->get_count(i1,j1);
				t_front_count += bfs_onepart(part_edge, cedge,
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
					t_front_count += bfs_onepart_col(part_edge,
							cedge, i1 + b_i, j1 + b_j );
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


void algo_t::calc_part_needed()
{
    g->read_part->reset();
    #pragma omp parallel for schedule (dynamic, 64)
    for (uint32_t i = 0; i < p_s; ++i) {
        vertex_t v1 = (i << bit_shift2);
        vertex_t v2 = ((i+1) << bit_shift2);
        for (; v1 < v2; ++v1) {
            if (vertex_active(v1)) {
                g->read_part->set_bit_atomic(i);
                break;
            }
        }
    }
}

#ifdef VERIFY_BFS
index_t bfs_t::verify()
{
    bitmap_t* seen_edge = new bitmap_t(g->vert_count);
    #pragma omp parallel num_threads(NUM_THDS)
    {
	#pragma omp for schedule (dynamic, 1)// nowait//faster 
	for (part_t i = 0; i < p; ++i) {
	    #ifdef HALF_GRID 
	    for (part_t j = i; j < p; ++j)
	    #else 
	    for (part_t j = 0; j < p; ++j)
	    #endif
	    {
		#ifdef COMPACT_GRID
		vertex_t offset0 = (i << bits_in_block_shift);
		vertex_t offset1 = (j << bits_in_block_shift);
		#endif
		edge_t* part_edge = g->_edges + g->edge_start.get(i, j);
		index_t cedge = g->edge_start.get_count(i, j);
		for (uint64_t k = 0 ; k < cedge; ++k) {
		    #ifdef COMPACT_GRID
		    vertex_t v0 = offset0 + part_edge[k].v0;
		    vertex_t v1 = offset1 + part_edge[k].v1;
		    #else 
		    vertex_t v0 = part_edge[k].v0;
		    vertex_t v1 = part_edge[k].v1;
		    #endif
		    if (v0 != v1) {
			if (parent[v0] == v1) {
			    seen_edge->set_bit_atomic(v0);
			}
			if (parent[v1] == v0) {
			    seen_edge->set_bit_atomic(v1);
			}
			depth_t lvldiff = depth[v0] - depth[v1];
			if (!((lvldiff == 1) || (lvldiff == 255) || (lvldiff == 0))) {
			    cout << (index_t)lvldiff << endl;
			    cout << "bfs verification failed" << endl;
			    exit(-1);
			}
		    }
		}
	    }
	}
	#pragma omp barrier
	if (omp_get_thread_num() == 0) {
	    cout << "verifying next phase" << endl;
	}
	vertex_t vert_count = g->vert_count;
	#pragma omp for
	for (vertex_t k = 0; k < vert_count; ++k) {
	    if (k != root) {
	    	if (parent[k] != vertex_t(-1) && !seen_edge->get_bit(k)) {
		    cout << "bfs verification failed" << endl;
		    exit(-1);
		}
	    }
	}
	
    }
    return 0;
}
#endif
