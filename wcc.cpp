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
#include "wcc.h"


cid_t invalid_cid = -1;
extern grid* g;

void wcc2_t::init(vertex_t a_vert_count)
{
    vert_count = a_vert_count;
    memset(front_count, 0, sizeof(vertex_t)* NUM_THDS);
    
    vert_cid = (cid_t*)calloc(sizeof(cid_t), vert_count);
    memset(vert_cid, invalid_cid, sizeof(cid_t)*vert_count);

    cid = (cid_t*)calloc(sizeof(cid_t), vert_count);
    memset(cid, invalid_cid, sizeof(cid_t)*vert_count);
    
    wcc_group = 0;
	//cout << invalid_cid << endl;
    iteration = 0;
}

void wcc2_t::algo_mem_part(segment* seg)
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
		
            if (iteration == 0) {
                for (part_t j1 = j2; j1 <= j_end1; ++j1) {
                    edge_t* part_edge = (edge_t*)(edges + 
                        (start_edge->get(i1, j1) << bytes_in_edge_shift));
                    index_t cedge = start_edge->get_count(i1,j1);
                    t_front_count += wcc_onepart(part_edge, cedge,
                                   i1 + b_i, j1 + b_j);
                       
                }
            } else {
                for (part_t j1 = j2; j1 <= j_end1; ++j1) {
                    edge_t* part_edge = (edge_t*)(edges + 
                        (start_edge->get(i1, j1) << bytes_in_edge_shift));
                    index_t cedge = start_edge->get_count(i1,j1);
                    t_front_count += wcc_onepart2(part_edge, cedge,
                                   i1 + b_i, j1 + b_j);
                       
                }
            }
	    }
	
    }
    front_count[omp_get_thread_num()] += t_front_count;
    //double end = mywtime();
    //cout << "mem_iteration time = " << end - start << endl;
    return;
}
index_t wcc2_t::wcc_onepart(edge_t* part_edge, index_t cedge, part_t i, part_t j)
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
    cid_t c0, c1, m;    
    int sw;

    for (uint64_t k = 0 ; k < cedge; ++k) {
        v0 = part_edge[k].v0;
        v1 = part_edge[k].v1;
  
		c0 = vert_cid0[v0];		
		c1 = vert_cid1[v1];
        sw = (c0 != invalid_cid) + ((c1 != invalid_cid) << 1);
        switch(sw) {
        case 0: //if ((c0 == invalid_cid) && (c1 == invalid_cid)) {
			m = __sync_fetch_and_add(&wcc_group, 1);
			vert_cid1[v1] = m;
			vert_cid0[v0] = m;
			map_cid(m, m);
			++t_front_count;
            break;
        case 1: //} else if(c1 == invalid_cid) {
			vert_cid1[v1] = c0;
			++t_front_count;
            break;
        case 2: // } else if (c0 == invalid_cid) {
			vert_cid0[v0] = c1;
			++t_front_count;
            break;
        case 3: // } else if (c0 != c1) {
			if (c0 < c1) {
                //if (cid[c1] == c1)
                vert_cid1[v1] = c0;
				map_cid(c1, c0);
			    ++t_front_count;
			} else if (c0 > c1) {
                //if (cid[c0] == c0)
                vert_cid0[v0] = c1;
				map_cid(c0, c1);
			    ++t_front_count;
			}
            break;
        default:
            assert(0);
		}
    }
    return t_front_count;
}

index_t wcc2_t::wcc_onepart2(edge_t* part_edge, index_t cedge, part_t i, part_t j)
{
    index_t t_front_count = 0;
    
    #ifdef COMPACT_GRID
    cid_t* vert_cid0 = vert_cid + ((index_t)i << bit_shift2);
    cid_t* vert_cid1 = vert_cid + ((index_t)j << bit_shift2);
    #else
    cid_t* vert_cid0 = vert_cid;
    cid_t* vert_cid1 = vert_cid;
    #endif

    cid_t c0, c1;    
    vertex_t v0, v1;
    for (uint64_t k = 0 ; k < cedge; ++k) {
  
        v0 = part_edge[k].v0;
        v1 = part_edge[k].v1;
		c0 = vert_cid0[part_edge[k].v0];		
		c1 = vert_cid1[part_edge[k].v1];
        
        if (c0 < c1) {
            //if (cid[c1] == c1)
            vert_cid1[v1] = c0;
            map_cid(c1, c0);
            ++t_front_count;
        } else if (c0 > c1) {
            //if (cid[c0] == c0)
            vert_cid0[v0] = c1;
            map_cid(c0, c1);
            ++t_front_count;
        }
    }
    return t_front_count;
}


int wcc2_t::iteration_finalize() 
{
    iteration = 1;
    
    for(int ithd = 1; ithd < NUM_THDS; ++ithd) {
        front_count[0] += front_count[ithd];
        front_count[ithd] = 0;
    }

    cout << "level: " 
         << " frontier count: " << front_count[0] <<endl
         << endl;
    
    /*
    if (front_count[0] == 0) {
        vector<cid_t> vv = sort_indexes(vert_cid);
        int empty = 0;
        int single = 0;
        int count = 1;
        int cont = 0;
        
        cid_t tt = vert_cid[vv[0]];
        
        for (vertex_t j = 1; j < vert_count; ++j) {
            cid_t i = vv[j];
            if (vert_cid[i] == invalid_cid) ++empty;
            else if (tt != vert_cid[i] ){
                ++count;
                if (cont == 1) {
                   //cout << tt << "= " << cont << " " << vv[j - 2] << " " 
                   //<< vv[j - 1] << endl;
                    ++single;
                }
                //else  //cout << tt << "= " << cont << endl;
                //if (cont == 0)  cout << tt << "== " << cont << " " << vv[j - 1] << endl;
                cont = 0;
                tt = vert_cid[i];
            } else {
                ++cont;
            }
        }
    }*/

    if (front_count[0] == 0) {
        //XXX Make it parallel
        std::sort(cid, cid + wcc_group);
        cid_t tt;
        
        cid_t count = 0;
        #pragma omp parallel for reduction(+:count) private(tt)
        for(cid_t i = 1; i < wcc_group; ++i) {
            tt = cid[i- 1];
            count += (tt != cid[i]);
        }
        ++count;
        cout << "WCC count = " << count << endl;
        cout << "wcc_group used <debug> =" << wcc_group << endl;
        return 1;
    }

    /*
    if (front_count[0] == 0) {
        std::sort(vert_cid, vert_cid + vert_count);
        int empty = 0;
        int single = 0;
        int count = 1;
        int cont = 0;
        int max_cont = 0;
        cid_t tt = vert_cid[0];

        for (vertex_t i = 1; i < vert_count; ++i) {
            if (vert_cid[i] == invalid_cid) ++empty;
            else if (tt != vert_cid[i] ){
                ++count;
                single += (cont == 1);
                max_cont = max(cont, max_cont);
                cont = 0;
                tt = vert_cid[i];
            } else {
                ++cont;
            }
        }
        
        cout << "Empty vertices = " << empty << endl;
        cout << "Single vertices = " << single << endl;
        cout << "Max size Component vertices = " << max_cont << endl;
        cout << "WCC count = " << count << endl;
        cout << "wcc_group used <debug> =" << wcc_group << endl;
		return 1;
	}
    */
	#pragma omp parallel for //num_threads (NUM_THDS) 
    for(cid_t i = 0; i < wcc_group; i++) {
		cid_t n = i;
		cid_t m = cid[n];
		
		while (m < n) {
			n = m;
			m = cid[n];
		}
		cid[i] = m;
    }

	#pragma omp parallel for //num_threads (NUM_THDS) 
	for (vertex_t i = 0; i < vert_count; ++i) {
        if (vert_cid[i] == invalid_cid) continue;
		if(cid[vert_cid[i]] < vert_cid[i]) {
			vert_cid[i] = cid[vert_cid[i]];
		}
	}

    front_count[0] = 0;
    //calc_part_needed();
   return 0; 
}
    
