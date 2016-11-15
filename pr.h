
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



#ifndef _PR_H_
#define _PR_H_
#include "gstore.h"

typedef float rank_t;

class pr_t {
public:
    rank_t tmax[NUM_THDS];
    rank_t* vert_rank;
    rank_t* vert_rank_prior;
    int iteration_count;
    vertex_t vert_count;
public:
    void init(vertex_t vert_count, int a_iteration_count);
    void pagerank_onepart(edge_t* part_edge, index_t cedge, part_t i, part_t j);
    void iteration_finalize(int current_iteration);
    void algo_mem_part(segment* seg);
};
#endif
