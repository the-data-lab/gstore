
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


#ifndef _WCC2_H_
#define _WCC2_H_

#include "algo.h"
extern grid* g;

typedef uint32_t cid_t;

/*
#include <vector>
using namespace std;
template<typename T> class CompareIndicesByAnotherVectorValues { 
    std::vector<T>* _values; 
    
    public: CompareIndicesByAnotherVectorValues(std::vector<T>* values) 
            : _values(values) {} 
    
    public: bool operator() 
            (const int& a, const int& b) const 
            { 
                return (_values)[a] > (_values)[b]; 
            } 
}; 

inline
vector<cid_t> sort_indexes(const cid_t* v) 
{
      // initialize original index locations
      vector<cid_t> idx(g->vert_count);
      for (size_t i = 0; i != g->vert_count; ++i) idx[i] = i;
      
      // sort indexes based on comparing values in v
      sort(idx.begin(), idx.end(),
              [&v](cid_t i1, cid_t i2) {return v[i1] < v[i2];});
       
      return idx;
}
*/
class wcc2_t : public algo_t{
public:
    cid_t		wcc_group;
    vertex_t    vert_count;
    cid_t*		cid;
	cid_t*		vert_cid;
	//bitmap_t	read_part;
    vertex_t front_count[NUM_THDS];
    int iteration ;
	
public:
    void init(vertex_t vert_count);
    int iteration_finalize();
    index_t wcc_onepart(edge_t* part_edge, index_t cedge, part_t i, part_t j);
    index_t wcc_onepart2(edge_t* part_edge, index_t cedge, part_t i, part_t j);
    void algo_mem_part(segment* seg);

	inline cid_t get_cid(vertex_t v) 
	{
		return vert_cid[v];
	}

	inline void set_cid(vertex_t v, cid_t c) 
	{
		vert_cid[v] = c;
	}

	inline void map_cid(cid_t c1, cid_t c)
	{
		cid[c1] = c;
	}
};

#endif
