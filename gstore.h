
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


#ifndef __GRID_H__
#define __GRID_H__
#include <stdint.h>
#include <string.h>
#include "list.h"
#include "bitmap.h"
#include <libaio.h>

using namespace std;

#define handle_error(msg) \
                   do { perror(msg); exit(EXIT_FAILURE); } while (0)

typedef uint64_t index_t;
typedef uint32_t part_t;
typedef uint32_t spart_t;
typedef uint8_t depth_t;
typedef Bitmap bitmap_t;
typedef int32_t degree_t;
typedef int32_t bdegree_t;
typedef int16_t sdegree_t;

//typedef uint64_t vertex_t;
//typedef uint64_t word_t;
//#define bytes_in_edge_shift 4 
typedef uint64_t vertex_t;
typedef uint16_t word_t;
typedef uint64_t part_id;


//tile count in 1 dimension
extern index_t p_s;

//physical group count in 1 dimension
extern index_t p;

#define COMPACT_GRID
//#define HALF_GRID

#define NUM_THDS 32 
#define bytes_in_edge_shift 2 


//For tile count in the graph
#define bit_shift2 16
#define part_mask2 (-1UL << bit_shift2)
#define part_mask2_2 ~(part_mask2)

//For conversion XXX
#define bit_shift0 25
#define part_mask0 (-1UL << bit_shift0)
#define part_mask0_2 ~(part_mask0)

//--------******* ----------//

//For tiles in a physical group
#define p_p 256 
#define bit_shift3 8 
#define part_mask3 (-1UL << bit_shift3)
#define part_mask3_2 ~(part_mask3)

//For physical group count in the graph
//bit_shift3+bit_shift3
#define bit_shift1 24 
#define part_mask1 (-1UL << bit_shift1)
#define part_mask1_2 ~(part_mask1) 


//2*bit_shift3
#define bit_shift4 16 
#define part_mask4 (-1UL << bit_shift4)
#define part_mask4_2 ~(part_mask4)


#ifdef HALF_GRID
#define myassert(i, j) assert( i <= j) 
inline index_t 
calc_index(index_t i, index_t j, index_t part_count) 
{
    return ((((i)*((part_count << 1) - i - 1)) >> 1) + j);
}

inline index_t
calc_index_opt(index_t i, index_t j, int bit_shift) 
{
    return ((((i)*((1 << (bit_shift + 1)) - i - 1)) >> 1) + j);
}

inline index_t 
calc_total_part(index_t i) 
{
    return (((i)*(i + 1)) >> 1) ;
}

//+ (part - i)*256*256 + i*257*128;
inline index_t 
beg_edge_offset1(index_t i) 
{
    return (((calc_index(i,i,p) - i) << bit_shift4) 
              + (((i)*(p_p + 1)) << (bit_shift3 - 1)));
}

//+ (part - i - 1)*256*256 + (i+1)*257*128;
inline index_t
beg_edge_offset2(index_t i, index_t j)
{
    return ((calc_index(i,j, p) - i - 1)<< bit_shift4) 
                + (((i+1)*(p_p + 1))<< (bit_shift3 - 1));
}


#else
#define myassert(i, j) 

inline index_t 
calc_index(index_t i, index_t j, index_t part_count) 
{
    return (i*part_count + j);
}

inline index_t
calc_index_opt(index_t i, index_t j, int bit_shift) 
{
	return  ((i << bit_shift)  + j);
}


inline index_t 
calc_total_part(index_t i) 
{
    return (i*i);
}

//+ (part - i)*256*256 + i*257*128;
inline index_t 
beg_edge_offset1(index_t i) 
{
    return (((calc_index(i,0,p)) << bit_shift4));
}

//+ (part - i - 1)*256*256 + (i+1)*257*128;
inline index_t
beg_edge_offset2(index_t i, index_t j)
{
    return ((calc_index(i,j, p))<< bit_shift4);
}

#endif

inline index_t
calc_index_f(index_t i, index_t j, index_t part_count) 
{
  return  ((i)*part_count  + j);
}

inline index_t
calc_index_f_opt(index_t i, index_t j, int bit_shift) 
{
	return  ((i << bit_shift)  + j);
}

#define ALIGN_MASK 0xFFFFFFFFFFFE00
#define UPPER_ALIGN(x) (((x) + 511) & ALIGN_MASK)
#define LOWER_ALIGN(x) ((x) & ALIGN_MASK)


inline
void get_s_ij(part_id start, spart_t& i, spart_t& j) 
{
	uint32_t y = (start & part_mask4_2);
	i = (y >> bit_shift3);
	j = (y & part_mask3_2);
	return;
}

inline
void get_ij(part_id start, part_t& b_i, part_t& b_j, spart_t& i, spart_t& j) 
{
	uint32_t y = (start & part_mask4_2);
	i = (y >> bit_shift3);
	j = (y & part_mask3_2);

	index_t x = (start >> bit_shift4);
	b_i = x/p;
	b_j = x%p;

	return;
}

class algo_t;

#ifdef GENEARATOR_USE_PACKED_TYPE
struct gedge_t {
private:
    uint32_t v0_low;
    uint32_t v1_low;
    uint32_t high; /* v1 in high half, v0 in low half */

public:    
    inline bool is_self_loop() {
        return ((v0_low == v1_low) && 
		(high & 0xFFFF) == (high >> 16));
    }
    
    bool operator == (gedge_t e) {
        return ((v0_low == e.v0_low) && 
		(v1_low == e.v1_low) &&
		(high & 0xFFFF) == (high >> 16));
    }

    inline vertex_t get_v0()
    {
	return (v0_low | ((int64_t)((int16_t)(high & 0xFFFF)) << 32));
    }
    inline vertex_t get_v1()
    {
	return (v1_low | ((int64_t)((int16_t)(high >> 16)) << 32));
    }
    inline void set_v0(vertex_t a_v0) { assert(0); }
    inline void set_v1(vertex_t a_v1) { assert(0); } 
};
#else 
struct gedge_t {
private:
    uint32_t v0;
    uint32_t v1;

public:
    inline bool is_self_loop() {
       return (v0 == v1);
    }
    
    bool operator == (gedge_t e) {
        return ((v0 == e.v0) && (v1 == e.v1));
    }

    inline vertex_t get_v0()
    {
	return v0;
    }
    inline vertex_t get_v1()
    {
	return v1;
    }
    inline void set_v0(vertex_t a_v0) { v0 = a_v0; }
    inline void set_v1(vertex_t a_v1) { v1 = a_v1; } 
};
#endif

struct cedge_t {
    uint32_t v0;
    uint32_t v1;
    
    inline cedge_t () 
    {
        v0 = 0;
        v1 = 0;
    }

    //We are discarding some upper bits knowingly.
    inline cedge_t(gedge_t e) 
    { 
        v0 = e.get_v0(); 
	v1 = e.get_v1(); 
    }

    inline cedge_t operator = (gedge_t e)
    {
       v0 = e.get_v0(); 
       v1 = e.get_v1();
       return *this; 
    }
    inline bool operator == (cedge_t e) {
        return ((v0 == e.v0) && (v1 == e.v1));
    }
};


#ifdef COMPACT_GRID
struct edge_t {
    word_t v0;
    word_t v1;
    inline edge_t operator = (gedge_t e)
    {
       v0 = e.get_v0(); 
       v1 = e.get_v1();
       return *this; 
    }
    inline edge_t operator = (cedge_t e)
    {
       v0 = e.v0; v1 = e.v1;
       return *this; 
    }
    inline edge_t () 
    {
        v0 = 0;
        v1 = 0;
    }
    inline bool operator == (edge_t e) {
        return ((v0 == e.v0) && (v1 == e.v1));
    }
};
#else
typedef struct cedge_t edge_t;
#endif

#ifdef HALF_GRID
template <class T, class V>
class matrix {
public:
    inline matrix() {}
    inline virtual ~matrix() {
        //free(val);
    }
    inline void free_mem() {
        free(val);
        val = 0;
    }
    inline void virtual init(T i) {
        #ifdef HALF_GRID
        index_t size = (i*(i + 1))/2;
        #else
        index_t size = i*i;
        #endif
        //+1 is for setinel node
        val = (V*)calloc(sizeof(V), size + 1);
        part_count = i;
    }
    inline virtual void init1(T i) {
        val = 0;
        part_count = i;
    }
    inline virtual V get(T i, T j) {
        myassert(i, j); 
        V index = calc_index(i,j, part_count); 
        return val[index];
    }
    inline virtual V get(V index) {
        return val[index];
    }
    inline virtual V get_index(T i, T j) {
		return calc_index(i,j, part_count);
    }	
    inline V get_index_opt(T i, T j, int bitshift) {
		return calc_index_opt(i,j, bitshift);
    }
    inline virtual V get_end(T i, T j) {
        V index = calc_index(i,j, part_count); 
        myassert(i, j); 
        return val[index + 1];
    }
    inline virtual V get_end(V i) {
        return val[i + 1];
    }
    inline virtual void set(T i, T j, V v) {
        V index = calc_index(i,j, part_count); 
        myassert(i, j); 
        val[index] = v;
    }
    inline virtual index_t get_count(T i, T j) {
        V index = calc_index(i,j,part_count); 
        myassert(i, j); 
        return val[index + 1] - val[index];
    }
    
    inline virtual void incr(T i, T j) {
        myassert(i,j); 
        V index = calc_index(i,j,part_count); 
        ++val[index];
    }

    inline virtual index_t edge_count() {
        #ifdef HALF_GRID
        index_t size = (part_count*(part_count + 1))/2;
        #else
        index_t size = part_count*part_count;
        #endif
        return val[size];
    }
    
    inline virtual V atomic_incr(T i, T j) {
        myassert(i,j); 
        V index = calc_index(i,j,part_count);
        V ret =__sync_fetch_and_add(val + index, 1); 
        return ret;
    }
    inline virtual V atomic_incr_opt(T i, T j, int bitshift) {
        myassert(i,j); 
        V index = calc_index_opt(i,j, bitshift);
        V ret =__sync_fetch_and_add(val + index, 1); 
        return ret;
    }
    
    inline virtual V atomic_incr(V index) {
        V ret =__sync_fetch_and_add(val + index, 1); 
        return ret;
    }

    /*
    inline uint64_t find(V offset, T&l, T&m) {
        uint64_t edge_offset = (offset >> 2);
        part_t i = 0;
        part_t j = 0;
        for(i = 0; i < p; ++i) {
            for(j = i; j < p; ++j) {
                if(edge_offset < val[index + 1]) {
                    if (j == i) {
                        m = p -1;
                        l = i -1;  
                    } else {
                        l = i;
                        m = j -1;
                    }
                    return (val[index] << 2);
                }
            }
        }
        l = p - 1;
        m = p - 1;
        return (val[size] << 2);
        //assert(0);
    }
    
    inline index_t get_row_count(T i) {
        T j = i;
        myassert(i, j); 
        return val[index + p - j] - val[index];
    }
    */

public:
    V*      val;
    T       part_count;
};
#endif

template <class T, class V>
class matrix_f 
#ifdef HALF_GRID
: public matrix<T,V>
#endif
{ 
public:
    
    inline void free_mem() {
        free(this->val);
        this->val = 0;
    }

    //XXX Full matrix version starts
    inline void init(T i) {
        this->part_count = i;
        index_t size = i*i;
        //+1 is for setinel node
        this->val = (V*)calloc(sizeof(V), size + 1);
    }
    inline void init1(T i) {
        this->part_count = i;
        //+1 is for setinel node
        this->val = 0;
    }
    inline V get(T i, T j) {
        V index = calc_index_f(i,j, this->part_count); 
        return this->val[index];
    }
    inline V get(V index) {
        return this->val[index];
    }
    inline V get_index(T i, T j) {
	    return calc_index_f(i,j, this->part_count);
    }	
    
    inline V get_index_opt(T i, T j, int bitshift) {
	    return calc_index_f_opt(i,j, bitshift);
    }	
    
    inline V get_end(T i, T j) {
        V index = calc_index_f(i,j, this->part_count); 
        return this->val[index + 1];
    }
    inline V get_end(V i) {
        return this->val[i + 1];
    }
    inline void set(T i, T j, V v) {
        V index = calc_index_f(i,j, this->part_count); 
        this->val[index] = v;
    }
    inline index_t get_count(T i, T j) {
        V index = calc_index_f(i,j, this->part_count); 
        return this->val[index + 1] - this->val[index];
    }
    

    inline void incr(T i, T j) {
        V index = calc_index_f(i,j, this->part_count); 
        ++this->val[index];
    }
    
    inline V atomic_incr(T i, T j) {
        V index = calc_index_f(i,j, this->part_count); 
        V ret =__sync_fetch_and_add(this->val + index, 1); 
        return ret;
    }
    inline V atomic_incr_opt(T i, T j, int  bitshift) {
        V index = calc_index_f_opt(i,j, bitshift); 
        V ret =__sync_fetch_and_add(this->val + index, 1); 
        return ret;
    }
    
    inline V atomic_incr(V index) {
        V ret =__sync_fetch_and_add(this->val + index, 1); 
        return ret;
    }

    inline index_t edge_count() {
        index_t size = this->part_count*this->part_count;
        return this->val[size];
    }
#ifndef HALF_GRID
public:
    V*      val;
    T part_count;
#endif
};

#ifndef HALF_GRID
#define matrix matrix_f
#endif


class part_meta_t {
public:
    size_t offset;
    part_id start;
    part_id end;
};

class last_read_t {
public:
    part_t b_i; //big i
    part_t b_j;
    part_t s_i; //small i 
    part_t s_j; 
};

class segment {
public:
    char* buf;
    index_t ctx_count;
    part_meta_t* meta;
};

 class wcc_t;
class grid {
public:
    grid();
    ~grid();
    
    void init(int argc, char* argv[]);

    //In-memory conversions
    void pre_grid(string edgefile, gedge_t* edges, index_t nedges);
    void proc_grid(string edgefile, string part_file);
    
    //Out-of-core conversion for converting big files/dir to G-Store format
    void proc_grid_big(string edgefile, string part_file, bool is_odir);
    
    //They generate intermediate files
    void pre_grid_file(string edgefile);
    void pre_grid_dir(string idir);

    //It reads the intermediate files to generate final one giant files
    void post_grid_file(string edgefile, string part_file, bool is_dir);
    //It reads the intermediate files to generate many small files
    void post_grid_dir(string edgefile, string part_file, bool is_dir);

    
	
	void pre_csr(string edgefile, gedge_t* edges, index_t nedges);
    void proc_csr(string edgefile, string part_file);
	
    void save_grid(string edgefile);
    void save_grid_big(string edgefile, bool is_odir);
    void save_start_files(string edgefile, bool is_odir);
    void save_degree_files(string edgefile);
    void compress_degree();

	void analyze_grid_size(string edgefile);
    
    void read_grid(string edgefile);
    void read_start(string edgefile);
    void read_degree(string edgefile);
    
    void read_grid_in_mem(string edgefile);
    void read_degree_in_mem(string edgefile);
    void read_start_in_mem(string edgefile);
    
    int  read_aio_init(string edgefile);
    
    //cache handling
    void swap_parts();
    void new_swap_parts1();
    void new_swap_parts2();
    void swap_help(segment* seg, segment* tmp);
    void cache_to_cache();

public:
    void do_algo(algo_t* algo);
    void do_algo_prop(algo_t& algo);
    void bfs();
    void bfs2();
    void bfs_mmap();
    void pagerank();
    void pagerank_mmap();
    void kcore(int kc);
    void kcore_mmap(int kc);
	void wcc2();
	void wcc2_mmap();
    void traverse();
    void wcc();
    void wcc_mmap();
    void wcc_init(wcc_t* algo);
	void tc();

public:


public:
    vertex_t vert_count;
	vertex_t bdegree_count;

    bitmap_t* read_part;
    bitmap_t* read_part_next;
     
    //2-D edge partitions
    edge_t*  _edges; 

    //Bigger partition
    //matrix<part_t, index_t> edge_start;
    
    //smaller partitions of whole graph
    index_t* _s_start_edge;
    
	degree_t* vert_degree;
	sdegree_t* svert_degree;
	bdegree_t* bvert_degree;


    //move it to cache-driver/io-driver
    //We will maintain two copies, for two memory partitions.
    segment* seg1;
    segment* seg2;

    segment* seg3;

    //Caches random algo
    segment* cached_pool1;
    segment* cached_pool2;

};


class part_cache_t {
public:
    index_t start;
    index_t end;

public:
    /*
    inline bool operator < (part_cache_t& obj)
    {
        return (start < obj.start);
    }*/
};

inline bool operator < (const part_cache_t & o1, const part_cache_t& o2)
{
    return (o1.start < o2.start);
}

typedef struct __aio_meta {
    struct io_event* events;
    struct iocb** cb_list;
    io_context_t  ctx;
    int busy;
} aio_meta_t;

class io_driver {
public:
    size_t read_aio_random(segment* seg);
    int wait_aio_completion();
    
private:
    aio_meta_t* aio_meta;

public:
    io_driver();
};

class cache_driver {
public:
    //What is already processed as part of cached from last iteration.
    part_cache_t* part_cached;
    index_t       ctx_count2;
    index_t       start_index;


public:
    inline cache_driver()
    {
        //What is already processed as part of cached from last iteration.
        part_cached = (part_cache_t*)calloc(sizeof(part_cache_t), 1024*1024*16);
        ctx_count2  = 0;
    }
    int 
    prep_read_aio_random(const last_read_t& last_read, 
		       size_t to_read, bitmap_t* read_part, 
		       part_meta_t* part_meta,
		       int& iteration_done, 
		       size_t& start_addr);
    
    int 
    prep_read_aio_col(const last_read_t& last_read, part_t i, bitmap_t* read_part, 
		    part_meta_t* part_meta, 
		    size_t to_read, matrix<spart_t, index_t>* start_edge,
		    size_t& start_addr,  int& k, bool& found_start);

    int prep_read_aio_seq(const last_read_t& last_read, size_t to_read, 
                  part_meta_t* part_meta,
		  int& iteration_done, size_t& start_addr);

    size_t fetch_mem_part(last_read_t& last_read, segment* seg,
			  bitmap_t* read_part, int iotype = 0);

    void prep_pread_cached(segment* seg1, index_t ctx_count); 
    index_t prep_pread_cached(segment* seg); 
                                        
    index_t search_cached_part(index_t search, bool& found);

};


#endif
