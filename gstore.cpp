
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
#include <libaio.h>
#include <math.h>
#include <asm/mman.h>
#include <dirent.h>
#include "wtime.h"
#include "gstore.h"
#include "pr.h"
#include "bfs.h"
#include "bfs2.h"
#include "kcore.h"
#include "wcc.h"
#include "traverse.h"
#define AIO_MAXIO 16384 
#define AIO_BATCHIO 256
//These are IO threads
#define IO_THDS 1

index_t p = 0;
index_t p_s = 0;
size_t mp_count = 1;
int f_edge = 0;
size_t sz_copied = 0;

grid* g;
cache_driver*   cache;
io_driver*      io;
	
index_t* s_count_edge_in;  


int arg = -1;
size_t total_size = (8192L << 20); // 2048 MB
char* main_buf = 0;
char* init_buf = 0;
size_t cache_size = 0;
//size_t memory = (1024L << 20);// 512 MB;
size_t memory = (256L << 20);// 512 MB;
size_t read_sz = (2L << 20);
uint8_t compressed = 0;
////////////////////////

void conv_to_text(string edgefile);
void text_to_bin(string textfile);
void remove_dup(string edgefile);
void end_ctx(part_meta_t* part_meta, int k, part_t b_i, part_t b_j, spart_t i, spart_t j);
bool end_cache_ctx(segment* seg, segment* tmp, size_t& start_addr, 
                   matrix<spart_t, index_t>* start_edge, char* buf, char* edges); 
void shallow_copy(segment* cached_pool1, segment* seg);

inline bool 
is_end(last_read_t& last_read) {
    return (last_read.b_i == p -1 && last_read.b_j == p - 1 
	    && last_read.s_i == p_p - 1 && last_read.s_j == p_p - 1);
}

inline bool 
is_start(last_read_t last_read) 
{
    return (last_read.b_i == 0 && 
            last_read.b_j == 0 && 
            last_read.s_i == 0 && 
            last_read.s_j);

}

grid::grid()
{
    _edges = 0; 
    vert_degree = 0;
	bvert_degree = 0;
	bdegree_count = 0;
    svert_degree = 0;
}

grid::~grid()
{
    free(vert_degree);
    vert_degree = 0;
	free(bvert_degree);
	bvert_degree = 0;
	free(svert_degree);
    svert_degree = 0;
}

void grid::init(int argc, char * argv[])
{
    int o;
    int job = 0;
    uint32_t scale;
    int c = 0;
    string edgefile;
    string part_file;
    while ((o = getopt (argc, argv, "s:o:hi:j:c:m:v:a:")) != -1) {
        switch(o) {
            case 's': //scale
                scale = atoi(optarg);
                vert_count = (1L << scale);
                scale -= bit_shift1;
                p = (1 << scale);
                p_s = (vert_count >> bit_shift2) + (0 != (vert_count & part_mask2_2));
                read_part = new bitmap_t(p_s);
                read_part_next = new bitmap_t(p_s);
                read_part_next->reset();
                break;
            case 'v'://vert count
                //vert_count = atoi(optarg);
				sscanf(optarg, "%ld", &vert_count);
                p = (vert_count >> bit_shift1) + (0 != (vert_count & part_mask1_2));
                //p_s = (vert_count >> bit_shift2) + (0 != (vert_count & part_mask2_2));
                p_s = p_p*p;
                read_part = new bitmap_t(p_s);
                read_part_next = new bitmap_t(p_s);
                read_part_next->reset();
                break;
            case 'i':
                edgefile = optarg;
                break;
            case 'o':
                part_file = optarg;
                break;
            case 'j':
                job = atoi(optarg);
                break;
            case 'h':
                cout << "Coming soon" << endl;
                return;
            case 'c':
                c = atoi(optarg);
                break;
            case 'm':
                total_size = ((size_t)atoi(optarg) << 20);
                break;
            case 'a':
                arg = atoi(optarg);
                break;
            default:
               cout << "Unknown argument" << endl ; 
        }
    }

    cout << "Total Physcial Group Count (in One Direction) = " << p << endl;
    cout << "Total Tiles Count (in One Direction) = " << p_s << endl;
    double start, end;
  
    switch(c) {
    case 1:
        start = mywtime();
		proc_grid(edgefile, part_file);
		end = mywtime();
        save_grid(part_file);
		cout << "GStore Time = " << end - start << endl;
        return ;
    case 2:
        proc_grid_big(edgefile, part_file, false);
        save_grid_big(part_file, false);
        return;
    case 3:
        part_file += "/gstore";
        proc_grid_big(edgefile, part_file, true);
        save_grid_big(part_file, true);
        return;
    case 4:
        text_to_bin(edgefile);
        return;
	case 5:
		start = mywtime();
		proc_csr(edgefile, part_file);
		end = mywtime();
		cout << "CSR Time = " << end - start << endl;
		return;
	case 6:
        start = mywtime();
		proc_grid(edgefile, part_file);
		end = mywtime();
		cout << "GStore Time = " << end - start << endl;
        //break;
        return ;

    case 7:
        remove_dup(edgefile);
        return;
	case 8:
		conv_to_text(edgefile);
		return;
    } 

   
    switch(job) {
    case 0:
            start = mywtime();
            read_aio_init(edgefile);
            bfs();
            end = mywtime();
            cout << "Half BFS time = " << end-start << endl;
            break;    
    case 10:
            start = mywtime();
            read_grid(edgefile);
            bfs_mmap();
            end = mywtime();
            cout << "Half BFS time = " << end-start << endl;
            break;    
    case 20:
            read_grid_in_mem(edgefile);
            start = mywtime();
            bfs_mmap();
            end = mywtime();
            cout << "Half BFS time = " << end-start << endl;
            break;    
    case 30:
            start = mywtime();
            read_aio_init(edgefile);
            bfs2();
            end = mywtime();
            cout << "Half BFS time = " << end-start << endl;
            break;    
    case 1:
            start = mywtime();
            read_aio_init(edgefile);
            read_degree_in_mem(edgefile);
            pagerank();
            end = mywtime();
            cout << "Page-rank time = " << end-start << endl;
            break;
    case 11:
            start = mywtime();
            read_grid(edgefile);
            read_degree_in_mem(edgefile);
            pagerank_mmap();
            end = mywtime();
            cout << "Page-rank_mmap time = " << end-start << endl;
            break;
    case 21:
            read_grid_in_mem(edgefile);
            read_degree_in_mem(edgefile);
            start = mywtime();
            pagerank_mmap();
            end = mywtime();
            cout << "Page-rank_mmap time = " << end-start << endl;
            break;

    case 2:
            start = mywtime();
            read_aio_init(edgefile);
            read_degree_in_mem(edgefile);
            kcore(arg);
            end = mywtime();
            cout << "Kcore time = " << end-start << endl;
            break;
    case 12:
            start = mywtime();
            read_grid(edgefile);
            read_degree_in_mem(edgefile);
            kcore_mmap(arg);
            end = mywtime();
            cout << "Kcore time = " << end-start << endl;
            break;
    case 22:
            read_grid_in_mem(edgefile);
            read_degree_in_mem(edgefile);
            start = mywtime();
            kcore_mmap(arg);
            end = mywtime();
            cout << "Kcore time = " << end-start << endl;
            break;
    case 3:
            start = mywtime();
            read_aio_init(edgefile);
            wcc();
            end = mywtime();
            cout << "WCC time = " << end-start << endl;
            break;
    case 13:
            start = mywtime();
            read_grid(edgefile);
            wcc_mmap();
            end = mywtime();
            cout << "WCC time = " << end-start << endl;
            break;
    case 23:
            read_grid_in_mem(edgefile);
            start = mywtime();
            wcc_mmap();
            end = mywtime();
            cout << "WCC time = " << end-start << endl;
            break;
    case 4:
            start = mywtime();
            read_aio_init(edgefile);
            read_degree_in_mem(edgefile);
            traverse();
            end = mywtime();
            cout << "Single color propagation time = " << end-start << endl;
            break;
	case 100:
			analyze_grid_size(edgefile);
			return;
    default:
            cout << "Wrong value for -j argument" << endl;
    }
}


void grid::pre_grid_dir(string idir)
{
    string edgefile;
    //Number of files in one row.
    index_t p_v = (vert_count >> bit_shift0);

    //Total number of files.
    index_t vcount = calc_total_part(p_v);
    cout << "Creating " << vcount << " intermediate file, p_v =" << p_v << endl;
    
    string* tfile = new string [vcount];
    FILE** tf = (FILE**) calloc(sizeof(FILE*), vcount); 
    cedge_t** buf = (cedge_t**)calloc(sizeof(cedge_t), vcount); 
    int64_t count = (1L<<30); 
    gedge_t* edges = (gedge_t*)malloc(sizeof(gedge_t)*count);
    gedge_t* edges2 = (gedge_t*)malloc(sizeof(gedge_t)*count);
    
    for( index_t ifile = 0; ifile < vcount; ++ifile) {
        char tmp[64] = {0};
        sprintf(tmp, ".%ld", ifile);
        tfile[ifile] = idir + "/.tmp" + tmp;
        tf[ifile] = fopen(tfile[ifile].c_str(), "wb");
        assert(tf[ifile] != NULL);
        buf[ifile] = (cedge_t*)malloc(sizeof(cedge_t)*count);
        cout << tfile[ifile] << endl;
    }


    //files
    #ifdef HALF_GRID
    matrix<spart_t, index_t> count_edge;
    count_edge.init(p_v);
    #else
    matrix_f<spart_t, index_t> count_edge;
    count_edge.init(p_v);
    #endif
    
    //tiles
    index_t total_s_part = calc_total_part(p_s);
    _s_start_edge = (index_t*) calloc(sizeof(index_t), total_s_part);  
    vert_degree = (degree_t*)calloc(sizeof(degree_t), vert_count);
    

    //Input Directory read
    struct dirent *ptr;
    DIR* dir;


    dir = opendir(idir.c_str());
    while (NULL != (ptr= readdir(dir))) {
        if (ptr->d_name[0] == '.') continue;
        edgefile = idir + "/" + string(ptr->d_name);
        cout << "Reading " << edgefile << endl;

        //read the file and convert it to first group
        struct stat st_edge;
        FILE* fid_edge = fopen(edgefile.c_str(), "rb");
        stat(edgefile.c_str(), &st_edge);
        assert(st_edge.st_size != 0);
        //index_t nedges = st_edge.st_size/sizeof(gedge_t);
        size_t to_read = st_edge.st_size/sizeof(gedge_t);
        size_t sz_read = 0;
        index_t ecount = 0;
        sz_read = fread(edges2, sizeof(gedge_t), count, fid_edge);
        assert(sz_read != 0);

        while(to_read != 0) {
            to_read -= sz_read; 
            ecount = sz_read; 
            swap(edges, edges2);
            #pragma omp parallel num_threads(NUM_THDS) shared(edges,edges2) 
            {
                if (0 == omp_get_thread_num() && (to_read != 0)) {
                    sz_read = fread(edges2, sizeof(gedge_t), count, fid_edge);
                    if (sz_read == 0) {
                        int err = ferror(fid_edge);
                        cout << err << endl;
                        exit(-1);
                    }
                    assert(sz_read != 0);
                }

                #ifdef HALF_GRID
                matrix<spart_t, index_t> start_edge_half;
                start_edge_half.part_count = p_p;
                #endif

                matrix_f<spart_t, index_t> start_edge_full;
                start_edge_full.part_count = p_p;
                matrix<spart_t, index_t>* s_start_edge;

                gedge_t  edge;
                part_t j, i;
                part_t b_i, b_j;
                spart_t s_i, s_j;
                vertex_t v0, v1, v2, v3;
                index_t offset;
                #pragma omp for schedule (dynamic, 4096) 
                for (index_t k = 0; k < ecount; ++k) {
                    edge = edges[k];
                    if (edge.is_self_loop()) continue;
                    v2 = edge.get_v0();
                    v3 = edge.get_v1();
                    
                    #ifdef HALF_GRID 
                    v0 = min(v2, v3);
                    v1 = max(v3, v2);
                    #else
                    v0 = v2;
                    v1 = v3;
                    #endif

                    i = (v0 >> bit_shift0);
                    j = (v1 >> bit_shift0);
                
                    index_t n = count_edge.get_index(i,j); 
                    index_t m = count_edge.atomic_incr(n);

                    buf[n][m].v0 = (v0 & part_mask0_2); 
                    buf[n][m].v1 = (v1 & part_mask0_2); 
                    
                    b_i = (v0 >> bit_shift1);
                    b_j = (v1 >> bit_shift1);
                
                    s_i = ((v0 >> bit_shift2) & part_mask3_2);
                    s_j = ((v1 >> bit_shift2) & part_mask3_2);
                    
                    #ifdef HALF_GRID
                    if (b_i == b_j) {
                        s_start_edge = &start_edge_half;
                        offset = beg_edge_offset1(b_i);
                    } else {
                        s_start_edge = &start_edge_full;
                        offset = beg_edge_offset2(b_i, b_j);
                    }
                    #else	
                    s_start_edge = &start_edge_full;
                    offset = beg_edge_offset2(b_i, b_j);
                    #endif
                    
                    s_start_edge->val = _s_start_edge + offset;
                    s_start_edge->atomic_incr(s_i, s_j);

                    __sync_fetch_and_add(vert_degree + edge.get_v0(), 1);
                    #ifdef HALF_GRID
                    __sync_fetch_and_add(vert_degree + edge.get_v1(), 1);
                    #endif
                }
            } 
            for (index_t ifile = 0; ifile < vcount; ++ifile) {
                fwrite(buf[ifile], sizeof(cedge_t), count_edge.val[ifile], tf[ifile]);
                count_edge.val[ifile] = 0;
            }
        }    
        fclose(fid_edge);
    }
    closedir(dir);

    //Read the file and convert it to first group.
    
    //Calculate the edge start
    index_t prefix_sum = 0;
    index_t curr_value = 0;
    for (index_t ipart = 0; ipart < total_s_part; ++ipart) {
        curr_value = _s_start_edge[ipart];
        _s_start_edge[ipart] = prefix_sum;
        prefix_sum += curr_value;
    }
    _s_start_edge[total_s_part] = prefix_sum;
    cout << "Total edges = " << prefix_sum << endl;
    
    for (index_t ifile = 0; ifile < vcount; ++ifile) {
        fclose(tf[ifile]);
    }
}

void grid::pre_grid_file(string edgefile)
{
    //Number of files in one row.
    index_t p_v = (vert_count >> bit_shift0);

    //Total number of files.
    index_t vcount = calc_total_part(p_v);
    cout << "Creating " << vcount << " intermediate files, p_v =" << p_v << endl;
    
    string* tfile = new string [vcount];
    FILE** tf = (FILE**) calloc(sizeof(FILE*), vcount); 
    cedge_t** buf = (cedge_t**)calloc(sizeof(cedge_t), vcount); 
    int64_t count = (1L<<30); 
    gedge_t* edges = (gedge_t*)malloc(sizeof(gedge_t)*count);
    gedge_t* edges2 = (gedge_t*)malloc(sizeof(gedge_t)*count);
    
    for( index_t ifile = 0; ifile < vcount; ++ifile) {
        char tmp[64] = {0};
        sprintf(tmp, ".%ld", ifile);
        tfile[ifile] = edgefile + tmp;
        tf[ifile] = fopen(tfile[ifile].c_str(), "wb");
        assert(tf[ifile] != NULL);
        buf[ifile] = (cedge_t*)malloc(sizeof(cedge_t)*count);
        cout << tfile[ifile] << endl;
    }


    //files
    #ifdef HALF_GRID
    matrix<spart_t, index_t> count_edge;
    count_edge.init(p_v);
    #else
    matrix_f<spart_t, index_t> count_edge;
    count_edge.init(p_v);
    #endif
    
    //Smaller partitions
    index_t total_s_part = calc_total_part(p_s);
    _s_start_edge = (index_t*) calloc(sizeof(index_t), total_s_part);  
    vert_degree = (degree_t*)calloc(sizeof(degree_t), vert_count);
    
    
    //Read the file and convert it to first group.
    struct stat st_edge;
    FILE* fid_edge = fopen(edgefile.c_str(), "rb");
    stat(edgefile.c_str(), &st_edge);
    assert(st_edge.st_size != 0);
    //index_t nedges = st_edge.st_size/sizeof(gedge_t);
    size_t to_read = st_edge.st_size/sizeof(gedge_t);
    size_t sz_read = 0;
    index_t ecount = 0;
    cout << "Reading " << edgefile << endl;
    sz_read = fread(edges2, sizeof(gedge_t), count, fid_edge);
    assert(sz_read != 0);

    while(to_read != 0) {
        to_read -= sz_read; 
        ecount = sz_read; 
        swap(edges, edges2);
        #pragma omp parallel num_threads(NUM_THDS) shared(edges,edges2) 
        {
            if (0 == omp_get_thread_num() && (to_read != 0)) {
                sz_read = fread(edges2, sizeof(gedge_t), count, fid_edge);
                if (sz_read == 0) {
                    int err = ferror(fid_edge);
                    cout << err << endl;
                    exit(-1);
                }
                assert(sz_read != 0);
            }

			#ifdef HALF_GRID
            matrix<spart_t, index_t> start_edge_half;
            start_edge_half.part_count = p_p;
			#endif

            matrix_f<spart_t, index_t> start_edge_full;
            start_edge_full.part_count = p_p;
            matrix<spart_t, index_t>* s_start_edge;

            gedge_t  edge;
            part_t j, i;
            part_t b_i, b_j;
            spart_t s_i, s_j;
            vertex_t v0, v1, v2, v3;
            index_t offset;
            #pragma omp for schedule (dynamic, 4096) 
            for (index_t k = 0; k < ecount; ++k) {
                edge = edges[k];
                if (edge.is_self_loop()) continue;
                v2 = edge.get_v0();
                v3 = edge.get_v1();
                
                #ifdef HALF_GRID 
                v0 = min(v2, v3);
                v1 = max(v3, v2);
                #else
                v0 = v2;
                v1 = v3;
                #endif

                i = (v0 >> bit_shift0);
                j = (v1 >> bit_shift0);
            
                index_t n = count_edge.get_index(i,j); 
                index_t m = count_edge.atomic_incr(n);

                buf[n][m].v0 = (v0 & part_mask0_2); 
                buf[n][m].v1 = (v1 & part_mask0_2); 
                
                b_i = (v0 >> bit_shift1);
                b_j = (v1 >> bit_shift1);
            
                s_i = ((v0 >> bit_shift2) & part_mask3_2);
                s_j = ((v1 >> bit_shift2) & part_mask3_2);
                
                #ifdef HALF_GRID
                if (b_i == b_j) {
                    s_start_edge = &start_edge_half;
                    offset = beg_edge_offset1(b_i);
                } else {
                    s_start_edge = &start_edge_full;
                    offset = beg_edge_offset2(b_i, b_j);
                }
                #else	
                s_start_edge = &start_edge_full;
                offset = beg_edge_offset2(b_i, b_j);
                #endif
                
                s_start_edge->val = _s_start_edge + offset;
                s_start_edge->atomic_incr(s_i, s_j);

                __sync_fetch_and_add(vert_degree + edge.get_v0(), 1);
                #ifdef HALF_GRID
                __sync_fetch_and_add(vert_degree + edge.get_v1(), 1);
                #endif
            }
        } 
        for (index_t ifile = 0; ifile < vcount; ++ifile) {
            fwrite(buf[ifile], sizeof(cedge_t), count_edge.val[ifile], tf[ifile]);
            count_edge.val[ifile] = 0;
        }
    }    
    
    //Calculate the edge start
    index_t prefix_sum = 0;
    index_t curr_value = 0;
    for (index_t ipart = 0; ipart < total_s_part; ++ipart) {
        curr_value = _s_start_edge[ipart];
        _s_start_edge[ipart] = prefix_sum;
        prefix_sum += curr_value;
    }
    _s_start_edge[total_s_part] = prefix_sum;
    cout << "Total edges = " << prefix_sum << endl;
    
    for (index_t ifile = 0; ifile < vcount; ++ifile) {
        fclose(tf[ifile]);
    }
    fclose(fid_edge);
}

void grid::post_grid_dir(string edgefile, string part_file, bool is_dir)
{
    //Number of files in one dimension 
    index_t p_v = (vert_count >> bit_shift0);
    //Total number of files.
    index_t vcount = calc_total_part(p_v);
   
    //Number of physical groups in a file 
    index_t num_part = p/p_v;
    cout << "Physical groups in one intermediate file (in one direction) = " 
         << num_part <<endl;
    index_t total_parts = num_part*num_part;

    //Buffer for each physical partition in a file
    edge_t** buf = (edge_t**)malloc(total_parts*sizeof(edge_t*));
    
    //Counting number of edges in each tile.
    matrix_f<part_t, index_t>* count_edge = new matrix_f<part_t, index_t>[total_parts];
    
    //Final edge count in each physical partition.
    index_t* ecount = (index_t*)calloc(sizeof(index_t), total_parts);
    
    char tmp[64] = {0};
    string basefile = edgefile;
    if (is_dir) {
        basefile += "/.tmp";
    }
    
    //Read file[0]
    string tfile = basefile + string(".0");
    FILE* tf = fopen(tfile.c_str(), "rb");
    assert(tf != NULL);
    struct stat st;
    stat(tfile.c_str(), &st);
    //Number of edges in this file.
    index_t count2 = st.st_size/sizeof(cedge_t);
    cedge_t* edges2 = (cedge_t*)malloc(st.st_size);
    cout << "Reading " << tfile << endl; 
    fread(edges2, sizeof(cedge_t), count2, tf);
    fclose(tf);

    cedge_t* edges = 0; 
    index_t count = 0;

    #pragma omp parallel num_threads(NUM_THDS) shared(edges2, edges, count2, count) 
    {
        cedge_t edge;
        vertex_t v0, v1;
        part_t b_i, b_j;
        spart_t s_i, s_j;
        index_t offset, n , m;
        
		#ifdef HALF_GRID
		matrix<spart_t, index_t> start_edge_half;
        start_edge_half.part_count = p_p;
		#endif

        matrix_f<spart_t, index_t> start_edge_full;
        start_edge_full.part_count = p_p;
        matrix<spart_t, index_t>* start_edge;
        
        //For each physical group within an intermediate file
        #pragma omp for
        for (index_t ipart = 0; ipart < total_parts; ++ipart) {
            count_edge[ipart].init(p_p);
            buf[ipart] = (edge_t*)malloc(sizeof(edge_t)*(4L<<30));
        }
        
        for (index_t a = 0; a < p_v; ++a) {
            index_t base_b_i = (a*num_part);
            index_t buf_index = 0;
            
            #ifdef HALF_GRID
            for (index_t b = a; b < p_v; ++b) 
            #else
            for (index_t b = 0; b < p_v; ++b) 
            #endif 
            {
                #pragma omp barrier
                if (0 == omp_get_thread_num()) {
                    swap(edges2, edges);
                    swap(count2, count);
                    free(edges2);
                    edges2 = 0;
                }
                #pragma omp barrier

                index_t base_b_j = (b*num_part);
                index_t ifile = calc_index(a, b, p_v) + 1;

                //Read next file. i.e. file[ifile] in one omp thread only
                if (0 == omp_get_thread_num() && (ifile < vcount)) {
                    sprintf(tmp, ".%ld", ifile);
                    tfile = basefile + tmp;
                    tf = fopen(tfile.c_str(), "rb");
                    assert(tf != NULL);
                    cout << "Reading " << tfile << endl; 
                    stat(tfile.c_str(), &st);
                    count2 = st.st_size/sizeof(cedge_t);
                    edges2 = (cedge_t*)malloc(st.st_size);
                    fread(edges2, sizeof(cedge_t), count2, tf);
                    fclose(tf);
                }
                
                #pragma omp for schedule(dynamic, 4096)
                for (index_t k = 0; k < count; ++k) {
                    edge = edges[k];
                    v0 = edge.v0 + (a << bit_shift0);
                    v1 = edge.v1 + (b << bit_shift0);
                    b_i = (v0 >> bit_shift1);
                    b_j = (v1 >> bit_shift1);
                    s_i = ((v0 >> bit_shift2) & part_mask3_2);
                    s_j = ((v1 >> bit_shift2) & part_mask3_2);

                    
                    #ifdef HALF_GRID
                    if (b_i == b_j) {
                        start_edge = &start_edge_half;
                        offset = beg_edge_offset1(b_i);
                    } else {
                        start_edge = &start_edge_full;
                        offset = beg_edge_offset2(b_i,b_j);
                    }
                    #else 
                    start_edge = &start_edge_full;
                    offset = beg_edge_offset2(b_i,b_j);
                    #endif
                    start_edge->val  = _s_start_edge + offset;

                    buf_index = calc_index_f(b_i - base_b_i, b_j - base_b_j, num_part);
                    n = start_edge->get(s_i, s_j) - _s_start_edge[offset];
                    m = count_edge[buf_index].atomic_incr(s_i, s_j);

                    //This assignment takes last 2 bytes automatically
                    buf[buf_index][n + m].v0 = v0; 
                    buf[buf_index][n + m].v1 = v1;

                    //edge count in each physical group 
                    __sync_fetch_and_add(ecount + buf_index, 1);
                }

                //Write the physical partitions now using a single thread. 
                if (0 == omp_get_thread_num()) {

                    //create the file to write. We incremented ifile earlier.
                    sprintf(tmp, ".%ld", ifile - 1);
                    string ofile = part_file + ".grid" + tmp;
                    FILE* f_sep = fopen(ofile.c_str(), "wb");
                    assert(f_sep != 0); 

                    //For each physical group in x-direction
                    for (index_t c = base_b_i; c < base_b_i + num_part; ++c) {

                        #ifdef HALF_GRID
                        if (base_b_i == base_b_j) {
                            //For each physical group in y-direction
                            for (index_t d = c; d < base_b_j + num_part; ++d) {
                                buf_index = calc_index_f(c - base_b_i,
                                                       d - base_b_j, num_part);
                                fwrite(buf[buf_index], sizeof(edge_t), 
                                       ecount[buf_index], f_sep);

                                memset(count_edge[buf_index].val, 0, 
                                       sizeof(index_t)*p_p*p_p);
                                ecount[buf_index] = 0;
                            }

                        } else {
                            #endif
                            //For each physical group in y-direction
                            for (index_t d = base_b_j; d < base_b_j + num_part; ++d) {
                                buf_index = calc_index_f(c - base_b_i, 
                                                         d - base_b_j, num_part);
                                
                                fwrite(buf[buf_index], sizeof(edge_t), 
                                      ecount[buf_index], f_sep);
                                //reset all the counts.
                                memset(count_edge[buf_index].val, 0, 
                                       sizeof(index_t)*p_p*p_p);
                                ecount[buf_index] = 0;
                            }
                        #ifdef HALF_GRID
                        }
                        #endif
                    }
                    fclose(f_sep);
                }
                #pragma omp barrier
            }
        }
    }
}

void grid::post_grid_file(string edgefile, string part_file, bool is_dir)
{
    //Number of files in one dimension 
    index_t p_v = (vert_count >> bit_shift0);
    //Total number of files.
    index_t vcount = calc_total_part(p_v);
   
    //The final edge file.
    string file = part_file + ".grid";
    /*
    index_t total_s_part = calc_total_part(p_s);
    int fd = open(file.c_str(), O_RDWR|O_CREAT, S_IRWXU);
    ftruncate(fd, _s_start_edge[total_s_part]*sizeof(edge_t));
    cout << "size of grid file = " << _s_start_edge[total_s_part]*sizeof(edge_t) << endl;
    close(fd);
    */
    FILE* f = fopen(file.c_str(), "wb");
    assert(f != 0); 
    
    //Number of big partitions in a file row and buffer for each one
    index_t num_part = p/p_v;
    cout << "num_part = " << num_part <<endl;
    index_t total_parts = num_part*num_part;
    edge_t** buf = (edge_t**)malloc(total_parts*sizeof(edge_t*));
    //Counting number of edges in each smaller partitions.
    matrix_f<part_t, index_t>* count_edge = new matrix_f<part_t, index_t>[total_parts];
    
    //final edge count in each partition.
    index_t* ecount = (index_t*)calloc(sizeof(index_t), total_parts);
    
    //Read file[0]
    char tmp[64] = {0};
    string basefile = edgefile;
    if (is_dir) {
        basefile += "/.tmp";
    }
    string tfile = basefile + string(".0");
    FILE* tf = fopen(tfile.c_str(), "rb");
    assert(tf != NULL);
    struct stat st;
    stat(tfile.c_str(), &st);
    //Number of edges in this file.
    index_t count2 = st.st_size/sizeof(cedge_t);
    cedge_t* edges2 = (cedge_t*)malloc(st.st_size);
    cout << "Reading " << tfile << endl; 
    fread(edges2, sizeof(cedge_t), count2, tf);
    fclose(tf);
    cedge_t* edges = 0; 
    index_t count = 0;

    #pragma omp parallel num_threads(NUM_THDS) shared(edges2, edges, count2, count) 
    {
        cedge_t edge;
        vertex_t v0, v1;
        part_t b_i, b_j;
        spart_t s_i, s_j;
        index_t offset, n , m;
        
		#ifdef HALF_GRID
		matrix<spart_t, index_t> start_edge_half;
        start_edge_half.part_count = p_p;
		#endif

        matrix_f<spart_t, index_t> start_edge_full;
        start_edge_full.part_count = p_p;
        matrix<spart_t, index_t>* start_edge;
        
        #pragma omp for
        for (index_t ipart = 0; ipart < total_parts; ++ipart) {
            count_edge[ipart].init(p_p);
            buf[ipart] = (edge_t*)malloc(sizeof(edge_t)*(4L<<30));
        }
        
        for (index_t a = 0; a < p_v; ++a) {
            index_t base_b_i = (a*num_part);
            index_t buf_index = 0;
            
            #ifdef HALF_GRID
            for (index_t b = a; b < p_v; ++b) 
            #else
            for (index_t b = 0; b < p_v; ++b) 
            #endif 
            {
                #pragma omp barrier
                if (0 == omp_get_thread_num()) {
                    swap(edges2, edges);
                    swap(count2, count);
                    free(edges2);
                    edges2 = 0;
                }
                #pragma omp barrier

                index_t base_b_j = (b*num_part);
                index_t ifile = calc_index(a, b, p_v) + 1;
                if (0 == omp_get_thread_num() && (ifile < vcount)) {
                    //Read next file. i.e. file[ifile]
                    sprintf(tmp, ".%ld", ifile);
                    tfile = basefile + tmp;
                    tf = fopen(tfile.c_str(), "rb");
                    assert(tf != NULL);
                    cout << "Reading " << tfile << endl; 
                    stat(tfile.c_str(), &st);
                    count2 = st.st_size/sizeof(cedge_t);
                    edges2 = (cedge_t*)malloc(st.st_size);
                    fread(edges2, sizeof(cedge_t), count2, tf);
                    fclose(tf);
                }
                
                #pragma omp for schedule(dynamic, 4096)
                for (index_t k = 0; k < count; ++k) {
                    edge = edges[k];
                    v0 = edge.v0 + (a << bit_shift0);
                    v1 = edge.v1 + (b << bit_shift0);
                    b_i = (v0 >> bit_shift1);
                    b_j = (v1 >> bit_shift1);
                    s_i = ((v0 >> bit_shift2) & part_mask3_2);
                    s_j = ((v1 >> bit_shift2) & part_mask3_2);

                    
                    #ifdef HALF_GRID
                    if (b_i == b_j) {
                        start_edge = &start_edge_half;
                        offset = beg_edge_offset1(b_i);
                    } else {
                        start_edge = &start_edge_full;
                        offset = beg_edge_offset2(b_i,b_j);
                    }
                    #else 
                    start_edge = &start_edge_full;
                    offset = beg_edge_offset2(b_i,b_j);
                    #endif
                    start_edge->val  = _s_start_edge + offset;

                    buf_index = calc_index_f(b_i - base_b_i, b_j - base_b_j, num_part);
                    n = start_edge->get(s_i, s_j) - _s_start_edge[offset];
                    m = count_edge[buf_index].atomic_incr(s_i, s_j);

                    //This assignment takes last 2 bytes automatically
                    buf[buf_index][n + m].v0 = v0; 
                    buf[buf_index][n + m].v1 = v1;

                    //edge count in each bigger partition
                    __sync_fetch_and_add(ecount + buf_index, 1);
                }
                //Write the big partitions now. 
                if ( 0 == omp_get_thread_num() ) { 
                    for (index_t c = base_b_i; c < base_b_i + num_part; ++c) {
                        #ifdef HALF_GRID
                        if (base_b_i == base_b_j) {
                            offset = beg_edge_offset1(c);

                            for (index_t d = c; d < base_b_j + num_part; ++d) {
                                buf_index = calc_index_f(c - base_b_i, 
                                                       d - base_b_j, num_part);
                                cout << "seeking at" 
                                     << _s_start_edge[offset]*sizeof(edge_t) << endl;

                                if ( -1 == fseek(f, _s_start_edge[offset]*sizeof(edge_t),
                                                 SEEK_SET)) {
                                    cout << "fseek failed" << endl;
                                }
                                fwrite(buf[buf_index], sizeof(edge_t), 
                                       ecount[buf_index], f);

                                memset(count_edge[buf_index].val, 0, 
                                       sizeof(index_t)*p_p*p_p);
                                ecount[buf_index] = 0;
                                offset = beg_edge_offset2(c, d + 1);
                            }
                        } else {
                            #endif
                            for (index_t d = base_b_j; d < base_b_j + num_part; ++d) {
                                buf_index = calc_index_f(c - base_b_i, 
                                                         d - base_b_j, num_part);
                                offset = beg_edge_offset2(c,d);
                                cout << "seeking at" 
                                     << _s_start_edge[offset]*sizeof(edge_t) << endl;
                                if ( -1 == fseek(f, _s_start_edge[offset]*sizeof(edge_t),
                                                 SEEK_SET)) {
                                    cout << "fseek failed" << endl;
                                }
                                fwrite(buf[buf_index], sizeof(edge_t), 
                                      ecount[buf_index], f);
                                //reset all the counts.
                                memset(count_edge[buf_index].val, 0, 
                                       sizeof(index_t)*p_p*p_p);
                                ecount[buf_index] = 0;
                            }
                        #ifdef HALF_GRID
                        }
                        #endif
                    }
                }
                #pragma omp barrier
            }
        }
    }
    fclose(f);
}

void grid::proc_grid_big(string edgefile, string part_file, bool is_odir) 
{
    struct stat sb;
    bool is_dir = false;
    
    if (stat(edgefile.c_str(), &sb) == 0) {
        if (S_ISDIR(sb.st_mode)) {
            is_dir = true;
            pre_grid_dir(edgefile);
        } else {
            pre_grid_file(edgefile);
        }
    } else {
        assert("input file/dir doesn't exist");
    }

    cout << "Intermediate files created" << endl;

    //Do you want a single file or multiple files to be generated.
    if (is_odir) {
        post_grid_dir(edgefile, part_file, is_dir);
    } else {
        post_grid_file(edgefile, part_file, is_dir);
    }

}

void grid::pre_csr(string edgefile, gedge_t* edges, index_t nedges)
{
    _s_start_edge = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
	
	#ifndef HALF_GRID
	s_count_edge_in = (index_t*) calloc(sizeof(index_t), vert_count + 1);  
	#endif

   	vert_degree = (degree_t*)calloc(sizeof(degree_t), vert_count);
    
	#pragma omp parallel num_threads(NUM_THDS) 
    {
        gedge_t  edge;
        vertex_t v0, v1;
        
        #pragma omp for
        for(index_t k = 0; k < nedges; ++k) {
            edge = edges[k];
            if (edge.is_self_loop()) continue;
	    
            v0 = edge.get_v0();
            v1 = edge.get_v1();
            
            __sync_fetch_and_add(_s_start_edge + v0, 1);
			__sync_fetch_and_add(vert_degree + v0, 1);
			
			#ifdef HALF_GRID
            __sync_fetch_and_add(_s_start_edge + v1, 1);
			__sync_fetch_and_add(vert_degree + v1, 1);
			#else 
            __sync_fetch_and_add(s_count_edge_in + v1, 1);
			#endif
        }
    }
    
    //Calculate the CSR beg_pos
    index_t prefix_sum = 0;
	index_t curr_value = 0;
	#ifndef HALF_GRID
    index_t prefix_sum_in = 0;
	index_t curr_value_in = 0;
	#endif
    for (index_t ipart = 0; ipart < vert_count; ++ipart) {
		curr_value = _s_start_edge[ipart];
        _s_start_edge[ipart] = prefix_sum;
        prefix_sum += curr_value;
		#ifndef HALF_GRID
		curr_value_in = s_count_edge_in[ipart];
        s_count_edge_in[ipart] = prefix_sum_in;
        prefix_sum_in += curr_value_in;
		#endif
    }

    _s_start_edge[vert_count] = prefix_sum;
	#ifndef HALF_GRID
    s_count_edge_in[vert_count] = prefix_sum_in;
    cout << "Total edges = " << prefix_sum_in << endl;
	#endif
    cout << "Total edges = " << prefix_sum << endl;

}
void grid::proc_csr(string edgefile, string part_file)
{
    //read the binary edge file
    int fid_edge = open(edgefile.c_str(), O_RDONLY);
    struct stat st_edge;
    fstat(fid_edge, &st_edge);
    assert(st_edge.st_size != 0);
    index_t nedges = st_edge.st_size/sizeof(gedge_t);
    gedge_t* edges;
    /*
    edges = (gedge_t*)mmap(0, st_edge.st_size, PROT_READ, 
                                    MAP_PRIVATE, fid_edge, 0);
    madvise(edges, st_edge.st_size, MADV_SEQUENTIAL);
    */
    edges = (gedge_t*) malloc(st_edge.st_size);
    FILE* f = fopen(edgefile.c_str(), "rb");
    fread(edges, sizeof(gedge_t), nedges, f);

    double start = mywtime();

    pre_csr(edgefile, edges, nedges);

    index_t* count_edge = (index_t*) calloc(sizeof(index_t), vert_count);  
    uint32_t* adj = (uint32_t*)calloc(sizeof(uint32_t), _s_start_edge[vert_count]);
    
	#ifndef HALF_GRID
	uint32_t* adj_in = (uint32_t*)calloc(sizeof(uint32_t), s_count_edge_in[vert_count]);
    index_t* count_edge_in = (index_t*) calloc(sizeof(index_t), vert_count);  
	#endif

    //---classify the edges in the grid
    #pragma omp parallel num_threads(NUM_THDS) 
    { 
        gedge_t  edge;
        vertex_t v0, v1;
        
        index_t n, m;
        
        #pragma omp for
		for(index_t k = 0; k < nedges; ++k) {
			edge = edges[k];
			if (edge.is_self_loop()) continue;
	    
			v0 = edge.get_v0();
			v1 = edge.get_v1();
			
			n = _s_start_edge[v0];
			m = __sync_fetch_and_add(count_edge + v0, 1);
			adj[n + m] = v1;

			#ifdef HALF_GRID
			n = _s_start_edge[v1];
			m = __sync_fetch_and_add(count_edge + v1, 1);
			adj[n + m] = v0;
			#else
			n = s_count_edge_in[v1];
			m = __sync_fetch_and_add(count_edge_in + v1, 1);
			adj_in[n + m] = v0;
			#endif	
		}
    }
   
    //munmap (edges, st_edge.st_size);
	//free(adj);
	//free(_s_start_edge);
	//#ifndef HALF_GRID
	//free(adj_in);	
	//#endif

    close(fid_edge);
    fclose(f); 
    double end = mywtime();
    cout << "CSR conversion Time = " << end -  start << endl;
    cout << "classifcation done" << endl; 
}

void grid::pre_grid(string edgefile, gedge_t* edges, index_t nedges)
{
    index_t total_s_part = calc_total_part(p_s);
    _s_start_edge = (index_t*) calloc(sizeof(index_t), total_s_part);  
    vert_degree = (degree_t*)calloc(sizeof(degree_t), vert_count);
    
    //Do a dry run to know the size of each partition
    #pragma omp parallel num_threads(NUM_THDS) 
    {
        gedge_t  edge;
        part_t j, i;
        spart_t s_i, s_j;
        vertex_t v0, v1;
        vertex_t v2, v3;
        index_t offset = 0;
		
		#ifdef HALF_GRID
        matrix<spart_t, index_t> start_edge_half;
        start_edge_half.part_count = p_p;
		#endif

		matrix_f<spart_t, index_t> start_edge_full;
        start_edge_full.part_count = p_p;
        matrix<spart_t, index_t>* start_edge;
        
        #pragma omp for
        for(index_t k = 0; k < nedges; ++k) {
            edge = edges[k];
            if (edge.is_self_loop()) continue;
	    
            v2 = edge.get_v0();
            v3 = edge.get_v1();
            
			#ifdef HALF_GRID 
			v0 = min(v2, v3);
            v1 = max(v3, v2);
			#else
			v0 = v2;
			v1 = v3;
			#endif

            i = (v0 >> bit_shift1);
            j = (v1 >> bit_shift1);
            s_i = ((v0 >> bit_shift2) & part_mask3_2);
            s_j = ((v1 >> bit_shift2) & part_mask3_2);
            
			#ifdef HALF_GRID
            if (i == j) {
                start_edge = &start_edge_half;
                offset = beg_edge_offset1(i);
                start_edge->val = _s_start_edge + offset;
            } else {
                start_edge = &start_edge_full;
                offset = beg_edge_offset2(i, j);
                start_edge->val  = _s_start_edge + offset;
            }
			#else	
            start_edge = &start_edge_full;
            offset = beg_edge_offset2(i, j);
            start_edge->val  = _s_start_edge + offset;
			#endif

            start_edge->atomic_incr(s_i, s_j);
			__sync_fetch_and_add(vert_degree + edge.get_v0(), 1);
			#ifdef HALF_GRID
			__sync_fetch_and_add(vert_degree + edge.get_v1(), 1);
			#endif
        }
    }
    
    //Calculate the edge start
    index_t prefix_sum = 0;
    index_t curr_value = 0; 
    for (index_t ipart = 0; ipart < total_s_part; ++ipart) {
        curr_value = _s_start_edge[ipart];
        _s_start_edge[ipart] = prefix_sum;
        prefix_sum += curr_value;
    }

    _s_start_edge[total_s_part] = prefix_sum;
    cout << "Total edges = " << prefix_sum << endl;

}

void grid::proc_grid(string edgefile, string part_file)
{
    //read the binary edge file
    int fid_edge = open(edgefile.c_str(), O_RDONLY);
    struct stat st_edge;
    fstat(fid_edge, &st_edge);
    assert(st_edge.st_size != 0);
    index_t nedges = st_edge.st_size/sizeof(gedge_t);
    gedge_t* edges; 
    //gedge_t* edges = (gedge_t*)mmap(0, st_edge.st_size, PROT_READ, 
    //                                MAP_PRIVATE, fid_edge, 0);
    //madvise(edges, st_edge.st_size, MADV_SEQUENTIAL);
    
    edges = (gedge_t*) malloc(st_edge.st_size);
    
    FILE* f = fopen(edgefile.c_str(), "rb");
    fread(edges, sizeof(gedge_t), nedges, f);
    
    double start = mywtime();
    
    pre_grid(edgefile, edges, nedges);


    index_t total_s_part = calc_total_part(p_s);
    index_t* s_count_edge = (index_t*) calloc(sizeof(index_t), total_s_part);  
    _edges = (edge_t*)calloc(sizeof(edge_t), _s_start_edge[total_s_part]);

    //---classify the edges in the grid
    
    //Do it mmap
    //open and truncate
    /*string file = part_file + ".grid";
    int f = open(file.c_str(), O_RDWR|O_CREAT, S_IRWXU);
    if (ENOENT == f ) {
        cout << "failed to create binary grid file" << endl;
        exit(-1);
    }
    
    ftruncate(f, prefix_sum*sizeof(edge_t));
    _edges = (edge_t*)mmap(0, prefix_sum*sizeof(edge_t),
                           PROT_READ|PROT_WRITE,                                  
                           MAP_SHARED,
                           f, 0);
    
    if (MAP_FAILED == _edges) {
            handle_error("pmap alloc");
    }
    */ 
    #pragma omp parallel num_threads(NUM_THDS) 
    { 
        gedge_t  edge;
        part_t j, i;
        spart_t s_i, s_j;
        vertex_t v0, v1;
        vertex_t v2, v3;
        index_t index;
       
		#ifdef HALF_GRID
        matrix<spart_t, index_t> start_edge_half;
        start_edge_half.part_count = p_p;
        matrix<spart_t, index_t> count_edge_half;
        count_edge_half.part_count = p_p;
		#endif

		matrix_f<spart_t, index_t> start_edge_full;
        start_edge_full.part_count = p_p;
        matrix<spart_t, index_t>* start_edge;
        
        matrix_f<spart_t, index_t> count_edge_full;
        count_edge_full.part_count = p_p;
        matrix<spart_t, index_t>* count_edge;

        index_t offset, n, m;
        
        #pragma omp for
		for(index_t k = 0; k < nedges; ++k) {
			edge = edges[k];
			if (edge.is_self_loop()) continue;
	    
			v2 = edge.get_v0();
			v3 = edge.get_v1();
			
			#ifdef HALF_GRID
			v0 = min(v2, v3);
			v1 = max(v3, v2);
			#else 
			v0 = v2;
			v1 = v3;
			#endif

            i = (v0 >> bit_shift1);
            j = (v1 >> bit_shift1);
			s_i = ((v0 >> bit_shift2) & part_mask3_2);
			s_j = ((v1 >> bit_shift2) & part_mask3_2);
			
			#ifdef HALF_GRID
			if (i == j) {
				start_edge = &start_edge_half;
				count_edge = &count_edge_half;
				offset = beg_edge_offset1(i);
			} else {
				start_edge = &start_edge_full;
				count_edge = &count_edge_full;
				offset = beg_edge_offset2(i,j);
			}
			#else 
			start_edge = &start_edge_full;
			count_edge = &count_edge_full;
			offset = beg_edge_offset2(i,j);
			#endif
			
			start_edge->val  = _s_start_edge + offset;
			count_edge->val  = s_count_edge + offset;
            index = start_edge->get_index(s_i, s_j);
			n = start_edge->get(index);
			m = count_edge->atomic_incr(index);
			_edges[n + m].v0 = v0; 
			_edges[n + m].v1 = v1; 
		}
    }
   
    //munmap (edges, st_edge.st_size);
    close(fid_edge); 
    //munmap(_edges, prefix_sum*sizeof(edge_t));
    fclose(f);
    double end = mywtime();
    cout << "classifcation done" << endl; 
    cout << "GStore Conversion Time = " << end - start << endl;
    
}


void grid::compress_degree() 
{
	sdegree_t big_degree = 0;
	
	#pragma omp parallel for reduction(+:big_degree)
	for (vertex_t i = 0; i < vert_count; ++i) {
		big_degree += (vert_degree[i] > 32767);
	}
    
	++big_degree;
	svert_degree = (sdegree_t*)calloc(sizeof(sdegree_t), vert_count);
	bvert_degree = (bdegree_t*)calloc(sizeof(bdegree_t), big_degree);

	bdegree_count = big_degree;
	cout << bdegree_count << endl;
	big_degree = 1;
	bvert_degree[0] = 6000;//put any number
	
    #pragma omp parallel for 
	for (vertex_t i = 0; i < vert_count; ++i) {
		if ((vert_degree[i]) > 32767) {
			sdegree_t m = __sync_fetch_and_add(&big_degree, 1);
			bvert_degree[m] = vert_degree[i];
			svert_degree[i] = -m;
		} else {
			svert_degree[i] = vert_degree[i];
		}
	}
}

//Called after processing cache 
//modify the cache for next iteration 
//Think only from BFS perspective
void grid::cache_to_cache()
{
    //double start = mywtime();
    //source segment
    segment* tmp = cached_pool1;
   
    //destination segment
    segment* seg = cached_pool2;
    seg->buf = tmp->buf;
    seg->ctx_count = 0;

    swap_help(seg, tmp);

    swap(cached_pool1, cached_pool2);
    //double end = mywtime();
    //cout << "cache swap = " << end - start << endl;
}

//Think only from pagerank perspective
void grid::new_swap_parts1()
{
    //double start = mywtime();
    if (seg1->buf + memory - init_buf  <= cache_size) {
		shallow_copy(cached_pool1, seg1);
	}
    if (seg1->ctx_count != 0) {
        //Do some deep copy
        //swap_help(cached_pool1, seg1);
    }
    
    //Swap things now.
    swap(seg1, seg2);
    //double end = mywtime();
    //cout << "swap = " << end - start << endl;
}

void grid::swap_parts()
{
    swap(seg1, seg2);
}


void grid::new_swap_parts2()
{
    segment* seg =  seg1;
    if (seg->buf + memory - init_buf <= cache_size) {
		shallow_copy(cached_pool1, seg);
    } else if (compressed == 0) {
        compressed = 1;
    }

    if (seg1->ctx_count != 0) {
        //Do some deep copy
        swap_help(cached_pool1, seg1);
    }
    swap(seg1, seg2);
}
#define edge_count(i_start, j_start, i_end, j_end)  (start_edge->get_end(i_end, j_end) - start_edge->get(i_start,  j_start))

//destination, source
void grid::swap_help(segment* seg, segment* tmp)
{
    size_t start_addr = 0;
    part_t big_i;
    part_t big_j;
	index_t b_i;
	index_t b_j;
	spart_t i;
	spart_t j;
	spart_t i_end;
	spart_t j_end;
   
	#ifdef HALF_GRID	
    matrix<spart_t, index_t> start_edge_half;
    start_edge_half.part_count = p_p;
	#endif
    matrix_f<spart_t, index_t> start_edge_full;
    start_edge_full.part_count = p_p;
    matrix<spart_t, index_t>* start_edge;
    
    //Find the start_addr
	if (seg->ctx_count) {
		get_ij(seg->meta[seg->ctx_count - 1].start, big_i, big_j, i, j);
		get_s_ij(seg->meta[seg->ctx_count - 1].end, i_end, j_end);
        
		#ifdef HALF_GRID
		if (big_i == big_j) {
            start_edge = &start_edge_half;
            start_edge->val = _s_start_edge + beg_edge_offset1(big_i);
        } else {
            start_edge = &start_edge_full;
            start_edge->val = _s_start_edge + beg_edge_offset2(big_i, big_j);
        }
		#else
		start_edge = &start_edge_full;
        start_edge->val = _s_start_edge + beg_edge_offset2(big_i, big_j);
		#endif

        start_addr = seg->meta[seg->ctx_count - 1].offset 
                        + ((start_edge->get(i, j) 
                            << bytes_in_edge_shift) & 0x1FF)
                        + (edge_count(i, j, i_end, j_end) 
                             << bytes_in_edge_shift);
    }
	
    part_meta_t* meta = tmp->meta;
	part_meta_t* cached_meta = seg->meta;

    for (index_t k = 0; k < tmp->ctx_count; ++k) {
		get_ij(meta[k].start, big_i, big_j, i, j);
		get_s_ij(meta[k].end, i_end, j_end);
	
		#ifdef HALF_GRID	
        spart_t j1_end = p_p - 1;
		if (big_i == big_j) {
            start_edge = &start_edge_half;
            start_edge->val = _s_start_edge + beg_edge_offset1(big_i);
        } else {
            start_edge = &start_edge_full;
            start_edge->val = _s_start_edge + beg_edge_offset2(big_i, big_j);
        }
		#else
        start_edge = &start_edge_full;
        start_edge->val = _s_start_edge + beg_edge_offset2(big_i, big_j);
		#endif

		index_t offset = (calc_index_f(big_i, big_j, p) << bit_shift4);
		b_i = (big_i << bit_shift3);
		b_j = (big_j << bit_shift3);
            
		spart_t i1 = i;
        spart_t j1 = j;
    
		// Align new offset. Add the offset from start edge of i,j 
		char* new_offset = tmp->buf 
			+ tmp->meta[k].offset
			+ ((start_edge->get(i,j) << bytes_in_edge_shift) & 0x1FF);
	
		char* edges = new_offset 
				  - (start_edge->get(i, j) << bytes_in_edge_shift);

		bool start = false;

		//processes only row caches.
		for (i1 = i; i1 <= i_end; ++i1) { 
                    
			//End the ctx
			if (!(read_part_next->get_bit(i1 + b_i) )) {
				if (start == true) {
					start = false;
					cached_meta[seg->ctx_count].end 
							= calc_index_f_opt(i1-1, p_p - 1, bit_shift3) + offset;
					
					if(!end_cache_ctx(seg, tmp, start_addr, start_edge,
								  seg->buf, edges)) {
						//Fetch only some columns
						return;
					}
				}
			}

			//start ctx
			if (read_part_next->get_bit(i1 + b_i)) { 
				if (start == false) {
					start = true;
					cached_meta[seg->ctx_count].start 
							= calc_index_f(i1, j1, p_p) + offset;
				}
			} 
			#ifdef HALF_GRID
			else {//if(read_part[i1])  // row to col cache
				assert(false == start);
				if (i1 == i_end) j1_end = j_end;

				for (index_t m = j1; m <= j1_end; ++m) {
					if (!(read_part_next->get_bit(m + b_j))) {
						//end ctx
						if (start == true) {
							start = false;
							cached_meta[seg->ctx_count].end 
								= calc_index_f(i1, m - 1, p_p) + offset;
							//cached_meta->meta[cached_meta->ctx_count].i_end = i1;
							//cached_meta->meta[cached_meta->ctx_count].j_end = m - 1;
							
							if(!end_cache_ctx(seg, tmp, start_addr, 
											  start_edge, seg->buf, edges)) {
							   return;
							} 
						}
					}
					//start ctx
					if (read_part_next->get_bit(m + b_j)) {
						if (start == false) {
							start = true;
							cached_meta[seg->ctx_count].start 
								= calc_index_f(i1, m, p_p) + offset;
							//cached_meta->meta[cached_meta->ctx_count].i = i1; 
							//cached_meta->meta[cached_meta->ctx_count].j = m;
						}
					} 
				}
			}
			#endif
			j1 = 0;
			
			#ifdef HALF_GRID
			if(b_i == b_j) {
				j1 = i1 + 1;
			}
			#endif	
		}
			
		if (start == true) {
			start = false;
			cached_meta[seg->ctx_count].end 
					= calc_index_f_opt(i_end, j_end, bit_shift3) + offset;
			
			if(!end_cache_ctx(seg, tmp, start_addr, start_edge,
							  seg->buf, edges)) {
				return;
			}
		}
	}
}



bool end_cache_ctx(segment* seg, segment* tmp, size_t& start_addr, 
                   matrix<spart_t, index_t>* start_edge, char* buf, char* edges) 
{
    size_t sz_to_read = 0;
	spart_t i, j, i_end, j_end;
	get_s_ij(seg->meta[seg->ctx_count].start, i, j);
	get_s_ij(seg->meta[seg->ctx_count].end, i_end, j_end);

    sz_to_read = edge_count(i, j, i_end, j_end) << bytes_in_edge_shift;
    
    //Copy some columns as we can't copy rows.
    if (start_addr + sz_to_read > cache_size) {
        return false;
    } 
    sz_copied  += sz_to_read;
    memmove(buf + start_addr, 
           edges + (start_edge->get(i, j) << bytes_in_edge_shift),
           sz_to_read); 

    //The negative number is to fool the any algo_mem_part().
    seg->meta[seg->ctx_count].offset = start_addr 
			- ((start_edge->get(i, j) << bytes_in_edge_shift) & 0x1FF);
    
    seg->ctx_count += 1;
    start_addr += sz_to_read;
    return true;
}

void shallow_copy(segment* cached_pool1, segment* seg)
{
	index_t ctx_count = cached_pool1->ctx_count;
	part_meta_t* cached_meta =  cached_pool1->meta + ctx_count;
	part_meta_t* meta = seg->meta;
	
	//copy the seg1 to cached_pool:Shallow copy
	for(index_t k = 0; k < seg->ctx_count; ++k) {
		
		cached_meta[k].start = meta[k].start;
		cached_meta[k].end = meta[k].end;
		cached_meta[k].offset = 
				meta[k].offset + seg->buf - cached_pool1->buf;
	}
	
	cached_pool1->ctx_count += seg->ctx_count;
	seg->buf = 0;
	seg->ctx_count = 0;

}

void grid::do_algo(algo_t* algo)
{
    double start = mywtime();
    last_read_t	  last_read;
    
    last_read.s_i   = -1;
    last_read.b_i   = -1;
    last_read.b_j   = -1;
    last_read.s_j   = p_p - 1;
    cache->start_index = 0;
    
    cache->fetch_mem_part(last_read, seg1, read_part);
    
    #pragma omp parallel num_threads(IO_THDS)
    {
        io->read_aio_random(seg1);
        io->wait_aio_completion();
    }
    last_read_t	  last_read1 = last_read; 
    last_read_t   last_read2 = last_read;
    
    bool          iteration_start = false;

    if (!is_end(last_read2) && !iteration_start) {
        cache->fetch_mem_part(last_read2, seg3, read_part);
    }
    
    while(true) {        
        //double st1 = mywtime();
        if (!is_end(last_read1) && !iteration_start) {
            seg3->buf = seg2->buf;
            swap(seg3,seg2);
            seg3->buf = 0;
            last_read = last_read2;
        }

        #pragma omp parallel num_threads(NUM_THDS) //collapse(2)
        {
            if (omp_get_thread_num() < IO_THDS) 
            {
               //End of iteration. Do nothing
               //Start of iteration. Do nothing
                if (!is_end(last_read1) && !iteration_start) {
                    //cache->fetch_mem_part(last_read, seg2, read_part);
                    io->read_aio_random(seg2);
                }
            }

            if ((omp_get_thread_num() == IO_THDS) && !iteration_start) {
                if (compressed == 1) {
                    cache_to_cache();
                    compressed = 2;
                }
	        }
            
            if (omp_get_thread_num() == IO_THDS + 1) {
                if (!is_end(last_read)) {
                    cache->fetch_mem_part(last_read2, seg3, read_part);
                }
            }
        
            if (2+IO_THDS == omp_get_thread_num()) {
                if (is_end(last_read1)) { 
                    //Single threaded block for now.
                    cache->prep_pread_cached(seg2, 0);
                    cache->prep_pread_cached(seg1, seg2->ctx_count);
                    cache->prep_pread_cached(cached_pool1, 
                                        seg2->ctx_count + seg1->ctx_count);
                    cache->ctx_count2 = seg2->ctx_count + 
                                seg1->ctx_count + cached_pool1->ctx_count;
                    std::sort(cache->part_cached, 
                              cache->part_cached + cache->ctx_count2);
                }
            }

            algo->algo_mem_part(seg1);
            if (iteration_start) {
                algo->algo_mem_part(cached_pool1);
            }

            //wait for aio completion
            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    io->wait_aio_completion();
                }
            }
        }
           
        //double st = mywtime();
        //cout << " part time = " << st - st1 << endl; 
        
        if (iteration_start) { 
            iteration_start = false;
        }

        if (is_end(last_read1)) {
            //one iteration done
            compressed = 0;
            last_read.s_i   = -1;
            last_read.b_i   = -1;
            last_read.b_j   = -1;
            last_read.s_j   = p_p - 1;
            
            last_read1.s_i   = 0;
            last_read1.b_i   = 0;
            last_read1.b_j   = 0;
            last_read1.s_j   = 0;

            last_read2 = last_read;
            
            iteration_start = true;
            double end = mywtime();
            cout << "iteration time = " << end -start << endl;
            start = mywtime();
            if( 0 != algo->iteration_finalize()) { 
                break;
            }
            //cout << "size copied = " << sz_copied << endl;
            sz_copied = 0;
            cache->start_index = 0;

        } else {
            last_read1 = last_read;
	        new_swap_parts2();
        }
        //double en = mywtime();
        //cout << "post time = " << en - st << endl;
    }
}

void 
grid::bfs2()
{
    double start = mywtime();
    last_read_t	  last_read;
    last_read.s_i   = -1;
    last_read.b_i   = -1;
    last_read.b_j   = -1;
    last_read.s_j   = p_p - 1;
    cache->start_index = 0;
    
    bfs2_t pr;
    int root = arg;
    if (arg == -1) root = 1;
    pr.init(vert_count, root, read_part); 
   
	index_t front_count = 0;
	bool io_skip = 1;	
    cache->fetch_mem_part(last_read, seg1, read_part, io_skip);
    #pragma omp parallel num_threads(IO_THDS)
    {
        io->read_aio_random(seg1);
        io->wait_aio_completion();
    }
    last_read_t	  last_read1 = last_read; 
    last_read_t   last_read2 = last_read;
    bool iteration_start = false;
    
    if (!is_end(last_read2) && !iteration_start) {
        cache->fetch_mem_part(last_read2, seg3, read_part, io_skip);
    }
        
    while(true) {
        if (!is_end(last_read1) && !iteration_start) {
            seg3->buf = seg2->buf;
            swap(seg3,seg2);
            seg3->buf = 0;
            last_read = last_read2;
        }

        #pragma omp parallel num_threads(NUM_THDS) 
        {
            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    //cache->fetch_mem_part(last_read, seg2, read_part, false);
                    io->read_aio_random(seg2);
                }
                
            }
			if (compressed == 1 && !iteration_start) {
				//pr.algo_mem_part2(cached_pool1);
				if (omp_get_thread_num() == IO_THDS) {
                    cache_to_cache();
                    compressed = 2;
                }
			}
            
            if (omp_get_thread_num() == IO_THDS+1) {
                if (!is_end(last_read)) {
                    cache->fetch_mem_part(last_read2, seg3, read_part, io_skip);
                }
            }

            if (2+IO_THDS == omp_get_thread_num()) {
                if (is_end(last_read1)) { 
                    cache->prep_pread_cached(seg2, 0);
                    cache->prep_pread_cached(seg1, seg2->ctx_count);
                    cache->prep_pread_cached(cached_pool1, 
                                      seg2->ctx_count + seg1->ctx_count);
                    cache->ctx_count2 = seg2->ctx_count 
                                + seg1->ctx_count + cached_pool1->ctx_count;
                    std::sort(cache->part_cached, 
                              cache->part_cached + cache->ctx_count2);
                }
            }
            
            pr.algo_mem_part(seg1);
			if (iteration_start) {
				pr.algo_mem_part(cached_pool1);
			} 
			//else  if (io_skip == 0 && compressed != 1) {
				//pr.algo_mem_part2(cached_pool1);
			//}
            //pr.algo_mem_part2(seg1);
            
            //Wait for aio
            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    io->wait_aio_completion();
                }
            }
        }
        
        if (iteration_start) { 
            iteration_start = false;
        }
        
		if (is_end(last_read1)) {
            //one iteration done
			last_read.s_i   = -1;
			last_read.b_i   = -1;
			last_read.b_j   = -1;
			last_read.s_j   = p_p - 1;
			
			last_read1.s_i   = 0;
			last_read1.b_i   = 0;
			last_read1.b_j   = 0;
			last_read1.s_j   = 0;
            
            last_read2 = last_read;

            iteration_start = true;
            double end = mywtime();
            cout << "iteration time = " << end -start << endl;
            front_count = pr.iteration_finalize();
			if(0 == front_count) break;
			#ifdef HALF_GRID
			else if (front_count <= 16*p) {
			   	io_skip = 0;
				//cout << "io_skip = " << (int)io_skip << endl;
			}
			#else
			else if (front_count <= 4*p_s) {
			   	io_skip = 0;
				//cout << "io_skip = " << (int)io_skip << endl;
			}
			#endif
            cache->start_index = 0;
            start = mywtime();
            //cout << "post iteration time = " << start - end << endl;
        } else {
			if (io_skip == 1) {
				new_swap_parts1();
			} else {
				new_swap_parts2();
			}
            last_read1 = last_read;
		}
    }
}

void grid::do_algo_prop(algo_t& pr)
{
    double start = mywtime();
    
    last_read_t	  last_read;
    last_read.s_i   = -1;
    last_read.b_i   = -1;
    last_read.b_j   = -1;
    last_read.s_j   = p_p - 1;
    cache->start_index = 0;

    
    cache->fetch_mem_part(last_read, seg1, read_part, 1);
    #pragma omp parallel num_threads(IO_THDS)
    {
        io->read_aio_random(seg1);
        io->wait_aio_completion();
    }

    last_read_t	  last_read1 = last_read; 
    last_read_t   last_read2 = last_read;
    bool iteration_start = false;
    
    if (!is_end(last_read2) && !iteration_start) {
        cache->fetch_mem_part(last_read2, seg3, read_part, 1);
    }
        
    while(true) {
        if (!is_end(last_read1) && !iteration_start) {
            seg3->buf = seg2->buf;
            swap(seg3,seg2);
            seg3->buf = 0;
            last_read = last_read2;
        }

        #pragma omp parallel num_threads(NUM_THDS) 
        {
            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    //cache->fetch_mem_part(last_read, seg2, read_part, false);
                    io->read_aio_random(seg2);
                }
                
            }
            
            if (omp_get_thread_num() == IO_THDS) {
                if (!is_end(last_read)) {
                    cache->fetch_mem_part(last_read2, seg3, read_part, 1);
                }
            }

            if (IO_THDS+1 == omp_get_thread_num()) {
                if (is_end(last_read1)) { 
                    cache->prep_pread_cached(seg2, 0);
                    cache->prep_pread_cached(seg1, seg2->ctx_count);
                    cache->prep_pread_cached(cached_pool1, 
                                      seg2->ctx_count + seg1->ctx_count);
                    cache->ctx_count2 = seg2->ctx_count 
                                + seg1->ctx_count + cached_pool1->ctx_count;
                    std::sort(cache->part_cached, 
                              cache->part_cached + cache->ctx_count2);
                }
            }
            
			if (iteration_start) {
				pr.algo_mem_part(cached_pool1);
				//pr.algo_mem_part(cached_pool1);
			} else if (is_end(last_read1)) {
				pr.algo_mem_part(cached_pool1);
            }
            pr.algo_mem_part(seg1);
			//pr.algo_mem_part(cached_pool1);
            //pr.algo_mem_part(seg1);
            
            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    io->wait_aio_completion();
                }
            }
        }
        
        if (iteration_start) { 
            iteration_start = false;
        }
        
		if (is_end(last_read1)) {
            //one iteration done
			last_read.s_i   = -1;
			last_read.b_i   = -1;
			last_read.b_j   = -1;
			last_read.s_j   = p_p - 1;
			
			last_read1.s_i   = 0;
			last_read1.b_i   = 0;
			last_read1.b_j   = 0;
			last_read1.s_j   = 0;
            
            last_read2 = last_read;

            iteration_start = true;
            double end = mywtime();
            cout << "iteration time = " << end -start << endl;
            if(0 == pr.iteration_finalize()) break;
            cache->start_index = 0;
            start = mywtime();
            //cout << "post iteration time = " << start - end << endl;
        } else {
			new_swap_parts1();
            last_read1 = last_read;
		}
    }
}

void 
grid::wcc()
{
    //double start = mywtime();
    wcc2_t pr;
    pr.init(vert_count); 
   
    do_algo_prop(pr);
}

void grid::traverse()
{
    traverse_t pr;
    vertex_t pivot = arg;
    if (arg == -1) pivot = 1;
    pr.init(vert_count, pivot, 0);

    do_algo_prop(pr);
}

//step-back bypassed here. Used for measuring motivation.
void
grid::pagerank()
{
    double start = mywtime();
    
    last_read_t	  last_read;
    last_read.s_i   = -1;
    last_read.b_i   = -1;
    last_read.b_j   = -1;
    last_read.s_j   = p_p - 1;
    cache->start_index = 0;

    int iteration = 0;
    int iteration_count = arg;
	if (arg == -1) {
		iteration_count = 5;
	}
    pr_t pr;
    pr.init(vert_count, iteration_count); 
    
    cache->fetch_mem_part(last_read, seg1, read_part, 1);
    #pragma omp parallel num_threads(IO_THDS)
    {
        io->read_aio_random(seg1);
        io->wait_aio_completion();
    }

    last_read_t	  last_read1 = last_read; 
    last_read_t   last_read2 = last_read;
    bool iteration_start = false;
    
    if (!is_end(last_read2) && !iteration_start) {
        cache->fetch_mem_part(last_read2, seg3, read_part, 1);
    }
        
    while(true) {
        if (!is_end(last_read1) && !iteration_start) {
            seg3->buf = seg2->buf;
            swap(seg3,seg2);
            seg3->buf = 0;
            last_read = last_read2;
        }

        #pragma omp parallel num_threads(NUM_THDS) //collapse(2)
        {
            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    //cache->fetch_mem_part(last_read, seg2, read_part, false);
                    io->read_aio_random(seg2);
                }
                
            }
            
            if (omp_get_thread_num() == IO_THDS) {
                if (!is_end(last_read)) {
                    cache->fetch_mem_part(last_read2, seg3, read_part, 1);
                }
            }

            if (1+IO_THDS == omp_get_thread_num()) {
                if (is_end(last_read1)) { 
                    cache->prep_pread_cached(seg2, 0);
                    cache->prep_pread_cached(seg1, seg2->ctx_count);
                    cache->prep_pread_cached(cached_pool1, 
                                      seg2->ctx_count + seg1->ctx_count);
                    cache->ctx_count2 = seg2->ctx_count 
                                + seg1->ctx_count + cached_pool1->ctx_count;
                    std::sort(cache->part_cached, 
                              cache->part_cached + cache->ctx_count2);
                }
            }
            
            pr.algo_mem_part(seg1);
			if (iteration_start) {
				pr.algo_mem_part(cached_pool1);
			}
            
            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    io->wait_aio_completion();
                }
            }
        }
        
        if (iteration_start) { 
            iteration_start = false;
        }
        
		if (is_end(last_read1)) {
            //one iteration done
			last_read.s_i   = -1;
			last_read.b_i   = -1;
			last_read.b_j   = -1;
			last_read.s_j   = p_p - 1;
			
			last_read1.s_i   = 0;
			last_read1.b_i   = 0;
			last_read1.b_j   = 0;
			last_read1.s_j   = 0;
            
            last_read2 = last_read;

            iteration_start = true;
            double start1 = mywtime();
			#pragma omp parallel num_threads(NUM_THDS)
            pr.iteration_finalize(iteration); 
            double end = mywtime();
            cout << "post-processing time = " << end - start1 << endl;
            cout << "iteration time = " << end -start << endl;
            start = mywtime();
            ++iteration;
            cache->start_index = 0;
            if (iteration == iteration_count) break;
        } else {
			new_swap_parts1();
            last_read1 = last_read;
		}
    }
}

void grid::bfs()
{
	double start1 = mywtime();
    //vertex_t root = 250972645;
    vertex_t root = arg;
    if (arg == -1) root = 0; 
    bfs_t    bf;
    bf.init(vert_count, root, read_part);

    do_algo(&bf);
	double end1 = mywtime();
	cout << end1 - start1 << endl;
    return;
}


/*
void grid::kcore(int kc)
{
	#ifndef HALF_GRID
	cout << "Kcore not applicable for directed graphs" << endl;
	return;
	#endif

	if (kc == -1) kc = 10;
    //Initialization
    kcore_t kcorep;
    kcorep.init(vert_count, kc, vert_degree, read_part);
    if( 0 == kcorep.iteration_finalize()) return;
    
    //Main body
    do_algo(&kcorep);

	//Final result
    vertex_t total = 0;
    #pragma omp parallel for reduction(+:total)
    for(vertex_t i = 0; i < vert_count; i++) {
		if(g->vert_degree[i] != 0) ++total;
    }
    cout << ">= " << kc << "vertices count = " << total << endl;
}
*/


void grid::kcore(int kc)
{
	#ifndef HALF_GRID
	cout << "Kcore not applicable for directed graphs" << endl;
	return;
	#endif

	if (kc == -1) kc = 10;
    //Initialization
    kcore_t kcorep;
    kcorep.init(vert_count, kc, vert_degree, read_part);
    if(0 == kcorep.iteration_finalize()) return;
    kcore_t* algo = &kcorep;

    double start = mywtime();
    last_read_t	  last_read;
    
    last_read.s_i   = -1;
    last_read.b_i   = -1;
    last_read.b_j   = -1;
    last_read.s_j   = p_p - 1;
    cache->start_index = 0;
    
    int io_skip = 0;
    cache->fetch_mem_part(last_read, seg1, read_part, io_skip);
    #pragma omp parallel num_threads(IO_THDS)
    {
        io->read_aio_random(seg1);
        io->wait_aio_completion();
    }

    last_read_t	  last_read1 = last_read; 
    last_read_t   last_read2 = last_read;
    
    bool          iteration_start = false;
    if (!is_end(last_read2) && !iteration_start) {
        cache->fetch_mem_part(last_read2, seg3, read_part, io_skip);
    }
    
    while(true) {        
        //double st1 = mywtime();
        if (!is_end(last_read1) && !iteration_start) {
            seg3->buf = seg2->buf;
            swap(seg3,seg2);
            seg3->buf = 0;
            last_read = last_read2;
        }

        #pragma omp parallel num_threads(NUM_THDS) //collapse(2)
        {
            if (omp_get_thread_num() < IO_THDS) 
            {
               //End of iteration. Do nothing
               //Start of iteration. Do nothing
                if (!is_end(last_read1) && !iteration_start) {
                    //cache->fetch_mem_part(last_read, seg2, read_part);
                    io->read_aio_random(seg2);
                }
            }

            if (omp_get_thread_num() == IO_THDS && !iteration_start) {
                if (compressed == 1) {
                    cache_to_cache();
                    compressed = 2;
                }
	        }
            
            if (omp_get_thread_num() == IO_THDS + 1) {
                if (!is_end(last_read)) {
                    cache->fetch_mem_part(last_read2, seg3, read_part, io_skip);
                }
            }
        
            if (2+IO_THDS == omp_get_thread_num()) {
                if (is_end(last_read1)) { 
                    //Single threaded block for now.
                    cache->prep_pread_cached(seg2, 0);
                    cache->prep_pread_cached(seg1, seg2->ctx_count);
                    cache->prep_pread_cached(cached_pool1, 
                                        seg2->ctx_count + seg1->ctx_count);
                    cache->ctx_count2 = seg2->ctx_count + 
                                seg1->ctx_count + cached_pool1->ctx_count;
                    std::sort(cache->part_cached, 
                              cache->part_cached + cache->ctx_count2);
                }
            }

            algo->algo_mem_part(seg1);
            if (iteration_start) {
                algo->algo_mem_part(cached_pool1);
            }

            if (omp_get_thread_num() < IO_THDS) 
            {
                if (!is_end(last_read1) && !iteration_start) {
                    io->wait_aio_completion();
                }
            }
        }
           
        //double st = mywtime();
        //cout << " part time = " << st - st1 << endl; 
        
        if (iteration_start) { 
            iteration_start = false;
        }

        if (is_end(last_read1)) {
            //one iteration done
            compressed = 0;
            last_read.s_i   = -1;
            last_read.b_i   = -1;
            last_read.b_j   = -1;
            last_read.s_j   = p_p - 1;
            
            last_read1.s_i   = 0;
            last_read1.b_i   = 0;
            last_read1.b_j   = 0;
            last_read1.s_j   = 0;

            last_read2 = last_read;
            
            iteration_start = true;
            double end = mywtime();
            cout << "iteration time = " << end -start << endl;
            start = mywtime();
            index_t front_count = algo->iteration_finalize();
			if(0 == front_count) break;
			//else if (front_count <= 16*p) {
			   	//io_skip = 0;
				//cout << "io_skip = " << (int)io_skip << endl;
			//}

            //cout << "size copied = " << sz_copied << endl;
            sz_copied = 0;
            cache->start_index = 0;

        } else {
			if (io_skip == 1) {
				new_swap_parts1();
			}
            else {
				new_swap_parts2();
			}
            last_read1 = last_read;
        }
    }
	//Final result
    vertex_t total = 0;
    #pragma omp parallel for reduction(+:total)
    for(vertex_t i = 0; i < vert_count; i++) {
		if(g->vert_degree[i] != 0) ++total;
    }
    cout << ">= " << kc << "vertices count = " << total << endl;
}

void
grid::wcc_mmap()
{
    wcc2_t pr;
	pr.init(vert_count);
    int iteration = 0;
       
    while (true)	 {
        double start = mywtime();
		#pragma omp parallel num_threads(NUM_THDS) //collapse(2)
        {
			vertex_t t_front_count = 0;
			matrix<spart_t, index_t> start_edge_half;
			start_edge_half.part_count = p_p;
			matrix_f<spart_t, index_t> start_edge_full;
			start_edge_full.part_count = p_p;
			matrix<spart_t, index_t>* start_edge;
			spart_t n2 = 0;
	        for (part_t i = 0; i < p; ++i) {
				index_t b_i = (i << bit_shift3);
                start_edge = &start_edge_half;
                start_edge->val = _s_start_edge + beg_edge_offset1(i);

				#ifdef HALF_GRID
                for (part_t j = i; j < p; ++j) 
				#else
				for (part_t j = 0; j < p; ++j) 
				#endif
				{
					index_t b_j = (j << bit_shift3);
                    #pragma omp for nowait schedule (dynamic, 1)
                    for (spart_t m = 0; m < p_p; ++m) {
						n2 = 0;
						#ifdef HALF_GRID
						if (i == j) {
                            n2 = m;
                        }
						#endif
                        if (iteration == 0){
                        //#pragma omp for nowait schedule (dynamic, 8)
                        for (spart_t n = n2; n < p_p; ++n) {
                            edge_t* part_edge = _edges + start_edge->get(m, n);
                            index_t cedge = start_edge->get_count(m, n); 
                            t_front_count += pr.wcc_onepart(part_edge, cedge, m + b_i, n + b_j);
                        }
                        } else {
                        //#pragma omp for nowait schedule (dynamic, 8)
                        for (spart_t n = n2; n < p_p; ++n) {
                            edge_t* part_edge = _edges + start_edge->get(m, n);
                            index_t cedge = start_edge->get_count(m, n); 
                            t_front_count += pr.wcc_onepart2(part_edge, cedge, m + b_i, n + b_j);
                        }
                        }
                    }

                    start_edge = &start_edge_full;
                    start_edge->val  = _s_start_edge + beg_edge_offset2(i, j + 1);
                }
            }
			pr.front_count[omp_get_thread_num()] += t_front_count;
        }
        double end = mywtime();
        cout << "iteration time = " << end - start << endl;
	    if (0 == pr.iteration_finalize()) {
            break;
        }
        ++iteration;
    }
}

void
grid::pagerank_mmap()
{
    int iteration_count = arg;
	if (arg == -1) {
		iteration_count = 5;
	}
    pr_t pr;
    pr.init(vert_count, iteration_count);
        
    #pragma omp parallel num_threads(NUM_THDS) //collapse(2)
    {
        matrix<spart_t, index_t> start_edge_half;
        start_edge_half.part_count = p_p;
        matrix_f<spart_t, index_t> start_edge_full;
        start_edge_full.part_count = p_p;
        matrix<spart_t, index_t>* start_edge;
	    spart_t n2 = 0;
	    for(int iteration =0; iteration < iteration_count; ++iteration) {
	        for (part_t i = 0; i < p; ++i) {
				index_t b_i = (i << bit_shift3);
                start_edge = &start_edge_half;
                start_edge->val = _s_start_edge + beg_edge_offset1(i);

				#ifdef HALF_GRID
                for (part_t j = i; j < p; ++j) 
				#else
				for (part_t j = 0; j < p; ++j) 
				#endif
				{
					index_t b_j = (j << bit_shift3);
                    #pragma omp for nowait schedule (dynamic, 1)
                    for (spart_t m = 0; m < p_p; ++m) {
						n2 = 0;
						#ifdef HALF_GRID
						if (i == j) {
                            n2 = m;
                        }
						#endif
                        for (spart_t n = n2; n < p_p; ++n) {
                            edge_t* part_edge = _edges + start_edge->get(m, n);
                            index_t cedge = start_edge->get_count(m, n); 
                            pr.pagerank_onepart(part_edge, cedge, m + b_i, n + b_j);
                        }
                    }

                    start_edge = &start_edge_full;
                    start_edge->val  = _s_start_edge + beg_edge_offset2(i, j + 1);
                }
            }

	        #pragma omp barrier
	
	        pr.iteration_finalize(iteration);
			#pragma omp barrier
        }
    }
}

void grid::bfs_mmap()
{
    //vertex_t root = 2982970;
    vertex_t root =1 ;//28482313;// 2982970;
    bfs_t bf;
    bf.init(vert_count, root, read_part);
    
    while (true) {
        double start = mywtime();
   
        #pragma omp parallel num_threads(NUM_THDS)
        {
			vertex_t  t_front_count = 0;
			matrix<spart_t, index_t> start_edge_half;
			start_edge_half.part_count = p_p;
			matrix_f<spart_t, index_t> start_edge_full;
			start_edge_full.part_count = p_p;
			matrix<spart_t, index_t>* start_edge;
			spart_t n2 = 0;

            for (part_t i = 0; i < p; ++i) {
				index_t b_i = (i << bit_shift3);
				start_edge = &start_edge_half;
				start_edge->val = _s_start_edge + beg_edge_offset1(i);

                #ifdef HALF_GRID
                for (part_t j = i; j < p; ++j) 
				#else
                for (part_t j = 0; j < p; ++j) 
				#endif
				{
					index_t b_j = (j << bit_shift3);
                    #pragma omp for schedule (dynamic, 1) nowait
                    for (part_t k = 0; k < p_p; ++k) {
						n2 = 0;
						#ifdef HALF_GRID
						if (i == j) {
							n2 = k;
						}
                        #endif
					    #ifndef HALF_GRID	
						if (read_part->get_bit(k + b_i)) {
                            for (part_t l = n2; l < p_p; ++l) {
                                edge_t* part_edge = _edges 
                                        + start_edge->get(k, l);
                                index_t cedge = start_edge->get_count(k, l);
                                t_front_count += bf.bfs_onepart(part_edge, 
									cedge, k + b_i, l+b_j);
                            }
                        }
						#else	
                        for (part_t l = n2; l < p_p; ++l) {
                            edge_t* part_edge = _edges 
                                    + start_edge->get(k, l);
                            index_t cedge = start_edge->get_count(k, l);
                            
                            if (read_part->get_bit(k + b_i) 
                                && read_part->get_bit(l+b_j)) {
                                t_front_count += bf.bfs_onepart(
                                        part_edge, cedge, k + b_i, l + b_j);

                            } else if (read_part->get_bit(l + b_j)) {
                                t_front_count += bf.bfs_onepart_col(
                                        part_edge, cedge, k + b_i, l + b_j);
                            } else if (read_part->get_bit(k + b_i)) {
                                t_front_count += bf.bfs_onepart_row(
                                        part_edge, cedge, k + b_i, l + b_j);
                            }
                        }
						#endif
                    }
					start_edge = &start_edge_full;
					start_edge->val  = _s_start_edge + beg_edge_offset2(i, j + 1);
                }
            }
            bf.front_count[omp_get_thread_num()] += t_front_count;
        } 
        double end = mywtime();
        cout << "iteration time = " << end - start << endl; 
        if(bf.iteration_finalize()) break;
    }
    return;
}


void grid::kcore_mmap(int kc)
{
	#ifndef HALF_GRID
	cout << "Kcore not applicable for directed graphs" << endl;
	return;
	#endif
    
	kcore_t kcorep;
	if (kc == -1) kc = 10;
    kcorep.init(vert_count, kc, vert_degree, read_part);
    if(0 == kcorep.iteration_finalize()) return;
    
    while(true) {  
        double start = mywtime();      
        #pragma omp parallel num_threads(NUM_THDS) //collapse(2)
        {
	    matrix<spart_t, index_t> start_edge_half;
	    start_edge_half.part_count = p_p;
	    matrix_f<spart_t, index_t> start_edge_full;
	    start_edge_full.part_count = p_p;
	    matrix<spart_t, index_t>* start_edge;
	    spart_t n2 = 0;

        for (part_t i = 0; i < p; ++i) {
			index_t b_i = (i << bit_shift3);
			start_edge = &start_edge_half;
			start_edge->val = _s_start_edge + beg_edge_offset1(i);
			#ifdef HALF_GRID    
            for (part_t j = i; j < p; ++j) 
			#else
            for (part_t j = 0; j < p; ++j) 
			#endif
			{
				index_t b_j = (j << bit_shift3);
				#pragma omp for schedule (dynamic, 1) nowait
				for (part_t k = 0; k < p_p; ++k) {
					n2 = 0;
					#ifdef HALF_GRID
					if (i == j) {
						n2 = k;
					}
					#endif
					#ifndef HALF_GRID
					if (read_part->get_bit(k + b_i)) {
						for (part_t l = n2; l < p_p; ++l) {
							edge_t* part_edge = _edges 
									+ start_edge->get(k, l);
							index_t cedge = start_edge->get_count(k, l);
							kcorep.kcore_onepart(part_edge, 
							cedge, k + b_i, l+b_j);
						}
					}
					#else
                    for (part_t l = n2; l < p_p; ++l) {
                        edge_t* part_edge = _edges 
                                + start_edge->get(k, l);
                        index_t cedge = start_edge->get_count(k, l);
                        if (read_part->get_bit(k+b_i) && read_part->get_bit(l+b_j)) {
                            kcorep.kcore_onepart( part_edge, cedge, 
                                                  k + b_i, l + b_j);
                        } else if (read_part->get_bit(k + b_i)) {
                            kcorep.kcore_onepart_row( part_edge, cedge, 
                                                      k + b_i, l + b_j);
                        } else if (read_part->get_bit(l + b_j)) {
                            kcorep.kcore_onepart_col( part_edge, cedge, 
                                                      k + b_i, l + b_j);
                        } 
                    }
					#endif
				}
				start_edge = &start_edge_full;
				start_edge->val  = _s_start_edge + beg_edge_offset2(i, j + 1);
			}
        }
    } 
    double  end = mywtime();
    cout << "iteration time = " << end - start << endl;
	if(0 == kcorep.iteration_finalize()) break;
    }
	
    vertex_t total = 0;
    #pragma omp parallel for reduction(+:total)
    for(vertex_t i = 0; i < vert_count; i++) {
		if(g->vert_degree[i] != 0) ++total;
    }
    

	cout << ">= " << kc << "vertices count = " << total << endl;
}


int 
main(int argc, char *argv[]) 
{
    g = new grid;
    cache = new cache_driver;
    io = new io_driver;
    g->init(argc, argv);

    return 0;
}

void grid::save_grid(string edgefile)
{
    cout << "save grid start" << endl;
    string file = edgefile + ".grid";
    FILE* f = fopen(file.c_str(), "wb");
    assert(f != 0);
    index_t size = calc_total_part(p_s);
    fwrite(_edges, sizeof(edge_t), _s_start_edge[size], f);
    fclose(f);
    save_start_files(edgefile, false);
    save_degree_files(edgefile);
} 

void grid::save_grid_big(string edgefile, bool is_odir)
{
    save_start_files(edgefile, is_odir);
    save_degree_files(edgefile);
}

void grid::save_start_files(string edgefile, bool is_odir)
{
    if (is_odir) {
        //Number of files in one dimension 
        index_t p_v = (vert_count >> bit_shift0);
        
        //Total number of files.
        index_t vcount = calc_total_part(p_v);
        index_t num_part = p/p_v;
        index_t total_parts = num_part*num_part;

        
        index_t* start_edge = (index_t*)calloc(sizeof(index_t), p_p*total_parts*p_p);
        
        for (index_t a = 0; a < p_v; ++a) {
            index_t base_b_i = (a*num_part);
            #ifdef HALF_GRID
            for (index_t b = a; b < p_v; ++b) 
            #else
            for (index_t b = 0; b < p_v; ++b) 
            #endif 
            {
                //The code is for each intermediate file
                index_t base_tile_index = 0; 
                index_t tile_count = 0, total_tile_count = 0;
                index_t edge_count = 0, prefix = 0;
                index_t base_b_j = (b*num_part);
                index_t ifile = calc_index(a, b, p_v);
                
                // for each row of physical groups 
                for (index_t c = base_b_i; c < base_b_i + num_part; ++c) {
                    #ifdef HALF_GRID
                    if (base_b_i == base_b_j) 
                        base_tile_index = beg_edge_offset1(c);
                    else
                    #endif
                        base_tile_index = beg_edge_offset2(c, base_b_j);
                    
                        
                    tile_count = beg_edge_offset2(c, base_b_j + num_part) - base_tile_index;
                    
                    //For each tile in all physical groups in one row
                    //The last iteration is for +1 that we follow, will be overwritten in 
                    //next execution of this for loop
                    for (index_t d = 0; d <= tile_count; ++d) {
                        edge_count = _s_start_edge[base_tile_index + d] - _s_start_edge[base_tile_index];
                        prefix += edge_count;
                        start_edge[total_tile_count + d] = prefix;
                    }
                    total_tile_count += tile_count;
                }
                //Write the start_edge for each file
                char tmp[32];
                sprintf(tmp, ".%ld", ifile);
                string file = edgefile + ".start" + tmp;
                FILE* f = fopen(file.c_str(), "wb");
                assert(f != 0);
                fwrite(start_edge, sizeof(index_t), total_tile_count, f);
                fclose(f);
                cout << "Total Edges in file" << ifile << " = "<< prefix << endl;
            }
        }
    }else {
        index_t size = calc_total_part(p_s);
        string file = edgefile + ".start";
        FILE* f = fopen(file.c_str(), "wb");
        assert(f != 0);
        fwrite(_s_start_edge, sizeof(index_t), size + 1, f);
        fclose(f);
        //cout << _s_start_edge[size  - 1] << endl;
        cout << "Total Edges again = "<< _s_start_edge[size] << endl;
    }
}

void grid::save_degree_files(string edgefile)
{
    string file = edgefile + ".degree";
    FILE* f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(vert_degree, sizeof(degree_t), vert_count, f);
    fclose(f);
	
    cout << "Compressing degree start" << endl;
	compress_degree();
	cout << "Compressing done" << endl;
    
	file = edgefile + ".sdegree";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(svert_degree, sizeof(sdegree_t), vert_count, f);
    fclose(f);
	
	file = edgefile + ".bdegree";
    f = fopen(file.c_str(), "wb");
    assert(f != 0);
    fwrite(bvert_degree, sizeof(bdegree_t), bdegree_count, f);
    cout << "saving done" << endl;
    fclose(f);
}

void grid::read_grid(string edgefile)
{
    string file = edgefile + ".grid";
    f_edge = open(file.c_str(), O_RDONLY);
    struct stat st_edge;
    fstat(f_edge, &st_edge);
    assert(st_edge.st_size != 0);
    
    //index_t nedges = st_edge.st_size/sizeof(edge_t);
    _edges = (edge_t*)mmap(0, st_edge.st_size, PROT_READ,
                                  MAP_PRIVATE, f_edge, 0);

    read_start(edgefile);
    
}

void grid::read_grid_in_mem(string edgefile)
{

    string file = edgefile + ".grid";
    struct stat st_edge;
    stat(file.c_str(), &st_edge);
    assert(st_edge.st_size != 0);
    
    FILE* f = fopen(file.c_str(), "rb");
    assert(f != 0);
    _edges = (edge_t*) malloc(st_edge.st_size);
    assert(_edges);
    fread(_edges, sizeof(edge_t), st_edge.st_size/sizeof(edge_t), f);
    fclose(f);
   
	read_start_in_mem(edgefile);
}

void grid::read_degree_in_mem(string edgefile)
{
    string file;
    FILE* f_degree;
    /*
    file = edgefile + ".sdegree";
    f_degree = fopen(file.c_str(), "rb");
    assert(f_degree != 0);
    svert_degree = (sdegree_t*) malloc(sizeof(sdegree_t)*vert_count);
    assert(svert_degree);
    fread(svert_degree, sizeof(sdegree_t), vert_count, f_degree);
    fclose(f_degree);
    
    file = edgefile + ".bdegree";
    struct stat st_edge;
    stat(file.c_str(), &st_edge);
    bdegree_count = st_edge.st_size/sizeof(bdegree_t);
	cout << "bdegree_count = " << bdegree_count << endl;
    

    f_degree = fopen(file.c_str(), "rb");
    assert(f_degree != 0);
    bvert_degree = (bdegree_t*) malloc(st_edge.st_size);
    assert(bvert_degree);
    fread(bvert_degree, sizeof(bdegree_t), st_edge.st_size/sizeof(bdegree_t), f_degree);
    fclose(f_degree);
	*/
    
    file = edgefile + ".degree";
    f_degree = fopen(file.c_str(), "rb");
    assert(f_degree != 0);
    vert_degree = (degree_t*)mmap(NULL, sizeof(degree_t)*vert_count, 
                           PROT_READ|PROT_WRITE,
                           MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, 0 , 0);
    if (vert_degree == MAP_FAILED) {
        vert_degree = (degree_t*) calloc(sizeof(degree_t), vert_count);
    }
    fread(vert_degree, sizeof(degree_t), vert_count, f_degree);
    fclose(f_degree);

}

void grid::read_degree(string edgefile)
{
	
	string file = edgefile + ".bdegree";
    int f_bdegree = open(file.c_str(), O_RDONLY);
    assert(f_bdegree != 0);
    struct stat st_degree;
    fstat(f_bdegree, &st_degree);
    assert(st_degree.st_size != 0);
    bvert_degree = (bdegree_t*)mmap(0, st_degree.st_size, PROT_READ,
                                  MAP_PRIVATE, f_bdegree, 0);

	
	file = edgefile + ".sdegree";
    f_bdegree = open(file.c_str(), O_RDONLY);
    assert(f_bdegree != 0);
    fstat(f_bdegree, &st_degree);
    assert(st_degree.st_size != 0);
    svert_degree = (sdegree_t*)mmap(0, st_degree.st_size, PROT_READ,
                                  MAP_PRIVATE, f_bdegree, 0);
	
/*	
    file = edgefile + ".degree";
    f_degree = fopen(file.c_str(), "rb");
    assert(f_degree != 0);
    vert_degree = (degree_t*) calloc(sizeof(degree_t), vert_count);
    fread(vert_degree, sizeof(degree_t), vert_count, f_degree);
    fclose(f_degree);
*/	
}

int grid::read_aio_init(string edgefile)
{
    //read_start_in_mem(edgefile);
    read_start(edgefile);

    string file = edgefile + ".grid";
    f_edge = open(file.c_str(), O_RDONLY|O_DIRECT);
    //f_edge = open(file.c_str(), O_RDONLY);
    
    struct stat st_edge;
    fstat(f_edge, &st_edge);
    assert(st_edge.st_size != 0);
    
    cache_size = total_size - (memory << 1);
    
    //count of number of edges
    index_t nedges = st_edge.st_size/sizeof(edge_t);
    
    mp_count = ceil((double)sizeof(edge_t)*nedges/(memory));
    cout << "edge_count = " << nedges << endl
         << "memory size = " << total_size << endl
         << "segment size = " << memory << endl
         << "Cache size = " << cache_size << endl
         << "memory partition count = " << mp_count << endl;

    seg1 = (segment*)malloc(sizeof(segment));
    seg2 = (segment*)malloc(sizeof(segment));
    seg3 = (segment*)malloc(sizeof(segment));
   
    //Allocate aligned memory 
   
    init_buf = (char*)mmap(NULL, total_size, PROT_READ|PROT_WRITE, 
                           MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, 
                           0 , 0);

    if (MAP_FAILED == init_buf) {
        if(posix_memalign((void**)&init_buf,2097152 , (total_size))) {
			perror("posix_memalign");
			return -1;
        }
    }

    if (mlock(init_buf, total_size) < 0) {
        handle_error("mlock Failed");
    }

    seg1->ctx_count = 0;
    seg2->ctx_count = 0;
    seg3->ctx_count = 0;
    
    seg1->meta = (part_meta_t*)calloc(sizeof(part_meta_t), 1024*1024);
    seg2->meta = (part_meta_t*)calloc(sizeof(part_meta_t), 1024*1024);
    seg3->meta = (part_meta_t*)calloc(sizeof(part_meta_t), 1024*1024);
    
    //seg1->buf = init_buf;
    //seg2->buf = init_buf + memory;
    seg1->buf = 0;
    seg2->buf = 0;
    seg3->buf = 0;
    
    cached_pool1 = (segment*)malloc(sizeof(segment));
    cached_pool2 = (segment*)malloc(sizeof(segment));
    
    cached_pool1->ctx_count = 0;
    cached_pool2->ctx_count = 0;
    
    cached_pool1->meta = (part_meta_t*)calloc(sizeof(part_meta_t), 1024*1024*16);
    cached_pool2->meta = (part_meta_t*)calloc(sizeof(part_meta_t), 1024*1024*16);
    
    cached_pool1->buf = init_buf;// + 2*memory;
    cached_pool2->buf = init_buf;// + 2*memory;

	main_buf = init_buf;
    
    return 0;
}


//This can be done now multi-threaded with some changes.
/*
 * @to_read:     Memory Size of the segment
 * @start_addr:  offset in the memory segement.
 * @part_meta:   Partition/tile metadata for one tile.
 * @last_read:   Start reading from this partition.
 *
 */
int cache_driver::prep_read_aio_seq(const last_read_t& last_read, size_t to_read, 
                  part_meta_t* part_meta, 
		  int& iteration_done, size_t& start_addr)
{
	#ifdef HALF_GRID
	matrix<spart_t, index_t> start_edge_half;
	start_edge_half.part_count = p_p;
	#endif

	matrix_f<spart_t, index_t> start_edge_full;
	start_edge_full.part_count = p_p;
	matrix<spart_t, index_t>* start_edge;
        
	#ifdef HALF_GRID
	if (last_read.b_i == last_read.b_j) {
        start_edge = &start_edge_half;
        start_edge->val = g->_s_start_edge 
					+ beg_edge_offset1(last_read.b_i);
    } else {
		start_edge = &start_edge_full;
		start_edge->val = g->_s_start_edge 
			+ beg_edge_offset2(last_read.b_i, last_read.b_j);
	}
	#else 
	start_edge = &start_edge_full;
	start_edge->val = g->_s_start_edge 
			+ beg_edge_offset2(last_read.b_i, last_read.b_j);
	#endif
    start_edge->part_count = p_p;

    int     k = 0;
    //size_t  
    //start_addr = 0;
    bool found_start = false;

    index_t start = 0;
    index_t b_i = (last_read.b_i << bit_shift3);
    index_t b_j = (last_read.b_j << bit_shift3);
    index_t offset = (calc_index_f(last_read.b_i, last_read.b_j, p) << bit_shift4);

    index_t i_start = last_read.s_i;
    index_t j_start = last_read.s_j;
    
    iteration_done = 1;

    //Also track the continuity,
    //When buffers separate, start address should be aligned.
    for (part_t i = i_start ; i < p_p; ++i) {
		
		//End the ctx if we started.
        /*if ( found_start && 
             ((start_edge->get_count(i, j_start) << bytes_in_edge_shift) 
					+ start_addr > read_sz)
           )*/ 
        //if(found_start && (i & 0x1F)) 
        if(found_start) 
        {
			found_start = false;
			part_id end = calc_index_f_opt(i - 1, p_p - 1, bit_shift3) + offset; 
			part_meta[k].end = end;
			start_addr = UPPER_ALIGN(start_addr);
			++k;
		}
        
        for (part_t j = j_start; j < p_p; ++j) {
            //caculate the start and end index
            start = calc_index_f_opt(i, j, bit_shift3) + offset; 
            bool found = false;
            //Search it in cached_part.
            start_index = search_cached_part(start, found);
            if (found) {
                if(found_start) {
                    found_start = false;
                    //We may have to end if we have started.
                    end_ctx(part_meta, k, last_read.b_i, last_read.b_j, i, j);
                    start_addr = UPPER_ALIGN(start_addr);
                    ++k;
                }
                continue;
            }

            //Did we start, no => start it now.
            if (found_start == false) {
                //If at least we can read one partition ahead.
                if ((start_edge->get_count(i, j) << bytes_in_edge_shift) 
                    + start_addr > to_read) {
                    return k;
                }

                found_start = true;
				part_meta[k].start = start;
                //part_meta[k].i = i;
                //part_meta[k].j = j;
                part_meta[k].offset = start_addr;
                start_addr += ((start_edge->get(i, j) 
								<< bytes_in_edge_shift) & 0x1FF);
            }

            if ((start_edge->get_count(i, j) << bytes_in_edge_shift) 
					+ start_addr <= to_read) {
                start_addr += (start_edge->get_count(i, j) << bytes_in_edge_shift);
            } else {
                //Memory is full now.
                //Update the end positions. and return.
                end_ctx(part_meta, k, last_read.b_i, last_read.b_j, i, j);
                ++k;
                start_addr = UPPER_ALIGN(start_addr);
                return k;
            }
        }
		j_start = 0;
		#ifdef HALF_GRID
        if(b_i == b_j) {
			j_start = i + 1;
		}	
		#endif
    }

    if (true == found_start) {
        found_start = false;
        part_id end = calc_index_f_opt(p_p - 1, p_p - 1, bit_shift3) + offset; 
        part_meta[k].end = end;
        start_addr = ((start_addr + 511) & ALIGN_MASK);
        ++k;
    }
    //XXX
    iteration_done = 2;
    return k;
}

void end_ctx(part_meta_t* part_meta, int k, part_t b_i, part_t b_j, spart_t i, spart_t j)
{
	part_id offset = ((calc_index_f(b_i, b_j, p)) << bit_shift4);
	part_id end = 0;
   
	if (j == 0) {
		end = calc_index_f_opt(i-1, p_p-1, bit_shift3) + offset;
	} else {
		end = calc_index_f_opt(i, j-1, bit_shift3) + offset;
	}
   
	#ifdef HALF_GRID	
	if(b_i == b_j) {
		if (j == i) {
			end = calc_index_f_opt(i-1, p_p-1, bit_shift3) + offset;
		} else {
			end = calc_index_f_opt(i, j-1, bit_shift3) + offset;
		}
    }
	#endif

	part_meta[k].end = end;
}

//This can be done now multi-threaded with some changes.
int cache_driver::prep_read_aio_random(const last_read_t& last_read, 
				       size_t to_read, bitmap_t* read_part, 
				       part_meta_t* part_meta,
				       int& iteration_done, 
				       size_t& start_addr)
{
	#ifdef HALF_GRID
	matrix<spart_t, index_t> start_edge_half;
	start_edge_half.part_count = p_p;
	#endif

	matrix_f<spart_t, index_t> start_edge_full;
	start_edge_full.part_count = p_p;
    matrix<spart_t, index_t>* start_edge;

	#ifdef HALF_GRID
    if (last_read.b_i == last_read.b_j) {
        start_edge = &start_edge_half;
	    start_edge->val = g->_s_start_edge + beg_edge_offset1(last_read.b_i);
    } else {
        start_edge = &start_edge_full;
	    start_edge->val = g->_s_start_edge 
		    + beg_edge_offset2(last_read.b_i, last_read.b_j);
    }
	#else
    start_edge = &start_edge_full;
	start_edge->val = g->_s_start_edge 
		    + beg_edge_offset2(last_read.b_i, last_read.b_j);
	#endif

    start_edge->part_count = p_p;

    int     k = 0;
    //start_addr = 0;
    bool found_start = false;

    index_t b_i = (last_read.b_i << bit_shift3);
    index_t b_j = (last_read.b_j << bit_shift3);
    part_id offset = (calc_index_f(last_read.b_i, last_read.b_j, p) << bit_shift4);
    
	index_t i_start = last_read.s_i;
    index_t j_start = last_read.s_j;
    index_t start = 0;

    //Also track the continuity,
    //When buffers separate, start address should be aligned.
    index_t i1 = b_i + i_start;
    for (part_t i = i_start; i < p_p ; ++i, ++i1) {
        bool done = false;
		
		//End the ctx if we started.
		if(found_start) {
			found_start = false;
			part_id end = calc_index_f_opt(i - 1, p_p - 1, bit_shift3) + offset; 
			part_meta[k].end = end;
			start_addr = UPPER_ALIGN(start_addr);
			++k;
		}
        
		if (read_part->get_bit(i1)) {
            
			for (part_t j = j_start; j < p_p; ++j) {
                //caculate the start and end index
                start = calc_index_f_opt(i, j, bit_shift3) + offset; 
		        index_t start1 = start_edge->get_index(i, j);
                
				//Search it in cached_part.
                bool found = false;
                start_index = search_cached_part(start, found);
                
				if (found) {
                    if(found_start) {
                        //Check if end of this partition and end of next partitions are 
                        //within 512 bytes.
                        index_t end1 = start_edge->get_end(start1);
                        index_t end2 = start_edge->get_end(start1 + 1);
                        if (LOWER_ALIGN(end2 << bytes_in_edge_shift) 
                            - UPPER_ALIGN(end1 << bytes_in_edge_shift)
                            < 512) {
                            start_addr += (start_edge->get_count(i, j) 
                                           << bytes_in_edge_shift);
                            continue;
                        }

                        found_start = false;
                        end_ctx(part_meta, k, last_read.b_i, last_read.b_j, i, j);
                        start_addr = UPPER_ALIGN(start_addr);
                        ++k;
                    }
                    continue;
                }

                //Did we start, no => start it now.
                if (found_start == false) {
                    //If at least we can read one partition ahead.
                    if ((start_edge->get_count(i, j) << bytes_in_edge_shift) 
                        + start_addr > to_read) {
			            iteration_done = 1;
                        return k;
                    }

                    found_start = true;
                    //part_meta[k].i = i;
                    //part_meta[k].j = j;
					part_meta[k].start = start;
                    part_meta[k].offset = start_addr;
                    start_addr += ((start_edge->get(i, j) 
                                    << bytes_in_edge_shift) & 0x1FF);
                }

                if ((start_edge->get_count(i, j) << bytes_in_edge_shift) 
						+ start_addr <= to_read) {
                    start_addr += (start_edge->get_count(i, j) 
                                   << bytes_in_edge_shift);
                } else {
                    //Memory is full now.
                    //Update the end positions. and return.
                    end_ctx(part_meta, k, last_read.b_i, last_read.b_j, i, j);
                    ++k;
		            iteration_done = 1;
                    return k;
                }
            }
            done = true;
        }

        if (false == done) {
            
            //Did we start? complete the end
            if (true == found_start) {
                found_start = false;
				//end_ctx    
				part_id end = calc_index_f_opt(i-1, p_p-1, bit_shift3) + offset; 
				part_meta[k].end = end;
				assert(part_meta[k].end >= part_meta[k].start);
				//part_meta[k].i_end = i -1;
                //part_meta[k].j_end = p_p -1;
		        //assert(part_meta[k].i <= part_meta[k].i_end);
		        //assert(part_meta[k].j <= part_meta[k].j_end);
				
				start_addr = ((start_addr + 511) & ALIGN_MASK);
                ++k;
            }
            #ifdef HALF_GRID 
            if (prep_read_aio_col(last_read, i, read_part, part_meta, to_read,
                                start_edge, start_addr, k, found_start)) {
		    
		        iteration_done = 1;
	       	    return k;
	        }
            #endif
        }
	        
		j_start = 0;
		#ifdef HALF_GRID
        if (b_i == b_j) {
	        j_start = i + 1;
	    }
		#endif
    }

    if (true == found_start) {
        found_start = false;
		
		//end_ctx
		part_id end = calc_index_f_opt(p_p-1, p_p-1, bit_shift3) + offset; 
		part_meta[k].end = end;
		assert(part_meta[k].end >= part_meta[k].start);
        /*
		part_meta[k].i_end = p_p -1;
        part_meta[k].j_end = p_p -1;
	    assert(part_meta[k].i <= part_meta[k].i_end);
	    assert(part_meta[k].j <= part_meta[k].j_end);
		*/
		start_addr = ((start_addr + 511) & ALIGN_MASK);
        ++k;
    }

    //We are still left with some memory. informa the caller.
    iteration_done = 2;
    return k;
}

int 
cache_driver::prep_read_aio_col(const last_read_t& last_read, part_t i, 
		   bitmap_t* read_part, part_meta_t* part_meta, 
                   size_t to_read, matrix<spart_t, index_t>* start_edge,
		   size_t& start_addr,  int& k, bool& found_start)
{
    index_t offset = (calc_index_f(last_read.b_i, last_read.b_j, p) << bit_shift4);
    
    index_t j1 = (last_read.b_j << bit_shift3);
    spart_t j = 0;
  
    if (last_read.b_i == last_read.b_j) {
		j1 += i + 1;
        j = i + 1;	
    } else {
		j = 0;
    }
    
    index_t start = 0;
    
    //fetch the columns this time, start from i+1 (obvious).
    for (; j < p_p; ++j, ++j1) {
        start = calc_index_f(i, j, p_p) + offset; 
		index_t start1 = start_edge->get_index(i, j);
        bool done = false; 
        
		if (read_part->get_bit(j1)) {
            //caculate the start index
            bool found = false;
            //Search it in cached_part.
            start_index = search_cached_part(start, found);
        
			if (found) {
                //We may have to end if we have started.
                if(true == found_start) {
                    //Check if end of this partition and end of next partitions are 
                    //within 512 bytes.
                    index_t end1 = start_edge->get_end(start1);
                    index_t end2 = start_edge->get_end(start1 + 1);
                    if (LOWER_ALIGN(end2 << bytes_in_edge_shift) 
                        - UPPER_ALIGN(end1 << bytes_in_edge_shift) < 512) {
                        start_addr += (start_edge->get_count(i, j) 
				       << bytes_in_edge_shift);
                        continue;
                    }

                    found_start = false;
                    
					//end_ctx
					part_id end = calc_index_f(i, j-1, p_p) + offset; 
					part_meta[k].end = end;
					assert(part_meta[k].end >= part_meta[k].start);
					/*
					part_meta[k].j_end = j - 1;
                    part_meta[k].i_end = i;
					assert(part_meta[k].i <= part_meta[k].i_end);
					assert(part_meta[k].j <= part_meta[k].j_end);
                    */
					start_addr = UPPER_ALIGN(start_addr);
                    ++k;
                }
                continue;
            }

            //did we start? no => start again.
            if (found_start == false) {
                if ((start_edge->get_count(i,j) << bytes_in_edge_shift) 
                    + start_addr > to_read) {
                    return k;
                }

                found_start = true;
				part_meta[k].start = start;
                /*
				part_meta[k].i = i;
                part_meta[k].j = j;
                */
                part_meta[k].offset = start_addr;
				start_addr += ((start_edge->get(i,j) << bytes_in_edge_shift) & 0x1FF);
            }

            //Can we fetch the next column with available memory?
            if ((start_edge->get_count(i, j) << bytes_in_edge_shift) 
					+ start_addr <= to_read) {
                start_addr += (start_edge->get_count(i,j) << bytes_in_edge_shift);

            } else {
                //Memory is full now. Finish and return. end_ctx
				part_id end = calc_index_f(i, j-1, p_p) + offset; 
                if (end <= part_meta[k].start) {
                    return k; 
                } else {
					part_meta[k].end = end;
					assert(part_meta[k].end >= part_meta[k].start);
                    //part_meta[k].i_end = i;
                    //part_meta[k].j_end = j - 1;
                }
				//assert(part_meta[k].i <= part_meta[k].i_end);
				//assert(part_meta[k].j <= part_meta[k].j_end);
                ++k; 
                return k;
            }
            done = true;
        }

        if (false == done) {
            
			//Did we start? complete the end
            if (true == found_start) {
                //Check if end of this partition and end of next partitions are 
                //within 512 bytes.
                index_t end1 = start_edge->get_end(start1);
                index_t end2 = start_edge->get_end(start1 + 1);
                if (LOWER_ALIGN(end2 << bytes_in_edge_shift) 
                    - UPPER_ALIGN(end1 << bytes_in_edge_shift) < 512) {
                    start_addr += (start_edge->get_count(i, j) << bytes_in_edge_shift);
                    continue;
                }

                found_start = false;
                //end_ctx
				part_id end = calc_index_f(i, j-1, p_p) + offset; 
				part_meta[k].end = end;
				assert(part_meta[k].end >= part_meta[k].start);
				
				/*part_meta[k].j_end = j - 1;
                part_meta[k].i_end = i;
				assert(part_meta[k].i <= part_meta[k].i_end);
				assert(part_meta[k].j <= part_meta[k].j_end);
				*/
                
				start_addr = UPPER_ALIGN(start_addr);
                ++k;
            }
        }
    }

    if (true == found_start) {
		found_start = false;
		//end_ctx
		part_id end = calc_index_f(i, p_p-1, p_p) + offset; 
		part_meta[k].end = end;
		assert(part_meta[k].end >= part_meta[k].start);
		/*part_meta[k].j_end = p_p - 1;
		part_meta[k].i_end = i;
		assert(part_meta[k].i <= part_meta[k].i_end);
		assert(part_meta[k].j <= part_meta[k].j_end);
		*/
		start_addr = UPPER_ALIGN(start_addr);
		++k;
    }
    return 0;
}


io_driver::io_driver()
{
    aio_meta = new aio_meta_t [IO_THDS];
    
    for (int j = 0; j < IO_THDS; ++j) {
        aio_meta[j].events = new struct io_event [AIO_MAXIO];
        aio_meta[j].cb_list = new struct iocb*[AIO_MAXIO];
        aio_meta[j].ctx = 0;
        for(index_t i = 0; i < AIO_MAXIO; ++i) {	
            aio_meta[j].cb_list[i] = new struct iocb;
        }
        if(io_setup(AIO_MAXIO, &aio_meta[j].ctx) < 0) {
            cout << AIO_MAXIO << endl;
            perror("io_setup");
            assert(0);
        }
        aio_meta[j].busy = 0;
    }
}

size_t io_driver::read_aio_random(segment* seg)
{
    index_t ctx_count = seg->ctx_count;
    if (0 == ctx_count) return 0;
    
    char** pbuf = &seg->buf;

    if (__sync_bool_compare_and_swap(pbuf, 0, main_buf))
    {
        //seg->buf = main_buf;
        main_buf += memory;
        assert(main_buf <= init_buf + total_size);
    }
/*
	cout << seg->meta[0].start << " ";
	cout << seg->meta[ctx_count - 1].end << " ";
	cout << ctx_count << endl;
*/  
    //double start = mywtime();
    
    size_t  start_address = 0;
    size_t  sz_to_read = 0;
    size_t disk_offset = 0; 
    int ret = 0;
    
    #ifdef HALF_GRID
    matrix<spart_t, index_t> start_edge_half;
	start_edge_half.part_count = p_p;
    #endif
	matrix_f<spart_t, index_t> start_edge_full;
	start_edge_full.part_count = p_p;
	matrix<spart_t, index_t>* start_edge;
	part_t b_i;
	part_t b_j;
	spart_t s_i;
	spart_t s_j;
	spart_t s_i_end;
	spart_t s_j_end;
  
    int tid          = omp_get_thread_num();
    vertex_t per_thd_work = ctx_count/IO_THDS;
    vertex_t my_beg       = tid * per_thd_work;
    vertex_t my_end       = (tid + 1) * per_thd_work;
    if (tid == IO_THDS - 1 ) my_end= ctx_count;
    
    index_t k1= 0;
    index_t k = 0;
    
    //cout << omp_get_thread_num() << endl;

    //for (index_t k = 0; k < ctx_count; ++k) 
    for (k = my_beg; k < my_end; ++k) 
    {
		get_ij(seg->meta[k].start, b_i, b_j, s_i, s_j);
		get_s_ij(seg->meta[k].end, s_i_end, s_j_end);

		
		#ifdef HALF_GRID
		if ((b_i == b_j)) {
            start_edge = &start_edge_half;
            start_edge->val = g->_s_start_edge 
					+ beg_edge_offset1(b_i);
		} else {
			start_edge = &start_edge_full;
			start_edge->val = g->_s_start_edge + beg_edge_offset2(b_i, b_j);
		}
		#else
		start_edge = &start_edge_full;
		start_edge->val = g->_s_start_edge + beg_edge_offset2(b_i, b_j);
		#endif

		disk_offset = ((start_edge->get(s_i, s_j) << bytes_in_edge_shift) 
					   & ALIGN_MASK);
			
		//XXX
        start_address = seg->meta[k].offset;
			
		sz_to_read  = (start_edge->get_end(s_i_end, s_j_end) 
						   << bytes_in_edge_shift)  
						  - disk_offset;

		sz_to_read = ((sz_to_read + 511) & ALIGN_MASK);
		assert(sz_to_read < (size_t) 0x0FFFFFFFFF);

		//cout << "sz_to_read" << sz_to_read << endl;
				
		io_prep_pread(aio_meta[tid].cb_list[k1++], 
						  f_edge, 
						  (char*)seg->buf + start_address, 
						  sz_to_read, 
						  disk_offset);
			
		//start_address += sz_to_read;
        
        if (0 == ((k1) % AIO_BATCHIO)) {
            ret = io_submit(aio_meta[tid].ctx, AIO_BATCHIO, 
                            aio_meta[tid].cb_list+ k1  - AIO_BATCHIO);
            
            if (ret != AIO_BATCHIO) {
                cout << ret  << endl;
                perror("io_submit");
                assert(0);
            }
            aio_meta[tid].busy += AIO_BATCHIO;
        }

        if (aio_meta[tid].busy == AIO_MAXIO) {
            ret = io_getevents(aio_meta[tid].ctx, AIO_BATCHIO, 
                               AIO_BATCHIO, aio_meta[tid].events, 0);

            if (ret != AIO_BATCHIO) {
                cout << ctx_count << " " << ret << endl;
                perror(" io_getevents");
                assert(0);
            }
            aio_meta[tid].busy -= AIO_BATCHIO;
            k1 = 0;
        }
    }
    int rem_aio = (k1 % AIO_BATCHIO);
    if (rem_aio)
    {
        ret = io_submit(aio_meta[tid].ctx, rem_aio, 
                        aio_meta[tid].cb_list + k1 - rem_aio);
        
        if (ret != rem_aio) {
            cout << ret <<" "  << rem_aio << endl;
            perror("io_submit");
            assert(0);
        }
        aio_meta[tid].busy += rem_aio;
    }

    //double end = mywtime();
    //assert(events[0].obj == cb_list[0]);
    return 0;
}

int io_driver::wait_aio_completion()
{
    uint32_t tid = omp_get_thread_num();
    int ret = 0;
    if (aio_meta[tid].busy) {
        ret = io_getevents(aio_meta[tid].ctx, aio_meta[tid].busy, 
                           aio_meta[tid].busy, aio_meta[tid].events, 0);

        if (ret != aio_meta[tid].busy) {
            cout << aio_meta[tid].busy << " " << ret << endl;
            perror(" io_getevents");
        }
    }

    aio_meta[tid].busy = 0;
    return 0;
}

void grid::read_start_in_mem(string edgefile)
{
    string file = edgefile + ".start";
    FILE* f_start = fopen(file.c_str(), "rb");
    
    struct stat st_edge;
    stat(file.c_str(), &st_edge);
    assert(st_edge.st_size != 0);
    
    _s_start_edge = (index_t*)malloc(st_edge.st_size);
    assert(_s_start_edge);
    fread(_s_start_edge, sizeof(index_t), st_edge.st_size/sizeof(index_t), f_start);
    fclose(f_start);
    //assert(edge_start.val[edge_start.size] == nedges);
}

void grid::read_start(string edgefile)
{
    string file = edgefile + ".start";
    int f_start = open(file.c_str(), O_RDONLY);
    
    struct stat st_edge;
    fstat(f_start, &st_edge);
    assert(st_edge.st_size != 0);
    
    _s_start_edge = (index_t*)mmap(0, st_edge.st_size, PROT_READ,
                                  MAP_PRIVATE, f_start, 0);

    if (MAP_FAILED == _s_start_edge) {
            handle_error("start file mmap");
    }
    assert(MAP_FAILED != _s_start_edge);
    //assert(edge_start.val[edge_start.size] == nedges);
}

inline bool 
next_part(last_read_t& last_read) 
{
    ++last_read.s_j; 
    last_read.s_j = (last_read.s_j & part_mask3_2);

    if (last_read.s_j == 0) {
		++last_read.s_i;
		last_read.s_i = (last_read.s_i & part_mask3_2); 
    }

    if ((last_read.s_i == 0) && (last_read.s_j == 0)) {
		++last_read.b_j;
		last_read.b_j = (last_read.b_j % p);
		
		if (last_read.b_j == 0) {
			++last_read.b_i;
			last_read.b_i = (last_read.b_i % p);
			#ifdef HALF_GRID
			last_read.b_j = last_read.b_i;
			#endif
		}
    }
	
	#ifdef HALF_GRID
    if (last_read.b_j == last_read.b_i) {
        if (last_read.s_j == 0) {
			last_read.s_j = last_read.s_i;
		}
    }
	#endif
    return true;
}
size_t cache_driver::fetch_mem_part(last_read_t& last_read, segment* seg, 
				    bitmap_t* read_part, int random /*= 0*/)
{
    //double start = mywtime();
    size_t sz_to_read = 0;
    seg->ctx_count = 0;
    next_part(last_read);
    
    //cout << last_read.b_i << " " << last_read.b_j << " " << last_read.s_i << " " << last_read.s_j << " ";
    
    while(true) {
		part_meta_t* s_meta = seg->meta + seg->ctx_count; 
		index_t s_ctx_count = 0;
		int iteration_done = 0;
		
		if (0 == random) {
			s_ctx_count = prep_read_aio_random(last_read, memory, read_part, s_meta,
				iteration_done, sz_to_read);
		} else {
			s_ctx_count = prep_read_aio_seq(last_read, memory, s_meta, 
						   iteration_done, sz_to_read);
		} 
		
		seg->ctx_count += s_ctx_count;

		if (iteration_done != 2) {
			if (s_ctx_count !=0) {
				get_s_ij(s_meta[s_ctx_count - 1].end, last_read.s_i, last_read.s_j);
			}
		} else {
			last_read.s_i = p_p - 1;
			last_read.s_j = p_p - 1;
		}

		//If we are done or buf got full, exit.
		if (is_end(last_read) || iteration_done == 1) {
			break; 
		}		
	
		//Fetch one physical group at a time.
		if (2 == random) {
			break;
		}	

		//Do we still have more space.
		if (iteration_done == 2) {
			next_part(last_read);
		}
    } 
    //double end = mywtime();
    return sz_to_read;
}


//Search from start and returns the upper bound of k and found.
index_t cache_driver::search_cached_part(index_t search, bool& found)
{
    found = false;
    //return -1;
    for (index_t k = start_index; k < ctx_count2; ++k ) {
        if (search < part_cached[k].start) return k;
        if (search <= part_cached[k].end) {
            found = true;
            return k;
        }
    }
    
    return -1;
}

//Call it at the begining of iteration.
void cache_driver::prep_pread_cached(segment* seg, index_t l)
{
    //#pragma omp for
    for (index_t k = 0; k < seg->ctx_count; ++k) {
		part_cached[l + k].start = seg->meta[k].start;
		part_cached[l + k].end = seg->meta[k].end;
	}
    return;
}
