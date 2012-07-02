/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*-----------------------------------------------------------------------------*/
/*                                                		               */
/*     File: queues_3d.h                         	     	               */
/*                                     		                               */
/*-----------------------------------------------------------------------------*/
/*       Data structures for FIFO, FILO queues and GRAPH                       */

#include   <math.h>

/* This value is enough for 128*128*128 images, but not 256*256*256 images */
#define MAX_VERTEX    32768

#define DYNAMIC_REG   200
#define MERGING_REG   100
#define MINIMUN_SIZE  5
enum
{
  OFF_SET = -32768
};
#define IS_OFF_SET(voxel_val)	/*((voxel_val) >= OFF_SET)*/((voxel_val) < 0)
#define UNCERTAINTY   2

#define WSHED         0
#define MASK          1
#define INIT          2
#define INIT_VALUE    0

/* definitions _selection in watershed_3d.c*/
#define CB_MERGE_WTSHED 47
#define MANUAL_MARKER_WSHED 48
#define AUTOMAT_MARKER_WSHED 49
#define MIN_MARKERS_WSHED 50
#define FLAT_ZONES_WSHED 51
#define DYN_MERGE_WTSHED 52
#define WATERFALLS 53
#define WATERSHED_3D 54
#define WATERSHED_SIMPLE 55


/*-------------   Structures for the queue nodes   -------------------*/

struct tpinfo {
       int  coord_x;
       int  coord_y;
       int  coord_z;
       long grey_level;
};
typedef struct tpinfo TPINFO;

struct node {
       TPINFO 	   infolist;
       struct node *next;
};
typedef struct node NODE;
typedef NODE *NODEPTR;



/*-------------     Structures for FIFO queues   -------------------- */

typedef struct {		
          NODEPTR head, tail;
} FIFO;



/* ----------------   Structures for the graph   -------------------- */

/*
    the graphs here are used for 3D watersheds

    I guess the vertices (this is the name for nodes in graphs) are valleys
    segmented by the watershed, and the arcs... err... seem to be a list of
    voxels in each region close to the other one.

    It seems that the voxels in arcs are actually regarded to be not part of
    the vertices themselve, i.e. when we merge a voxel from arcs into the
    vertex, we increment num_pixel...
*/

#define NEW_ARCS 1

#define HASH_STAT 1

struct tpvertex{
       int 		region_num;
       long		size_vertex;
       struct tpvertex	*next_vertex;
};
typedef struct tpvertex TPVERTEX;
typedef TPVERTEX	*VERTEX;

#if NEW_ARCS
typedef struct NeighbourNode
{
  int region_num;               /* id of the neighbour vertex */
  NODEPTR arc;                  /* list of voxels on the boundary of the
                                current vertex with the neighbour vertex */
  struct NeighbourNode *next;	/* next node in list */
  struct NeighbourNode *prev;	/* previous node in list */

  struct NeighbourNode *inverse;/* points to the node which connect the
                                same vertices in inverse order, i.e. if node
                                represents arc(i->j), node.inverse represents
                                arc(j->i) */
  /* NOTA : the code is designed so that
  find_neighbour_node(graph, j, i)->inverse is a valid alternative to
  find_neighbour_node(graph, i, j), i.e. if find_neighbour_node(graph, i, j)
  is non-NULL, find_neighbour_node(graph, j, i) is non-NULL, too. */
} NeighbourNode;
#endif

struct tpinfovertex {
	long 		num_pixels;
	long		grey_region;		/* mean of the grey value in the vertex ? */
	double		sqr_grey_region;	/* standard deviation of the grey value in the vertex ? */
	long		grey_min;		/* grey value at the minimum in the vertex ? */
	long		dynamic;		/* ??? */
	VERTEX	 	merged_vertex;		/* list of vertices we want to merge with the current one, as found in various merging algorithms ? */
	int	 	vertex_merged;		/* number of the vertex this vertex is to be merged with ? */
	NODEPTR		neighbourg_vertex;	/* list of adjacent vertices ? (used only in dynamic merging) */
	NODEPTR		voxels_vertex;		/* list of voxels in this vertex */
#if NEW_ARCS
	NeighbourNode	*real_neighbour_vertices;	/* list of adjacent vertices */
#endif
};

typedef struct tpinfovertex TPINFO_VERTEX;

#if !NEW_ARCS

struct tpgraph{
       NODEPTR 		arcs[MAX_VERTEX][MAX_VERTEX];
       TPINFO_VERTEX	vertex_field[MAX_VERTEX];
};

#else

/* we use a hash table to try to accelerate searching for an arc list or
determining whether 2 vertices are close */
/*
  the hash function is :
    f(x,y) = x[1-HASH_BITS] ^ rotate(y[1-HASH_BITS], HASH_SHIFT)
  where a[1-HASH_BITS] is the word made of the HASH_BITS LSBits of a,
  and rotate(x,y) is x bit-rotated y times (since HASH_SHIFT = HASH_BITS/2,
  I did not bother to define the rotation direction)

  The table will be 2^HASH_BITS * sizeof(hash_node) byte-long.
  (sizeof(hash_node) is typically 12 on 32-bit machines, 24 on 64-bit ones)
*/
#define HASH_SHIFT 10
#define HASH_BITS (2 * HASH_SHIFT)
#define HASH_SIZE (1 << HASH_BITS)
#define HASH_MASK (HASH_SIZE - 1)

typedef struct hash_node
{
	int vertex;
        int neighbour_vertex;
        NeighbourNode *the_node;
} hash_node;

/* this implementation assumes that neighbour_vertex is always positive
(i.e. it shifts neighbour_vertex right without bothering to convert it to
the unsigned type of the same size or masking the end result) */
/* we may want to fix this if we ever get rid of the OFF_SET stuff */
#define hash_function(vertex, neighbour_vertex)				\
  ((vertex ^ (((neighbour_vertex) << HASH_SHIFT)			\
                | ((neighbour_vertex) >> HASH_SHIFT))) & HASH_MASK)

#if HASH_STAT
void reset_hash_stats(void);
void print_hash_stats(void);
#endif

struct tpgraph{
       TPINFO_VERTEX	vertex_field[MAX_VERTEX];
       hash_node	arc_hash_table[HASH_SIZE];
};

#endif

typedef struct tpgraph TPGRAPH;
typedef TPGRAPH *GRAPH;

#if !NEW_ARCS

#define find_arc(the_graph, vertex, neighbour_vertex)			\
	((the_graph)->arcs[(vertex)][(neighbour_vertex)])

#else

extern NeighbourNode *find_neighbour_node(GRAPH graph, int vertex, int neighbour_vertex);
NeighbourNode *find_or_create_neighbour_node(GRAPH graph, int vertex, int neighbour_vertex);
NODEPTR find_arc(GRAPH graph, int vertex, int neighbour_vertex);

#endif

#if !NEW_ARCS

#define are_vertex_connected(the_graph, vertex, neighbour_vertex)	\
	(find_arc(the_graph, vertex, neighbour_vertex) != NULL)

#else

#define are_vertex_connected(the_graph, vertex, neighbour_vertex)	\
	(find_neighbour_node(the_graph, vertex, neighbour_vertex) != NULL)

#endif

#if !NEW_ARCS

#define add_arc(the_graph, vertex, neighbour_vertex, data)		\
	arc_graph(&((the_graph)->arcs[(vertex)][(neighbour_vertex)]), (data));

#else

#if 1

void add_arc(GRAPH the_graph, int vertex, int neighbour_vertex, TPINFO data);

#else

/* This macro does not work, because cc does not accept nested preprocessor
arguments.  This results into one of the weirdest error message I have ever
read : 'unexpected symbol: "32768"' */

#define add_arc(the_graph, vertex, neighbour_vertex, data)		\
{									\
  NeighbourNode *r;							\
									\
  r = find_or_create_neighbour_node((the_graph), (vertex), (neighbour_vertex)); \
  arc_graph(& (r->arc), (data));					\
}

#endif

#endif

/*------------- Structures used by the manual merging  -------------*/

typedef struct {
       int		merger_region;
       int		num_objects;	
       GRAPH 		graph_field;
       grphic3d		*imres_input, *imres_output, *matrix_regions;       
       VERTEX		vertex_field;
}TPMERGING;
typedef TPMERGING	*MERGING_PTR;


/*------------- Structures used by the manual marker selection -----------*/

typedef struct {
       int       current_region;
       NODEPTR	 selected;
       FIFO 	 OQ[256];  /* ordered queue*/
       GRAPH	 graph_markers;
       grphic3d	 *im_gradient, *img_regions, *img_input;
}TPMANUAL_MARKERS;

typedef TPMANUAL_MARKERS *MANUAL_MARKERS_PTR;


/*------------- Structures used by the manual merging&marking  -------------*/

typedef struct {
       int		merger_region;
       int		num_objects;	
       GRAPH 		graph_field;
       grphic3d		*imres_input, *imres_output, *matrix_regions, 
       			*imregionsOQ;       
       VERTEX		vertex_field;
       FIFO 	 	merge_mark_OQ[256];  /* ordered queue*/       
}TPMERGING_MARKING;
typedef TPMERGING_MARKING	*MERGING_MARKING_PTR;


/* ------------- Structures for array_graph  -------------*/

typedef struct
{	int 	  index;
	double    contrast;
}ARRAY_GRAPH;


typedef struct
{	int 	  indexation;
	long      frequency;
	VERTEX	  regions_field;
}ARRAY_VERTEX;

/*Structure usually used by the histogram*/
typedef struct
{	int 	grey_level;
	double  frequency;
} HISTO_ARRAY;


 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern void create_p_fifo_3d(FIFO *fifo) ;
extern void insert_p_fifo_3d(FIFO *fifo, TPINFO info) ;
extern void tidy_insert_fifo_3d(FIFO *fifo, TPINFO info) ;
extern int fifo_empty_3d(FIFO fifo) ;
extern void free_p_fifo_3d(FIFO *fifo) ;
extern TPINFO fifo_first_3d(FIFO *fifo) ;
extern void initialize_graph_3d(GRAPH *graph, int min) ;
extern void free_graph(GRAPH *graph) ;
extern void reset2_graph_3d(GRAPH *graph) ;
extern void reset2_graph_3d_vertex_list(GRAPH *graph, VERTEX *vertex_list, grphic3d *o, grphic3d *imres) ;
extern void reset2_graph_histogram_3d(GRAPH *graph, VERTEX *vertex_list, ARRAY_VERTEX *histogram_regions) ;
extern void vertex_graph(GRAPH *graph, int region, long int grey) ;
extern void merge_graph_3d(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc) ;
extern void merge_graph_3d_histogram(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc, ARRAY_VERTEX *histogram_regions) ;
extern void merge_graph_3d_vertex_list(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc, VERTEX vertex_list) ;
extern void arc_graph(NODEPTR *list, TPINFO info) ;
extern int graph_empty_3d(GRAPH *graph, int region) ;
extern void initialize_array_graph(GRAPH graph, ARRAY_GRAPH *A, int condition) ;
extern void create_list_3d(NODEPTR *list) ;
extern void insert_list_3d(NODEPTR *list, TPINFO info) ;
extern TPINFO delete_first_list_3d(NODEPTR *list) ;
extern void free_memory_list_3d(NODEPTR *list) ;
extern unsigned long calculate_mean_3d(NODEPTR list) ;
extern int list_empty_3d(NODEPTR list) ;
extern void tidy_insert_list(NODEPTR *list, TPINFO info) ;
extern void create_vertex_list(VERTEX *vertex_list) ;
extern void free_vertex_list(VERTEX *vertex_list) ;
extern void insert_list_vertex_list(VERTEX *vertex_list, int region, long int size, VERTEX next) ;
extern void insert_vertex_list(VERTEX *vertex_list, int region, long int size) ;
extern void tidy_insert_vertex_list(VERTEX *vertex_list, int region, long int size) ;
extern void initialize_vertex_list_graph(GRAPH *graph, VERTEX *vertex_list, grphic3d *o, int *label, grphic3d *im1resc) ;
extern void initialize_array_vertex_graph(GRAPH *graph, grphic3d *o, int *label, grphic3d *im1resc, ARRAY_VERTEX *histogram_regions) ;
extern void sort_array3d(HISTO_ARRAY *A) ;
extern void sort_array_vertex(ARRAY_VERTEX *A) ;
extern int look_for_merger(GRAPH graph, int region) ;
extern void delete_vertex_histogram(GRAPH graph, ARRAY_VERTEX *histogram_regions, int index, int region) ;
extern void sort_array_graph(ARRAY_GRAPH *A) ;
extern void copy_array_vertex(ARRAY_VERTEX *A, ARRAY_VERTEX *A_indexed) ;
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
