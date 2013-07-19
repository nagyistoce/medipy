/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!----------------------------------------------------------------------------*
**                                                       		                   *
**      \file: queues_3d.c                              	     	               *
**      \brief Description : gestion de listes          		                   *
**                                                                             *
** --------------------------------------------------------------------------- */

#include <config.h>
#include <stdio.h>
#include <math.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_lang.h"
#include "noyau/imx_3d.h"

#include "outils/queues_3d.h"


#ifndef __DOXYGEN_SHOULD_SKIP_THIS__
/* Function prototypes */

/* FIFOs */
void create_p_fifo_3d(FIFO *fifo);
void insert_p_fifo_3d(FIFO *fifo, TPINFO info);
int fifo_empty_3d(FIFO fifo);
void tidy_insert_fifo_3d(FIFO *fifo, TPINFO info);
void free_p_fifo_3d(FIFO *fifo);
TPINFO fifo_first_3d(FIFO *fifo);

/* GRAPHs */
void initialize_array_graph(GRAPH graph, ARRAY_GRAPH *A, int condition);
void free_graph(GRAPH *graph);
void reset2_graph_3d(GRAPH *graph);
void vertex_graph(GRAPH *graph, int region, long int grey);
void merge_graph_3d(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc);
void merge_graph_3d_vertex_list(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc, VERTEX vertex_list);
void arc_graph(NODEPTR *list, TPINFO info);
int graph_empty_3d(GRAPH *graph, int region);
void sort_array_graph(ARRAY_GRAPH *A);

/* FILOs */
void create_list_3d(NODEPTR *list);
void insert_list_3d(NODEPTR *list, TPINFO info);
TPINFO delete_first_list_3d(NODEPTR *list);
void free_memory_list_3d(NODEPTR *list);
unsigned long calculate_mean_3d(NODEPTR list);
int list_empty_3d(NODEPTR list);
void tidy_insert_list(NODEPTR *list, TPINFO info);

/* FILOs' vertex_node */
void create_vertex_list(VERTEX *vertex_list);
void insert_list_vertex_list(VERTEX *vertex_list, int region, long int size, VERTEX next);
void insert_vertex_list(VERTEX *vertex_list, int region, long int size);
void tidy_insert_vertex_list(VERTEX *vertex_list, int region, long int size);


void initialize_vertex_list_graph(GRAPH *graph, VERTEX *vertex_list, grphic3d *o, int *label, grphic3d *im1resc);

void sort_array_vertex(ARRAY_VERTEX *A);
void delete_vertex_histogram(GRAPH graph, ARRAY_VERTEX *histogram_regions, int index, int region);
void sort_array_graph(ARRAY_GRAPH *A);
void copy_array_vertex(ARRAY_VERTEX *A, ARRAY_VERTEX *A_indexed);

void delete_vertex_histogram(GRAPH graph, ARRAY_VERTEX *histogram_regions, int index, int region);
#endif /*__DOXYGEN_SHOULD_SKIP_THIS__*/


/* Function definition */

#if NEW_ARCS

#if HASH_STAT

static long search_count;
static long collision_count;

void reset_hash_stats(void)
{
  search_count = collision_count = 0; 
}

void print_hash_stats(void)
{
  fprintf(stderr, "searches : %ld collisions : %ld\n", search_count, collision_count);
}

#endif

/*
  Find the neighbour list node for neighbour_vertex in the
  real_neighbour_vertices list of vertex in GRAPH.

  Return a reference to this node.  Return NULL if node does not exist.
*/
NeighbourNode *find_neighbour_node(GRAPH graph, int vertex, int neighbour_vertex)
{
  hash_node *p;
  NeighbourNode *reply;

  p = & graph->arc_hash_table[hash_function(vertex, neighbour_vertex)];

#if HASH_STAT
  if ((p->vertex != -1) && (p->neighbour_vertex != -1))
  {
    search_count++;

    if ((p->vertex != vertex) || (p->neighbour_vertex != neighbour_vertex))
      collision_count++;
  }
#endif

  if ((p->vertex == vertex) && (p->neighbour_vertex == neighbour_vertex))
  {
    reply = p->the_node;
  }
  else
  {
    reply = graph->vertex_field[vertex].real_neighbour_vertices;

    while (reply && (reply->region_num != neighbour_vertex))
      reply = reply->next;

    p->vertex = vertex;
    p->neighbour_vertex = neighbour_vertex;
    p->the_node = reply;
  }

  return reply;
}

/*
  Find the neighbour list node for neighbour_vertex in the
  real_neighbour_vertices list of vertex in GRAPH.  Create a new node if it
  does not exist.

  Return a reference to this node.  Return NULL if node does not exist, and
  malloc failed.

  Only legal way to alloc a neighbour node.
*/
NeighbourNode *find_or_create_neighbour_node(GRAPH graph, int vertex, int neighbour_vertex)
{
  hash_node *p;
  NeighbourNode *reply;

  p = & graph->arc_hash_table[hash_function(vertex, neighbour_vertex)];

#if HASH_STAT
  if ((p->vertex != -1) && (p->neighbour_vertex != -1))
  {
    search_count++;

    if ((p->vertex != vertex) || (p->neighbour_vertex != neighbour_vertex))
      collision_count++;
  }
#endif

  if ((p->vertex == vertex) && (p->neighbour_vertex == neighbour_vertex))
  {
    reply = p->the_node;
  }
  else
  {
    reply = graph->vertex_field[vertex].real_neighbour_vertices;

    while (reply && (reply->region_num != neighbour_vertex))
      reply = reply->next;

    p->vertex = vertex;
    p->neighbour_vertex = neighbour_vertex;
    p->the_node = reply;
  }

  if (reply == NULL)
  {
    /* we need to allocate both the current node and the inverse node */
    reply = malloc(sizeof(NeighbourNode));
    if (! reply)
    {
      fprintf(stderr, "malloc failed (file %s, line %d)!\n", __FILE__, __LINE__);
      return NULL;
    }
    reply->inverse = malloc(sizeof(NeighbourNode));
    if (! reply->inverse)
    {
      free(reply);
      fprintf(stderr, "malloc failed (file %s, line %d)!\n", __FILE__, __LINE__);
      return NULL;
    }
    reply->region_num = neighbour_vertex;
    reply->arc = NULL;
    reply->next = graph->vertex_field[vertex].real_neighbour_vertices;
    reply->prev = NULL;
    reply->inverse->region_num = vertex;
    reply->inverse->arc = NULL;
    reply->inverse->next = graph->vertex_field[neighbour_vertex].real_neighbour_vertices;
    reply->inverse->prev = NULL;
    reply->inverse->inverse = reply;

    if (graph->vertex_field[vertex].real_neighbour_vertices)
      graph->vertex_field[vertex].real_neighbour_vertices->prev = reply;
    graph->vertex_field[vertex].real_neighbour_vertices = reply;

    if (graph->vertex_field[neighbour_vertex].real_neighbour_vertices)
      graph->vertex_field[neighbour_vertex].real_neighbour_vertices->prev = reply->inverse;
    graph->vertex_field[neighbour_vertex].real_neighbour_vertices = reply->inverse;

    /* always put the new node in hash table, correcting possible outdated data */
    p->vertex = vertex;
    p->neighbour_vertex = neighbour_vertex;
    p->the_node = reply;

    /* do not put the inverse node in hash table, since it was not directly
    requested, but check for outdated data and correct it */
    p = & graph->arc_hash_table[hash_function(neighbour_vertex, vertex)];
    if ((p->vertex == neighbour_vertex) && (p->neighbour_vertex == vertex))
    {
      p->the_node = reply->inverse;
    }
  }

  return reply;
}

/*
  Free the neighbour list node the_node linking neighbour_vertex to vertex in
  GRAPH.
*/

/*
  Does one half of the job - frees one node without freeing the inverse node
*/
void free_neighbour_node_internal(GRAPH graph, int vertex, int neighbour_vertex,
                          		NeighbourNode *the_node)
{
  hash_node *p;
  NODEPTR q;

  p = & graph->arc_hash_table[hash_function(vertex, neighbour_vertex)];
  if ((p->vertex == vertex) && (p->neighbour_vertex == neighbour_vertex))
  {
    p->the_node = NULL;
  }

  if (the_node->prev)
    the_node->prev->next = the_node->next;
  else
    graph->vertex_field[vertex].real_neighbour_vertices = the_node->next;

  while ((q = the_node->arc) != NULL)
  {
    the_node->arc = q->next;
    free(q);
  }

  free(the_node);
}

/*
  Calls free_neighbour_node_internal for the_node and its inverse
*/
void free_neighbour_node(GRAPH graph, int vertex, int neighbour_vertex,
                          NeighbourNode *the_node)
{
  free_neighbour_node_internal(graph, neighbour_vertex, vertex, the_node->inverse);
  free_neighbour_node_internal(graph, vertex, neighbour_vertex, the_node);
}

/*
  Add a new voxel to an arc.  Alloc the arc if needed.
*/
void add_arc(GRAPH the_graph, int vertex, int neighbour_vertex, TPINFO data)
{
  NeighbourNode *r;

  r = find_or_create_neighbour_node((the_graph), (vertex), (neighbour_vertex));
  /*if (! r)
    return 1;*/
  arc_graph(& (r->arc), (data));
}


/*
  The voxel list associated with an arc.
*/
NODEPTR find_arc(GRAPH graph, int vertex, int neighbour_vertex)
{
  NeighbourNode *r;

  r = find_neighbour_node(graph, vertex, neighbour_vertex);

  return r ? (r->arc) : NULL;
}

/*
  Remove a voxel from an arc.  The ownership of the storage is passed to the
  caller, i.e. you should either free the returned node or keep the reference
  somewhere (presumably in another vixel list) for later use.
*/
NODEPTR remove_arc(GRAPH the_graph, int vertex, int neighbour_vertex, NeighbourNode *the_node)
{
  NODEPTR reply;

  reply = the_node->arc;

  if (reply)
  {
    the_node->arc = reply->next;
    /*reply->next = NULL;*/	/* cleaner this way, but actually not indispensable */

    if ((! the_node->arc) && (! the_node->inverse->arc))
      free_neighbour_node(the_graph, vertex, neighbour_vertex, the_node);
  }

  return reply;
}

#endif


/************* FIFO queues' functions  ******************** */
/* FIFO  - first input first output linked_list library     */

/*-----------------   create_p_fifo_3d	--------------------*/
/*  Initialise the fifo					    */
/*----------------------------------------------------------*/

void create_p_fifo_3d(FIFO *fifo)
{
  (*fifo).head = NULL;
  (*fifo).tail = NULL;  
}

/*------------------ insert_p_fifo_3d ----------------------*/
/*  Insert a new node in the fifo			    */
/*----------------------------------------------------------*/

void insert_p_fifo_3d(FIFO *fifo, TPINFO info)
{
    NODEPTR newone;	

    newone = (NODEPTR) malloc(sizeof(NODE)); 
    newone->infolist.coord_x = info.coord_x;
    newone->infolist.coord_y = info.coord_y;
    newone->infolist.coord_z = info.coord_z;
    newone->infolist.grey_level = info.grey_level;
    newone->next = NULL;

    if  (fifo_empty_3d(*fifo)){
	(*fifo).head = newone;
	(*fifo).tail = newone;
    } else {
	(*fifo).tail->next =  newone;
	(*fifo).tail = newone;
    }
}




/*----------------------------------------------------------------------*/
/*    tidy_insert_fifo							*/
/* insert a node in fifolist according to the info.grey_level		*/
/*----------------------------------------------------------------------*/

void tidy_insert_fifo_3d(FIFO *fifo, TPINFO info)
{
    NODEPTR newone;	
    int end;
    NODEPTR p;
  
    newone = (NODEPTR) malloc(sizeof(NODE)); 
    newone->infolist.coord_x = info.coord_x;
    newone->infolist.coord_y = info.coord_y;
    newone->infolist.coord_z = info.coord_z;
    newone->infolist.grey_level = info.grey_level;
    newone->next = NULL;

    if  (fifo_empty_3d(*fifo)){
	(*fifo).head = newone;
	(*fifo).tail = newone;
    } else {
      p= (*fifo).head;		
      if (p->infolist.grey_level >= newone->infolist.grey_level){
 	  newone->next = p;
	  (*fifo).head = newone;  
	  end=1; 
      	}
       else{
         end=0;
	 while (end != 1){
	   if (p->next == NULL) 
	      end=1;
	     else
              if (p->next->infolist.grey_level >= newone->infolist.grey_level)
	          end=1;
		else
	  	  p=p->next;
	 }
	 newone->next = p->next;
	 p->next=newone;
       }
     }     
}


/*------------------------------------------------------------*/
/*  fifo_empty_3d					      */
/*  returns true if fifo is empty			      */
/*------------------------------------------------------------*/

int fifo_empty_3d(FIFO fifo)
{
   return ((fifo.head == NULL) && (fifo.tail == NULL));
}
 
/*-------------------------------------------------------------*/
/*  free_p_fifo_3d					       */
/*  Release the memory allocated by the fifo		       */
/*-------------------------------------------------------------*/

void free_p_fifo_3d(FIFO *fifo)
{
    NODEPTR p;
    
  if  (fifo_empty_3d(*fifo) == FALSE){
    while ((*fifo).head != (*fifo).tail) {
	p= (*fifo).head;
	(*fifo).head = (*fifo).head->next;
	free(p);
    }
    p=(*fifo).head;
   (*fifo).head=(*fifo).tail->next;
   (*fifo).tail=(*fifo).tail->next;
   free(p);
  } 
}


/*-------------------------------------------------------------*/
/*  fifo_first_3d					       */
/*  removes a pixel from the head of the fifo stack			       */
/*-------------------------------------------------------------*/

TPINFO fifo_first_3d(FIFO *fifo)
{

   TPINFO p_info;
   NODEPTR p;
		
    if  (fifo_empty_3d(*fifo)){
      fprintf(stderr,"error\n");exit(1);}   
   else{
      p= (*fifo).head;
      p_info = (*fifo).head->infolist;
      if ((*fifo).head == (*fifo).tail){
         (*fifo).head = NULL; 
	 (*fifo).tail = NULL;
       }
       else
        (*fifo).head =(*fifo).head->next;
    free(p);	 
    return p_info; 
   }
}


/************* graph functions  ******************************* */


/*--------------------------------------------------------------*/
/*  initialize_graph_3d						*/
/*--------------------------------------------------------------*/
void initialize_graph_3d(GRAPH *graph, int min)
{
  int i;
  GRAPH aux;

   (*graph)=NULL;
   aux = (GRAPH) malloc(sizeof(TPGRAPH));   
   for (i= 0; i< MAX_VERTEX; i++){
    aux->vertex_field[i].num_pixels = 0;
    aux->vertex_field[i].grey_region = 0;
    aux->vertex_field[i].sqr_grey_region = 0.0;
    aux->vertex_field[i].grey_min = min;
    if (i==(INIT+1))
    	aux->vertex_field[i].dynamic = 255;
      else	
    	aux->vertex_field[i].dynamic = -1;
    aux->vertex_field[i].merged_vertex = NULL;    
    aux->vertex_field[i].vertex_merged = 0;
    aux->vertex_field[i].neighbourg_vertex =NULL;
    aux->vertex_field[i].voxels_vertex =NULL;
#if !NEW_ARCS
    for (j= 0; j< MAX_VERTEX; j++) 
	aux->arcs[i][j]=NULL;	
#else
    aux->vertex_field[i].real_neighbour_vertices = NULL;
#endif
   }

#if NEW_ARCS
  for (i=0; i<HASH_SIZE; i++)
  {
    aux->arc_hash_table[i].vertex = -1;
    aux->arc_hash_table[i].neighbour_vertex = -1;
    aux->arc_hash_table[i].the_node = NULL;
  }
#endif
 
 (*graph)=aux;
}

/*--------------------------------------------------------------*/
/*  free_graph							*/
/*--------------------------------------------------------------*/
void free_graph(GRAPH *graph)
{
  int i;
  NODEPTR p;  
  VERTEX vertex;
      

  for (i= 0; i< MAX_VERTEX; i++){
   if ((*graph)->vertex_field[i].merged_vertex != NULL){
    while( (*graph)->vertex_field[i].merged_vertex != NULL){
      vertex = (*graph)->vertex_field[i].merged_vertex;        
      (*graph)->vertex_field[i].merged_vertex= (*graph)->vertex_field[i].merged_vertex->next_vertex;
      free(vertex); 
     }	
    }
   }
  

#if !NEW_ARCS
  for (i= 0; i< MAX_VERTEX; i++)
   for (j= 0; j< MAX_VERTEX; j++){
    if ((*graph)->arcs[i][j] != NULL){
     while((*graph)->arcs[i][j] != NULL){     
       p=(*graph)->arcs[i][j];
       (*graph)->arcs[i][j]= (*graph)->arcs[i][j]->next; 
       free(p);
      }
     }
   }
#endif
   
  for (i= 0; i< MAX_VERTEX; i++){
    if ((*graph)->vertex_field[i].voxels_vertex != NULL){
    while( (*graph)->vertex_field[i].voxels_vertex != NULL){
      p = (*graph)->vertex_field[i].voxels_vertex;        
      (*graph)->vertex_field[i].voxels_vertex=(*graph)->vertex_field[i].voxels_vertex->next;
      free(p); 
     }	
    }
   }


  for (i= 0; i< MAX_VERTEX; i++){
    if ((*graph)->vertex_field[i].neighbourg_vertex != NULL){
    while( (*graph)->vertex_field[i].neighbourg_vertex != NULL){
      p = (*graph)->vertex_field[i].neighbourg_vertex;        
      (*graph)->vertex_field[i].neighbourg_vertex=(*graph)->vertex_field[i].neighbourg_vertex->next;
      free(p); 
     }	
    }
   }

#if NEW_ARCS
  for (i= 0; i< MAX_VERTEX; i++)
  {
    NeighbourNode *q;
    while( (q = (*graph)->vertex_field[i].real_neighbour_vertices) != NULL)
    {
      (*graph)->vertex_field[i].real_neighbour_vertices=q->next;
      while ((p = q->arc) != NULL)
      {
        q->arc = p->next;
        free(p);
      }
      free(q);
    }
  }
#endif
}


/*------------  reset2_graph_3d ----------------------------- */    

void reset2_graph_3d(GRAPH *graph)
{
  int i;
  long total_grey, num_merged;
  VERTEX  vertex;  
  double total_sqr_grey;
    
   for (i= 0; i< MAX_VERTEX; i++){
    total_grey= (*graph)->vertex_field[i].grey_region * (*graph)->vertex_field[i].num_pixels;   
    total_sqr_grey =(double)((*graph)->vertex_field[i].sqr_grey_region * (*graph)->vertex_field[i].num_pixels);   
    num_merged = (*graph)->vertex_field[i].num_pixels;
    if ((*graph)->vertex_field[i].merged_vertex != NULL){   
      vertex = (*graph)->vertex_field[i].merged_vertex; 
      while (vertex != NULL){
        total_grey= total_grey +((*graph)->vertex_field[vertex->region_num].grey_region
				* vertex->size_vertex);
	total_sqr_grey = (double)(total_sqr_grey +((*graph)->vertex_field[vertex->region_num].sqr_grey_region
				* vertex->size_vertex));
        num_merged= num_merged + vertex->size_vertex;
	vertex = vertex->next_vertex; 
      }
     (*graph)->vertex_field[i].grey_region = (long int)((double)total_grey/num_merged);
     (*graph)->vertex_field[i].sqr_grey_region = (double)(total_sqr_grey/num_merged);
     (*graph)->vertex_field[i].num_pixels = num_merged;
     
      while((*graph)->vertex_field[i].merged_vertex != NULL){
        vertex = (*graph)->vertex_field[i].merged_vertex;
	(*graph)->vertex_field[vertex->region_num].grey_region = (*graph)->vertex_field[i].grey_region;
	(*graph)->vertex_field[vertex->region_num].num_pixels = num_merged;
	(*graph)->vertex_field[vertex->region_num].vertex_merged = i;
	(*graph)->vertex_field[i].merged_vertex= (*graph)->vertex_field[i].merged_vertex->next_vertex;
        free(vertex); 
       }	      
   }  

  }       
}


/*------------  reset2_graph_3d_vertex_list ----------------------------- */    

void reset2_graph_3d_vertex_list(GRAPH *graph, VERTEX *vertex_list, grphic3d *o, grphic3d *imres)
{
  long total_grey, num_merged;
  VERTEX  vertex, aux_list;  
  double total_sqr_grey;
  /*NODEPTR p;*/
  

  while (*vertex_list != NULL){
    aux_list= *vertex_list;
    total_grey= (*graph)->vertex_field[aux_list->region_num].grey_region * (*graph)->vertex_field[aux_list->region_num].num_pixels;   
    total_sqr_grey =(double)((*graph)->vertex_field[aux_list->region_num].sqr_grey_region * (*graph)->vertex_field[aux_list->region_num].num_pixels);   
    num_merged = (*graph)->vertex_field[aux_list->region_num].num_pixels;
    if ((*graph)->vertex_field[aux_list->region_num].merged_vertex != NULL){   
      vertex = (*graph)->vertex_field[aux_list->region_num].merged_vertex; 
      while (vertex != NULL){
        total_grey= total_grey +((*graph)->vertex_field[vertex->region_num].grey_region
				* vertex->size_vertex);
	total_sqr_grey = (double)(total_sqr_grey +((*graph)->vertex_field[vertex->region_num].sqr_grey_region
				* vertex->size_vertex));
        num_merged= num_merged + vertex->size_vertex;
	vertex = vertex->next_vertex; 
      }
     (*graph)->vertex_field[aux_list->region_num].grey_region = (long int)((double)total_grey/num_merged);
     (*graph)->vertex_field[aux_list->region_num].sqr_grey_region = (double)(total_sqr_grey/num_merged);
     (*graph)->vertex_field[aux_list->region_num].num_pixels = num_merged;
     
      while((*graph)->vertex_field[aux_list->region_num].merged_vertex != NULL){
        vertex = (*graph)->vertex_field[aux_list->region_num].merged_vertex;
	(*graph)->vertex_field[vertex->region_num].grey_region = (*graph)->vertex_field[aux_list->region_num].grey_region;
	(*graph)->vertex_field[vertex->region_num].num_pixels = num_merged;/* deberia de ser 0 ?? */
	(*graph)->vertex_field[vertex->region_num].vertex_merged = aux_list->region_num;
	
	/* ii)
	while((*graph)->vertex_field[vertex->region_num].voxels_vertex != NULL){
	 p= (*graph)->vertex_field[vertex->region_num].voxels_vertex;
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= aux_list->region_num;
	 imres->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= 
	 				(*graph)->vertex_field[aux_list->region_num].grey_region;
	 insert_list_3d(&((*graph)->vertex_field[aux_list->region_num].voxels_vertex), p->infolist); 
	 (*graph)->vertex_field[vertex->region_num].voxels_vertex= p->next;
	 free(p);	
	 }
	 (ii*/
	 
	(*graph)->vertex_field[aux_list->region_num].merged_vertex= (*graph)->vertex_field[aux_list->region_num].merged_vertex->next_vertex;
        free(vertex); 
       }	      
   }

   *vertex_list= (*vertex_list)->next_vertex;
   free(aux_list);
  }       
}



/*------------  reset2_graph_histogram_3d ----------------------------- */    

void reset2_graph_histogram_3d(GRAPH *graph, VERTEX *vertex_list, ARRAY_VERTEX *histogram_regions)
{
  long total_grey, num_merged;
  VERTEX  vertex, aux_list;  
  double total_sqr_grey;  
  

  while (*vertex_list != NULL){
    aux_list= *vertex_list;
    total_grey= (*graph)->vertex_field[aux_list->region_num].grey_region * (*graph)->vertex_field[aux_list->region_num].num_pixels;   
    total_sqr_grey =(double)((*graph)->vertex_field[aux_list->region_num].sqr_grey_region * (*graph)->vertex_field[aux_list->region_num].num_pixels);   
    num_merged = (*graph)->vertex_field[aux_list->region_num].num_pixels;
    if ((*graph)->vertex_field[aux_list->region_num].merged_vertex != NULL){   
      vertex = (*graph)->vertex_field[aux_list->region_num].merged_vertex; 
      while (vertex != NULL){
        total_grey= total_grey +((*graph)->vertex_field[vertex->region_num].grey_region
				* vertex->size_vertex);
	total_sqr_grey = (double)(total_sqr_grey +((*graph)->vertex_field[vertex->region_num].sqr_grey_region
				* vertex->size_vertex));
        num_merged= num_merged + vertex->size_vertex;
	vertex = vertex->next_vertex; 
      }

      /* in ->size_vertex is stored the index_histogram where vertex is found*/ 
      
     delete_vertex_histogram(*graph, histogram_regions, aux_list->size_vertex, aux_list->region_num);

     (*graph)->vertex_field[aux_list->region_num].grey_region = (long int)((double)total_grey/num_merged);
     (*graph)->vertex_field[aux_list->region_num].sqr_grey_region = (double)(total_sqr_grey/num_merged);
     (*graph)->vertex_field[aux_list->region_num].num_pixels = num_merged;

      histogram_regions[(*graph)->vertex_field[aux_list->region_num].grey_region].frequency += num_merged; 
      tidy_insert_vertex_list(&(histogram_regions[(*graph)->vertex_field[aux_list->region_num].grey_region].regions_field),aux_list->region_num, num_merged);      
     
      while((*graph)->vertex_field[aux_list->region_num].merged_vertex != NULL){
        vertex = (*graph)->vertex_field[aux_list->region_num].merged_vertex;
	(*graph)->vertex_field[vertex->region_num].grey_region = (*graph)->vertex_field[aux_list->region_num].grey_region;
	(*graph)->vertex_field[vertex->region_num].num_pixels = num_merged;/* deberia de ser 0 ?? */
	(*graph)->vertex_field[vertex->region_num].vertex_merged = aux_list->region_num;
	(*graph)->vertex_field[aux_list->region_num].merged_vertex= (*graph)->vertex_field[aux_list->region_num].merged_vertex->next_vertex;
        free(vertex); 
       }	      
   }  
   
   *vertex_list= (*vertex_list)->next_vertex;
   free(aux_list);
  }       
}



/*----------------------------------------------------------------------*/
/*  vertex_graph							*/
/*----------------------------------------------------------------------*/

void vertex_graph(GRAPH *graph, int region, long int grey)
{
   (*graph)->vertex_field[region].num_pixels++;

   (*graph)->vertex_field[region].grey_region  = (*graph)->vertex_field[region].grey_region + grey;
   
   (*graph)->vertex_field[region].sqr_grey_region = (double)((*graph)->vertex_field[region].sqr_grey_region + (grey*grey));
}



/*--------------  merge_graph_3d --------------------------------------------*/    
/*- merging regiona into regionb, regiona with less num_pixels than regionb -*/
/*---------------------------------------------------------------------------*/

/* -> dynamic merge watersheds ? */

#if !NEW_ARCS

#define insert_arc(the_graph, vertex, neighbour_vertex, the_arc)	\
{									\
  the_arc->next = (*graph)->arcs[vertex][neighbour_vertex];		\
  (*graph)->arcs[vertex][neighbour_vertex] = the_arc;			\
}

#else

#if 0

#define insert_arc(the_graph, vertex, neighbour_vertex, the_arc)	\
{									\
  NeighbourNode *q;							\
									\
  q = find_or_create_neighbour_node(the_graph, vertex, neighbour_vertex); \
									\
  the_arc->next = q->arc;						\
									\
  q->arc = the_arc;							\
}

#else

static void insert_arc(GRAPH the_graph, int vertex, int neighbour_vertex, NODEPTR the_arc)
{
  NeighbourNode *q;

  q = find_or_create_neighbour_node(the_graph, vertex, neighbour_vertex);

  the_arc->next = q->arc;

  q->arc = the_arc;
}

#endif

#endif


void merge_graph_3d(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc)
{
#if !NEW_ARCS
  int i/*, reg_aux*/;
#endif
  long sizea, sizeb;
  long real_grey_merger,real_grey_merged;
  NODEPTR p;

#if NEW_ARCS
  NeighbourNode *q;

  NeighbourNode *r;
#endif

 if ( (are_vertex_connected(*graph, regiona, regionb) || are_vertex_connected(*graph, regionb, regiona)) 
 	&& (regiona > (INIT+1)) && (regionb > (INIT+1))){		

#if NEW_ARCS
  for (r = (*graph)->vertex_field[regiona].real_neighbour_vertices;
        r != NULL; r = r->next)
  {
    int i = r->region_num;
#else
  for (i= 0; i< MAX_VERTEX ; i++)
   if (i != regiona){
#endif
    if (are_vertex_connected(*graph, regiona, i)){
     if(i==regionb){
      /* pixels of regiona close to regionb are added to regionb */
#if NEW_ARCS
      q = r;
      while ((p = remove_arc(*graph, regiona, i, q)) != NULL)
      {
        (*graph)->vertex_field[regionb].num_pixels++;
        o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] = regionb;
        free(p);
      }
      //free_neighbour_node(*graph, regiona, i, q);
#else
      while((*graph)->arcs[regiona][i] != NULL){
	p=(*graph)->arcs[regiona][i];
        (*graph)->vertex_field[regionb].num_pixels++;
     	o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= regionb;
     	(*graph)->arcs[regiona][i]= (*graph)->arcs[regiona][i]->next; 
     	free(p);
      }
#endif
    }   
   else{
     /* pixels of regiona close to other regions are either added to regionb
     or to the other region */
#if NEW_ARCS
      q = r;
      while ((p = q->arc) != NULL)
      {
        q->arc = p->next;
#else
     while ((*graph)->arcs[regiona][i] != NULL){
      p=(*graph)->arcs[regiona][i]; 
      (*graph)->arcs[regiona][i]=(*graph)->arcs[regiona][i]->next;
#endif

      real_grey_merger=(*graph)->vertex_field[regionb].grey_region;
      real_grey_merged=(*graph)->vertex_field[i].grey_region;

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]-real_grey_merger) < 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]-real_grey_merged)){
         insert_arc(*graph, regionb, i, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb;
       }

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) > 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, i, regionb, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + i;	 
       }
     
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) == 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){    
      if ((*graph)->vertex_field[regionb].grey_min > (*graph)->vertex_field[i].grey_min){
          insert_arc(*graph, regionb, i, p); 
	  o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb;	  
        }
       else{
           insert_arc(*graph, i, regionb, p); 
	   o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + i;
         }
       }
                  
    }
   }   
  }  
  
  if (are_vertex_connected(*graph, i, regiona)){
   if(i==regionb){
      /* pixels of regionb close to regiona are added to regionb (-- what ???) */
#if NEW_ARCS
      q = r->inverse;
      while ((p = remove_arc(*graph, i, regiona, q)) != NULL)
      {
        (*graph)->vertex_field[regionb].num_pixels++;
        o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] = regionb;
        free(p);
      }
      //free_neighbour_node(*graph, i, regiona, q);
#else
      while((*graph)->arcs[i][regiona] != NULL){
	p=(*graph)->arcs[i][regiona];
        (*graph)->vertex_field[regionb].num_pixels++;
	
     	o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= regionb;
     	(*graph)->arcs[i][regiona]= (*graph)->arcs[i][regiona]->next;  
     	free(p);
      }
#endif
    }   
   else{ 
      /* pixels of other regions close to regiona are either added to regionb
      or to the other region */
#if NEW_ARCS
      q = r->inverse;
      while ((p = q->arc) != NULL)
      {
        q->arc = p->next;
#else
      while ((*graph)->arcs[i][regiona] != NULL){
      p=(*graph)->arcs[i][regiona]; 
      (*graph)->arcs[i][regiona]=(*graph)->arcs[i][regiona]->next;
#endif

      real_grey_merger=(*graph)->vertex_field[regionb].grey_region;
      real_grey_merged=(*graph)->vertex_field[i].grey_region;

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) < 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, regionb, i, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb; 
       }

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) > 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, i, regionb, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + i; 
       }
     
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) == 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){    
      if ((*graph)->vertex_field[regionb].grey_min > (*graph)->vertex_field[i].grey_min){
          insert_arc(*graph, regionb, i, p);
	  o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb; 
        }
       else{
           insert_arc(*graph, i, regionb, p);
	   o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + i; 
         }
       }
                   
     }
    } 
   }   
  } 
 
  sizea=(*graph)->vertex_field[regiona].num_pixels;
  sizeb=(*graph)->vertex_field[regionb].num_pixels;

  /* now the boundaries have been fixed, add the vertice to be merged into a
  to the list of vertice to be merged into b */
  insert_list_vertex_list(&((*graph)->vertex_field[regionb].merged_vertex), regiona, sizea, (*graph)->vertex_field[regiona].merged_vertex);
  (*graph)->vertex_field[regiona].merged_vertex =NULL; 

 }
}




/*--------------  merge_graph_3d_histogram ----------------------------------*/    
/*- merging regiona into regionb, regiona with less num_pixels than regionb -*/
/*---------------------------------------------------------------------------*/

/* -> CB merge watersheds (merging_catchment_basins_higher_probability) */

void merge_graph_3d_histogram(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc, ARRAY_VERTEX *histogram_regions)
{
  int j/*, reg_aux*/;
  long sizea, sizeb;
  NODEPTR p;
  VERTEX vertex;
  long real_grey_merger,real_grey_merged;  

#if NEW_ARCS
  NeighbourNode *q;

  NeighbourNode *r;
#endif

      
 if ( (are_vertex_connected(*graph, regiona, regionb) || are_vertex_connected(*graph, regionb, regiona)) 
 	&& (regiona > (INIT+1)) && (regionb > (INIT+1))){		
  
 for (j=0; j<= 255;j++){
  if (histogram_regions[j].regions_field != NULL){
   vertex= histogram_regions[j].regions_field;
   while (vertex!= NULL){
   if (vertex->region_num != regiona){
#if NEW_ARCS
    if ((r = find_neighbour_node(*graph, regiona, vertex->region_num)) != NULL){
#else
    if (are_vertex_connected(*graph, regiona, vertex->region_num)){
#endif
     if(vertex->region_num==regionb){
#if NEW_ARCS
      q = r;
      while ((p = remove_arc(*graph, regiona, vertex->region_num, q)) != NULL)
      {
        (*graph)->vertex_field[regionb].num_pixels++;
        o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] = regionb;
        free(p);
      }
      //free_neighbour_node(*graph, regiona, vertex->region_num, q);
#else
      while((*graph)->arcs[regiona][vertex->region_num] != NULL){
	p=(*graph)->arcs[regiona][vertex->region_num];
        (*graph)->vertex_field[regionb].num_pixels++;
     	o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= regionb;
     	(*graph)->arcs[regiona][vertex->region_num]= (*graph)->arcs[regiona][vertex->region_num]->next; 
     	free(p);
      }
#endif
    }   
   else{
#if NEW_ARCS
      q = r;
      while ((p = q->arc) != NULL)
      {
        q->arc = p->next;
#else
     while ((*graph)->arcs[regiona][vertex->region_num] != NULL){
      p=(*graph)->arcs[regiona][vertex->region_num]; 
      (*graph)->arcs[regiona][vertex->region_num]=(*graph)->arcs[regiona][vertex->region_num]->next;
#endif

      real_grey_merger=(*graph)->vertex_field[regionb].grey_region;
      real_grey_merged=(*graph)->vertex_field[vertex->region_num].grey_region;
      
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) < 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, regionb, vertex->region_num, p); 
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb;
       }

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) > 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]  - real_grey_merged)){
         insert_arc(*graph, vertex->region_num, regionb, p); 
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num;	 
       }
     
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) == 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){    
      if ((*graph)->vertex_field[regionb].grey_min > (*graph)->vertex_field[vertex->region_num].grey_min){
          insert_arc(*graph, regionb, vertex->region_num, p); 
	  o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb;	  
        }
       else{
           insert_arc(*graph, vertex->region_num, regionb, p); 
	   o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num;
         }
       }
                  
    }
   }   
  }  
  
  if (are_vertex_connected(*graph, vertex->region_num, regiona)){
   if(vertex->region_num==regionb){
#if NEW_ARCS
      q = r->inverse;
      while ((p = remove_arc(*graph, vertex->region_num, regiona, q)) != NULL)
      {
        (*graph)->vertex_field[regionb].num_pixels++;
        o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] = regionb;
        free(p);
      }
      //free_neighbour_node(*graph, vertex->region_num, regiona, q);
#else
      while((*graph)->arcs[vertex->region_num][regiona] != NULL){
	p=(*graph)->arcs[vertex->region_num][regiona];
        (*graph)->vertex_field[regionb].num_pixels++;
	
     	o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= regionb;
     	(*graph)->arcs[vertex->region_num][regiona]= (*graph)->arcs[vertex->region_num][regiona]->next;  
     	free(p);
      }
#endif
    }   
   else{ 
#if NEW_ARCS
      q = r->inverse;
      while ((p = q->arc) != NULL)
      {
        q->arc = p->next;
#else
      while ((*graph)->arcs[vertex->region_num][regiona] != NULL){
      p=(*graph)->arcs[vertex->region_num][regiona]; 
      (*graph)->arcs[vertex->region_num][regiona]=(*graph)->arcs[vertex->region_num][regiona]->next;
#endif

      real_grey_merger=(*graph)->vertex_field[regionb].grey_region;
      real_grey_merged=(*graph)->vertex_field[vertex->region_num].grey_region;

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) < 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, regionb, vertex->region_num, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb; 
       }

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) > 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]  - real_grey_merged)){
         insert_arc(*graph, vertex->region_num, regionb, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num; 
       }
     
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) == 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){    
      if ((*graph)->vertex_field[regionb].grey_min > (*graph)->vertex_field[vertex->region_num].grey_min){
          insert_arc(*graph, regionb, vertex->region_num, p);
	  o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb; 
        }
       else{
           insert_arc(*graph, vertex->region_num, regionb, p);
	   o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num; 
         }
        }
                   
       }
      } 
     }   
    } 
     vertex= vertex->next_vertex; 
    }/* while*/
   }
  } 
    sizea=(*graph)->vertex_field[regiona].num_pixels;
    sizeb=(*graph)->vertex_field[regionb].num_pixels;
    insert_list_vertex_list(&((*graph)->vertex_field[regionb].merged_vertex), regiona, sizea, (*graph)->vertex_field[regiona].merged_vertex);
    (*graph)->vertex_field[regiona].merged_vertex =NULL; 
  }   
}

    
/*--------------  merge_graph_3d_vertex_list ---------------------------------*/    
/*- merging regiona into regionb, regiona with less num_pixels than regionb -*/
/*---------------------------------------------------------------------------*/

/* -> CB merge watersheds (merging_catchment_basins_vertex) */
/* -> manual merge ??? */

void merge_graph_3d_vertex_list(GRAPH *graph, int regiona, int regionb, grphic3d *o, grphic3d *im1resc, VERTEX vertex_list)
{
  long sizea, sizeb;
  NODEPTR p;
  VERTEX vertex;
  long real_grey_merger,real_grey_merged;

#if NEW_ARCS
  NeighbourNode *q;

  NeighbourNode *r;
#endif

      
 if ( (are_vertex_connected(*graph, regiona, regionb) || are_vertex_connected(*graph, regionb, regiona)) 
 	&& (regiona > (INIT+1)) && (regionb > (INIT+1))){		  
  if (vertex_list != NULL){
   vertex= vertex_list;
   while (vertex!= NULL){
   if (vertex->region_num != regiona){
#if NEW_ARCS
    if ((r = find_neighbour_node(*graph, regiona, vertex->region_num)) != NULL){
#else
    if (are_vertex_connected(*graph, regiona, vertex->region_num)){
#endif
     if(vertex->region_num==regionb){
#if NEW_ARCS
      q = r;
      while ((p = remove_arc(*graph, regiona, vertex->region_num, q)) != NULL)
      {
        (*graph)->vertex_field[regionb].num_pixels++;
        o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] = regionb;
        free(p);
      }
      //free_neighbour_node(*graph, regiona, vertex->region_num, q);
#else
      while((*graph)->arcs[regiona][vertex->region_num] != NULL){
	p=(*graph)->arcs[regiona][vertex->region_num];
        (*graph)->vertex_field[regionb].num_pixels++;
     	o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= regionb;
     	(*graph)->arcs[regiona][vertex->region_num]= (*graph)->arcs[regiona][vertex->region_num]->next; 
     	free(p);
      }
#endif
    }   
   else{
#if NEW_ARCS
      q = r;
      while ((p = q->arc) != NULL)
      {
        q->arc = p->next;
#else
     while ((*graph)->arcs[regiona][vertex->region_num] != NULL){
      p=(*graph)->arcs[regiona][vertex->region_num]; 
      (*graph)->arcs[regiona][vertex->region_num]=(*graph)->arcs[regiona][vertex->region_num]->next;
#endif

      real_grey_merger=(*graph)->vertex_field[regionb].grey_region;
      real_grey_merged=(*graph)->vertex_field[vertex->region_num].grey_region;
      
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) < 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, regionb, vertex->region_num, p); 
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb;
       }

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) > 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, vertex->region_num, regionb, p); 
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num;	 
       }
     
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] -real_grey_merger) == 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]-real_grey_merged)){    
      if ((*graph)->vertex_field[regionb].grey_min > (*graph)->vertex_field[vertex->region_num].grey_min){
          insert_arc(*graph, regionb, vertex->region_num, p); 
	  o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb;	  
        }
       else{
           insert_arc(*graph, vertex->region_num, regionb, p); 
	   o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num;
         }
       }
                  
    }
   }   
  }  
  
  if (are_vertex_connected(*graph, vertex->region_num, regiona)){
   if(vertex->region_num==regionb){
#if NEW_ARCS
      q = r->inverse;
      while ((p = remove_arc(*graph, vertex->region_num, regiona, q)) != NULL)
      {
        (*graph)->vertex_field[regionb].num_pixels++;
        o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] = regionb;
        free(p);
      }
      //free_neighbour_node(*graph, vertex->region_num, regiona, q);
#else
      while((*graph)->arcs[vertex->region_num][regiona] != NULL){
	p=(*graph)->arcs[vertex->region_num][regiona];
        (*graph)->vertex_field[regionb].num_pixels++;
	
     	o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= regionb;
     	(*graph)->arcs[vertex->region_num][regiona]= (*graph)->arcs[vertex->region_num][regiona]->next;  
     	free(p);
      }
#endif
    }   
   else{ 
#if NEW_ARCS
      q = r->inverse;
      while ((p = q->arc) != NULL)
      {
        q->arc = p->next;
#else
      while ((*graph)->arcs[vertex->region_num][regiona] != NULL){
      p=(*graph)->arcs[vertex->region_num][regiona]; 
      (*graph)->arcs[vertex->region_num][regiona]=(*graph)->arcs[vertex->region_num][regiona]->next;
#endif

      real_grey_merger=(*graph)->vertex_field[regionb].grey_region;
      real_grey_merged=(*graph)->vertex_field[vertex->region_num].grey_region;

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) < 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){
         insert_arc(*graph, regionb, vertex->region_num, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb; 
       }

      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) > 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]  - real_grey_merged)){
         insert_arc(*graph, vertex->region_num, regionb, p);
	 o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num; 
       }
     
      if (labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merger) == 
	 labs(im1resc->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z] - real_grey_merged)){    
      if ((*graph)->vertex_field[regionb].grey_min > (*graph)->vertex_field[vertex->region_num].grey_min){
          insert_arc(*graph, regionb, vertex->region_num, p);
	  o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + regionb; 
        }
       else{
           insert_arc(*graph, vertex->region_num, regionb, p);
	   o->mri[p->infolist.coord_x][p->infolist.coord_y][p->infolist.coord_z]= OFF_SET + vertex->region_num; 
         }
        }
                   
       }
      } 
     }   
    } 
     vertex= vertex->next_vertex; 
    }/* while*/
   }
  
    sizea=(*graph)->vertex_field[regiona].num_pixels;
    sizeb=(*graph)->vertex_field[regionb].num_pixels;
    insert_list_vertex_list(&((*graph)->vertex_field[regionb].merged_vertex), regiona, sizea, (*graph)->vertex_field[regiona].merged_vertex);
    (*graph)->vertex_field[regiona].merged_vertex =NULL; 
  }   
}
   
      
/*----------------------------------------------------------------------*/
/*  arc_graph								*/
/*----------------------------------------------------------------------*/

/*
    insert a copy of 'info' into 'list'.  Supposed to be used only with the list
    of a graph.
*/

void arc_graph(NODEPTR *list, TPINFO info)
{
    NODEPTR newone;  

    newone = (NODEPTR) malloc(sizeof(NODE)); 
    newone->infolist.coord_x = info.coord_x;
    newone->infolist.coord_y = info.coord_y;
    newone->infolist.coord_z = info.coord_z;
    newone->infolist.grey_level = info.grey_level;
    newone->next = NULL;

    if (*list == NULL) { *list = newone; }
    else { newone->next = *list;
     	   *list = newone; }
}


/*----------------------------------------------------------------------*/
/*  graph_empty_3d							*/
/*  returns true if fifo is empty					*/
/*----------------------------------------------------------------------*/

int graph_empty_3d(GRAPH *graph, int region)
{
    return ((*graph)->vertex_field[region].num_pixels == 0) ;
}

 
/*----------------------------------------------------------------------*/
/*   initialize_array_graph						*/
/*----------------------------------------------------------------------*/

void initialize_array_graph(GRAPH graph, ARRAY_GRAPH *A, int condition)
{
    int i/*, j*/;
    /*ARRAY_GRAPH aux;*/

       for (i= 0; i< MAX_VERTEX; i++){
	 if (condition == DYNAMIC_REG)
	    A[i].contrast=(double)(graph->vertex_field[i].dynamic);
	 if (condition == MERGING_REG)	   
	   A[i].contrast=(double)(graph->vertex_field[i].num_pixels);
/*	   A[i].contrast=(double)(graph->vertex_field[i].grey_region);*/
	 A[i].index=i;
	}

}
 


/************* FILO queues' (list) functions  ******************** */
/* FILO queue's (list) = first input last output		  */

/*----------------------------------------------------------------------*/
/*  create_list_3d							*/
/*!	\brief Initialise the list							*/
/*!	\param list : la liste a initialiser */
/*! \remark list vaut NULL a la sortie  */
/*----------------------------------------------------------------------*/

void create_list_3d(NODEPTR *list)
{
    *list = NULL;
}


/*----------------------------------------------------------------------*/
/*  insert_list_3d								*/
/*! \brief Insert a new node in the list					*/
/*!	\param list : la liste	*/
/*! \param info : l'element a inserer */
/*! \remark l'element est integralement recopie dans la liste, il n'est pas simplement insere*/
/*----------------------------------------------------------------------*/

void insert_list_3d(NODEPTR *list, TPINFO info)
{
    NODEPTR newone;
  
    newone = (NODEPTR) malloc(sizeof(NODE)); 
    newone->infolist.coord_x = info.coord_x;
    newone->infolist.coord_y = info.coord_y;
    newone->infolist.coord_z = info.coord_z;
    newone->infolist.grey_level = info.grey_level;
    newone->next = NULL;

    if (*list == NULL) { *list = newone; }
    else { newone->next = *list;
     	   *list = newone; }
}


/*----------------------------------------------------------------------*/
/*  delete_first_list_3d							*/
/*! \brief  Delete and return the first node in the list			*/
/*! \param list : the list	(E/S)										*/
/*! \return the ancient first node								*/
/*----------------------------------------------------------------------*/

TPINFO delete_first_list_3d(NODEPTR *list)
{
    NODEPTR top;
    TPINFO deleted;

    if (*list != NULL) {
	top = *list;
	deleted.coord_x = top->infolist.coord_x;
	deleted.coord_y = top->infolist.coord_y;
	deleted.coord_z = top->infolist.coord_z;
	deleted.grey_level = top->infolist.grey_level;
	*list = (*list)->next;
	free(top);
    }
    return(deleted);
}

/*----------------------------------------------------------------------*/
/*  free_memory_list_3d							*/
/*!  Release the memory allocated by the list				*/
/*!  \param list : the list  							*/
/*----------------------------------------------------------------------*/

void free_memory_list_3d(NODEPTR *list)
{
    NODEPTR p;

    while (*list!=NULL) {
	p = *list;
	*list = (*list)->next;
	free(p);
    }
}



/*----------------------------------------------------------------------*/
/*  calculate_mean_3d							*/
/*!  Calculate the mean value of the pixels within the list		*/
/*  \param list : the list
	\retval the mean   */
/*----------------------------------------------------------------------*/
unsigned long calculate_mean_3d(NODEPTR list)
{
    NODEPTR p;
    long sum;
    int cont;
    
    sum= 0; 
    cont=0;
    p = list;
    if ( p!=NULL ) {	
     cont=1;
     while (p != NULL) {
	    sum+= p->infolist.grey_level;
	    cont++;
	    p = p->next;
	 }
      cont--;	 
      if (cont > 0)
       sum= (long)((double)((sum/cont)+0.5));
     } 

  return(sum);
}

/*----------------------------------------------------------------------*/
/*  list_empty_3d							*/
/*  Returns True if the list is empty, otherwise returns False		*/
/*----------------------------------------------------------------------*/

int list_empty_3d(NODEPTR list)
{
   return(list==NULL);
}


/*----------------------------------------------------------------------*/
/*    tidy_insert_list							*/
/* insert a node in list according to the info.grey_level		*/
/*----------------------------------------------------------------------*/

void tidy_insert_list(NODEPTR *list, TPINFO info)
{
    int /*i, j, k,*/ end;
    NODEPTR newone,p;
  

    newone = (NODEPTR) malloc(sizeof(NODE)); 
    newone->infolist.coord_x = info.coord_x;
    newone->infolist.coord_y = info.coord_y;
    newone->infolist.coord_z = info.coord_z;
    newone->infolist.grey_level = info.grey_level;
    newone->next = NULL;

    if (*list == NULL)
      *list=newone;
     else{
      p= *list;		
      if (p->infolist.grey_level > newone->infolist.grey_level){
 	  newone->next = p;
	  *list = newone;  
	  end=1; 
      	}
       else{
         end=0;
	 while (end != 1){
	   if (p->next == NULL) 
	      end=1;
	     else
              if (p->next->infolist.grey_level > newone->infolist.grey_level)
	          end=1;
		else
	  	  p=p->next;
	 }
	 newone->next = p->next;
	 p->next=newone;
       }
     }     
}



/************* vertex_list functions  ************************************/
/* a FILO queue of vertex						*/

/*----------------------------------------------------------------------*/
/*  create_vertex_list							*/
/*  Initialise the list							*/
/*----------------------------------------------------------------------*/
/*
  create an empty list
*/
void create_vertex_list(VERTEX *vertex_list)
{
    *vertex_list = NULL;
}


/*----------------------------------------------------------------------*/
/*  free_vertex_list							*/
/*  Release the memory allocated by the list				*/
/*----------------------------------------------------------------------*/

void free_vertex_list(VERTEX *vertex_list)
{
   VERTEX p;
   
   if (*vertex_list!=NULL) {
    while (*vertex_list!=NULL) {
	p = *vertex_list;
	*vertex_list = (*vertex_list)->next_vertex;
	free(p);
    }
  }  
}


/*----------------------------------------------------------------------*/
/*  insert_list_vertex_list						*/
/*  Insert a new vertex_node and its next address in a vertex_list	*/
/*----------------------------------------------------------------------*/
/*
  If next is NULL, behave like insert_vertex_list.
  Else, create a new node, append it to 'vertex_list', then append 'next'
  to the resulting list.
*/
void insert_list_vertex_list(VERTEX *vertex_list, int region, long int size, VERTEX next)
{
   VERTEX vertex, aux;  

   vertex = (VERTEX) malloc(sizeof(TPVERTEX));        
   vertex->region_num  = region;
   vertex->size_vertex = size;
   vertex->next_vertex = next; 

    if (*vertex_list == NULL) { *vertex_list = vertex; }
    else {
      if (next == NULL){    
	  vertex->next_vertex = *vertex_list;
	}
      else{
        aux=vertex;
        while (aux->next_vertex != NULL)  /* look for the last vertex*/
         aux=aux->next_vertex;
        aux->next_vertex = *vertex_list;
       }   
      *vertex_list = vertex;    
     }
}

/*----------------------------------------------------------------------*/
/*  insert_vertex_list						*/
/*  Insert a new vertex_node in a vertex_list				*/
/*----------------------------------------------------------------------*/
void insert_vertex_list(VERTEX *vertex_list, int region, long int size)
{
   VERTEX vertex/*, aux*/;  

   vertex = (VERTEX) malloc(sizeof(TPVERTEX));        
   vertex->region_num  = region;
   vertex->size_vertex = size;
   vertex->next_vertex = NULL; 

    if (*vertex_list == NULL) { *vertex_list = vertex; }
    else {
	vertex->next_vertex = *vertex_list;	
	*vertex_list = vertex;    
     }
}

/*----------------------------------------------------------------------*/
/*    tidy_insert_vertex_list							*/
/*----------------------------------------------------------------------*/
/*
  Insert a new node in a vertex_list, keeping the list sorted by growing
  region size.

  This must be the worst way to implement this, since the execution
  time is O(n^2).  We had much better append nodes in whatever order suit us
  (O(n)), then perform a fusion sort (O(n*log(n))).
*/
void tidy_insert_vertex_list(VERTEX *vertex_list, int region, long int size)
{
   int /*i, j, k,*/ end;
   VERTEX vertex, aux;  

   vertex = (VERTEX) malloc(sizeof(TPVERTEX));        
   vertex->region_num  = region;
   vertex->size_vertex = size;
   vertex->next_vertex = NULL; 

    if (*vertex_list == NULL) { 
    	*vertex_list = vertex; }
     else{
      aux= *vertex_list;		
      if (aux->size_vertex > vertex->size_vertex){
 	  vertex->next_vertex = aux;
	  *vertex_list = vertex;  
	  end= TRUE; 
      	}
       else{
         end= FALSE;
	 while (end != TRUE){
	   if (aux->next_vertex == NULL) 
	      end= TRUE;
	     else 
              if (aux->next_vertex->size_vertex > vertex->size_vertex)
	          end= TRUE;
		else
	  	  aux= aux->next_vertex;
	 }
	 vertex->next_vertex= aux->next_vertex;
	 aux->next_vertex= vertex;
       }
     }      
}


/*************   queues & graphs   ************************************/
 
/*--------------------------------------------------------------------*/
/*   initialize_vertex_list_graph				      */
/* build a vertex_list from the graph				      */
/*--------------------------------------------------------------------*/
/*
  build a list of every vertice, sorted by growing size.
*/
void initialize_vertex_list_graph(GRAPH *graph, VERTEX *vertex_list, grphic3d *o, int *label, grphic3d *im1resc)
{
    int i, region_vertex, /*merger,*/ current_label, number_mergings;
    /*VERTEX aux;*/ 
    long size;  

    number_mergings= 0;
    current_label= *label;
    
       for (i= (INIT+1); i<= current_label; i++){
         /* if ((*graph)->vertex_field[i].num_pixels < MINIMUN_SIZE){
	      merger= look_for_merger((*graph), i);
	      merge_graph_3d(&(*graph), i, merger, o, im1resc);
              printf("initialize_vertex_list_graph:  merger= %d,  merged= %d \n", merger, i);
	      /reset2_graph_3d(&(*graph));	/
	      (*graph)->vertex_field[i].grey_region = (*graph)->vertex_field[merger].grey_region;
      	      (*graph)->vertex_field[i].vertex_merged = merger;
    	      (*graph)->vertex_field[i].num_pixels = 0;

	      number_mergings++;
	    }
	   else{*/
	    size= (*graph)->vertex_field[i].num_pixels;
	    region_vertex= i;
	    tidy_insert_vertex_list(&(*vertex_list), region_vertex, size);
	  /* }        */
	}
  
  *label=current_label - number_mergings;	
}

 
/*----------------------------------------------------------------------*/
/*   initialize_array_vertex_graph						*/
/*----------------------------------------------------------------------*/
void initialize_array_vertex_graph(GRAPH *graph, grphic3d *o, int *label, grphic3d *im1resc, ARRAY_VERTEX *histogram_regions)
{
    int i, region_vertex, /*merger,*/ current_label, number_mergings;
    /*VERTEX aux;*/ 
    long size;  

    number_mergings= 0;
    current_label= *label;

    
       for (i= (INIT+2); i<= current_label; i++){
    /*      if ((*graph)->vertex_field[i].num_pixels < MINIMUN_SIZE){	  
	      merger= look_for_merger((*graph), i);
	      merge_graph_3d(&(*graph), i, merger, o, im1resc);
       printf("initialize: merger= %d, merged= %d de %d pixels \n", merger, i, (*graph)->vertex_field[i].num_pixels);
       	     histogram_regions[(*graph)->vertex_field[merger].grey_region].frequency+= (*graph)->vertex_field[i].num_pixels;
	      (*graph)->vertex_field[i].grey_region = (*graph)->vertex_field[merger].grey_region;
      	      (*graph)->vertex_field[i].vertex_merged = merger;
    	      (*graph)->vertex_field[i].num_pixels = 0;

	      number_mergings++;
	    }
	   else{*/
	    size= (*graph)->vertex_field[i].num_pixels;
	    region_vertex= i;
/*	    insert_vertex_list(&(histogram_regions[(*graph)->vertex_field[i].grey_region].regions_field), region_vertex,
size);*/
	    tidy_insert_vertex_list(&(histogram_regions[(*graph)->vertex_field[i].grey_region].regions_field), region_vertex, size);	    
	    histogram_regions[(*graph)->vertex_field[i].grey_region].frequency+= (*graph)->vertex_field[i].num_pixels;
	   /*}        */
	}
  
  *label=current_label - number_mergings;	
}


 
/*----------------------------------------------------------------------*/
/*   sort_ array3d							*/
/*! \brief  Sort an array in ascending order					*/
/*! \param   A : the array  (E/S)*/
/*----------------------------------------------------------------------*/

void sort_array3d(HISTO_ARRAY *A)
{
    int i, j;
    HISTO_ARRAY aux;   
  
    for (i=1; i<=255; i++)
	for (j=255; j>=i; j--)
	    if (A[j].frequency < A[j-1].frequency)
		{
		    aux.frequency = A[j-1].frequency;
		    aux.grey_level = A[j-1].grey_level;
		    A[j-1].frequency=A[j].frequency;
		    A[j-1].grey_level=A[j].grey_level;
		    A[j].frequency = aux.frequency;
		    A[j].grey_level = aux.grey_level;
		}
}

  
 
/*----------------------------------------------------------------------*/
/*   sort_array_vertex							*/
/*   Sort an array in ascending order					*/
/*----------------------------------------------------------------------*/

void sort_array_vertex(ARRAY_VERTEX *A)
{
    int i, j;
    ARRAY_VERTEX  aux;   
  
    for (i=1; i< 256; i++)
	for (j= (256-1); j>=i; j--)
	    if (A[j].frequency < A[j-1].frequency)
		{
		    aux.frequency = A[j-1].frequency;
		    aux.indexation = A[j-1].indexation;
		    aux.regions_field= A[j-1].regions_field;
		    A[j-1].frequency=A[j].frequency;
		    A[j-1].indexation=A[j].indexation;
		    A[j-1].regions_field= A[j].regions_field;
		    A[j].frequency = aux.frequency;
		    A[j].indexation= aux.indexation;
		    A[j].regions_field= aux.regions_field;		    
		}
}

/*----------------------------------------------------------------------*/
/*  look_for_merger							*/
/*----------------------------------------------------------------------*/
int look_for_merger(GRAPH graph, int region)
{
   int i, region_merger, found;
  
   region_merger= 0;
   found= FALSE; 
   i= INIT+2;
    
   /* merger= the first neighbourg(region) found with size >= MINIMUN_SIZE*/
   while ((i<= MAX_VERTEX) &&(found == FALSE)){
    if (are_vertex_connected(graph, region, i) || are_vertex_connected(graph, i, region) )
     if (graph->vertex_field[i].num_pixels >= MINIMUN_SIZE){
       found = TRUE;
       region_merger=i;
      }
    i++;  
    }
        
   return (region_merger);      
}



/*----------------------------------------------------------------------*/
/*  delete_vertex_histogram							*/
/*----------------------------------------------------------------------*/
    
void delete_vertex_histogram(GRAPH graph, ARRAY_VERTEX *histogram_regions, int index, int region)
{
    int found;
    VERTEX vertex, deleted; 

   found= FALSE;    
   
   if (histogram_regions[index].regions_field != NULL){
    vertex= histogram_regions[index].regions_field;
    if (vertex->region_num == region){
      found= TRUE;
      histogram_regions[index].regions_field= vertex->next_vertex;
      free(vertex);
     }
    else{
       while (found != TRUE){
	   if (vertex->next_vertex == NULL) 
	      found= TRUE; /* error the vertex does not exit, it is assumed this case never happen*/
	     else 
              if (vertex->next_vertex->region_num == region)
	          found= TRUE;
		else
	  	  vertex= vertex->next_vertex;
	 }
	 if (vertex->next_vertex != NULL){ /* else error*/	    	    
	    deleted = vertex->next_vertex;
	    vertex->next_vertex= deleted->next_vertex;
	    deleted->next_vertex=NULL;
/*	    if (index == 125)
	      printf("es \n"); */
	    histogram_regions[index].frequency -= deleted->size_vertex; 
	    free(deleted);	    
	   }
       }
     }   
}

/*----------------------------------------------------------------------*/
/*   sort_ array_graph							*/
/*   Sort an array in ascending order					*/
/*----------------------------------------------------------------------*/

void sort_array_graph(ARRAY_GRAPH *A)
{
    int i, j;
    ARRAY_GRAPH aux;   
  
    for (i=1; i< MAX_VERTEX; i++)
	for (j= (MAX_VERTEX-1); j>=i; j--)
	    if (A[j].contrast < A[j-1].contrast)
		{
		    aux.contrast = A[j-1].contrast;
		    aux.index = A[j-1].index;
		    A[j-1].contrast=A[j].contrast;
		    A[j-1].index=A[j].index;
		    A[j].contrast = aux.contrast;
		    A[j].index= aux.index;
		}
}

 
/*----------------------------------------------------------------------*/
/*   copy_sort_array							*/
/*----------------------------------------------------------------------*/

void copy_array_vertex(ARRAY_VERTEX *A, ARRAY_VERTEX *A_indexed)
{
    int i;
  
     for (i= 0; i< 256; i++){	
	A_indexed[i].frequency=A[i].frequency;
	A_indexed[i].indexation= A[i].indexation;
	A_indexed[i].regions_field= NULL /* A[i].regions_field*/;		    
      }

}



