/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------
***	
***	file:		imx_stack.h
***
***	project:	Imagix 1.01 
***			
***
***	description:	Header dor imx_stack
***			/ fct about stack
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***---------------------------------------------------------------------*/


#ifndef	__IMX_STACK_INCL
#define	__IMX_STACK_INCL


struct stack {
  unsigned 	long nr;
  struct stack *next;
  struct stack *prev;
  int 		ptr_size;
  void	       *ptr;
} ;

extern	unsigned long	depth		(struct stack *stack);
extern	struct stack   *init_stack	(void);
extern	struct stack   *push		(struct stack *stack, void *elt, int elt_size);
extern	struct stack   *pull		(struct stack *stack);
extern	struct stack   *top		(struct stack *stack);
extern  struct stack   *top2    (struct stack *stack);

/* listes */

#endif

 
/* AUTOMATICALLY EXTRACTED PROTOTYPES (MB 6/12/2000) */
#ifdef __cplusplus 
extern "C" 
{
#endif /*__cplusplus*/ 
extern struct stack   *init_stack(void) ;
extern struct stack   *push(struct stack *stack, void *elt, int elt_size) ;
extern struct stack   *top(struct stack *stack) ;
extern struct stack   *pull(struct stack *stack) ;
#ifndef __cplusplus 
extern struct stack   *delete(struct stack *stack) ;
#endif
#ifdef __cplusplus 
}
#endif /*__cplusplus*/ 
