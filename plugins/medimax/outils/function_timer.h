/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*
*** timevar.h
*** Login : <wiest@localhost.localdomain>
*** Started on  Thu Mar 11 08:38:09 2004 Nicolas Wiest
*** $Id: function_timer.h,v 1.1 2004/03/11 08:29:30 wiest Exp $
*** 	
***	Copyright (c) 2004, ULP-IPB Strasbourg.
***	All rights are reserved.
***
*/

#ifndef TIMEVAR_H_
#define TIMEVAR_H_

# ifdef TIME_FUNCTIONS


/*
** Using Gettimeofday on linux allows us to have a sub millisecond ganularity
*/

#  include <sys/time.h>

#  define START_TIMER \
static struct timeval _tstart, _tend; \
static struct timezone tz; \
double t1, t2; \
  gettimeofday(&_tstart, &tz);

#  define END_TIMER \
  gettimeofday(&_tend,&tz); \
  t1 =  (double)_tstart.tv_sec + (double)_tstart.tv_usec/(1000*1000); \
  t2 =  (double)_tend.tv_sec + (double)_tend.tv_usec/(1000*1000); \
  printf("Timed Function: %s run in: %f\n", __FUNCTION__, t2 - t1);
# else

#  define START_TIMER 
#  define END_TIMER 

#endif /* !END_TIMER */

#endif /* !TIMEVAR_H_ */
