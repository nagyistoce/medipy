#ifndef _CONFIG_H_
#define _CONFIG_H_ 

#define DMALLOC_DISABLE 
#define NO_GUI 
#define NO_VERIF 
#define COMPILE_FOR_MEDIPY

#ifdef WIN32
#include <float.h>
#define finite(x) _finite(x)
#define isinf(x) (!_isnan(x) && !_finite(x))
#define isnan(x) _isnan(x)
#define bzero(p, size) (void)memset((p), 0, (size))
#include "windowsOnly/misc.h"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#endif

#endif 
