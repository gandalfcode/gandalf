// ============================================================================
// Debug.h
// Contains macro definitions of use debug functions.
// If the appropriate debug level is set in the Makefile, then the 
// macros are defined to contain code that outputs messages to screen.
// If not, then functions are undefined and no additional code is included.
// ============================================================================

#if defined(DEBUG1)
#define debug1(x)   cout << x << endl
#else
#define debug1(x)
#endif

#if defined(DEBUG2)
#define debug2(x)   cout << x << endl
#else
#define debug2(x)
#endif
