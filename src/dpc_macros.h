/* Does not work on Darwin
#ifdef DBL_DECIMAL_DIG
  #define OP_DBL_Digs (DBL_DECIMAL_DIG)
#else  
  #ifdef DECIMAL_DIG
    #define OP_DBL_Digs (DECIMAL_DIG)
  #else  
    #define OP_DBL_Digs (DBL_DIG + 3)
  #endif
#endif
*/
// Digits to print to file
#define MYE 16

#define ERR_MEM printf("ERROR: Not enough memory. File: %s (line %d)\n.",__FILE__,__LINE__); \
  *dpc_status = 2; \
  return;


#define ERR_MEM_LP printf("ERROR: Not enough memory. File: %s (line %d)\n.",__FILE__,__LINE__); \
  *dpc_status = 2; \
  return(-1.0);
