// Time-stamp: <2009-01-17 15:01:29 hamada>

/*
 * Copyright (C) 2007 
 *      Tsuyoshi Hamada <hamada@progrape.jp>
 *      All rights reserved.
 * This code is released under version 2 of the GNU GPL.
 */

#define IDIM  (4)
#define JDIM  (4)
#define FDIM  (3)
//#undef SM_MAX_BYTE
//#define SM_MAX_BYTE (16384-32)
//#define NJ_SHMEM 256
#define NJ_SHMEM 128  // ** Technic ** 
#define NSP 8
#define NVSP 16
//#define NVSP 24 // GTX260
#define NPIPE (NSP*NVSP)
#define   KIRIAGE(x,y)     (((x) % (y)) ?  ((x/y)+1) : (x/y))
#define   MAX(x,y)     (((x) > (y)) ?  (x) : (y))
#define   MIN(x,y)     (((x) < (y)) ?  (x) : (y))
