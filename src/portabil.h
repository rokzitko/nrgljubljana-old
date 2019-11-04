// portabil.h - Portability code. Low-level stuff belongs here.
// Copyright (C) 2005-2019 Rok Zitko

#ifndef _portabil_h_
#define _portabil_h_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <complex>
#include <limits>

using namespace std;

// Autotools definitions, generated by "configure" script
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Square raises x to the power of two
#ifndef HAVE_SQR
#define sqr(x) ((x)*(x))
#endif

// Backtrace. The code needs to be linked using the -rdynamic flag on
// systems using GNU ld to make the function names available to the
// program.
#ifdef HAVE_EXECINFO_H
 #include <execinfo.h>

 #include <cstdio>
 #include <cstdlib>
#endif
void print_trace(void)
{
#ifdef HAVE_EXECINFO_H
   const unsigned int max_trace_len = 100;

   void *array[max_trace_len];
   size_t size;
   char **strings;
   size_t i;
        
   size = backtrace (array, max_trace_len);
   strings = backtrace_symbols (array, size);
        
   printf ("Obtained %zd stack frames.\n", size);
   
   for (i = 0; i < size; i++)
     printf ("%s\n", strings[i]);
   
   free (strings);
#endif
}

#define assert_isfinite(x) finite_test_fnc(x, __FILE__, __LINE__)

bool my_isfinite(double x)
{
   return isfinite(x);
}

typedef complex<double> cmpl;

bool my_isfinite(cmpl z)
{
   return isfinite(z.real()) && isfinite(z.imag());
}

inline int isfinite(cmpl z)
{
   return (isfinite(z.real()) && isfinite(z.imag()) ? 1 : 0);
}

template <typename T>
T finite_test_fnc(T x, const char *file, int line)
{
   if (!my_isfinite(x)) {
      cout << "Finiteness assertion failed." << endl;
      cout << "File " << file << ", line " << line << "." << endl;
      cout << "Exiting." << endl; 
      cout << "#### EXITING DUE TO FAILED ASSERTION." << endl;
      print_trace();
      exit(1);
   }
   return x;
}

#define my_assert(x) do { \
   if (!(x)) { \
      cout << "Assertion failed." << endl; \
      cout << "File " << __FILE__ << ", line " << __LINE__ << "." << endl; \
      const char sAssert[] = #x; \
      cout << sAssert << endl; \
      cout << "Exiting." << endl; \
      cout << "#### EXITING DUE TO FAILED ASSERTION." << endl; \
      print_trace(); \
      exit(1); \
   } \
} while(0)

// Assert a==b
#define my_assert_equal(a, b) do { \
  if (!(a==b)) { \
      cout << "Assertion of equality failed." << endl; \
      cout << "File " << __FILE__ << ", line " << __LINE__ << "." << endl; \
      const char sAssert_a[] = #a; \
      const char sAssert_b[] = #b; \
      cout << sAssert_a << "=" << a << endl; \
      cout << sAssert_b << "=" << b << endl; \
      cout << "Exiting." << endl; \
      cout << "#### EXITING DUE TO FAILED ASSERTION." << endl; \
      print_trace(); \
      exit(1); \
   } \
} while(0)

#define my_warning(...) do { \
   cout << "Warning." << endl; \
   cout << "File " << __FILE__ << ", line " << __LINE__ << "." << endl; \
   char buf[256]; \
   sprintf(buf, __VA_ARGS__); \
   cout << buf << endl; \
   cout << "Continuing." << endl; \
} while(0)

#define my_error(...)  do { \
   cout << "Error." << endl; \
   cout << "File " << __FILE__ << ", line " << __LINE__ << "." << endl; \
   char buf[256]; \
   sprintf(buf, __VA_ARGS__); \
   cout << buf << endl; \
   cout << "Exiting" << endl; \
   cout << "#### EXITING DUE TO ERROR." << endl; \
   print_trace(); \
   exit(1); \
} while(0)

#define my_assert_not_reached() do { \
   cout << "Should never be reached." << endl; \
   cout << "File " << __FILE__ << ", line " << __LINE__ << "." << endl; \
   cout << "Exiting." << endl; \
   cout << "#### EXITING DUE TO FAILED ASSERTION." << endl; \
   print_trace(); \
   exit(1); \
} while(0)

// Workaround
#ifndef __TIMESTAMP__
#define __TIMESTAMP__ ""
#endif

#ifndef ARRAYLENGTH
# define ARRAYLENGTH(x) (sizeof(x) / sizeof(*(x)))
#endif

// TO DO: std::clamp in C++17
#ifndef CLIP
# define CLIP(x, xmin, xmax) min(max(x, xmin), xmax)
#endif

// *** Timing

#include <time.h>                                                                  
#include <sys/time.h>

// Return time in seconds                                                          
double gettime (void)                                                              
{                                                                                  
#ifdef HAVE_LIBRT                                                                  
// Low-level timing code pinched from glib-2.18.2                                  
// Other choices: CLOCK_MONOTONIC, CLOCK_PROCESS_CPUTIME_ID                        
#define USE_CLOCK CLOCK_REALTIME                                                   
#define NSEC_PER_SEC 1000000000                                                    
                                                                                   
   struct timespec tv;                                                             
                                                                                   
   clock_gettime(USE_CLOCK, &tv);                                                  
                                                                                   
   return tv.tv_sec + (double)tv.tv_nsec / NSEC_PER_SEC;                           
#else                                                                              
   struct timeval tv;
   
   gettimeofday(&tv, NULL);
   const double MILLION = 1.0e6;
   return (double)(tv.tv_sec+tv.tv_usec/MILLION);
#endif                                                                             
}               

// *** Memory usage

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach_init.h>
#include <mach/task.h>

// Source: http://blog.kuriositaet.de/?p=257
int getmem (unsigned int *rss, unsigned int *vs)
{
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS != task_info(mach_task_self(),
       TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
    {
        return -1;
    }
    *rss = t_info.resident_size;
    *vs  = t_info.virtual_size;
    return 0;
}

int memoryused() {
  unsigned int rss;
  unsigned int vs;
  
  getmem(&rss, &vs);
  
  return rss/1024; // result in kB block of memory!
}
#define HAS_MEMORY_USAGE
#endif

#if defined(__linux)
#include <string.h>
#include <stdio.h>
// Returns the number of kB blocks of memory used by this process (self)
int memoryused() {                                                         
   int used = 0;
   
   // See /usr/src/linux/Documentation/filesystems/proc.txt
   if (FILE *info = fopen("/proc/self/status", "r")) {                            
     char buf[512];                                                               
     int val;                                                                     
                                                                                   
     while (fgets(buf, sizeof(buf), info)) {                                      
       if (!strncmp(buf, "VmSize:", 7) && sscanf(buf+7, "%d", &val))              
         used = val;                                                              
     }                                                                            
     fclose(info);                                                                
   }
   
   return used;
}
#define HAS_MEMORY_USAGE                                                           
#endif                                                                             
                                                                                   
// Fall-back
#if !defined(HAS_MEMORY_USAGE)
int memoryused() { return 0; }
#endif

// More considerate way of opening output files: 
// - check if we really want to open the file (non-empty filename)
// - check for open errors
ofstream safeopen(string filename, bool binary = false)
{
   my_assert(filename != "");
   ios::openmode flag = ios::out;
   if (binary)
      flag |= ios::binary;
   ofstream F(filename.c_str(), flag);
   if (!F)
      my_error("safeopen: Can't open %s", filename.c_str());
   return F;
}

#endif