#ifndef MBO_SYS_H
#define MBO_SYS_H

#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define MBO_API __attribute__ ((dllexport))
    #else
      #define MBO_API __declspec(dllexport)
    #endif
  #else
    #ifdef __GNUC__
      #define MBO_API __attribute__ ((dllimport))
    #else
      #define MBO_API __declspec(dllimport)
    #endif
  #endif
  #define MBO_Internal
#else
  #if __GNUC__ >= 4
    #define MBO_API __attribute__ ((visibility ("default")))
    #define MBO_Internal  __attribute__ ((visibility ("hidden")))
  #else
    #define MBO_API
    #define MBO_Internal
  #endif
#endif

#endif

