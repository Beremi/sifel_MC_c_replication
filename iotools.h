#ifndef SIFEL_COMPAT_IOTOOLS_H
#define SIFEL_COMPAT_IOTOOLS_H

#include <stdarg.h>
#include <stdio.h>

// Minimal standalone replacement for the SIFEL XFILE input wrapper used by
// Tomas Koudelka's mohrc_ugn files.
struct XFILE
{
  FILE *file;
};

inline int xfscanf(XFILE *in, const char *fmt, ...)
{
  if ((in == NULL) || (in->file == NULL))
    return EOF;

  va_list args;
  va_start(args, fmt);
  const int ret = vfscanf(in->file, fmt, args);
  va_end(args);
  return ret;
}

// SIFEL stack-allocation macros. In this standalone repo the existing vector
// and matrix constructors allocate normally, so the macros collapse to sizes.
#define ASTCKVEC(n) (n)
#define ASTCKMAT(m, n) (m), (n)
#define RSTCKVEC(n, v) (n), (v)

#endif
