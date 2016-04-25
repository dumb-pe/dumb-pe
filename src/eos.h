#ifndef __EOS_CUBIC_H__
#define __EOS_CUBIC_H__

#include "libguile.h"

struct CubicEosModel;
typedef struct CubicEosModel CubicEosModel;

static SCM eos_cubic_model;

void eos_cubic_init (void);

SCM eos_cubic_make (SCM type, SCM nc, SCM rest);

SCM eos_cubic_print (SCM port, SCM model);

#endif  /* __EOS_CUBIC_H__ */
