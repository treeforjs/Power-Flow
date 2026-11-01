// Stand-in driver library: forwards all entry points, times kernel launches (see README.txt).
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "stubs.h"

typedef int CUresult;
typedef void* CUfunction;
typedef void* CUstream;
typedef void* CUevent;

#define X(n) __attribute__((visibility("hidden"))) void* kt_real_##n;
KT_FORWARDED(X)
#undef X

static void* real_lib;
static CUresult (*r_launch)(CUfunction, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned, CUstream, void**, void**);
static CUresult (*r_stream_sync)(CUstream);
static CUresult (*r_ctx_sync)(void);
static CUresult (*r_event_create)(CUevent*, unsigned);
static CUresult (*r_event_record)(CUevent, CUstream);
static CUresult (*r_event_elapsed)(float*, CUevent, CUevent);
static CUresult (*r_event_destroy)(CUevent);
static CUresult (*r_func_name)(const char**, CUfunction);

typedef struct { CUevent a, b; int k; } Pending;
typedef struct { char name[96]; double ms; long launches; double threads; } Row;
static Pending* pend; static size_t npend, cpend;
static Row rows[512]; static int nrows;
static int depth;

__attribute__((constructor)) static void kt_init(void) {
    real_lib = dlopen("/usr/lib/libcuda.so.615.78.08", RTLD_NOW | RTLD_LOCAL);
    if (!real_lib) { fprintf(stderr, "kernel timer: %s\n", dlerror()); abort(); }
#define X(n) kt_real_##n = dlsym(real_lib, #n);
    KT_FORWARDED(X)
#undef X
    r_launch = dlsym(real_lib, "cuLaunchKernel");
    r_stream_sync = dlsym(real_lib, "cuStreamSynchronize");
    r_ctx_sync = dlsym(real_lib, "cuCtxSynchronize");
    r_event_create = dlsym(real_lib, "cuEventCreate");
    r_event_record = dlsym(real_lib, "cuEventRecord");
    r_event_elapsed = dlsym(real_lib, "cuEventElapsedTime");
    r_event_destroy = dlsym(real_lib, "cuEventDestroy_v2");
    r_func_name = dlsym(real_lib, "cuFuncGetName");
}

static int row_of(CUfunction f) {
    const char* nm = "?";
    if (r_func_name) r_func_name(&nm, f);
    for (int i = 0; i < nrows; ++i) if (!strncmp(rows[i].name, nm, sizeof rows[i].name - 1)) return i;
    if (nrows == 512) return 511;
    strncpy(rows[nrows].name, nm, sizeof rows[nrows].name - 1);
    return nrows++;
}

static void drain(void) {
    for (size_t i = 0; i < npend; ++i) {
        float ms = 0;
        if (r_event_elapsed(&ms, pend[i].a, pend[i].b) == 0) rows[pend[i].k].ms += ms;
        r_event_destroy(pend[i].a); r_event_destroy(pend[i].b);
    }
    npend = 0;
}

CUresult cuLaunchKernel(CUfunction f, unsigned gx, unsigned gy, unsigned gz, unsigned bx, unsigned by, unsigned bz,
                        unsigned shmem, CUstream s, void** params, void** extra) {
    CUevent a = 0, b = 0;
    int timed = depth == 0 && r_event_create(&a, 0) == 0 && r_event_create(&b, 0) == 0;
    if (timed) r_event_record(a, s);
    CUresult r = r_launch(f, gx, gy, gz, bx, by, bz, shmem, s, params, extra);
    if (timed) {
        r_event_record(b, s);
        int k = row_of(f);
        rows[k].launches += 1; rows[k].threads += (double)gx * gy * gz * bx * by * bz;
        if (npend == cpend) { cpend = cpend ? 2 * cpend : 4096; pend = realloc(pend, cpend * sizeof *pend); }
        pend[npend++] = (Pending){a, b, k};
    }
    return r;
}

CUresult cuStreamSynchronize(CUstream s) { CUresult r = r_stream_sync(s); if (r == 0) drain(); return r; }
CUresult cuCtxSynchronize(void) { CUresult r = r_ctx_sync(); if (r == 0) drain(); return r; }

__attribute__((destructor)) static void kt_fini(void) {
    const char* out = getenv("HALL_KTIMER_OUT");
    FILE* f = out ? fopen(out, "w") : NULL;
    if (!f) return;
    fprintf(f, "{\"kernels\": [\n");
    for (int i = 0; i < nrows; ++i)
        fprintf(f, " {\"name\": \"%s\", \"launches\": %ld, \"device_ms\": %.6f, \"threads\": %.0f}%s\n", rows[i].name,
                rows[i].launches, rows[i].ms, rows[i].threads, i + 1 < nrows ? "," : "");
    fprintf(f, "], \"undrained_launches\": %zu}\n", npend);
    fclose(f);
}
