/* Minimal MPI placement probe: one JSON line per rank (rank, size, affinity, nice,
   CUDA_VISIBLE_DEVICES, thread variables) plus an Allreduce to show the ranks communicate. */
#define _GNU_SOURCE
#include <errno.h>
#include <mpi.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>

static void json_env(const char *key, int last)
{
    const char *value = getenv(key);
    if (value) printf("\"%s\":\"%s\"%s", key, value, last ? "" : ",");
    else printf("\"%s\":null%s", key, last ? "" : ",");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, size = 0, sum = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Allreduce(&rank, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    cpu_set_t set;
    CPU_ZERO(&set);
    sched_getaffinity(0, sizeof set, &set);
    char affinity[4096] = "";
    size_t used = 0;
    for (int cpu = 0; cpu < CPU_SETSIZE && used + 16 < sizeof affinity; ++cpu)
        if (CPU_ISSET(cpu, &set)) used += (size_t)snprintf(affinity + used, sizeof affinity - used, used ? ",%d" : "%d", cpu);
    errno = 0;
    int nice_value = getpriority(PRIO_PROCESS, 0);
    char host[256] = "";
    gethostname(host, sizeof host - 1);
    for (int turn = 0; turn < size; ++turn) {
        if (turn == rank) {
            printf("{\"probe\":\"mpi_probe\",\"host\":\"%s\",\"rank\":%d,\"size\":%d,\"allreduce_rank_sum\":%d,"
                   "\"affinity\":[%s],\"running_cpu\":%d,\"nice\":%d,\"env\":{",
                   host, rank, size, sum, affinity, sched_getcpu(), nice_value);
            json_env("CUDA_VISIBLE_DEVICES", 0);
            json_env("PF_LAUNCH_RANK", 0);
            json_env("PF_CORE", 0);
            json_env("OMP_NUM_THREADS", 0);
            json_env("OPENBLAS_NUM_THREADS", 0);
            json_env("MKL_NUM_THREADS", 1);
            printf("}}\n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
