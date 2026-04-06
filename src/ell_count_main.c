#include "common.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#define MAX_RANK 20
#define MAX_EXP 20
#define COUNTS_N (12 + 6 * (MAX_RANK + 1) + 4 * (MAX_EXP + 1))

#define C_TOTAL 0
#define C_INERT 1
#define C_RAMIFIED 2
#define C_SPLIT 3
#define C_INERT_TRIVIAL 4
#define C_INERT_NONTRIVIAL 5
#define C_RAMIFIED_TRIVIAL 6
#define C_RAMIFIED_NONTRIVIAL 7
#define C_INERT_RANK_INC 8
#define C_INERT_RANK_SAME 9
#define C_RAMIFIED_RANK_INC 10
#define C_RAMIFIED_RANK_SAME 11
#define C_INERT_BY_RANK 12
#define C_INERT_BY_RANK_INC (12 + MAX_RANK + 1)
#define C_INERT_BY_RANK_SAME (12 + 2 * (MAX_RANK + 1))
#define C_RAMIFIED_BY_RANK (12 + 3 * (MAX_RANK + 1))
#define C_RAMIFIED_BY_RANK_INC (12 + 4 * (MAX_RANK + 1))
#define C_RAMIFIED_BY_RANK_SAME (12 + 5 * (MAX_RANK + 1))
#define C_INERT_CYCLIC (12 + 6 * (MAX_RANK + 1))
#define C_INERT_CYCLIC_RANK_SAME (12 + 6 * (MAX_RANK + 1) + (MAX_EXP + 1))
#define C_RAMIFIED_CYCLIC (12 + 6 * (MAX_RANK + 1) + 2 * (MAX_EXP + 1))
#define C_RAMIFIED_CYCLIC_RANK_SAME (12 + 6 * (MAX_RANK + 1) + 3 * (MAX_EXP + 1))

/* Response buffer sent from worker to coordinator:
 *   resp[0]       = work_idx (cast to uint64_t)
 *   resp[1..N]    = counts[COUNTS_N]
 */
#define RESP_N (1 + COUNTS_N)

/* ell-adic valuation: largest v such that ell^v divides n. */
static uint32_t ell_val(uint64_t n, uint64_t ell) {
    if (n == 0) {
        return 0;
    }
    uint32_t v = 0;
    while (n % ell == 0) {
        n /= ell;
        v++;
    }
    return v;
}

/* Validation helpers */

#define CHECK(cond, ...)                                                                           \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            fprintf(stderr, "ERROR: " __VA_ARGS__);                                                \
            fprintf(stderr, " at boundary=%ld\n", boundary);                                       \
            MPI_Abort(MPI_COMM_WORLD, 1);                                                          \
        }                                                                                          \
    } while (0)

typedef struct {
    int total, trivial, nontrivial;
    int rank_inc, rank_same;
    int by_rank, by_rank_inc, by_rank_same;
    int cyclic, cyclic_rank_same;
    const char *label;
} split_indices_t;

static const split_indices_t INERT = {
    C_INERT,
    C_INERT_TRIVIAL,
    C_INERT_NONTRIVIAL,
    C_INERT_RANK_INC,
    C_INERT_RANK_SAME,
    C_INERT_BY_RANK,
    C_INERT_BY_RANK_INC,
    C_INERT_BY_RANK_SAME,
    C_INERT_CYCLIC,
    C_INERT_CYCLIC_RANK_SAME,
    "inert",
};

static const split_indices_t RAMIFIED = {
    C_RAMIFIED,
    C_RAMIFIED_TRIVIAL,
    C_RAMIFIED_NONTRIVIAL,
    C_RAMIFIED_RANK_INC,
    C_RAMIFIED_RANK_SAME,
    C_RAMIFIED_BY_RANK,
    C_RAMIFIED_BY_RANK_INC,
    C_RAMIFIED_BY_RANK_SAME,
    C_RAMIFIED_CYCLIC,
    C_RAMIFIED_CYCLIC_RANK_SAME,
    "ramified",
};

static void validate_counts(const uint64_t *c, long boundary, const split_indices_t *s) {
    CHECK(c[s->trivial] + c[s->nontrivial] == c[s->total], "%s_trivial+%s_non_trivial != %s",
          s->label, s->label, s->label);
    CHECK(c[s->by_rank + 0] == c[s->trivial], "%s_rank_0 != %s_trivial", s->label, s->label);

    uint64_t rank_pos = 0;
    for (int r = 1; r <= MAX_RANK; r++) {
        rank_pos += c[s->by_rank + r];
    }
    CHECK(rank_pos == c[s->nontrivial], "sum(%s_rank_r, r>=1) != %s_non_trivial", s->label,
          s->label);

    CHECK(c[s->rank_inc] + c[s->rank_same] == c[s->total], "%s_rank_inc+%s_rank_same != %s",
          s->label, s->label, s->label);

    uint64_t sum_inc = 0, sum_same = 0;
    for (int r = 0; r <= MAX_RANK; r++) {
        CHECK(c[s->by_rank_inc + r] + c[s->by_rank_same + r] == c[s->by_rank + r],
              "%s_rank_%d_inc+%s_rank_%d_same != %s_rank_%d", s->label, r, s->label, r, s->label,
              r);
        sum_inc += c[s->by_rank_inc + r];
        sum_same += c[s->by_rank_same + r];
    }
    CHECK(sum_inc == c[s->rank_inc], "sum(%s_rank_r_inc) != %s_rank_inc", s->label, s->label);
    CHECK(sum_same == c[s->rank_same], "sum(%s_rank_r_same) != %s_rank_same", s->label, s->label);
    CHECK(c[s->by_rank_same + 0] == 0, "%s_rank_0_same != 0", s->label);

    /* Cyclic validation */
    CHECK(c[s->cyclic + 0] == c[s->trivial], "%s_cyclic_0 != %s_trivial", s->label, s->label);

    uint64_t sum_cyc = 0, sum_cyc_rs = 0;
    for (int e = 0; e <= MAX_EXP; e++) {
        sum_cyc += c[s->cyclic + e];
        sum_cyc_rs += c[s->cyclic_rank_same + e];
        CHECK(c[s->cyclic_rank_same + e] <= c[s->cyclic + e],
              "%s_cyclic_%d_rank_same > %s_cyclic_%d", s->label, e, s->label, e);
    }
    CHECK(sum_cyc == c[s->by_rank + 0] + c[s->by_rank + 1],
          "sum(%s_cyclic) != %s_rank_0 + %s_rank_1", s->label, s->label, s->label);
    CHECK(sum_cyc_rs == c[s->by_rank_same + 1], "sum(%s_cyclic_rank_same) != %s_rank_1_same",
          s->label, s->label);
    CHECK(c[s->cyclic_rank_same + 0] == 0, "%s_cyclic_0_rank_same != 0", s->label);
}

/*
 * Process one (congruence class, file index) pair and write aggregate
 * counts into `counts`.
 *
 * Reads two gzipped files in lockstep:
 *   fundamental: {folder}/cl{a}mod{m}/cl{a}mod{m}.{file_idx}.gz
 *                line format: dist h c1 c2 ... ct
 *   ell:         {folder}/cl{a}mod{m}l{ell}/cl{a}mod{m}l{ell}.{file_idx}.gz
 *                line format: dist kron c1 c2 ... ct
 */
static void process_file(const char *folder, int a, int m, long ell, long file_idx, long files,
                         long d_max, uint64_t counts[COUNTS_N]) {
    memset(counts, 0, COUNTS_N * sizeof(uint64_t));

    long d_total = d_max / (files * (long)m);
    long starting_d = file_idx * d_total * (long)m + (long)a;

    char fund_cmd[2048], ell_cmd[2048];
    snprintf(fund_cmd, sizeof(fund_cmd), "gunzip -c '%s/cl%dmod%d/cl%dmod%d.%ld.gz'", folder, a, m,
             a, m, file_idx);
    snprintf(ell_cmd, sizeof(ell_cmd), "gunzip -c '%s/cl%dmod%dl%ld/cl%dmod%dl%ld.%ld.gz'", folder,
             a, m, ell, a, m, ell, file_idx);

    FILE *fund_f = popen(fund_cmd, "r");
    if (fund_f == NULL) {
        fprintf(stderr, "Error opening fundamental file: a=%d m=%d idx=%ld\n", a, m, file_idx);
        return;
    }
    FILE *ell_f = popen(ell_cmd, "r");
    if (ell_f == NULL) {
        fprintf(stderr, "Error opening ell file: a=%d m=%d l=%ld idx=%ld\n", a, m, ell, file_idx);
        pclose(fund_f);
        return;
    }

    char fund_line[LINE_BUFFER_SIZE], ell_line[LINE_BUFFER_SIZE];
    long fund_d = starting_d;
    long ell_d = starting_d;

    while (fgets(fund_line, LINE_BUFFER_SIZE, fund_f) != NULL &&
           fgets(ell_line, LINE_BUFFER_SIZE, ell_f) != NULL) {
        /* Parse fundamental entry: dist h c1 c2 ... ct */
        char *tok = strtok(fund_line, " \t\n");
        if (tok == NULL) {
            continue;
        }
        long fund_dist = atol(tok);

        tok = strtok(NULL, " \t\n"); /* skip class number h */
        if (tok == NULL) {
            continue;
        }

        int fund_ell_rank = 0;
        uint32_t fund_ell_exp = 0;
        while ((tok = strtok(NULL, " \t\n")) != NULL) {
            uint64_t c = strtoull(tok, NULL, 10);
            if (c > 0) {
                uint32_t v = ell_val(c, (uint64_t)ell);
                if (v > 0) {
                    fund_ell_rank++;
                    if (v > fund_ell_exp) {
                        fund_ell_exp = v;
                    }
                }
            }
        }
        fund_d += fund_dist * (long)m;

        /* Parse ell entry: dist kron c1 c2 ... ct */
        tok = strtok(ell_line, " \t\n");
        if (tok == NULL) {
            continue;
        }
        long ell_dist = atol(tok);

        tok = strtok(NULL, " \t\n");
        if (tok == NULL) {
            continue;
        }
        int kron = atoi(tok);

        int ell_ell_rank = 0;
        while ((tok = strtok(NULL, " \t\n")) != NULL) {
            uint64_t c = strtoull(tok, NULL, 10);
            if (c > 0 && c % (uint64_t)ell == 0) {
                ell_ell_rank++;
            }
        }
        ell_d += ell_dist * (long)m;

        if (fund_d != ell_d) {
            fprintf(stderr,
                    "Discriminant mismatch: fund=%ld ell=%ld"
                    " (a=%d m=%d idx=%ld)\n",
                    fund_d, ell_d, a, m, file_idx);
            break;
        }

        if (fund_ell_rank > MAX_RANK) {
            fprintf(stderr,
                    "WARNING: fund_ell_rank=%d > MAX_RANK=%d"
                    " at d=%ld (a=%d m=%d idx=%ld)\n",
                    fund_ell_rank, MAX_RANK, fund_d, a, m, file_idx);
        }
        if (fund_ell_exp > MAX_EXP) {
            fprintf(stderr,
                    "WARNING: fund_ell_exp=%u > MAX_EXP=%d"
                    " at d=%ld (a=%d m=%d idx=%ld)\n",
                    fund_ell_exp, MAX_EXP, fund_d, a, m, file_idx);
        }

        int trivial = (fund_ell_rank == 0);
        int r = (fund_ell_rank <= MAX_RANK) ? fund_ell_rank : MAX_RANK;
        int is_cyclic = (fund_ell_rank <= 1);
        int cyc_exp = (int)((fund_ell_exp <= MAX_EXP) ? fund_ell_exp : MAX_EXP);

        counts[C_TOTAL]++;

        if (kron == -1) {
            counts[C_INERT]++;
            if (trivial) {
                counts[C_INERT_TRIVIAL]++;
            } else {
                counts[C_INERT_NONTRIVIAL]++;
            }
            if (ell_ell_rank > fund_ell_rank) {
                counts[C_INERT_RANK_INC]++;
                counts[C_INERT_BY_RANK_INC + r]++;
            } else {
                counts[C_INERT_RANK_SAME]++;
                counts[C_INERT_BY_RANK_SAME + r]++;
                if (is_cyclic) {
                    counts[C_INERT_CYCLIC_RANK_SAME + cyc_exp]++;
                }
            }
            if (is_cyclic) {
                counts[C_INERT_CYCLIC + cyc_exp]++;
            }
            counts[C_INERT_BY_RANK + r]++;
        } else if (kron == 0) {
            counts[C_RAMIFIED]++;
            if (ell != 3) {
                if (trivial) {
                    counts[C_RAMIFIED_TRIVIAL]++;
                } else {
                    counts[C_RAMIFIED_NONTRIVIAL]++;
                }
                if (ell_ell_rank > fund_ell_rank) {
                    counts[C_RAMIFIED_RANK_INC]++;
                    counts[C_RAMIFIED_BY_RANK_INC + r]++;
                } else {
                    counts[C_RAMIFIED_RANK_SAME]++;
                    counts[C_RAMIFIED_BY_RANK_SAME + r]++;
                    if (is_cyclic) {
                        counts[C_RAMIFIED_CYCLIC_RANK_SAME + cyc_exp]++;
                    }
                }
                if (is_cyclic) {
                    counts[C_RAMIFIED_CYCLIC + cyc_exp]++;
                }
                counts[C_RAMIFIED_BY_RANK + r]++;
            }
        } else {
            /* kron == 1: split */
            counts[C_SPLIT]++;
        }
    }

    pclose(fund_f);
    pclose(ell_f);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    if (argc != 5) {
        fprintf(stderr, "Format: mpirun -np [#procs] ./ell_count"
                        " [d_max] [files] [ell] [folder]\n");
        MPI_Finalize();
        exit(1);
    }

    const long d_max = atol(argv[1]);
    const long files = atol(argv[2]);
    const long ell = atol(argv[3]);
    const char *folder = argv[4];
    const long step = d_max / files;
    const int N = NUM_CONGRUENCES * (int)files; /* total work items */

    int myrank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (myrank == 0) {
        /* Coordinator */
        int W = num_procs - 1;         /* number of workers */
        int initial = (W < N) ? W : N; /* initial batch size */

        uint64_t *by_index = (uint64_t *)calloc(files * COUNTS_N, sizeof(uint64_t));
        if (by_index == NULL) {
            fprintf(stderr, "Memory allocation failed.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fprintf(stderr, "ell_count: l=%ld, D_max=%ld, files=%ld, step=%ld\n", ell, d_max, files,
                step);
        fprintf(stderr, "Processing %d file pairs ...\n", N);
        fflush(stderr);

        /* Send initial batch */
        for (int i = 0; i < initial; i++) {
            MPI_Send(&i, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
        }

        uint64_t resp[RESP_N];
        MPI_Status status;
        int done = 0;

        /* Main distribution loop */
        for (int i = initial; i < N; i++) {
            MPI_Recv(resp, RESP_N, MPI_UINT64_T, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            int src = status.MPI_SOURCE;
            int file_idx = (int)(resp[0]) % (int)files;

            uint64_t *slot = &by_index[(long)file_idx * COUNTS_N];
            for (int k = 0; k < COUNTS_N; k++) {
                slot[k] += resp[1 + k];
            }

            done++;
            fprintf(stderr, "\r  %d/%d", done, N);
            fflush(stderr);

            MPI_Send(&i, 1, MPI_INT, src, 0, MPI_COMM_WORLD);
        }

        /* Drain remaining outstanding responses */
        for (int i = 0; i < initial; i++) {
            MPI_Recv(resp, RESP_N, MPI_UINT64_T, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            int file_idx = (int)(resp[0]) % (int)files;

            uint64_t *slot = &by_index[(long)file_idx * COUNTS_N];
            for (int k = 0; k < COUNTS_N; k++) {
                slot[k] += resp[1 + k];
            }

            done++;
            fprintf(stderr, "\r  %d/%d", done, N);
            fflush(stderr);
        }
        fprintf(stderr, "\n");
        fflush(stderr);

        /* Terminate all workers */
        int sentinel = -1;
        for (int i = 1; i <= W; i++) {
            MPI_Send(&sentinel, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        /* Compute cumulative sums across file indices */
        uint64_t *cumul = (uint64_t *)calloc(files * COUNTS_N, sizeof(uint64_t));
        if (cumul == NULL) {
            fprintf(stderr, "Memory allocation failed.\n");
            free(by_index);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        {
            uint64_t running[COUNTS_N];
            memset(running, 0, sizeof(running));
            for (long fi = 0; fi < files; fi++) {
                uint64_t *src_slot = &by_index[fi * COUNTS_N];
                for (int k = 0; k < COUNTS_N; k++) {
                    running[k] += src_slot[k];
                }
                memcpy(&cumul[fi * COUNTS_N], running, COUNTS_N * sizeof(uint64_t));
            }
        }

        /* Find max non-zero rank/exponent columns to display */
        int max_inert_rank = 0, max_ram_rank = 0;
        int max_inert_exp = 0, max_ram_exp = 0;
        for (long fi = 0; fi < files; fi++) {
            uint64_t *c = &cumul[fi * COUNTS_N];
            for (int r = MAX_RANK; r > max_inert_rank; r--) {
                if (c[C_INERT_BY_RANK + r] > 0) {
                    max_inert_rank = r;
                    break;
                }
            }
            for (int e = MAX_EXP; e > max_inert_exp; e--) {
                if (c[C_INERT_CYCLIC + e] > 0) {
                    max_inert_exp = e;
                    break;
                }
            }
            if (ell != 3) {
                for (int r = MAX_RANK; r > max_ram_rank; r--) {
                    if (c[C_RAMIFIED_BY_RANK + r] > 0) {
                        max_ram_rank = r;
                        break;
                    }
                }
                for (int e = MAX_EXP; e > max_ram_exp; e--) {
                    if (c[C_RAMIFIED_CYCLIC + e] > 0) {
                        max_ram_exp = e;
                        break;
                    }
                }
            }
        }

        /* CSV header */
        printf("boundary,total,inert,ramified,split"
               ",inert_trivial,inert_non_trivial"
               ",inert_rank_inc,inert_rank_same");
        for (int r = 0; r <= max_inert_rank; r++) {
            printf(",inert_rank_%d,inert_rank_%d_inc,inert_rank_%d_same", r, r, r);
        }
        for (int e = 0; e <= max_inert_exp; e++) {
            printf(",inert_cyclic_%d,inert_cyclic_%d_rank_same", e, e);
        }
        if (ell != 3) {
            printf(",ramified_trivial,ramified_non_trivial"
                   ",ramified_rank_inc,ramified_rank_same");
            for (int r = 0; r <= max_ram_rank; r++) {
                printf(",ramified_rank_%d,ramified_rank_%d_inc"
                       ",ramified_rank_%d_same",
                       r, r, r);
            }
            for (int e = 0; e <= max_ram_exp; e++) {
                printf(",ramified_cyclic_%d,ramified_cyclic_%d_rank_same", e, e);
            }
        }
        printf("\n");

        /* CSV data rows */
        for (long fi = 0; fi < files; fi++) {
            uint64_t *c = &cumul[fi * COUNTS_N];
            long boundary = (fi + 1) * step;

            CHECK(c[C_INERT] + c[C_RAMIFIED] + c[C_SPLIT] == c[C_TOTAL],
                  "inert+ramified+split != total");
            validate_counts(c, boundary, &INERT);
            if (ell != 3) {
                validate_counts(c, boundary, &RAMIFIED);
            }

            printf("%ld,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64
                   ",%" PRIu64 ",%" PRIu64,
                   boundary, c[C_TOTAL], c[C_INERT], c[C_RAMIFIED], c[C_SPLIT], c[C_INERT_TRIVIAL],
                   c[C_INERT_NONTRIVIAL], c[C_INERT_RANK_INC], c[C_INERT_RANK_SAME]);
            for (int r = 0; r <= max_inert_rank; r++) {
                printf(",%" PRIu64 ",%" PRIu64 ",%" PRIu64, c[C_INERT_BY_RANK + r],
                       c[C_INERT_BY_RANK_INC + r], c[C_INERT_BY_RANK_SAME + r]);
            }
            for (int e = 0; e <= max_inert_exp; e++) {
                printf(",%" PRIu64 ",%" PRIu64, c[C_INERT_CYCLIC + e],
                       c[C_INERT_CYCLIC_RANK_SAME + e]);
            }
            if (ell != 3) {
                printf(",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64, c[C_RAMIFIED_TRIVIAL],
                       c[C_RAMIFIED_NONTRIVIAL], c[C_RAMIFIED_RANK_INC], c[C_RAMIFIED_RANK_SAME]);
                for (int r = 0; r <= max_ram_rank; r++) {
                    printf(",%" PRIu64 ",%" PRIu64 ",%" PRIu64, c[C_RAMIFIED_BY_RANK + r],
                           c[C_RAMIFIED_BY_RANK_INC + r], c[C_RAMIFIED_BY_RANK_SAME + r]);
                }
                for (int e = 0; e <= max_ram_exp; e++) {
                    printf(",%" PRIu64 ",%" PRIu64, c[C_RAMIFIED_CYCLIC + e],
                           c[C_RAMIFIED_CYCLIC_RANK_SAME + e]);
                }
            }
            printf("\n");
        }
        fflush(stdout);

        free(by_index);
        free(cumul);
    } else {
        /* Worker */
        int work_idx;
        uint64_t resp[RESP_N];
        MPI_Status status;

        MPI_Recv(&work_idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        while (work_idx != -1) {
            int class_idx = work_idx / (int)files;
            int file_idx = work_idx % (int)files;
            int a = CONGRUENCES[class_idx][0];
            int m = CONGRUENCES[class_idx][1];

            uint64_t counts[COUNTS_N];
            process_file(folder, a, m, ell, (long)file_idx, files, d_max, counts);

            resp[0] = (uint64_t)work_idx;
            memcpy(&resp[1], counts, COUNTS_N * sizeof(uint64_t));

            MPI_Send(resp, RESP_N, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&work_idx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        }
    }

    MPI_Finalize();
    return 0;
}
