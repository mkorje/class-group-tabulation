#include "clgrp_ell.h"

#include "common.h"
#include <libclgrp/clgrp.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#ifdef WITH_PARI
#include <pari/pari.h>

/* Verify class group structure against PARI/GP. Aborts on mismatch. */
void pari_verify(int *result, const long D) {
    GEN g = stoi(D);
    GEN v = quadclassunit0(g, 0, 0, 0);
    GEN clgrp = compo(v, 2);

    for (int i = 1; i < lg(clgrp); i++) {
        if (result[i - 1] != gtolong(compo(clgrp, i))) {
            char err[100];
            sprintf(err, "D=%ld, Cl[%d] = %d, actual %ld\n", D, i - 1, result[i - 1],
                    gtolong(compo(clgrp, i)));
            fprintf(stderr, "%s", err);
            exit(1);
        }
    }

    cgiv(clgrp);
    cgiv(v);
    cgiv(g);
}
#endif

#define MAX_INVARIANTS 20

int verify_input_files_exist(const char *folder, int a, int m, long files) {
    char name[512];

    for (long i = 0; i < files; i++) {
        snprintf(name, sizeof(name), "%s/cl%dmod%d/cl%dmod%d.%ld.gz", folder, a, m, a, m, i);
        if (access(name, F_OK) == -1) {
            fprintf(stderr, "Missing input file: %s\n", name);
            return 0;
        }
    }

    return 1;
}

void process_clgrp_file(const int index, const long D_total, const char *folder, const int a,
                        const int m, const long ell, const int *spf) {
    char cmd[1024], output_name[512], output_dir[512];
    char line[LINE_BUFFER_SIZE];
    FILE *inf, *outf;

    struct timeval begin, end;

    /* Check if compressed output file already exists and is valid */
    snprintf(output_name, sizeof(output_name), "%s/cl%dmod%dl%ld/cl%dmod%dl%ld.%d.gz", folder, a, m,
             ell, a, m, ell, index);
    if (access(output_name, F_OK) != -1) {
        snprintf(cmd, sizeof(cmd), "gzip -t '%s'", output_name);
        if (system(cmd) == 0) {
            printf("Output file %s already exists, skipping.\n", output_name);
            return;
        }
        /* Corrupt gz from interrupted run; remove and reprocess */
        fprintf(stderr, "Removing corrupt output file %s\n", output_name);
        remove(output_name);
    }

    /* Remove partial output from interrupted run */
    snprintf(output_name, sizeof(output_name), "%s/cl%dmod%dl%ld/cl%dmod%dl%ld.%d", folder, a, m,
             ell, a, m, ell, index);
    if (access(output_name, F_OK) != -1) {
        remove(output_name);
    }

    /* Open input file via gunzip */
    snprintf(cmd, sizeof(cmd), "gunzip -c '%s/cl%dmod%d/cl%dmod%d.%d.gz'", folder, a, m, a, m,
             index);
    inf = popen(cmd, "r");
    if (inf == NULL) {
        fprintf(stderr, "Failed to open input file for index %d\n", index);
        return;
    }

    const long D_max = (index + 1) * D_total * m;
    long h_max = h_upper_bound(-D_max * ell * ell * ell * ell);
    int table_size = next_prime((((long)sqrt(h_max)) << 1) - 1);
    if (table_size == -1) {
        fprintf(stderr, "Not enough primes in liboptarith/primes.h\n");
        exit(1);
    }

    htab_t R, Q;
    htab_init(&R, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);
    htab_init(&Q, table_size, sizeof(form_t), &hash_form_t, &eq_form_t, &del_form_t);

    snprintf(output_dir, sizeof(output_dir), "%s/cl%dmod%dl%ld", folder, a, m, ell);
    mkdir(output_dir, 0744);
    snprintf(output_name, sizeof(output_name), "%s/cl%dmod%dl%ld/cl%dmod%dl%ld.%d", folder, a, m,
             ell, a, m, ell, index);
    outf = fopen(output_name, "w");
    if (outf == NULL) {
        fprintf(stderr, "Failed to create output file %s\n", output_name);
        pclose(inf);
        return;
    }

    long D = (long)index * D_total * m + a;
    int dist, h;
    int result[MAX_INVARIANTS];
    int output_rank;
    int kron;
    char output_line[LINE_BUFFER_SIZE];
    long D_sub;
    int init_pow, h_temp;

    gettimeofday(&begin, NULL);

    while (fgets(line, LINE_BUFFER_SIZE, inf) != NULL) {
        /* Parse input line: dist h c1 c2 ... ct */
        char *tok = strtok(line, " \t\n");
        if (tok == NULL) {
            continue;
        }

        dist = atoi(tok);

        tok = strtok(NULL, " \t\n");
        if (tok == NULL) {
            continue;
        }
        h = atoi(tok);

        D += (long)dist * m;

        kron = kronecker_symbol(-D, ell);

        D_sub = D;
        if (kron == 0) {
            h *= ell;
            D_sub *= ell * ell;
        } else if (kron == -1) {
            h *= (ell + 1) * ell;
            D_sub *= ell * ell * ell * ell;
        } else if (kron == 1) {
            h *= (ell - 1) * ell;
            D_sub *= ell * ell * ell * ell;
        }

        if (D == 4) {
            h /= 2;
        } else if (D == 3) {
            h /= 3;
        }

        if (kron == 1 || (kron == 0 && ell == 3)) {
            snprintf(output_line, sizeof(output_line), "%d\t%d\t0\n", dist, kron);
        } else {
            /* Compute class structure of order of index ell^2 */
            init_pow = 1;
            h_temp = h;
            {
                int h_rem = h;
                while (h_rem > 1) {
                    int p = spf[h_rem];
                    h_temp /= p;
                    do {
                        h_rem /= p;
                    } while (h_rem % p == 0);
                    if (h_temp % p != 0) {
                        init_pow *= p;
                    }
                }
            }

            h /= init_pow;
            output_rank = compute_group_bjt(result, -D_sub, init_pow, h, &R, &Q);

            h = result[0] * init_pow;
            result[1] *= init_pow;

            /* Format output line: dist kron c1 c2 ... ct */
            int off = snprintf(output_line, sizeof(output_line), "%d\t%d\t", dist, kron);
            for (int r = 1; r < output_rank; r++) {
                off += snprintf(output_line + off, sizeof(output_line) - off, "%d ", result[r]);
            }
            snprintf(output_line + off, sizeof(output_line) - off, "%d\n", result[output_rank]);

#ifdef WITH_PARI
            pari_verify(result + 1, -D_sub);
#endif
        }

        fputs(output_line, outf);
    }

    pclose(inf);
    fclose(outf);

    /* Compress output file */
    snprintf(cmd, sizeof(cmd), "gzip '%s'", output_name);
    if (system(cmd) != 0) {
        fprintf(stderr, "Warning: gzip failed for %s\n", output_name);
    }

    htab_clear(&R);
    htab_clear(&Q);

    gettimeofday(&end, NULL);
    long exec_us = (end.tv_sec - begin.tv_sec) * 1000000L + (end.tv_usec - begin.tv_usec);
    printf("index=%d, a=%d, m=%d, ell=%ld, took %.3f\n", index, a, m, ell, exec_us / 1e6);
    fflush(stdout);
}
