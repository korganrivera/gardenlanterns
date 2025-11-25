// lantern_ga.c
// Allele-based Genetic Algorithm for lantern placement (C)
// Compile: gcc -O3 -std=c99 lantern_ga.c -o lantern_ga -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/////// Parameters (edit and recompile) //////////////////////////////////////
#define POWER_LM          700.0      // lumens per lantern
#define LANTERN_HEIGHT_FT 6.0        // mount height in feet
#define TARGET_LUX        1.0        // target minimum lux
#define GRID_NX           100        // grid resolution X (coarse->fast, increase for accuracy)
#define GRID_NY           100        // grid resolution Y
#define POP_SIZE          48
#define GENERATIONS       200
#define MAX_SLOTS         10         // maximum loci (set to number of lanterns to test)
#define TOURNAMENT_K      3
#define CROSSOVER_RATE    0.90
#define MUTATION_RATE_DOM 0.12
#define MUTATION_RATE_REC 0.06
#define ELITISM           4
#define RANDOM_SEED       42
#define SAVE_CSV_HEATMAP  1          // 1 -> writes heatmap CSV "heatmap.csv"
///////////////////////////////////////////////////////////////////////////////

const double PI = 3.14159265358979323846;

typedef struct {
    double x, y;
    char type; // 'D' (dominant) or 'R' (recessive)
    int is_off; // 1 if OFF (value unused)
} Allele;

typedef struct {
    Allele a1, a2;
} Locus;

typedef struct {
    Locus loci[MAX_SLOTS]; // number of active loci determined by N_slots
} Chromosome;

// polygon: user-provided quadrilateral (feet)
double verts_ft[4][2] = {
    {56.509,   6.20},
    {156.641, 95.380},
    {125.857,170.050},
    {56.931, 162.723}
};

// conversions and global derived variables
const double FT_TO_M = 0.3048;
double verts_m[4][2];
double min_x, max_x, min_y, max_y;
double lantern_height_m;

// grid points inside polygon
double *grid_x = NULL, *grid_y = NULL;
int grid_count = 0;

// small RNG wrapper (LCG)
static unsigned long rng_state;
void rng_seed(unsigned long s) { rng_state = s ? s : 1; }
unsigned long rng_next() {
    rng_state = rng_state * 1664525UL + 1013904223UL;
    return rng_state;
}
double rand01() { return (double)(rng_next() & 0xFFFFFF) / (double)0x1000000; }

// point-in-polygon (ray casting)
int point_in_poly(double x, double y, double poly[][2], int n) {
    int inside = 0;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double xi = poly[i][0], yi = poly[i][1];
        double xj = poly[j][0], yj = poly[j][1];
        int intersect = ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi + 1e-18) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
}

// sample uniform random position inside polygon bounding box
void random_position_in_poly(double *rx, double *ry) {
    for (int attempt = 0; attempt < 200; ++attempt) {
        double x = min_x + rand01() * (max_x - min_x);
        double y = min_y + rand01() * (max_y - min_y);
        if (point_in_poly(x, y, verts_m, 4)) {
            *rx = x; *ry = y; return;
        }
    }
    // fallback: linear scan pick random interior grid point if sampling fails
    if (grid_count > 0) {
        int idx = (int)(rand01() * grid_count);
        *rx = grid_x[idx]; *ry = grid_y[idx];
    } else {
        *rx = 0; *ry = 0;
    }
}

// lumen model: lux on point from distance r (meters) using isotropic inverse-square
double illum_from_distance(double power, double r) {
    if (r < 1e-6) r = 1e-6;
    return power / (4.0 * PI * r * r);
}

// express locus per diploid dominance rules
int express_locus(const Locus *l, double *outx, double *outy) {
    Allele a = l->a1, b = l->a2;
    if (a.is_off && b.is_off) return 0;
    if (a.type == 'D' || b.type == 'D') {
        if (a.type == 'D' && !a.is_off && b.type != 'D') {
            *outx = a.x; *outy = a.y; return 1;
        }
        if (b.type == 'D' && !b.is_off && a.type != 'D') {
            *outx = b.x; *outy = b.y; return 1;
        }
        double sx = 0, sy = 0; int cnt = 0;
        if (!a.is_off) { sx += a.x; sy += a.y; cnt++; }
        if (!b.is_off) { sx += b.x; sy += b.y; cnt++; }
        if (cnt == 0) return 0;
        *outx = sx / cnt; *outy = sy / cnt; return 1;
    } else {
        if (!a.is_off) { *outx = a.x; *outy = a.y; return 1; }
        if (!b.is_off) { *outx = b.x; *outy = b.y; return 1; }
        return 0;
    }
}

// compute illuminance at all interior grid points for phenotype positions
void compute_illum_for_phenotype(double *posX, double *posY, int npos, double *out_illum) {
    for (int i = 0; i < grid_count; ++i) out_illum[i] = 0.0;
    for (int j = 0; j < npos; ++j) {
        double lx = posX[j], ly = posY[j];
        for (int i = 0; i < grid_count; ++i) {
            double dx = grid_x[i] - lx;
            double dy = grid_y[i] - ly;
            double horiz = sqrt(dx*dx + dy*dy);
            double r = sqrt(horiz*horiz + lantern_height_m*lantern_height_m);
            out_illum[i] += illum_from_distance(POWER_LM, r);
        }
    }
}

// fitness: maximize minimum lux across grid; mild penalty for number of lanterns
double fitness_of_chromosome(const Chromosome *c, int N_slots, double *out_minlux, int *out_lanterns, double *tmp_illum) {
    double posX[MAX_SLOTS], posY[MAX_SLOTS];
    int npos = 0;
    for (int i = 0; i < N_slots; ++i) {
        double px, py;
        if (express_locus(&c->loci[i], &px, &py)) {
            posX[npos] = px; posY[npos] = py; npos++;
        }
    }
    if (npos == 0) {
        *out_minlux = 0.0; *out_lanterns = 0; return -1e9;
    }
    compute_illum_for_phenotype(posX, posY, npos, tmp_illum);
    double minlux = tmp_illum[0];
    for (int i = 0; i < grid_count; ++i) if (tmp_illum[i] < minlux) minlux = tmp_illum[i];
    *out_minlux = minlux;
    *out_lanterns = npos;
    double fitness = minlux - 0.02 * npos;
    return fitness;
}

// initialize a random allele
Allele random_allele() {
    Allele a;
    if (rand01() < 0.25) {
        a.is_off = 1; a.x = a.y = 0; a.type = 'R';
    } else {
        double rx, ry;
        random_position_in_poly(&rx, &ry);
        a.x = rx; a.y = ry; a.is_off = 0;
        a.type = (rand01() < 0.35) ? 'D' : 'R';
    }
    return a;
}

// init random individual
void init_individual(Chromosome *c, int N_slots) {
    for (int i = 0; i < N_slots; ++i) {
        c->loci[i].a1 = random_allele();
        c->loci[i].a2 = random_allele();
    }
}
void copy_chromosome(const Chromosome *src, Chromosome *dst, int N_slots) { memcpy(dst, src, sizeof(Chromosome)); }

// tournament selection
void tournament_select(Chromosome *population, double *scores, int population_size, Chromosome *out) {
    int best = (int)(rand01() * population_size);
    for (int i = 1; i < TOURNAMENT_K; ++i) {
        int cand = (int)(rand01() * population_size);
        if (scores[cand] > scores[best]) best = cand;
    }
    copy_chromosome(&population[best], out, MAX_SLOTS);
}

// crossover and mutation helpers
void crossover_chromosome(const Chromosome *p1, const Chromosome *p2, Chromosome *child, int N_slots) {
    for (int i = 0; i < N_slots; ++i) {
        if (rand01() < 0.5) {
            if (rand01() < 0.5) child->loci[i] = p1->loci[i];
            else child->loci[i] = p2->loci[i];
        } else {
            if (rand01() < 0.5) child->loci[i].a1 = p1->loci[i].a1;
            else child->loci[i].a1 = p2->loci[i].a1;
            if (rand01() < 0.5) child->loci[i].a2 = p1->loci[i].a2;
            else child->loci[i].a2 = p2->loci[i].a2;
        }
    }
}

double gauss_rand() {
    double u = rand01(), v = rand01();
    if (u < 1e-12) u = 1e-12;
    return sqrt(-2.0 * log(u)) * cos(2.0 * PI * v);
}

void mutate_allele(Allele *a) {
    if (a->is_off) {
        if (rand01() < 0.03) {
            double rx, ry; random_position_in_poly(&rx, &ry);
            a->x = rx; a->y = ry; a->is_off = 0;
        }
    } else {
        if (rand01() < 0.7) {
            double move_sigma = 2.0;
            for (int trial = 0; trial < 12; ++trial) {
                double nx = a->x + gauss_rand() * move_sigma;
                double ny = a->y + gauss_rand() * move_sigma;
                if (nx < min_x || nx > max_x || ny < min_y || ny > max_y) continue;
                if (point_in_poly(nx, ny, verts_m, 4)) { a->x = nx; a->y = ny; break; }
            }
        }
        if (rand01() < 0.03) { a->is_off = 1; a->x = a->y = 0; }
    }
    if (rand01() < 0.01) a->type = (a->type == 'D') ? 'R' : 'D';
}
void mutate_chromosome(Chromosome *c, int N_slots) {
    for (int i = 0; i < N_slots; ++i) {
        double r = (c->loci[i].a1.type == 'D') ? MUTATION_RATE_DOM : MUTATION_RATE_REC;
        if (rand01() < r) mutate_allele(&c->loci[i].a1);
        r = (c->loci[i].a2.type == 'D') ? MUTATION_RATE_DOM : MUTATION_RATE_REC;
        if (rand01() < r) mutate_allele(&c->loci[i].a2);
    }
}

// run GA for N_slots
void run_ga_for_N(int N_slots, Chromosome *best_out, double *best_minlux_out, double *best_illum_out) {
    Chromosome *pop = malloc(sizeof(Chromosome) * POP_SIZE);
    Chromosome *newpop = malloc(sizeof(Chromosome) * POP_SIZE);
    double *scores = malloc(sizeof(double) * POP_SIZE);
    double *tmp_illum = malloc(sizeof(double) * grid_count);

    for (int i = 0; i < POP_SIZE; ++i) init_individual(&pop[i], N_slots);

    double best_score = -1e12, best_minlux = 0.0;
    Chromosome best_chrom;
    double *best_illum = malloc(sizeof(double) * grid_count);

    for (int gen = 0; gen < GENERATIONS; ++gen) {
        for (int i = 0; i < POP_SIZE; ++i) {
            double minlux; int lanterns;
            double fit = fitness_of_chromosome(&pop[i], N_slots, &minlux, &lanterns, tmp_illum);
            scores[i] = fit;
            if (fit > best_score) {
                best_score = fit; best_minlux = minlux;
                copy_chromosome(&pop[i], &best_chrom, N_slots);
                memcpy(best_illum, tmp_illum, sizeof(double) * grid_count);
            }
        }
        if ((gen % 50) == 0) printf("gen %d: best_score %.5f best_minlux %.4f\n", gen, best_score, best_minlux);

        int idxs[POP_SIZE];
        for (int i = 0; i < POP_SIZE; ++i) idxs[i] = i;
        for (int i = 0; i < ELITISM; ++i) {
            int besti = i;
            for (int j = i+1; j < POP_SIZE; ++j) if (scores[j] > scores[besti]) besti = j;
            int tmpi = idxs[i]; idxs[i] = idxs[besti]; idxs[besti] = tmpi;
        }
        for (int k = 0; k < ELITISM; ++k) copy_chromosome(&pop[idxs[k]], &newpop[k], N_slots);

        int cur = ELITISM;
        while (cur < POP_SIZE) {
            Chromosome parent1, parent2, child;
            tournament_select(pop, scores, POP_SIZE, &parent1);
            tournament_select(pop, scores, POP_SIZE, &parent2);
            if (rand01() < CROSSOVER_RATE) crossover_chromosome(&parent1, &parent2, &child, N_slots);
            else copy_chromosome(&parent1, &child, N_slots);
            mutate_chromosome(&child, N_slots);
            copy_chromosome(&child, &newpop[cur], N_slots);
            cur++;
        }
        memcpy(pop, newpop, sizeof(Chromosome) * POP_SIZE);
    }

    copy_chromosome(&best_chrom, best_out, N_slots);
    *best_minlux_out = best_minlux;
    memcpy(best_illum_out, best_illum, sizeof(double) * grid_count);

    free(pop); free(newpop); free(scores); free(tmp_illum); free(best_illum);
}

// prepare grid points inside polygon
void build_grid() {
    min_x =  1e9; min_y = 1e9; max_x = -1e9; max_y = -1e9;
    for (int i = 0; i < 4; ++i) {
        verts_m[i][0] = verts_ft[i][0] * FT_TO_M;
        verts_m[i][1] = verts_ft[i][1] * FT_TO_M;
        if (verts_m[i][0] < min_x) min_x = verts_m[i][0];
        if (verts_m[i][0] > max_x) max_x = verts_m[i][0];
        if (verts_m[i][1] < min_y) min_y = verts_m[i][1];
        if (verts_m[i][1] > max_y) max_y = verts_m[i][1];
    }

    int capacity = GRID_NX * GRID_NY;
    grid_x = (double*)malloc(sizeof(double) * capacity);
    grid_y = (double*)malloc(sizeof(double) * capacity);
    int count = 0;
    for (int iy = 0; iy < GRID_NY; ++iy) {
        double ty = (double)iy / (GRID_NY - 1);
        double y = min_y + ty * (max_y - min_y);
        for (int ix = 0; ix < GRID_NX; ++ix) {
            double tx = (double)ix / (GRID_NX - 1);
            double x = min_x + tx * (max_x - min_x);
            if (point_in_poly(x, y, verts_m, 4)) {
                grid_x[count] = x; grid_y[count] = y; count++;
            }
        }
    }
    grid_count = count;
    printf("Grid built: %d interior sample points (NX=%d NY=%d)\n", grid_count, GRID_NX, GRID_NY);
}

int main(int argc, char **argv) {
    rng_seed(RANDOM_SEED);
    lantern_height_m = LANTERN_HEIGHT_FT * FT_TO_M;

    build_grid();
    if (grid_count <= 0) {
        fprintf(stderr, "No interior grid points found. Exiting.\n");
        return 1;
    }

    Chromosome best_chrom;
    double best_minlux;
    double *best_illum = (double*)malloc(sizeof(double) * grid_count);

    int N_slots = 5;
    printf("Enter number of lantern slots to optimize (1..%d), or 0 to use default %d: ", MAX_SLOTS, N_slots);
    int inputN = 0;
    if (scanf("%d", &inputN) == 1) {
        if (inputN >= 1 && inputN <= MAX_SLOTS) N_slots = inputN;
    }

    printf("Running GA for N_slots=%d ... (POP=%d, GENS=%d)\n", N_slots, POP_SIZE, GENERATIONS);
    run_ga_for_N(N_slots, &best_chrom, &best_minlux, best_illum);

    double posX[MAX_SLOTS], posY[MAX_SLOTS];
    int npos = 0;
    for (int i = 0; i < N_slots; ++i) {
        double px, py;
        if (express_locus(&best_chrom.loci[i], &px, &py)) {
            posX[npos] = px; posY[npos] = py; npos++;
        }
    }

    int covered = 0;
    for (int i = 0; i < grid_count; ++i) if (best_illum[i] >= TARGET_LUX) covered++;
    double coverage_pct = 100.0 * covered / grid_count;

    printf("\nBest found (min lux): %.4f lux | coverage >= %.2f lux: %.2f%% | lanterns used: %d\n",
        best_minlux, TARGET_LUX, coverage_pct, npos);

    printf("Lantern positions (feet):\n");
    for (int i = 0; i < npos; ++i) {
        printf("  %2d: x = %8.2f ft, y = %8.2f ft\n", i+1, posX[i] / FT_TO_M, posY[i] / FT_TO_M);
    }

#if SAVE_CSV_HEATMAP
    {
        FILE *f = fopen("heatmap.csv", "w");
        if (f) {
            fprintf(f, "x_ft,y_ft,lux\n");
            for (int i = 0; i < grid_count; ++i) {
                fprintf(f, "%.6f,%.6f,%.6f\n", grid_x[i] / FT_TO_M, grid_y[i] / FT_TO_M, best_illum[i]);
            }
            fclose(f);
            printf("Wrote heatmap CSV to heatmap.csv (%d points)\n", grid_count);
        }
    }
#endif

    free(grid_x); free(grid_y); free(best_illum);
    return 0;
}

