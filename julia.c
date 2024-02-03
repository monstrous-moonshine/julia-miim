#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <sys/stat.h>

#define IMG_RES      800
#define PRE_ITER_CNT 200
#define IIM_ITER_CNT 2000000000
#define QUEUE_SIZE   (64 * 1024 * 1024)
#define MAX_HIT      1000
#define ORBIT_ITER_CNT 10000000
// Antialiasing
#define ALIAS_LEN 2
#define ALIAS_PTS (ALIAS_LEN * ALIAS_LEN)
#define RAW_PIX_CNT (IMG_RES * IMG_RES * ALIAS_PTS)

// Notes on some of the parameters above:
//
// Arguably the most important parameter is MAX_HIT, which is the number of
// times a pixel must be hit to stop being considered for further iteration.
// This cutoff is the modification that gives the process the name "modified"
// inverse iteration method.
//
// Getting the "lobes" to touch is very difficult for parabolic parameters.
// Even with MAX_HIT = 1000, there are large gaps; but only so much computing
// power can be lavished on it, and there's only so much that brute force can
// accomplish. To get good results for these cases, a better method is required,
// like distance estimation. Execution time increases rapidly with increasing
// MAX_HIT.
//
// The QUEUE_SIZE should be big enough that points encountered in the breadth
// first search are not lost unnecessarily because of space constraint, as much
// as possible.
//
// For IIM_ITER_COUNT, bigger values are obviously better. For a particular
// resolution, MAX_HIT, QUEUE_SIZE, and possibly starting value, I suspect the
// process saturates at some point. The IIM_ITER_COUNT should be big enough to
// not cut off the process too long before that. Otherwise, execution time does
// not seem to increase dramatically with this parameter, so it can be set to a
// quite high value.

#define MIN(a, b) ((a) < (b) ? a : (b))
#define ARRAY_LEN(a) ((sizeof a) / (sizeof a[0]))

static void write_pgm(uint8_t *image, const char *name) {
    FILE *fp = fopen(name, "w");
    if (!fp) {
        perror("fopen");
        return;
    }
    fprintf(fp, "P5\n%u %u\n255\n", IMG_RES, IMG_RES);
    fwrite(image, 1, IMG_RES * IMG_RES, fp);
    fclose(fp);
}

static void antialias(uint32_t *hit_counts, uint8_t *image) {
    for (int i = 0; i < IMG_RES; i++) {
        for (int j = 0; j < IMG_RES; j++) {
            int idx = i * IMG_RES * ALIAS_PTS + j * ALIAS_LEN;
            // "unroll" a 2 * 2 loop
            int sum = (hit_counts[idx] > 0) + (hit_counts[idx + 1] > 0);
            idx += IMG_RES * ALIAS_LEN;
            sum += (hit_counts[idx] > 0) + (hit_counts[idx + 1] > 0);
            // done
            int value = MIN(sum * 64, 255);
            image[i * IMG_RES + j] = value;
        }
    }
}

static int random_bit() {
    // store: reserve of random bits
    // bit_count: number of random bits remaining in store
    static int store = 0, bit_count = 0;
    if (bit_count == 0) {
        store = rand();
        bit_count = 32;
    }
    int bit = store & 1;
    store >>= 1;
    bit_count--;
    return bit;
}

static double complex pre_iterate(double complex c) {
    double complex z = 0;
    for (int i = 0; i < PRE_ITER_CNT; i++) {
        z = sqrt(z - c);
        z = random_bit() ? -z : z;
    }
    return z;
}

// ---------- Simple circular queue implementation ----------

static struct {
    int head, tail;
    double complex values[QUEUE_SIZE];
} iim_queue = {.head = 0, .tail = 0};

static void queue_push(double complex z) {
    iim_queue.values[iim_queue.head++] = z;
    iim_queue.head %= QUEUE_SIZE;
}

static double complex queue_pop() {
    double complex value = iim_queue.values[iim_queue.tail++];
    iim_queue.tail %= QUEUE_SIZE;
    return value;
}

static int queue_empty() {
    return iim_queue.head == iim_queue.tail;
}

static int queue_full() {
    int dist = iim_queue.tail - iim_queue.head;
    return dist == 1 || dist == QUEUE_SIZE - 1;
}

// ----------------------------------------------------------

static int c2idx(double complex z) {
    double x = creal(z), y = cimag(z);
    int col = (x + 2) / 4 * IMG_RES * ALIAS_LEN;
    int row = (2 - y) / 4 * IMG_RES * ALIAS_LEN;
    int idx = (row * IMG_RES * ALIAS_LEN) + col;
    return idx;
}

static void handle_value(uint32_t *hit_counts, double complex z) {
    int idx = c2idx(z);
    if (!queue_full() && hit_counts[idx] < MAX_HIT)
        queue_push(z);
}

static int julia_miim(uint32_t *hit_counts, double complex c) {
    queue_push(pre_iterate(c));
    int i;
    for (i = 0; i < IIM_ITER_CNT; i++) {
        if (queue_empty()) break;
        double complex p = queue_pop();
        int idx = c2idx(p);
        hit_counts[idx]++;
        double complex z = sqrt(p - c);
        handle_value(hit_counts, z);
        handle_value(hit_counts, -z);
    }
    return i;
}

static void draw_orbit(uint32_t *hit_counts, double complex c) {
    double complex z = 0;
    for (int i = 0; i < ORBIT_ITER_CNT; i++) {
        int idx = c2idx(z);
        hit_counts[idx]++;
        z = z * z + c;
    }
}

int main() {
    int ret = mkdir("images", 0777);
    if (ret == -1 && errno != EEXIST) {
        perror("mkdir");
        return -1;
    }
    uint32_t *hit_counts = malloc(RAW_PIX_CNT * sizeof *hit_counts);
    uint8_t *image = malloc(IMG_RES * IMG_RES * sizeof *image);
    if (!hit_counts || !image) {
        perror("malloc");
        return -1;
    }
    srand(42);
    double farey[] = {0./1,
                                              // 1./5,
                                        1./4,
                                              // 2./7,
                                  1./3,
                                              // 3./8,
                                        2./5,
                                              // 3./7,
                            1./2,
                                              // 4./7,
                                        3./5,
                                              // 5./8, 
                                  2./3,
                                              // 5./7,
                                        3./4,
                                              // 4./5,
                      1./1};
    for (unsigned i = 0; i < ARRAY_LEN(farey); i++) {
        char name[32];
        printf("\r%2d/%lu", i, ARRAY_LEN(farey));
        fflush(stdout);
        memset(hit_counts, 0, RAW_PIX_CNT * sizeof *hit_counts);
        // Scan period-2 lobe boundary
        // double complex param = -1 + 0.25 * exp(2 * M_PI * I * farey[i]);
        // Or, scan main cardioid boundary
        double complex z = 0.5 * exp(2 * M_PI * I * farey[i]);
        double complex param = z * (1 - z);
        julia_miim(hit_counts, param);
        draw_orbit(hit_counts, param);
        antialias(hit_counts, image);
        sprintf(name, "images/out%02d.pgm", i);
        write_pgm(image, name);
    }
    printf("\n");
    free(hit_counts);
    free(image);
}
