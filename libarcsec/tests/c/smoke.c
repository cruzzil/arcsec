/*
 * End-to-end test of libarcsec from C: every entry point, on a synthetic field.
 *
 * Built and run by tests/c_api.rs, which makes the field first:
 *
 *   smoke <catalog_dir> <image.f32> <image.fits> <width> <height> <ra> <dec> <scale>
 *
 * image.f32 is the raw image (native float, row 0 = FITS row 1), image.fits the same
 * pixels as a FITS file with RA/DEC/FOCALLEN/XPIXSZ keywords, the catalogue
 * directory holds a d50 database of the sky around (ra, dec), and scale is the
 * pixel scale in arcseconds. Exits 0 when every check passes. Written in C99 so it
 * builds with MSVC as well as gcc and clang, and so that the sanitizer and
 * valgrind runs see only plain C on this side.
 */

#include "arcsec.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond)                                                                    \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            fprintf(stderr, "FAIL %s:%d: %s (last error: %s)\n", __FILE__, __LINE__,   \
                    #cond, arcsec_last_error());                                       \
            failures++;                                                                \
        }                                                                              \
    } while (0)

static float *read_raw(const char *path, size_t n) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    float *buf = (float *)malloc(n * sizeof(float));
    if (buf && fread(buf, sizeof(float), n, f) != n) {
        free(buf);
        buf = NULL;
    }
    fclose(f);
    return buf;
}

/* The callback can run on any of the solver's threads; this test keeps it to one
 * (threads = 1, and it is unregistered before the multi-threaded solves), so a
 * plain int will do. */
static int log_lines = 0;
static void on_log(void *user, int level, const char *message) {
    (void)user;
    if (level >= ARCSEC_LOG_ERROR && level <= ARCSEC_LOG_TRACE && message && message[0])
        log_lines++;
}

static int progress_calls = 0;
static void on_progress(void *user, double fraction, const char *stage) {
    (void)user;
    if (stage && stage[0] && fraction <= 1.0) progress_calls++;
}

static int cancel_always(void *user) {
    (*(int *)user)++;
    return 1;
}

static double angle_diff(double a, double b) {
    double d = fmod(fabs(a - b), 360.0);
    return d > 180.0 ? 360.0 - d : d;
}

int main(int argc, char **argv) {
    if (argc != 9) {
        fprintf(stderr, "usage: %s catalog image.f32 image.fits width height ra dec scale\n",
                argv[0]);
        return 2;
    }
    const char *catalog = argv[1];
    const char *raw_path = argv[2];
    const char *fits_path = argv[3];
    const uint32_t width = (uint32_t)atoi(argv[4]);
    const uint32_t height = (uint32_t)atoi(argv[5]);
    const double ra = atof(argv[6]), dec = atof(argv[7]), scale = atof(argv[8]);
    const size_t npix = (size_t)width * height;

    /* ── Version ────────────────────────────────────────────────────────────── */
    CHECK(arcsec_abi_version() == ARCSEC_ABI_VERSION);
    CHECK(arcsec_version() != NULL && strlen(arcsec_version()) >= 5);
    printf("libarcsec %s, ABI %u\n", arcsec_version(), arcsec_abi_version());
    CHECK(strcmp(arcsec_status_string(ARCSEC_NO_SOLUTION), "no solution") == 0);
    CHECK(arcsec_last_error() != NULL);

    arcsec_set_log_callback(on_log, NULL, ARCSEC_LOG_INFO);

    float *pixels = read_raw(raw_path, npix);
    CHECK(pixels != NULL);
    if (!pixels) return 1;

    arcsec_solver *solver = arcsec_solver_new();
    CHECK(solver != NULL);

    /* ── Solve the buffer ───────────────────────────────────────────────────── */
    arcsec_image image;
    memset(&image, 0, sizeof image);
    image.struct_size = sizeof image;
    image.data = pixels;
    image.pixel_type = ARCSEC_PIXEL_F32;
    image.width = width;
    image.height = height;
    image.channels = 1;

    arcsec_solve_options opts;
    CHECK(arcsec_solve_options_init(&opts) == ARCSEC_OK);
    CHECK(opts.struct_size == sizeof opts);
    CHECK(opts.max_stars == 500 && opts.search_radius_deg == 180.0);
    opts.has_hint = 1;
    opts.ra_deg = ra + 0.2;
    opts.dec_deg = dec - 0.1;
    opts.search_radius_deg = 3.0;
    opts.pixel_scale_arcsec = scale;
    opts.catalog_dir = catalog;
    opts.sip_order = 3;
    opts.threads = 1;
    opts.progress = on_progress;

    arcsec_result *result = NULL;
    arcsec_status st = arcsec_solve(solver, &image, &opts, &result);
    CHECK(st == ARCSEC_OK);
    CHECK(result != NULL);
    arcsec_wcs wcs;
    memset(&wcs, 0, sizeof wcs);
    if (result) {
        wcs.struct_size = sizeof wcs;
        CHECK(arcsec_result_wcs(result, &wcs) == ARCSEC_OK);
        printf("solved: CRVAL %.5f %.5f, CRPIX %.1f %.1f, %.3f\"/px, %u stars, rms %.2f\", "
               "SIP order %d, database %s\n",
               wcs.crval1, wcs.crval2, wcs.crpix1, wcs.crpix2, wcs.pixel_scale_arcsec,
               wcs.matched_stars, wcs.rms_arcsec, wcs.sip_order,
               arcsec_result_database(result));
        CHECK(angle_diff(wcs.crval1, ra) < 0.01 && fabs(wcs.crval2 - dec) < 0.01);
        CHECK(fabs(wcs.pixel_scale_arcsec - scale) < 0.01);
        CHECK(wcs.matched_stars >= 10);
        CHECK(wcs.cdelt1 < 0.0 && wcs.cdelt2 > 0.0 && !wcs.mirrored);
        CHECK(fabs(wcs.crota1 - wcs.crota2) < 0.1);
        CHECK(wcs.sip_order == 0 || wcs.sip_order == 3);
        CHECK(strcmp(arcsec_result_database(result), "d50") == 0);

        arcsec_solve_info info;
        memset(&info, 0, sizeof info);
        info.struct_size = sizeof info;
        CHECK(arcsec_result_info(result, &info) == ARCSEC_OK);
        CHECK(info.binning == 1 && info.elapsed_seconds > 0.0);

        size_t n = arcsec_result_matched_stars(result, NULL, 0);
        CHECK(n == wcs.matched_stars);
        arcsec_matched_star *pairs = (arcsec_matched_star *)calloc(n, sizeof *pairs);
        CHECK(arcsec_result_matched_stars(result, pairs, n) == n);
        for (size_t i = 0; i < n; i++) {
            double x, y;
            CHECK(arcsec_result_sky_to_pixel(result, pairs[i].ra_deg, pairs[i].dec_deg, &x,
                                             &y) == ARCSEC_OK);
            CHECK(fabs(x - pairs[i].x) < 2.0 && fabs(y - pairs[i].y) < 2.0);
        }
        free(pairs);

        size_t hlen = arcsec_result_fits_header(result, NULL, 0);
        char *cards = (char *)malloc(hlen + 1);
        CHECK(arcsec_result_fits_header(result, cards, hlen + 1) == hlen);
        CHECK(hlen % 80 == 0 && strlen(cards) == hlen);
        CHECK(strncmp(cards, "CTYPE1  = 'RA---TAN", 19) == 0);
        CHECK(strstr(cards, "CRVAL1  =") != NULL && strncmp(cards + hlen - 80, "END ", 4) == 0);
        free(cards);

        double r2, d2;
        CHECK(arcsec_result_pixel_to_sky(result, wcs.crpix1, wcs.crpix2, &r2, &d2) ==
              ARCSEC_OK);
        CHECK(angle_diff(r2, wcs.crval1) < 1e-6 && fabs(d2 - wcs.crval2) < 1e-6);
    }
    arcsec_result_free(result);
    result = NULL;
    CHECK(log_lines > 0);
    CHECK(progress_calls > 0);
    arcsec_set_log_callback(NULL, NULL, ARCSEC_LOG_OFF);
    opts.progress = NULL;
    opts.threads = 0; /* every core from here on */

    /* ── The same image as 16-bit planes, top-down, with the hint in a header ── */
    uint16_t *u16 = (uint16_t *)malloc(npix * 2 * sizeof(uint16_t));
    for (uint32_t y = 0; y < height; y++)
        for (uint32_t x = 0; x < width; x++) {
            float v = pixels[(size_t)(height - 1 - y) * width + x];
            uint16_t s = (uint16_t)(v < 0 ? 0 : v > 65535 ? 65535 : v + 0.5f);
            u16[(size_t)y * width + x] = s;        /* plane 0 */
            u16[npix + (size_t)y * width + x] = s; /* plane 1: the same */
        }
    const void *planes[2] = {u16, u16 + npix};
    arcsec_image image16 = image;
    image16.data = NULL;
    image16.planes = planes;
    image16.pixel_type = ARCSEC_PIXEL_U16;
    image16.channels = 2;
    image16.flags = ARCSEC_IMAGE_TOP_DOWN;
    char header[5 * 80 + 1];
    snprintf(header, sizeof header, "%-80s%-80s%-80s%-80s%-80s", "RA      = 84.45",
             "DEC     = -5.25", "FOCALLEN= 206.265", "XPIXSZ  = 5.0", "END");
    arcsec_solve_options o16 = opts;
    o16.has_hint = 0;
    o16.pixel_scale_arcsec = 0.0;
    o16.fits_header = header;
    o16.sip_order = 1; /* linear, as Siril means it */
    st = arcsec_solve(solver, &image16, &o16, &result);
    CHECK(st == ARCSEC_OK);
    if (result) {
        arcsec_wcs w16;
        memset(&w16, 0, sizeof w16);
        w16.struct_size = sizeof w16;
        CHECK(arcsec_result_wcs(result, &w16) == ARCSEC_OK);
        CHECK(angle_diff(w16.crval1, wcs.crval1) < 1e-4 && fabs(w16.crval2 - wcs.crval2) < 1e-4);
        CHECK(w16.sip_order == 0);
    }
    arcsec_result_free(result);
    result = NULL;
    free(u16);

    /* ── Solve the FITS file, hint and scale from its header ─────────────────── */
    arcsec_solve_options ofile;
    arcsec_solve_options_init(&ofile);
    ofile.catalog_dir = catalog;
    ofile.search_radius_deg = 3.0;
    st = arcsec_solve_file(solver, fits_path, &ofile, &result);
    CHECK(st == ARCSEC_OK);
    arcsec_result_free(result);
    result = NULL;

    /* ── Analysis ───────────────────────────────────────────────────────────── */
    arcsec_analysis analysis;
    memset(&analysis, 0, sizeof analysis);
    analysis.struct_size = sizeof analysis;
    arcsec_star stars[8];
    CHECK(arcsec_analyse(&image, 0.0, 0, &analysis, stars, 8) == ARCSEC_OK);
    printf("analysis: %llu stars, median HFD %.2f px\n",
           (unsigned long long)analysis.star_count, analysis.hfd_median);
    CHECK(analysis.star_count > 20 && analysis.hfd_median > 1.0);
    CHECK(stars[0].hfd > 0.0 && stars[0].x >= 1.0);

    /* ── Catalogue queries ──────────────────────────────────────────────────── */
    char name[16];
    CHECK(arcsec_has_star_database(catalog) == 1);
    CHECK(arcsec_has_blind_index(catalog) == 0);
    CHECK(arcsec_select_database(catalog, 0.5, name, sizeof name) == 3);
    CHECK(strcmp(name, "d50") == 0);
    size_t need = arcsec_default_catalog_dir(NULL, 0);
    CHECK(need > 0);
    char *dir = (char *)malloc(need + 1);
    CHECK(arcsec_default_catalog_dir(dir, need + 1) == need && strlen(dir) == need);
    free(dir);

    /* ── Failures: codes, messages, and no result ───────────────────────────── */
    arcsec_solve_options bad = opts;
    bad.catalog_dir = "/nonexistent/arcsec/catalogs";
    CHECK(arcsec_solve(solver, &image, &bad, &result) == ARCSEC_DATABASE_NOT_FOUND);
    CHECK(result == NULL && strlen(arcsec_last_error()) > 0);
    printf("expected failure: %s\n", arcsec_last_error());

    CHECK(arcsec_solve_file(solver, "/nonexistent/arcsec/image.fits", &opts, &result) ==
          ARCSEC_FILE_ERROR);
    CHECK(arcsec_solve(NULL, &image, &opts, &result) == ARCSEC_INVALID_ARGUMENT);
    CHECK(arcsec_solve(solver, NULL, &opts, &result) == ARCSEC_INVALID_ARGUMENT);
    CHECK(arcsec_solve(solver, &image, &opts, NULL) == ARCSEC_INVALID_ARGUMENT);
    arcsec_solve_options uninitialised;
    memset(&uninitialised, 0, sizeof uninitialised);
    CHECK(arcsec_solve(solver, &image, &uninitialised, &result) == ARCSEC_INVALID_ARGUMENT);
    arcsec_image bad_image = image;
    bad_image.pixel_type = 42;
    CHECK(arcsec_solve(solver, &bad_image, &opts, &result) == ARCSEC_INVALID_ARGUMENT);
    bad_image = image;
    bad_image.row_stride = 3; /* less than a row */
    CHECK(arcsec_solve(solver, &bad_image, &opts, &result) == ARCSEC_INVALID_ARGUMENT);
    arcsec_wcs small;
    small.struct_size = 8;
    CHECK(arcsec_result_wcs(NULL, &small) == ARCSEC_INVALID_ARGUMENT);

    /* Searched only far from the field: no solution. */
    arcsec_solve_options far = opts;
    far.ra_deg = 200.0;
    far.dec_deg = 40.0;
    far.search_radius_deg = 0.0;
    far.auto_index = 0;
    CHECK(arcsec_solve(solver, &image, &far, &result) == ARCSEC_NO_SOLUTION);

    /* Cancelled by the callback. */
    int polls = 0;
    arcsec_solve_options cancelled = opts;
    cancelled.cancel = cancel_always;
    cancelled.cancel_user = &polls;
    cancelled.threads = 1;
    CHECK(arcsec_solve(solver, &image, &cancelled, &result) == ARCSEC_CANCELLED);
    CHECK(polls >= 1 && result == NULL);
    arcsec_solver_cancel(solver); /* nothing running: a no-op */

    /* NULL is safe to free. */
    arcsec_result_free(NULL);
    arcsec_solver_free(NULL);

    arcsec_solver_free(solver);
    free(pixels);

    if (failures) {
        fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    printf("all C API checks passed\n");
    return 0;
}
