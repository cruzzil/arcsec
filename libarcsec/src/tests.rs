//! Tests of the C entry points, called from Rust as C would call them.

use alloc::ffi::CString;
use alloc::sync::Arc;
use core::ffi::{CStr, c_char, c_int, c_void};
use core::ptr;
use core::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use arcsec_core::ImageBuffer;
use arcsec_core::test_support::{TempDir, TruthWcs, synthetic_field};

use crate::arcsec_status::{
    ARCSEC_BUSY, ARCSEC_CANCELLED, ARCSEC_DATABASE_NOT_FOUND, ARCSEC_FILE_ERROR,
    ARCSEC_INSUFFICIENT_STARS, ARCSEC_INTERNAL_ERROR, ARCSEC_INVALID_ARGUMENT, ARCSEC_OK,
};
use crate::*;

fn deg(d: f64) -> f64 {
    d.to_radians()
}

fn last_error() -> String {
    // SAFETY: arcsec_last_error never returns NULL.
    unsafe { CStr::from_ptr(arcsec_last_error()) }
        .to_string_lossy()
        .into_owned()
}

/// A synthetic 400×320 field at 5″/px around (84.3°, -5.2°), with its `d50`
/// database.
struct Field {
    dir: TempDir,
    img: ImageBuffer,
    truth: TruthWcs,
}

fn field() -> Field {
    let truth = TruthWcs::new(deg(84.3), deg(-5.2), 5.0, 23.0, false, 400, 320);
    let dir = TempDir::new("capi");
    let img = synthetic_field(dir.path(), "d50", &truth, 130, 1);
    Field { dir, img, truth }
}

fn c_path(p: &std::path::Path) -> CString {
    CString::new(p.to_str().unwrap()).unwrap()
}

fn options(catalog: &CString) -> arcsec_solve_options {
    let mut o = core::mem::MaybeUninit::<arcsec_solve_options>::uninit();
    // SAFETY: a writable struct.
    assert_eq!(
        unsafe { arcsec_solve_options_init(o.as_mut_ptr()) },
        ARCSEC_OK
    );
    // SAFETY: initialised just above.
    let mut o = unsafe { o.assume_init() };
    o.has_hint = 1;
    o.ra_deg = 84.3 + 0.2;
    o.dec_deg = -5.2 - 0.1;
    o.search_radius_deg = 2.0;
    o.pixel_scale_arcsec = 5.0;
    o.catalog_dir = catalog.as_ptr();
    o.threads = 2;
    o
}

fn f32_image(img: &ImageBuffer) -> arcsec_image {
    arcsec_image {
        struct_size: core::mem::size_of::<arcsec_image>(),
        data: img.data.as_ptr().cast(),
        planes: ptr::null(),
        pixel_type: arcsec_pixel_type::ARCSEC_PIXEL_F32 as u32,
        width: u32::try_from(img.width).unwrap(),
        height: u32::try_from(img.height).unwrap(),
        channels: 1,
        row_stride: 0,
        plane_stride: 0,
        flags: 0,
        reserved: 0,
    }
}

/// Solve and return the result, asserting the status.
fn solve(
    solver: *mut arcsec_solver,
    image: &arcsec_image,
    opts: &arcsec_solve_options,
    want: arcsec_status,
) -> *mut arcsec_result {
    let mut out: *mut arcsec_result = ptr::null_mut();
    // SAFETY: valid solver, image, options and out-pointer.
    let st = unsafe { arcsec_solve(solver, image, opts, &raw mut out) };
    assert_eq!(st, want, "{}", last_error());
    assert_eq!(out.is_null(), want != ARCSEC_OK);
    out
}

fn wcs_of(r: *const arcsec_result) -> arcsec_wcs {
    let mut w = core::mem::MaybeUninit::<arcsec_wcs>::zeroed();
    // SAFETY: zeroed is a valid arcsec_wcs.
    let mut w = unsafe { w.assume_init_mut() }.to_owned();
    w.struct_size = core::mem::size_of::<arcsec_wcs>();
    // SAFETY: live result, sized output.
    assert_eq!(unsafe { arcsec_result_wcs(r, &raw mut w) }, ARCSEC_OK);
    w
}

#[test]
fn version_and_status_strings() {
    // SAFETY: static NUL-terminated strings.
    let v = unsafe { CStr::from_ptr(arcsec_version()) }
        .to_str()
        .unwrap();
    assert_eq!(v, env!("CARGO_PKG_VERSION"));
    assert_eq!(arcsec_abi_version(), ARCSEC_ABI_VERSION);
    for (code, text) in [
        (0, "ok"),
        (1, "no solution"),
        (101, "cancelled"),
        (-5, "unknown status"),
    ] {
        // SAFETY: static NUL-terminated strings.
        let s = unsafe { CStr::from_ptr(arcsec_status_string(code)) };
        assert_eq!(s.to_str().unwrap(), text);
    }
}

#[test]
fn the_soname_follows_the_abi_version() {
    assert_eq!(env!("ARCSEC_BUILD_ABI"), ARCSEC_ABI_VERSION.to_string());
}

#[test]
fn a_synthetic_field_solves_through_the_c_api() {
    let f = field();
    let catalog = c_path(f.dir.path());
    let opts = options(&catalog);
    let solver = arcsec_solver_new();
    assert!(!solver.is_null());

    let r = solve(solver, &f32_image(&f.img), &opts, ARCSEC_OK);
    let w = wcs_of(r);
    let (ra, dec) = f.truth.pixel_to_sky(
        (f.img.width as f64 - 1.0) / 2.0,
        (f.img.height as f64 - 1.0) / 2.0,
    );
    // CRPIX is the image centre, so CRVAL is the truth there.
    assert!(
        (w.crval1 - ra.to_degrees()).abs() < 2e-4,
        "{} vs {}",
        w.crval1,
        ra.to_degrees()
    );
    assert!((w.crval2 - dec.to_degrees()).abs() < 2e-4);
    assert!((w.pixel_scale_arcsec - 5.0).abs() < 0.01);
    assert!(w.matched_stars >= 10 && w.mirrored == 0);
    assert!(w.cdelt1 < 0.0);
    assert!((w.crota2 - 23.0).abs() < 0.1, "crota2 {}", w.crota2);

    // The database was chosen by field size.
    // SAFETY: live result; the string lives as long as it.
    let db = unsafe { CStr::from_ptr(arcsec_result_database(r)) };
    assert_eq!(db.to_str().unwrap(), "d50");

    // Pixel <-> sky through the solution.
    let (mut ra2, mut dec2, mut x, mut y) = (0.0, 0.0, 0.0, 0.0);
    // SAFETY: live result and writable outputs.
    unsafe {
        assert_eq!(
            arcsec_result_pixel_to_sky(r, w.crpix1, w.crpix2, &raw mut ra2, &raw mut dec2),
            ARCSEC_OK
        );
        assert_eq!(
            arcsec_result_sky_to_pixel(r, ra2, dec2, &raw mut x, &raw mut y),
            ARCSEC_OK
        );
    }
    assert!((ra2 - w.crval1).abs() < 1e-9 && (dec2 - w.crval2).abs() < 1e-9);
    assert!((x - w.crpix1).abs() < 1e-6 && (y - w.crpix2).abs() < 1e-6);

    // The FITS cards, read back as a header.
    // SAFETY: live result; NULL buffer with length 0 sizes it.
    let need = unsafe { arcsec_result_fits_header(r, ptr::null_mut(), 0) };
    assert_eq!(need % 80, 0);
    let mut cards = vec![0 as c_char; need + 1];
    // SAFETY: room for need + 1 bytes.
    assert_eq!(
        unsafe { arcsec_result_fits_header(r, cards.as_mut_ptr(), cards.len()) },
        need
    );
    // SAFETY: NUL-terminated by the call.
    let text = unsafe { CStr::from_ptr(cards.as_ptr()) }
        .to_str()
        .unwrap()
        .to_owned();
    assert!(text.starts_with("CTYPE1  = 'RA---TAN'"), "{text}");
    assert!(text.ends_with(&format!("{:<80}", "END")));
    let h = arcsec_io::header::HeaderCards::parse(&text);
    assert!((h.number("CRVAL1").unwrap() - w.crval1).abs() < 1e-9);
    assert!((h.number("CD2_1").unwrap() - w.cd2_1).abs() < 1e-15);
    assert_eq!(h.raw("PLTSOLVD"), Some("T"));

    // Matched stars: the count, then the pairs.
    // SAFETY: live result; NULL output with capacity 0.
    let n = unsafe { arcsec_result_matched_stars(r, ptr::null_mut(), 0) };
    assert_eq!(n, w.matched_stars as usize);
    let mut pairs = vec![arcsec_matched_star::default(); n + 3];
    // SAFETY: room for n + 3.
    let n2 = unsafe { arcsec_result_matched_stars(r, pairs.as_mut_ptr(), pairs.len()) };
    assert_eq!(n2, n);
    assert!(pairs[..n].iter().all(|p| p.x >= 0.5 && p.ra_deg > 80.0));

    let mut info = arcsec_solve_info {
        struct_size: core::mem::size_of::<arcsec_solve_info>(),
        ..crate::util::Versioned::defaults()
    };
    // SAFETY: live result, sized output.
    assert_eq!(unsafe { arcsec_result_info(r, &raw mut info) }, ARCSEC_OK);
    assert_eq!(info.binning, 1);
    assert!(info.search_distance_deg >= 0.0 && info.has_index_estimate == 0);
    assert!(info.elapsed_seconds > 0.0 && info.mag_limit > 10.0);
    // SAFETY: freeing what the solve returned, once; NULL is ignored.
    unsafe {
        arcsec_result_free(r);
        arcsec_result_free(ptr::null_mut());
    }

    // The same pixels as u16, upside down and flagged as such, with the hint and
    // scale in a FITS header instead: the same solution.
    let rows: Vec<u16> = (0..f.img.height)
        .rev()
        .flat_map(|y| f.img.data[y * f.img.width..(y + 1) * f.img.width].iter())
        .map(|&v| v.round().clamp(0.0, 65535.0) as u16)
        .collect();
    let mut img16 = f32_image(&f.img);
    img16.data = rows.as_ptr().cast();
    img16.pixel_type = arcsec_pixel_type::ARCSEC_PIXEL_U16 as u32;
    img16.flags = ARCSEC_IMAGE_TOP_DOWN;
    let header = CString::new(format!(
        "{:<80}{:<80}{:<80}{:<80}",
        "RA      = 84.5", "DEC     = -5.3", "FOCALLEN= 206.265", "XPIXSZ  = 5.0"
    ))
    .unwrap();
    let mut o2 = opts;
    o2.has_hint = 0;
    o2.pixel_scale_arcsec = 0.0;
    o2.fits_header = header.as_ptr();
    o2.sip_order = 3;
    let r2 = solve(solver, &img16, &o2, ARCSEC_OK);
    let w2 = wcs_of(r2);
    assert!((w2.crval1 - w.crval1).abs() < 5e-5 && (w2.crval2 - w.crval2).abs() < 5e-5);
    assert!(
        (w2.cd1_1 - w.cd1_1).abs() < 1e-6 && (w2.cd2_1 - w.cd2_1).abs() < 1e-6,
        "{} {} / {} {}",
        w2.cd1_1,
        w.cd1_1,
        w2.cd2_1,
        w.cd2_1
    );
    // An undistorted field: SIP is either absent or negligible.
    assert!(w2.sip_order == 0 || w2.sip_order == 3);
    // SAFETY: freeing once.
    unsafe {
        arcsec_result_free(r2);
        arcsec_solver_free(solver);
    }
}

#[test]
fn a_fits_file_solves_with_its_own_header() {
    let f = field();
    let catalog = c_path(f.dir.path());
    let path = f.dir.path().join("field.fits");
    std::fs::write(
        &path,
        arcsec_core::test_support::fits_f32_bytes(
            &f.img,
            &[
                ("RA", "84.45"),
                ("DEC", "-5.25"),
                ("FOCALLEN", "206.265"),
                ("XPIXSZ", "5.0"),
            ],
        ),
    )
    .unwrap();
    let mut opts = options(&catalog);
    opts.has_hint = 0;
    opts.pixel_scale_arcsec = 0.0;
    let solver = arcsec_solver_new();
    let cpath = c_path(&path);
    let mut out = ptr::null_mut();
    // SAFETY: valid arguments.
    let st = unsafe { arcsec_solve_file(solver, cpath.as_ptr(), &raw const opts, &raw mut out) };
    assert_eq!(st, ARCSEC_OK, "{}", last_error());
    assert!(wcs_of(out).matched_stars >= 10);

    // A missing file and a non-image are file errors.
    let missing = c_path(&f.dir.path().join("nope.fits"));
    let mut out2 = ptr::null_mut();
    // SAFETY: valid arguments.
    let st = unsafe { arcsec_solve_file(solver, missing.as_ptr(), &raw const opts, &raw mut out2) };
    assert_eq!(st, ARCSEC_FILE_ERROR);
    assert!(out2.is_null());
    let junk = f.dir.path().join("junk.dat");
    std::fs::write(&junk, b"not an image").unwrap();
    let junk = c_path(&junk);
    // SAFETY: valid arguments.
    let st = unsafe { arcsec_solve_file(solver, junk.as_ptr(), &raw const opts, &raw mut out2) };
    assert_eq!(st, ARCSEC_FILE_ERROR);
    assert!(last_error().contains("not a FITS"), "{}", last_error());
    // SAFETY: freeing once.
    unsafe {
        arcsec_result_free(out);
        arcsec_solver_free(solver);
    }
}

#[test]
fn bad_arguments_are_errors_not_crashes() {
    let px = [0f32; 16];
    let img = arcsec_image {
        data: px.as_ptr().cast(),
        width: 4,
        height: 4,
        ..f32_image(&ImageBuffer::new(4, 4))
    };
    let solver = arcsec_solver_new();
    let mut out = ptr::null_mut();
    // SAFETY: every pointer is NULL or valid; that is the point.
    unsafe {
        assert_eq!(
            arcsec_solve(ptr::null(), &raw const img, ptr::null(), &raw mut out),
            ARCSEC_INVALID_ARGUMENT
        );
        assert!(last_error().contains("solver is NULL"));
        assert_eq!(
            arcsec_solve(solver, ptr::null(), ptr::null(), &raw mut out),
            ARCSEC_INVALID_ARGUMENT
        );
        assert_eq!(
            arcsec_solve(solver, &raw const img, ptr::null(), ptr::null_mut()),
            ARCSEC_INVALID_ARGUMENT
        );
        let mut opts = core::mem::zeroed::<arcsec_solve_options>();
        assert_eq!(
            arcsec_solve(solver, &raw const img, &raw const opts, &raw mut out),
            ARCSEC_INVALID_ARGUMENT,
            "an uninitialised options struct"
        );
        assert!(last_error().contains("struct_size"), "{}", last_error());
        arcsec_solve_options_init(&raw mut opts);
        // A blank image in a directory with no database: the database is the
        // first thing missing.
        let nowhere = c"/nonexistent/arcsec";
        opts.catalog_dir = nowhere.as_ptr();
        assert_eq!(
            arcsec_solve(solver, &raw const img, &raw const opts, &raw mut out),
            ARCSEC_DATABASE_NOT_FOUND
        );
        assert!(out.is_null());
        // Paths that do not exist answer "no".
        assert_eq!(arcsec_has_star_database(nowhere.as_ptr()), 0);
        assert_eq!(arcsec_has_blind_index(nowhere.as_ptr()), 0);
        // Result accessors on NULL.
        assert_eq!(
            arcsec_result_wcs(ptr::null(), ptr::null_mut()),
            ARCSEC_INVALID_ARGUMENT
        );
        assert!(arcsec_result_database(ptr::null()).is_null());
        assert_eq!(
            arcsec_result_matched_stars(ptr::null(), ptr::null_mut(), 9),
            0
        );
        arcsec_solver_cancel(ptr::null());
        arcsec_solver_free(ptr::null_mut());
        arcsec_solver_free(solver);
    }
}

#[test]
fn too_few_stars_and_no_solution_have_their_codes() {
    let f = field();
    let catalog = c_path(f.dir.path());
    let opts = options(&catalog);
    let solver = arcsec_solver_new();
    // A blank frame: nothing to detect.
    let blank = ImageBuffer {
        data: vec![1000.0; 400 * 320],
        width: 400,
        height: 320,
    };
    solve(solver, &f32_image(&blank), &opts, ARCSEC_INSUFFICIENT_STARS);
    // The right image, searched only far from where it is.
    let mut far = opts;
    far.ra_deg = 200.0;
    far.dec_deg = 40.0;
    far.search_radius_deg = 0.0;
    far.auto_index = 0;
    solve(
        solver,
        &f32_image(&f.img),
        &far,
        arcsec_status::ARCSEC_NO_SOLUTION,
    );
    assert_eq!(last_error(), "no solution found");
    // SAFETY: freeing once.
    unsafe { arcsec_solver_free(solver) };
}

unsafe extern "C" fn cancel_now(user: *mut c_void) -> c_int {
    // SAFETY: the tests pass a pointer to an AtomicUsize.
    unsafe { &*user.cast::<AtomicUsize>() }.fetch_add(1, Ordering::Relaxed);
    1
}

#[test]
fn the_cancel_callback_stops_a_solve() {
    let f = field();
    let catalog = c_path(f.dir.path());
    let mut opts = options(&catalog);
    let polls = AtomicUsize::new(0);
    opts.cancel = Some(cancel_now);
    opts.cancel_user = (&raw const polls).cast_mut().cast();
    let solver = arcsec_solver_new();
    solve(solver, &f32_image(&f.img), &opts, ARCSEC_CANCELLED);
    assert!(polls.load(Ordering::Relaxed) >= 1);
    // The solver is reusable afterwards.
    opts.cancel = None;
    let r = solve(solver, &f32_image(&f.img), &opts, ARCSEC_OK);
    // SAFETY: freeing once.
    unsafe {
        arcsec_result_free(r);
        arcsec_solver_free(solver);
    }
}

/// A raw pointer to share with another thread in a test.
#[derive(Clone, Copy)]
struct Shared(*mut arcsec_solver);
// SAFETY: arcsec_solver is internally synchronised; the tests keep it alive.
unsafe impl Send for Shared {}

/// Options (raw pointers inside) sent to another thread in a test.
struct SendOptions(arcsec_solve_options);
// SAFETY: the strings and user data they point to outlive the thread (joined).
unsafe impl Send for SendOptions {}

impl SendOptions {
    /// Taken by a method so a closure captures the whole wrapper, not the field.
    fn get(self) -> arcsec_solve_options {
        self.0
    }
}

unsafe extern "C" fn slow_poll(user: *mut c_void) -> c_int {
    // Signal that the solve is under way, then let it continue.
    // SAFETY: the test passes a pointer to an AtomicBool.
    unsafe { &*user.cast::<AtomicBool>() }.store(true, Ordering::Relaxed);
    std::thread::sleep(core::time::Duration::from_millis(2));
    0
}

#[test]
fn another_thread_can_cancel_and_a_busy_solver_says_so() {
    let f = field();
    let catalog = c_path(f.dir.path());
    let mut opts = options(&catalog);
    // Somewhere the field is not, over a wide area: a long search.
    opts.ra_deg = 200.0;
    opts.dec_deg = 40.0;
    opts.search_radius_deg = 60.0;
    opts.auto_index = 0;
    opts.threads = 1;
    let started = Arc::new(AtomicBool::new(false));
    opts.cancel = Some(slow_poll);
    opts.cancel_user = Arc::as_ptr(&started).cast_mut().cast();
    let solver = Shared(arcsec_solver_new());
    let img = f.img.clone();
    let sent = SendOptions(opts);
    let worker = std::thread::spawn(move || {
        let s = solver;
        let opts = sent.get();
        let mut out = ptr::null_mut();
        let image = f32_image(&img);
        // SAFETY: valid arguments; the solver outlives the thread.
        let st = unsafe { arcsec_solve(s.0, &raw const image, &raw const opts, &raw mut out) };
        (st, out.is_null())
    });
    while !started.load(Ordering::Relaxed) {
        std::thread::yield_now();
    }
    // A second solve on the same handle is refused, not queued.
    solve(solver.0, &f32_image(&f.img), &opts, ARCSEC_BUSY);
    // SAFETY: a live solver.
    unsafe { arcsec_solver_cancel(solver.0) };
    let (st, null) = worker.join().unwrap();
    assert_eq!(st, ARCSEC_CANCELLED);
    assert!(null);
    // SAFETY: the solve has returned.
    unsafe { arcsec_solver_free(solver.0) };
}

/// Progress reports seen: (stage, fraction) pairs.
type Seen = std::sync::Mutex<Vec<(String, f64)>>;

unsafe extern "C" fn record_progress(user: *mut c_void, fraction: f64, stage: *const c_char) {
    // SAFETY: the test passes a pointer to a Seen; arcsec passes a C string.
    let (seen, stage) = unsafe { (&*user.cast::<Seen>(), CStr::from_ptr(stage)) };
    if let Ok(mut v) = seen.lock() {
        v.push((stage.to_string_lossy().into_owned(), fraction));
    }
}

#[test]
fn progress_is_reported_by_stage() {
    let f = field();
    let catalog = c_path(f.dir.path());
    let mut opts = options(&catalog);
    // Start far enough off that the spiral has to move.
    opts.ra_deg = 84.3 + 1.2;
    opts.dec_deg = -5.2 - 0.8;
    let seen: Seen = std::sync::Mutex::new(Vec::new());
    opts.progress = Some(record_progress);
    opts.progress_user = (&raw const seen).cast_mut().cast();
    let solver = arcsec_solver_new();
    let r = solve(solver, &f32_image(&f.img), &opts, ARCSEC_OK);
    let seen = seen.into_inner().unwrap();
    assert_eq!(seen.first().map(|s| s.0.as_str()), Some("detecting stars"));
    let search: Vec<f64> = seen
        .iter()
        .filter(|s| s.0 == "searching")
        .map(|s| s.1)
        .collect();
    assert!(!search.is_empty(), "{seen:?}");
    assert!(search.iter().all(|&f| f > 0.0 && f < 1.0), "{search:?}");
    // SAFETY: freeing once.
    unsafe {
        arcsec_result_free(r);
        arcsec_solver_free(solver);
    }
}

#[test]
fn analysis_counts_and_measures_stars() {
    let f = field();
    let mut a = arcsec_analysis {
        struct_size: core::mem::size_of::<arcsec_analysis>(),
        ..crate::util::Versioned::defaults()
    };
    let mut stars = vec![arcsec_star::default(); 5];
    let img = f32_image(&f.img);
    // SAFETY: valid image, sized output, room for 5 stars.
    let st = unsafe {
        arcsec_analyse(
            &raw const img,
            0.0,
            0,
            &raw mut a,
            stars.as_mut_ptr(),
            stars.len(),
        )
    };
    assert_eq!(st, ARCSEC_OK, "{}", last_error());
    assert!(a.star_count > 20, "{}", a.star_count);
    assert!(
        a.hfd_median > 1.0 && a.hfd_median < 10.0,
        "{}",
        a.hfd_median
    );
    assert!(a.noise > 0.0 && a.background > 900.0);
    assert!(
        stars
            .iter()
            .all(|s| s.hfd > 0.0 && s.x >= 1.0 && s.snr > 30.0)
    );
    // A missing output is refused.
    // SAFETY: as above, with a NULL output.
    let st = unsafe { arcsec_analyse(&raw const img, 0.0, 0, ptr::null_mut(), ptr::null_mut(), 0) };
    assert_eq!(st, ARCSEC_INVALID_ARGUMENT);
}

#[test]
fn catalogue_queries() {
    let f = field();
    let dir = c_path(f.dir.path());
    let mut buf = [0 as c_char; 8];
    // SAFETY: valid strings and buffers.
    unsafe {
        assert_eq!(arcsec_has_star_database(dir.as_ptr()), 1);
        assert_eq!(arcsec_has_blind_index(dir.as_ptr()), 0);
        let n = arcsec_select_database(dir.as_ptr(), 0.5, buf.as_mut_ptr(), buf.len());
        assert_eq!(n, 3);
        assert_eq!(CStr::from_ptr(buf.as_ptr()).to_str().unwrap(), "d50");
        // Too small a buffer: cut, terminated, and the full length returned.
        let n = arcsec_select_database(dir.as_ptr(), 0.5, buf.as_mut_ptr(), 2);
        assert_eq!(n, 3);
        assert_eq!(CStr::from_ptr(buf.as_ptr()).to_str().unwrap(), "d");
        assert_eq!(
            arcsec_select_database(dir.as_ptr(), -1.0, buf.as_mut_ptr(), buf.len()),
            0
        );
        let need = arcsec_default_catalog_dir(ptr::null_mut(), 0);
        assert!(need > 0);
        let mut path = vec![0 as c_char; need + 1];
        assert_eq!(
            arcsec_default_catalog_dir(path.as_mut_ptr(), path.len()),
            need
        );
        let p = CStr::from_ptr(path.as_ptr()).to_string_lossy().into_owned();
        assert_eq!(
            std::path::Path::new(&p),
            arcsec_core::auto::default_catalog_dir()
        );
        assert!(arcsec_default_database_dir(ptr::null_mut(), 0) > 0);
    }
}

static LOG_LINES: AtomicUsize = AtomicUsize::new(0);

/// Counts well-formed lines. (No assertions here: a panic cannot leave an
/// `extern "C"` function, it would abort the test run.)
unsafe extern "C" fn count_lines(_user: *mut c_void, level: c_int, msg: *const c_char) {
    // SAFETY: arcsec passes a NUL-terminated string.
    let text = unsafe { CStr::from_ptr(msg) }.to_bytes();
    if (1..=5).contains(&level) && text == b"a line from arcsec" {
        LOG_LINES.fetch_add(1, Ordering::Relaxed);
    }
}

#[test]
fn the_log_callback_hears_progress_and_panics_become_errors() {
    // SAFETY: a thread-safe callback.
    unsafe { arcsec_set_log_callback(Some(count_lines), ptr::null_mut(), 3) };
    log::info!("a line from arcsec");
    assert!(LOG_LINES.load(Ordering::Relaxed) >= 1);

    // A panic inside an entry point is caught, reported and turned into a code.
    let st = crate::error::guard(|| panic!("deliberate test panic"));
    assert_eq!(st, ARCSEC_INTERNAL_ERROR);
    assert!(
        last_error().contains("deliberate test panic"),
        "{}",
        last_error()
    );
    // SAFETY: turning the callback off.
    unsafe { arcsec_set_log_callback(None, ptr::null_mut(), 0) };
}

#[test]
fn header_is_current() {
    let dir = env!("CARGO_MANIFEST_DIR");
    let config = cbindgen::Config::from_file(format!("{dir}/cbindgen.toml")).unwrap();
    let bindings = cbindgen::Builder::new()
        .with_crate(dir)
        .with_config(config)
        .generate()
        .expect("cbindgen could not generate the header");
    let mut generated = Vec::new();
    bindings.write(&mut generated);
    let generated = String::from_utf8(generated).unwrap();
    let path = std::path::Path::new(dir).join("include").join("arcsec.h");
    if std::env::var_os("ARCSEC_BLESS").is_some() {
        std::fs::write(&path, &generated).unwrap();
        return;
    }
    let committed = std::fs::read_to_string(&path)
        .unwrap_or_default()
        .replace("\r\n", "\n");
    assert!(
        committed == generated,
        "include/arcsec.h is out of date: regenerate it with \
         `ARCSEC_BLESS=1 cargo test -p libarcsec header_is_current` and commit it"
    );
}
