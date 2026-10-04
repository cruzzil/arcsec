# arcsec-catalogue

Catalogue management for the [arcsec](https://github.com/cruzzil/arcsec) plate solver:
the star databases (ASTAP's D05 to W08, and the V05/V50 photometric ones) and blind
indexes it solves with.

- **What exists and what is installed:** the registry of installable catalogues, how to
  recognise each one's files on disk, and where they live (`$ARCSEC_CATALOG_DIR`, else
  the platform's data directory).
- **Which database for a field:** the choice a solve makes when it is not told,
  re-exported from `arcsec-core`.
- **Installing:** resumable HTTPS downloads and zip and `.deb` extraction that cannot
  write outside the catalogue directory, then removing and verifying.
- **The blind index:** which index to build for the databases installed, what it will
  cost in disk, time and memory, whether the one installed is stale, and building it.

Nothing prints or exits: results and errors are typed, and long operations report
through a callback and honour `arcsec_core::cancel`.

```rust
use arcsec_catalogue::{default_dir, registry, select_db_for_fov};

let dir = default_dir();
for e in registry::REGISTRY {
    let state = if registry::is_installed(&dir, e) { "installed" } else { "-" };
    println!("{:<10} {:<10} {state}", e.id, e.purpose.label());
}
println!("for a 1.5° field: {:?}", select_db_for_fov(&dir, 1.5));
```

Downloading is behind the `download` feature, on by default. A program that only needs
to know where catalogues are, what is installed and which to use can turn it off
(`default-features = false`) and leaves out the HTTP, TLS and archive code.

It is shared by the `arcsec` command line (`arcsec catalog ...`) and is there for other
front-ends. Most people want the command-line tool (`cargo install arcsec`). Its API
follows those users and may change between minor versions.

Licensed under the MIT licence.
