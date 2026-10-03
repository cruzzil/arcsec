//! Sizes and durations as people say them, as the `arcsec catalog` command prints
//! them. Used by this crate's error messages and [`crate::index::Estimate::summary`].

/// Bytes, in decimal units, to one decimal place: `901.3 MB`.
#[must_use]
pub fn human_bytes(bytes: u64) -> String {
    const U: [&str; 5] = ["B", "kB", "MB", "GB", "TB"];
    let mut v = bytes as f64;
    let mut i = 0;
    while v >= 1000.0 && i < U.len() - 1 {
        v /= 1000.0;
        i += 1;
    }
    if i == 0 {
        format!("{bytes} {}", U[i])
    } else {
        format!("{v:.1} {}", U[i])
    }
}

/// A duration rounded the way a person would say it: `under 10 s`, `~20 s`,
/// `~13 min`, `~1 h 30 min`.
#[must_use]
pub fn duration(secs: f64) -> String {
    if secs < 10.0 {
        "under 10 s".to_string()
    } else if secs < 55.0 {
        format!("~{} s", ((secs / 5.0).round() * 5.0) as u64)
    } else if secs < 3600.0 {
        format!("~{} min", ((secs / 60.0).round() as u64).max(1))
    } else {
        let m = (secs / 60.0).round() as u64;
        format!("~{} h {} min", m / 60, m % 60)
    }
}

/// `0.6 GB`, or `90.0 MB` below 0.1 GB.
#[must_use]
pub fn gigabytes(bytes: u64) -> String {
    if bytes < 100_000_000 {
        human_bytes(bytes)
    } else {
        format!("{:.1} GB", bytes as f64 / 1e9)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn human_sizes_read_sensibly() {
        assert_eq!(human_bytes(500), "500 B");
        assert_eq!(human_bytes(1_500_000), "1.5 MB");
        assert_eq!(human_bytes(901_300_000), "901.3 MB");
        assert_eq!(human_bytes(1_213_400_000), "1.2 GB");
    }

    #[test]
    fn durations_read_naturally() {
        assert_eq!(duration(3.0), "under 10 s");
        assert_eq!(duration(22.0), "~20 s");
        assert_eq!(duration(107.0), "~2 min");
        assert_eq!(duration(782.0), "~13 min");
        assert_eq!(duration(5400.0), "~1 h 30 min");
        assert_eq!(gigabytes(1_276_000_000), "1.3 GB");
        assert_eq!(gigabytes(90_000_000), "90.0 MB");
    }
}
