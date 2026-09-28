// The catalogues `arcsec catalog` can install, and its `recommend` rule.
//
// This mirrors REGISTRY in arcsec/src/catalog_cmd/registry.rs and cmd_recommend in
// arcsec/src/catalog_cmd/mod.rs. Keep the two in step: if a catalogue, size or field
// range changes there, change it here too, and re-check a few fields against
// `arcsec catalog recommend --fov <deg>`.

export type Purpose = 'solving' | 'photometry' | 'blind';

export interface Catalogue {
	id: string;
	purpose: Purpose;
	/** Approximate download size in bytes. */
	bytes: number;
	/** Field-of-view range the catalogue is built for, in degrees (inclusive). */
	fov: [number, number];
	desc: string;
}

export const CATALOGUES: Catalogue[] = [
	{ id: 'd05', purpose: 'solving', bytes: 102_200_000, fov: [0.6, 6.0], desc: 'Gaia DR3 to 500 stars/deg². Smallest useful solving database.' },
	{ id: 'd20', purpose: 'solving', bytes: 399_600_000, fov: [0.3, 6.0], desc: 'Gaia DR3 to 2000 stars/deg².' },
	{ id: 'd50', purpose: 'solving', bytes: 901_300_000, fov: [0.2, 6.0], desc: 'Gaia DR3 to 5000 stars/deg². The usual choice for solving.' },
	{ id: 'd80', purpose: 'solving', bytes: 1_213_400_000, fov: [0.15, 6.0], desc: 'Gaia DR3 to 8000 stars/deg². Densest; needed below ~0.2° fields.' },
	{ id: 'g05', purpose: 'solving', bytes: 101_600_000, fov: [3.0, 20.0], desc: 'Wide fields, 3°–20°. The D-series stops at 6°.' },
	{ id: 'w08', purpose: 'solving', bytes: 330_000, fov: [20.0, 80.0], desc: 'Very wide fields, 20°–80°, to magnitude 8. Tiny.' },
	{ id: 'v05', purpose: 'photometry', bytes: 116_900_000, fov: [0.6, 6.0], desc: 'Johnson-V magnitudes plus Gaia BP-RP colour, 500 stars/deg².' },
	{ id: 'v50', purpose: 'photometry', bytes: 1_011_000_000, fov: [0.2, 6.0], desc: 'Johnson-V plus BP-RP colour, 5000 stars/deg². Deeper photometry.' },
	{ id: 'anet-4100', purpose: 'blind', bytes: 355_500_000, fov: [0.7, 180.0], desc: 'Tycho-2 blind indexes, scales 07–19 (fields ~0.7° and wider).' },
	{ id: 'anet-5200', purpose: 'blind', bytes: 8_800_000_000, fov: [0.1, 2.0], desc: 'Gaia LITE blind indexes 5200/5201/5202, 48 HEALPix each. Large.' },
];

/**
 * The smallest download of `purpose` whose field range covers `fovDeg`, bounds
 * inclusive - exactly `pick` in cmd_recommend. Undefined if none does.
 */
export function pick(purpose: Purpose, fovDeg: number): Catalogue | undefined {
	let best: Catalogue | undefined;
	for (const c of CATALOGUES) {
		if (c.purpose !== purpose || fovDeg < c.fov[0] || fovDeg > c.fov[1]) continue;
		// min_by_key keeps the first of equal keys; so does a strict "<".
		if (!best || c.bytes < best.bytes) best = c;
	}
	return best;
}

export interface Recommendation {
	solving?: Catalogue;
	photometry?: Catalogue;
	blind?: Catalogue;
	/** The `arcsec catalog install ...` line: solving, plus photometry if asked for. */
	install?: string;
}

export function recommend(fovDeg: number, wantPhotometry: boolean): Recommendation {
	const solving = pick('solving', fovDeg);
	const photometry = wantPhotometry ? pick('photometry', fovDeg) : undefined;
	const blind = pick('blind', fovDeg);
	const ids = [solving, photometry].filter((c): c is Catalogue => !!c).map((c) => c.id);
	return {
		solving,
		photometry,
		blind,
		install: ids.length ? `arcsec catalog install ${ids.join(' ')}` : undefined,
	};
}

/** Image scale in arcseconds per (binned) pixel - image_io::pixel_scale_from. */
export function pixelScale(focalMm: number, pixelUm: number, binning = 1): number {
	return ((pixelUm * binning) / focalMm) * 206.265;
}

/**
 * Field of view in degrees along each axis. Binning changes the pixel scale but not
 * the field, so the sensor size is given unbinned.
 */
export function fieldOfView(focalMm: number, pixelUm: number, widthPx: number, heightPx: number) {
	const scale = pixelScale(focalMm, pixelUm, 1);
	return { width: (scale * widthPx) / 3600, height: (scale * heightPx) / 3600 };
}

/** Download size in the decimal units `arcsec catalog list` prints ("102.2 MB"). */
export function human(bytes: number): string {
	const units = ['B', 'kB', 'MB', 'GB', 'TB'];
	let v = bytes;
	let i = 0;
	while (v >= 1000 && i < units.length - 1) {
		v /= 1000;
		i++;
	}
	return i === 0 ? `${bytes} B` : `${v.toFixed(1)} ${units[i]}`;
}

/** Field range as `arcsec catalog list` prints it ("0.15°–6°"). */
export function fovRange(c: Catalogue): string {
	return `${c.fov[0]}°–${c.fov[1]}°`;
}
