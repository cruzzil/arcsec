// The latest arcsec release, fetched from the GitHub API once per build.
//
// Asset names carry the version (arcsec-v0.1.2-x86_64-unknown-linux-gnu.tar.gz), so
// the download links cannot be written by hand. The site is rebuilt whenever a
// release is published (see .github/workflows/site.yml), which keeps them current.
//
// A failed fetch never fails the build: the page falls back to linking the releases
// page. Set GITHUB_TOKEN to avoid the API's anonymous rate limit (CI does), or
// ARCSEC_SITE_OFFLINE=1 to skip the request entirely.

export const REPO = 'cruzzil/arcsec';
export const REPO_URL = `https://github.com/${REPO}`;
export const RELEASES_URL = `${REPO_URL}/releases`;
export const LATEST_URL = `${RELEASES_URL}/latest`;

export type Os = 'windows' | 'macos' | 'linux';
export type Arch = 'x64' | 'arm64';

export interface Platform {
	/** Rust target triple, as it appears in the asset name. */
	target: string;
	os: Os;
	arch: Arch;
	/** Human name, e.g. "Linux (x86-64)". */
	label: string;
	/** Short qualifier shown under the button. */
	note: string;
}

/** The platforms the Release workflow builds, in display order. */
export const PLATFORMS: Platform[] = [
	{ target: 'x86_64-pc-windows-msvc', os: 'windows', arch: 'x64', label: 'Windows (x86-64)', note: '.zip' },
	{ target: 'aarch64-apple-darwin', os: 'macos', arch: 'arm64', label: 'macOS (Apple silicon)', note: '.tar.gz' },
	{ target: 'x86_64-unknown-linux-gnu', os: 'linux', arch: 'x64', label: 'Linux (x86-64)', note: '.tar.gz, glibc 2.35+' },
	{ target: 'aarch64-unknown-linux-gnu', os: 'linux', arch: 'arm64', label: 'Linux (arm64)', note: '.tar.gz, glibc 2.35+' },
];

export interface Asset extends Platform {
	name: string;
	url: string;
	size: number;
	/** The asset's own .sha256 file, if published. */
	sha256Url?: string;
}

export interface Release {
	ok: boolean;
	/** "0.1.2", without the leading v. */
	version?: string;
	tag?: string;
	date?: string;
	pageUrl: string;
	checksumsUrl?: string;
	assets: Asset[];
}

interface GhAsset {
	name: string;
	size: number;
	browser_download_url: string;
}

let cached: Promise<Release> | undefined;

export function getLatestRelease(): Promise<Release> {
	cached ??= fetchLatest();
	return cached;
}

async function fetchLatest(): Promise<Release> {
	const fallback: Release = { ok: false, pageUrl: LATEST_URL, assets: [] };
	if (process.env.ARCSEC_SITE_OFFLINE) return fallback;

	const headers: Record<string, string> = {
		Accept: 'application/vnd.github+json',
		'X-GitHub-Api-Version': '2022-11-28',
		'User-Agent': 'arcsec-site-build',
	};
	const token = process.env.GITHUB_TOKEN;
	if (token) headers.Authorization = `Bearer ${token}`;

	try {
		const res = await fetch(`https://api.github.com/repos/${REPO}/releases/latest`, {
			headers,
			signal: AbortSignal.timeout(15_000),
		});
		if (!res.ok) throw new Error(`HTTP ${res.status} ${res.statusText}`);
		const json = (await res.json()) as {
			tag_name: string;
			html_url: string;
			published_at: string;
			assets: GhAsset[];
		};
		const byName = new Map(json.assets.map((a) => [a.name, a]));
		const assets: Asset[] = [];
		for (const p of PLATFORMS) {
			const a = json.assets.find(
				(a) => a.name.includes(`-${p.target}.`) && /\.(tar\.gz|zip)$/.test(a.name),
			);
			if (!a) continue;
			assets.push({
				...p,
				name: a.name,
				url: a.browser_download_url,
				size: a.size,
				sha256Url: byName.get(`${a.name}.sha256`)?.browser_download_url,
			});
		}
		const release: Release = {
			ok: assets.length > 0,
			tag: json.tag_name,
			version: json.tag_name.replace(/^v/, ''),
			date: json.published_at?.slice(0, 10),
			pageUrl: json.html_url || LATEST_URL,
			checksumsUrl: byName.get('SHA256SUMS')?.browser_download_url,
			assets,
		};
		console.log(`[release] ${release.tag}: ${assets.length} platform archives`);
		return release;
	} catch (e) {
		console.warn(`[release] could not fetch the latest release (${String(e)}); linking the releases page instead`);
		return fallback;
	}
}
