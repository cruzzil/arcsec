#!/usr/bin/env node
// Check every internal link in the built site (dist/) resolves, including #fragments,
// and that none escapes the configured base path. Run after `npm run build`:
//
//   npm run check-links
//
// External links are not fetched. Exits 1 if anything is broken.

import { existsSync, readFileSync, readdirSync, statSync } from 'node:fs';
import { join, relative, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const dist = fileURLToPath(new URL('../dist/', import.meta.url));
const config = readFileSync(fileURLToPath(new URL('../astro.config.mjs', import.meta.url)), 'utf8');
const base = (config.match(/const base = '([^']*)'/)?.[1] ?? '/').replace(/\/$/, '');
const origin = 'https://site.invalid';

function* htmlFiles(dir) {
	for (const name of readdirSync(dir)) {
		const p = join(dir, name);
		if (statSync(p).isDirectory()) yield* htmlFiles(p);
		else if (name.endsWith('.html')) yield p;
	}
}

/** The URL path a dist file is served at. */
function pagePath(file) {
	const rel = relative(dist, file).split(sep).join('/');
	return `${base}/${rel.replace(/(^|\/)index\.html$/, '$1')}`;
}

/** The dist file a URL path is served from, or undefined. */
function fileFor(path) {
	if (path !== base && !path.startsWith(`${base}/`)) return undefined;
	const rel = decodeURIComponent(path.slice(base.length)).replace(/^\//, '');
	for (const candidate of [rel, join(rel, 'index.html')]) {
		const p = join(dist, candidate);
		if (existsSync(p) && statSync(p).isFile()) return p;
	}
	return undefined;
}

const ids = new Map();
function idsOf(file) {
	if (!ids.has(file)) {
		const html = readFileSync(file, 'utf8');
		ids.set(file, new Set([...html.matchAll(/\sid="([^"]+)"/g)].map((m) => m[1])));
	}
	return ids.get(file);
}

let checked = 0;
const broken = [];
for (const file of htmlFiles(dist)) {
	const html = readFileSync(file, 'utf8');
	const page = pagePath(file);
	const refs = [
		...[...html.matchAll(/\s(?:href|src)="([^"]*)"/g)].map((m) => m[1]),
		...[...html.matchAll(/\ssrcset="([^"]*)"/g)].flatMap((m) => m[1].split(',').map((s) => s.trim().split(/\s+/)[0])),
	];
	for (const raw of refs) {
		const ref = raw.replaceAll('&amp;', '&');
		if (!ref || /^(https?:|mailto:|data:|javascript:)/.test(ref) || ref.startsWith('//')) continue;
		const url = new URL(ref, origin + page);
		if (url.origin !== origin) continue;
		checked++;
		const target = fileFor(url.pathname);
		if (!target) {
			broken.push(`${page}: ${ref} (no such page or file)`);
			continue;
		}
		const frag = decodeURIComponent(url.hash.slice(1));
		if (frag && target.endsWith('.html') && !idsOf(target).has(frag)) {
			broken.push(`${page}: ${ref} (no element with id "${frag}")`);
		}
	}
}

if (broken.length) {
	console.error(`${broken.length} broken internal link(s):\n  ${broken.join('\n  ')}`);
	process.exit(1);
}
console.log(`All ${checked} internal links resolve under ${base || '/'}.`);
