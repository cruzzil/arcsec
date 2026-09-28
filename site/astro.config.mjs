// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// Hosting: GitHub Pages project site at https://cruzzil.github.io/arcsec/.
//
// To move to a custom domain (say https://arcsec.example.org/):
//   1. set `site` to that origin and `base` to '/' (or delete `base`);
//   2. add public/CNAME containing the bare domain name;
//   3. set the same domain under Settings > Pages in the repository.
// Nothing else needs to change: the sidebar is built from slugs, components link
// through src/lib/url.ts, and links in the Markdown pages are relative ("../nina/").
const site = 'https://cruzzil.github.io';
const base = '/arcsec';

const repo = 'https://github.com/cruzzil/arcsec';

export default defineConfig({
	site,
	base,
	trailingSlash: 'always',
	vite: {
		build: {
			rolldownOptions: {
				// Astro marks MDX pages with a "use astro:head-inject" directive that
				// Rolldown warns it may not preserve. It is Astro's own marker and harmless;
				// drop just that warning so real ones stand out.
				onLog(level, log, handler) {
					if (log.code === 'MODULE_LEVEL_DIRECTIVE' && log.message.includes('astro:head-inject')) return;
					handler(level, log);
				},
			},
		},
	},
	integrations: [
		starlight({
			title: 'arcsec',
			description:
				'An accurate astrometric plate solver and a drop-in replacement for ASTAP: same command line, same star databases, same output files. Works as the "ASTAP" solver in N.I.N.A.',
			logo: {
				src: './src/assets/logo.svg',
				alt: '',
			},
			favicon: '/favicon.svg',
			// src/pages/404.astro replaces Starlight's 404 page, with links that respect `base`.
			disable404Route: true,
			social: [{ icon: 'github', label: 'arcsec on GitHub', href: repo }],
			editLink: {
				baseUrl: `${repo}/edit/main/site/`,
			},
			customCss: ['./src/styles/custom.css'],
			components: {
				Hero: './src/components/Hero.astro',
			},
			head: [
				{ tag: 'meta', attrs: { name: 'theme-color', content: '#0b1422' } },
			],
			sidebar: [
				{
					label: 'Start here',
					items: [
						{ label: 'Getting started', slug: 'getting-started' },
						{ label: 'Use with N.I.N.A.', slug: 'nina' },
					],
				},
				{
					label: 'Catalogues',
					items: [
						{ label: 'Catalogue picker', slug: 'catalogues/picker' },
						{ label: 'Catalogue guide', slug: 'catalogues' },
					],
				},
				{
					label: 'Reference',
					items: [
						{ label: 'Command line', slug: 'reference/cli' },
						{ label: 'FAQ', slug: 'faq' },
					],
				},
				{
					label: 'Project',
					items: [
						{ label: 'GitHub', link: repo, attrs: { rel: 'noopener' } },
						{ label: 'Releases', link: `${repo}/releases` },
						{ label: 'Changelog', link: `${repo}/blob/main/CHANGELOG.md` },
						{ label: 'Benchmark results', link: `${repo}/blob/main/docs/test-images.md#6-results--103-images-arcsec-vs-astap` },
						{ label: 'How it works', link: `${repo}/blob/main/docs/plate-solving.md` },
						{ label: 'Contributing', link: `${repo}/blob/main/CONTRIBUTING.md` },
						{ label: 'crates.io', link: 'https://crates.io/crates/arcsec' },
						{ label: 'Library docs (docs.rs)', link: 'https://docs.rs/arcsec-core' },
					],
				},
			],
		}),
	],
});
