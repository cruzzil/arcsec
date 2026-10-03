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

// The logo. Change this one line to switch: 'asterism' (four stars joined, the star
// pattern arcsec matches), 'angle' (an angle and its arc), or 'wordmark' (the name
// alone, followed by the arcsecond sign ″). Files: src/assets/logos/, public/favicons/.
const LOGO = 'asterism';

/** @type {Record<string, { logo?: { light: string, dark: string, alt: string }, css: string[] }>} */
const logos = {
	asterism: {
		logo: { light: './src/assets/logos/asterism-light.svg', dark: './src/assets/logos/asterism-dark.svg', alt: '' },
		css: [],
	},
	angle: {
		logo: { light: './src/assets/logos/angle-light.svg', dark: './src/assets/logos/angle-dark.svg', alt: '' },
		css: [],
	},
	wordmark: { css: ['./src/styles/wordmark.css'] },
};

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
				'Free plate solving for astrophotographers. arcsec works out where your telescope is pointing from the stars in your image, and works in N.I.N.A. in place of ASTAP.',
			logo: logos[LOGO].logo,
			favicon: `/favicons/${LOGO}.svg`,
			// src/pages/404.astro replaces Starlight's 404 page, with links that respect `base`.
			disable404Route: true,
			social: [{ icon: 'github', label: 'arcsec on GitHub', href: repo }],
			editLink: {
				baseUrl: `${repo}/edit/main/site/`,
			},
			customCss: ['./src/styles/custom.css', ...logos[LOGO].css],
			components: {
				Hero: './src/components/Hero.astro',
				// Adds text links to the main pages beside the GitHub icon in the header.
				SocialIcons: './src/components/HeaderLinks.astro',
			},
			head: [
				{ tag: 'meta', attrs: { name: 'theme-color', content: '#0b1422' } },
			],
			sidebar: [
				{ label: 'Download and install', slug: 'getting-started' },
				{ label: 'Use with N.I.N.A.', slug: 'nina' },
				{ label: 'Use with Siril', slug: 'siril' },
				{ label: 'Which catalogue do I need?', slug: 'which-catalogue' },
				{ label: 'Catalogue guide', slug: 'catalogues' },
				{ label: 'FAQ', slug: 'faq' },
				{
					label: 'Reference',
					items: [
						{ label: 'Command line', slug: 'reference/cli' },
						{ label: 'For developers', slug: 'developers' },
						{ label: 'C library', slug: 'c-library' },
					],
				},
				{
					label: 'Project',
					items: [
						{ label: 'GitHub', link: repo },
						{ label: 'Releases', link: `${repo}/releases` },
						{ label: 'Changelog', link: `${repo}/blob/main/CHANGELOG.md` },
						{ label: 'Report a problem', link: `${repo}/issues` },
					],
				},
			],
		}),
	],
});
