// Typesets $…$ and $$…$$ with KaTeX at build time, for the "How it works" pages.
// Sätteri parses the math (its `math` feature); this turns each node into KaTeX's
// HTML, so the browser needs only the stylesheet (katex.min.css, in customCss).
import katex from 'katex';

/** @param {string} tex @param {boolean} displayMode */
function render(tex, displayMode) {
	// throwOnError: a typo in a formula fails the build instead of shipping red text.
	const html = katex.renderToString(tex, { displayMode, throwOnError: true, output: 'html' });
	// An mdast `html` node is spliced in where the math was. Raw content would be
	// re-parsed as Markdown, which makes inline math a block of its own.
	return { type: /** @type {const} */ ('html'), value: html };
}

export const katexPlugin = {
	name: 'arcsec-katex',
	/** @param {{ value: string }} node */
	math: (node) => render(node.value, true),
	/** @param {{ value: string }} node */
	inlineMath: (node) => render(node.value, false),
};
