// Loads the web page the way a browser does: index.html plus the local
// stylesheets and classic scripts it references, in document order.
// POLYTOOL_HTML=/path/to/index.html checks a staged copy (its sibling files are read).
const fs = require('node:fs');
const path = require('node:path');

const isGenerated = url => Boolean(url) && /(^|\/)pkg\//.test(url);
const isLocal = url => !/^(?:[a-z][a-z0-9+.-]*:)?\/\//i.test(url) && !/^(?:data|blob):/i.test(url);

function loadPage(htmlPath = process.env.POLYTOOL_HTML || path.join(__dirname, '..', 'index.html')) {
  const dir = path.dirname(htmlPath);
  const markup = fs.readFileSync(htmlPath, 'utf8');
  const stylesheets = [];
  const scripts = [];
  const tag = /<link\b([^>]*)>|<script\b([^>]*)>([\s\S]*?)<\/script>/gi;
  for (const match of markup.matchAll(tag)) {
    if (match[1] !== undefined) {
      if (!/rel=["']stylesheet["']/i.test(match[1])) continue;
      const href = /href=["']([^"']+)["']/i.exec(match[1]);
      if (href && isLocal(href[1])) stylesheets.push({url: href[1], file: path.join(dir, href[1].split('?')[0])});
    } else {
      const src = /\bsrc=["']([^"']+)["']/i.exec(match[2]);
      if (src && isLocal(src[1])) scripts.push({url: src[1], file: path.join(dir, src[1].split('?')[0])});
      else if (!src) scripts.push({url: null, file: null, text: match[3]});
    }
  }
  for (const item of stylesheets) item.text = fs.readFileSync(item.file, 'utf8');
  // pkg/ holds generated wasm glue that may not be built; tests load it separately.
  for (const item of scripts) {
    if (item.file) item.text = isGenerated(item.url) ? '' : fs.readFileSync(item.file, 'utf8');
  }
  return {
    dir, markup, stylesheets, scripts,
    // Application source in load order, excluding the generated wasm glue (pkg/).
    js: scripts.filter(item => !isGenerated(item.url)).map(item => item.text).join(''),
    css: stylesheets.map(item => item.text).join(''),
  };
}

module.exports = {loadPage};
