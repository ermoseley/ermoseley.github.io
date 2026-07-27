/* =========================================================================
   blog.js — renders the post list from blog/posts.json.

   Adding a post: write blog/posts/<slug>.html (copy an existing one), then
   add an entry to the top of the "posts" array in blog/posts.json.
   ========================================================================= */

(function (global) {
  'use strict';

  function fmt(iso) {
    const d = new Date(iso + 'T12:00:00');
    if (isNaN(d)) return iso;
    return d.toLocaleDateString('en-GB', { day: '2-digit', month: 'short', year: 'numeric' });
  }

  function esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  function card(p, base) {
    const href = base + 'posts/' + p.file;
    return '<a class="post" href="' + href + '">' +
      '<span class="when">' + esc(fmt(p.date)) + (p.draft ? ' · draft' : '') + '</span>' +
      '<h3>' + esc(p.title) + '</h3>' +
      '<p>' + esc(p.summary) + '</p>' +
      (p.tags && p.tags.length ? '<div class="tags">' + p.tags.map(esc).join(' · ') + '</div>' : '') +
      '</a>';
  }

  function empty(msg) {
    return '<div class="post"><span class="when">nothing yet</span>' +
      '<h3>' + esc(msg) + '</h3>' +
      '<p>Add an entry to <span class="mono">blog/posts.json</span> and a matching file in ' +
      '<span class="mono">blog/posts/</span>.</p></div>';
  }

  function render(el, base, limit) {
    if (!el) return;
    base = base || '';
    fetch(base + 'posts.json', { cache: 'no-cache' })
      .then(function (r) {
        if (!r.ok) throw new Error(r.status);
        return r.json();
      })
      .then(function (data) {
        let posts = (data && data.posts) || [];
        posts = posts.slice().sort(function (a, b) { return a.date < b.date ? 1 : -1; });
        if (limit) posts = posts.slice(0, limit);
        el.innerHTML = posts.length
          ? posts.map(function (p) { return card(p, base); }).join('')
          : empty('No posts yet.');
      })
      .catch(function () {
        // file:// blocks fetch; say so rather than showing a broken spinner
        el.innerHTML = empty('Post list unavailable — serve over HTTP.');
      });
  }

  global.Blog = { render: render };
})(window);
