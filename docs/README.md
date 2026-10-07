# Animated site

`index.html` loads `data.js` (the numerical series) and `app.js` (the animations).
Everything is static: no build step, no dependencies, no network calls.

The data is regenerated from the package by

```
python scripts/build_site_data.py          # all seven series
python scripts/build_site_data.py --only pi,quasi
```

`data.js` is a JavaScript assignment rather than a JSON file so that the page also
works when opened directly from disk: a `fetch` of a local JSON file is blocked by
the browser origin rules, while a `<script>` is not.

To preview locally:

```
python -m http.server 8000 --directory docs
```
