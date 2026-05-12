# ermoseley.github.io

Personal academic site built with [Quarto](https://quarto.org/).

## Local preview

Install [Quarto](https://quarto.org/docs/get-started/), then from this directory:

```bash
quarto preview
```

Build static output to `_site/`:

```bash
quarto render
```

## Deploy

GitHub Actions ([`.github/workflows/deploy.yml`](.github/workflows/deploy.yml)) runs `quarto render` and publishes `_site/` to GitHub Pages. You can also trigger a deploy from the Actions tab (**Run workflow**).

## Content

| File / directory | Role |
|------------------|------|
| [`_quarto.yml`](_quarto.yml) | Site + navbar + theme |
| [`index.qmd`](index.qmd) | Home |
| [`publications.qmd`](publications.qmd) | Full bibliography from [`references.bib`](references.bib) |
| [`projects/`](projects/) | Project write-ups; [`projects.qmd`](projects.qmd) is the grid listing |
| [`cv.qmd`](cv.qmd) | CV summary (extend or add `assets/pdf/moseley_cv.pdf`) |
| [`assets/`](assets/) | Images, [`assets/json/resume.json`](assets/json/resume.json) |

The empty [`.nojekyll`](.nojekyll) file is copied into `_site` so GitHub Pages does not run Jekyll on the published HTML.

## License

See [LICENSE](LICENSE) (inherited from the prior al-folio template where applicable).
