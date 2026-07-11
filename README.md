# SPARTA DSMC Drag Simulation — BEP Handover

A DSMC (Direct Simulation Monte Carlo) study of the aerodynamic drag on spherical
tin micro-particles travelling through a low-pressure hydrogen buffer gas, run with
the [SPARTA](https://sparta.github.io/) code. TU/e — ASML Bachelor End Project (BEP),
Irene Lozano Gonzalez.

**Live documentation:** https://irenelozano.github.io/Final-BEP-website/

## What's in this repository

| Path | What it is |
|------|------------|
| `index.html` | The documentation website — renders the guide below. |
| `HANDOVER_GUIDE.txt` | The complete step-by-step handover guide (Markdown source). |
| `code/` | All project code: simulation scripts, physical model, and analysis. |
| `README.md` | This file. |

## The `code/` folder

| Folder | Contents | Guide section |
|--------|----------|---------------|
| `pipeline/` | SPARTA input templates (`in.drag`, …) and the `run_*.sh` scripts that drive the simulations. | §5, §6 |
| `global_assets/` | Physical model definitions — H₂ gas species, VHS collision parameters, sphere surface mesh. | §4 |
| `matlab postprocessing/` | MATLAB scripts that validate the results and generate the report figures. | §7 |

For details on every file, how to run each simulation, and how to reproduce the
results, read the handover guide — on the live site above or in `HANDOVER_GUIDE.txt`.

## How the website works

`index.html` is a single static page. On load it fetches `HANDOVER_GUIDE.txt`,
renders the Markdown, and builds the sidebar navigation. The two files must stay in
the same folder and keep their exact names. The page must be served over HTTP
(GitHub Pages) — opening `index.html` directly from disk (`file://`) will not load
the guide.

## Download everything

Use **Code → Download ZIP** on the repository home page to get the guide and all
code in a single archive.
