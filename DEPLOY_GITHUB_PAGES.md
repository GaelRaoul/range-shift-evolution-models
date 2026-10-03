# GitHub Pages deployment checklist

If the whole scientific-code repository is used, the simplest approach is to publish the `interactive/` directory with a GitHub Pages workflow or to copy its contents to the branch/folder chosen for Pages.

For a dedicated simulator repository:

1. Copy all files from `interactive/` to the repository root.
2. Commit and push to `main`.
3. Open **Settings → Pages**.
4. Under **Build and deployment**, select **Deploy from a branch**.
5. Select branch **main** and folder **/(root)**.
6. Save.
7. Open the URL shown by GitHub Pages and run short A/B comparisons covering all three models.
8. Once available, put the article/archive/source URLs in `config.js`.

Recommended smoke test:

- test a current Firefox and Chromium/Chrome browser;
- select each of the three model types in Simulation A or B;
- verify the heatmap and the `N`, `Z`, `Zopt`, and `V` displays;
- verify that the two range-marker positions agree across the four plots for a simulation;
- test more than one `VLE` value in the sexual model;
- for the stochastic browser model, verify that repeating the same user-selected seed reproduces the same interactive realization.

The browser seed is only an interactive-simulator control; the archived research stochastic script does not set an explicit seed.
