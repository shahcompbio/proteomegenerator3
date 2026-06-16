# kentsislab/proteomegenerator3: Documentation

The kentsislab/proteomegenerator3 documentation is split into the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.

## Pipeline diagram

The metromap shown in the main [README](../README.md) is rendered from
[`metromap.mmd`](metromap.mmd) with [nf-metro](https://github.com/pinin4fjords/nf-metro).

To regenerate the `metromap.svg` and `metromap.html` renders after editing the
`.mmd` source:

```bash
pip install nf-metro
make -C docs
```

> **Note:** the `--theme light` and `--x-spacing 80` options are render-time
> flags and are **not** stored in `metromap.mmd`, so they must be passed on
> every render. The [`Makefile`](Makefile) captures them.

