# Makefile for regenerating the ProteomeGenerator3 pipeline metromap.
#
# The metromap is authored in docs/metromap.mmd (nf-metro syntax). The SVG and
# HTML renders are produced with nf-metro. Note that --theme and --x-spacing are
# render-time flags (they are NOT stored in the .mmd source), so they must be
# passed on every render — this Makefile captures them.
#
# Requirements:
#   pip install nf-metro
#
# Usage:
#   make -C docs            # regenerate both SVG and HTML
#   make -C docs metromap.svg
#   make -C docs metromap.html

MMD       := metromap.mmd
SVG       := metromap.svg
HTML      := metromap.html
THEME     := light
X_SPACING := 80

.PHONY: all clean
all: $(SVG) $(HTML)

$(SVG): $(MMD)
	nf-metro render $< -o $@ --theme $(THEME) --x-spacing $(X_SPACING)

$(HTML): $(MMD)
	nf-metro render $< -o $@ --format html --theme $(THEME) --x-spacing $(X_SPACING)

clean:
	rm -f $(SVG) $(HTML)
