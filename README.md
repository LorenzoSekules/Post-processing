# Post-processing of reflector radiation patterns

This repo has the far-field output files of two ideal offset parabolic reflectors. Both are fed by a Gaussian RHCP feed, so the reflected field is LHCP:

| Folder     | Reflector                       | Files                                      |
|------------|---------------------------------|--------------------------------------------|
| `S_band/`  | D = 30 m, F = 30 m, f = 2.2 GHz | `*_uv` (u,v grid), `*_tp` (theta-phi cuts) |
| `Ka_band/` | D = 10 m, F = 10 m, f = 30 GHz  | `*_uv` (u,v grid), `*_tp` (theta-phi cuts) |

`output_files_description_and_reflectors_optics.pptx` describes the file format and the reflector optics.

## MATLAB script `reflector_pattern_analysis.m`

Open MATLAB in this folder and press **Run**. You can also run one section at a time with *Run Section* (Ctrl+Enter); every intermediate result stays in the workspace. No toolbox is needed.

| Section | Content |
|---------|---------|
| 0 | user configuration (files, frequency, geometry, options) |
| 1 | reading of the pattern files |
| 2 | **[A] theoretical expected values** (uniform aperture, reflector with a -11 dB taper, CP beam squint) |
| 3 | **[B] performance from the raw file data** (file samples only) |
| 4 | **[C] performance from the interpolated data** (cubic spline of the complex field) |
| 5 | figures |
| 6 | report (screen, `results/performance_summary.txt`, `results/performance_results.mat`) |

Symbols: **G** is the co-polar (LHCP) gain and **G_x** the cross-polar (RHCP) gain. **eta** = G_max / (pi D/lambda)^2 is the aperture efficiency. **HPBW** is the -3 dB beamwidth. **FML** is the first minor lobe, and **FMLL** its level relative to G_max. **XPD** = G_max - G_x,max is the cross-polar discrimination. **AR** is the axial ratio.

### Interpolation

The script interpolates the complex far field E = 10^(G/20) exp(j phase), built from both the amplitude and the phase in the files, and then computes G = 20 log10|E|. The real and imaginary parts are interpolated with cubic splines. The far field is a smooth function of the direction. G in dB is not: it goes to -infinity at the nulls, and a spline applied to dB values oscillates there.

The method was checked on the S-band cut, subsampled to 2–2.5 samples per HPBW (the Ka-band case):

* complex-field spline: HPBW error 0.1–0.4 %, G_max error at most 0.03 dB;
* spline on dB values: HPBW error up to 28 %, G_max error up to 0.4 dB.

### Figures (`results/`)

| File | Content |
|------|---------|
| `<band>_uv_contour_full.png` | contours of G and G_x over the whole u,v grid |
| `<band>_uv_contour_zoom.png` | main beam: raw samples (left) and interpolation (right), G (top) and G_x (bottom); the plots mark the -3 dB contour, the main-beam region and the G_max and G_x,max points |
| `<band>_theta_cuts_full.png` | theta cuts from -90 to 90 deg in the 4 phi planes |
| `<band>_theta_cuts_zoom.png` | main beam of every plane, with samples, spline, G_max, -3 dB points, first nulls, first minor lobes with their FMLL, and G_x,max |
| `<band>_theta_cuts_overlay.png` | all the planes overlaid, to check the beam symmetry |
| `<band>_axial_ratio.png` | AR vs theta inside the main lobe |
| `<band>_uv_vs_tp_comparison.png` | the same plane from the two files |
| `comparison_S_vs_Ka_normalized.png` | both reflectors plotted against (pi D/lambda) sin(theta) |

### Notes

* The u-v grid step is 0.002, which is 0.115 deg. At Ka band the HPBW is 0.067 deg, so the grid has fewer than one sample per beamwidth there. The Ka-band u-v values are therefore not reliable, whether raw or interpolated, and the report says so. Use the theta-phi file for Ka band.
* Some values in the files are glued together, for example `0.89E-01-0.56`, because they are written in fixed-width Fortran fields. The reader handles this.
* The figures in `results/` were generated with GNU Octave as a reference. They include everything except `<band>_uv_contour_full.png`, which takes too long to draw in Octave. In MATLAB it takes a few seconds, and running the script there regenerates everything.
