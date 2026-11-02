# Post-processing of reflector radiation patterns

This repo has the far-field output files of two ideal offset parabolic reflectors. Both are fed by a Gaussian RHCP feed, so the reflected field is LHCP:

| Folder     | Reflector                                 | Files                                   |
|------------|-------------------------------------------|-----------------------------------------|
| `S_band/`  | D = 30 m, F = 30 m, f = 2.2 GHz           | `*_uv` (u,v grid), `*_tp` (theta-phi cuts) |
| `Ka_band/` | D = 10 m, F = 10 m, f = 30 GHz            | `*_uv` (u,v grid), `*_tp` (theta-phi cuts) |

`output_files_description_and_reflectors_optics.pptx` describes the file format and the reflector optics.

## MATLAB script

`reflector_pattern_analysis.m` reads the four files. It plots the patterns and computes the performance metrics. To run it, open MATLAB in this folder and press **Run** (or type `reflector_pattern_analysis`). It needs no toolbox, and it also runs in GNU Octave, although Octave is much slower for the filled contour plots.

### What it computes

* peak gain (co-polar LHCP). It gives the raw sample and a parabolic interpolation. It also gives the aperture efficiency relative to `(pi D / lambda)^2`.
* peak location as (u,v) and (theta,phi). This shows the beam squint caused by circular polarisation in an offset reflector.
* HPBW in the phi = 0/45/90/135 deg planes, from both files.
* first null and first side-lobe level.
* cross-polar (RHCP) peak and its location. The peak is computed inside the main lobe and over the whole file.
* XPD, the worst XPD inside the -3 dB beam, and the axial ratio at the peak.
* reference values from theory (uniform aperture directivity and HPBW, CP beam squint).

### Output (`results/`)

* `results/<band>/<band>_uv_contour_full.png`: contour of co- and cross-polar gain over the whole u,v grid.
* `results/<band>/<band>_uv_contour_zoom.png`: the same plot zoomed on the main beam, with the -3 dB contour.
* `results/<band>/<band>_theta_cuts_full.png`: theta cuts from -90 to 90 deg in four phi planes.
* `results/<band>/<band>_theta_cuts_zoom.png`: the same cuts zoomed on the beam, with the HPBW and the peak marked.
* `results/<band>/<band>_theta_cuts_overlay.png`: all the planes overlaid, to check the beam symmetry.
* `results/<band>/<band>_uv_vs_tp_comparison.png`: the same plane taken from the uv file and from the theta-phi file.
* `results/comparison_S_vs_Ka_normalized.png`: both reflectors plotted against `(pi D/lambda) sin(theta)`.
* `results/performance_summary.txt`: the table of all the metrics.
* `results/performance_results.mat`: all the numbers as a MATLAB struct.

### Notes

* The uv grid step is 0.002, which is about 0.115 deg. At Ka band the HPBW is about 0.07 deg, so the uv grid has fewer than 1 sample per beamwidth there. Use the theta-phi file for the Ka-band gain, HPBW and cross-pol. The script prints a warning about this.
* Some values in the files are glued together, for example `0.89E-01-0.56`, because they are written in fixed-width Fortran fields. The reader handles this.
