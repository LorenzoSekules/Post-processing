function results = reflector_pattern_analysis()
%REFLECTOR_PATTERN_ANALYSIS  Post-processing of reflector radiation patterns.
%
%   results = reflector_pattern_analysis()
%
%   Reads the far-field radiation pattern files of the two ideal offset
%   parabolic reflectors (S-band, D = 30 m and Ka-band, D = 10 m), plots the
%   radiation patterns and extracts the main antenna performance metrics.
%
%   Just press "Run" (F5) in MATLAB, or type  reflector_pattern_analysis
%   in the Command Window while this file's folder is the current folder.
%
%   ------------------------------------------------------------------------
%   INPUT FILES (two per reflector, in the folders S_band/ and Ka_band/)
%   ------------------------------------------------------------------------
%   *_uv : pattern on a regular (u,v) grid, u = v = [-0.4, 0.4], 401 points
%   *_tp : pattern on a (theta,phi) grid: theta = [0, 90] deg, 3601 points,
%          phi = [-180, 180] deg, 9 points (one cut every 45 deg)
%
%   Both files are expressed in the antenna boresight reference frame
%   (z axis = boresight), so theta = 0 <=> u = v = 0, with the usual
%   direction cosines (Balanis, "Antenna Theory", ch. 12-15):
%
%         u = sin(theta) cos(phi)          v = sin(theta) sin(phi)
%
%   File format (see slide 2 of the ppt note):
%     line 1      : title (which field components / which coordinates)
%     line 2      : N1  N2  -> number of values of 1st and 2nd coordinate
%     next lines  : the N1 values of the 1st coordinate (u or theta),
%                   5 values per line
%     then, for every value of the 2nd coordinate (v or phi), 4 blocks:
%        "<2nd coord> = value  E-l amplitude"  -> N1 values in dB
%        "<2nd coord> = value  E-l phase"      -> N1 values in deg
%        "<2nd coord> = value  E-r amplitude"  -> N1 values in dB
%        "<2nd coord> = value  E-r phase"      -> N1 values in deg
%     E-l = left-hand circular (LHCP) component, E-r = right-hand circular
%     (RHCP) component. The feed is RHCP, the reflection inverts the sense
%     of rotation, so the reflector radiates LHCP:
%        CO-POLAR    component = E-l (LHCP)
%        CROSS-POLAR component = E-r (RHCP)
%     The amplitudes are absolute gains in dBi (e.g. 55.5 dBi on boresight
%     for the 30 m dish at 2.2 GHz, as expected from (pi*D/lambda)^2).
%     The value -100 dB is used by the software as a "zero field" floor.
%
%   ------------------------------------------------------------------------
%   PERFORMANCE METRICS COMPUTED (for each reflector and each file type)
%   ------------------------------------------------------------------------
%   * Peak gain (co-polar)          [dBi]  raw sample and interpolated value
%   * Peak location                 (u,v) and (theta,phi)
%   * Half-power beamwidth (HPBW)   [deg]  in the phi = 0/45/90/135 planes
%   * First null and first side-lobe level (SLL) per plane
%   * Cross-polar peak              [dBi]  global and inside the main lobe
%   * Cross-polar peak location     (u,v) and (theta,phi)
%   * Cross-polar discrimination    XPD = co-polar peak - cross-polar peak
%   * Axial ratio at the beam peak  [dB]
%   * Reference values from theory: maximum directivity of a uniformly
%     illuminated circular aperture (pi*D/lambda)^2 -> aperture efficiency,
%     uniform-aperture HPBW, beam squint of an offset reflector fed by a
%     circularly polarised feed (Duan & Rahmat-Samii, 1991).
%
%   ------------------------------------------------------------------------
%   OUTPUT
%   ------------------------------------------------------------------------
%   * figures (PNG) in  results/<band>/  and  results/
%   * a text report     results/performance_summary.txt
%   * the struct array "results" with all the numbers (also saved in
%     results/performance_results.mat)
%
%   The code uses only base MATLAB functions (no toolbox required). It has
%   also been tested with GNU Octave 8.

% =========================================================================
% 0) USER CONFIGURATION
% =========================================================================
% Folder that contains this file: data folders and results are relative to it
rootDir = fileparts(mfilename('fullpath'));

cfg.outDir      = fullfile(rootDir, 'results');  % where figures/report go
cfg.saveFigures = true;      % save every figure as PNG
cfg.closeFigs   = false;     % close figures after saving (useful in batch)
cfg.dpi         = 150;       % resolution of the PNG files
cfg.coComp      = 'El';      % co-polar component  (E-l = LHCP)
cfg.xpComp      = 'Er';      % cross-polar component (E-r = RHCP)
cfg.planesDeg   = [0 45 90 135];   % phi planes used for the 1-D cuts
cfg.floorDb     = -99;       % values <= this are the "-100 dB" floor -> ignored
cfg.c0          = 299792458; % speed of light [m/s]

% Description of the two antennas (data from the ppt note / supervisor mail)
%   freqHz       : operating frequency
%   D, F         : projected aperture diameter and focal length [m]
%   offsetDeg    : offset angle of the feed (aoff in the ppt)
%   thetaZoomDeg : half-width of the zoomed theta cuts  [deg]
%   uvZoom       : half-width of the zoomed uv contour   [-]
ant = struct([]);
ant(1).name   = 'S-band';
ant(1).tag    = 'S_band';
ant(1).folder = fullfile(rootDir, 'S_band');
ant(1).uvFile = 'radiation_pattern_30m_ideal_reflector_S_band_LHCP_uv';
ant(1).tpFile = 'radiation_pattern_30m_ideal_reflector_S_band_LHCP_tp';
ant(1).freqHz = 2.2e9;
ant(1).D      = 30;
ant(1).F      = 30;
ant(1).offsetDeg    = 53.44;
ant(1).thetaZoomDeg = 1.5;
ant(1).uvZoom       = 0.025;

ant(2).name   = 'Ka-band';
ant(2).tag    = 'Ka_band';
ant(2).folder = fullfile(rootDir, 'Ka_band');
ant(2).uvFile = 'radiation_pattern_10m_ideal_reflector_Ka_band_LHCP_uv';
ant(2).tpFile = 'radiation_pattern_10m_ideal_reflector_Ka_band_LHCP_tp';
ant(2).freqHz = 30e9;
ant(2).D      = 10;
ant(2).F      = 10;
ant(2).offsetDeg    = 34.18;
ant(2).thetaZoomDeg = 0.5;
ant(2).uvZoom       = 0.02;

if ~exist(cfg.outDir, 'dir'); mkdir(cfg.outDir); end

% =========================================================================
% 1) ANALYSIS OF EACH REFLECTOR
% =========================================================================
results = struct([]);
for k = 1:numel(ant)
    fprintf('\n=============== %s reflector ===============\n', ant(k).name);
    r = analyze_antenna(ant(k), cfg);
    if isempty(results); results = r; else; results(k) = r; end
end

% =========================================================================
% 2) COMPARISON BETWEEN THE TWO REFLECTORS
% =========================================================================
plot_normalized_comparison(results, cfg);

% =========================================================================
% 3) REPORT (command window + text file + .mat file)
% =========================================================================
write_report(results, 1);                                   % 1 = screen
fid = fopen(fullfile(cfg.outDir, 'performance_summary.txt'), 'w');
write_report(results, fid);
fclose(fid);
save(fullfile(cfg.outDir, 'performance_results.mat'), 'results');
fprintf('\nFigures and report saved in: %s\n', cfg.outDir);
end


% #########################################################################
%                           LOCAL FUNCTIONS
% #########################################################################

function r = analyze_antenna(a, cfg)
%ANALYZE_ANTENNA  Read the two files of one reflector, compute the metrics
%and draw all the figures of that reflector.

outDir = fullfile(cfg.outDir, a.tag);
if ~exist(outDir, 'dir'); mkdir(outDir); end

% -------------------------------------------------------------------------
% 1.1 Read the files
% -------------------------------------------------------------------------
fprintf('Reading %s ...\n', a.uvFile);
uv = read_pattern_file(fullfile(a.folder, a.uvFile));
fprintf('Reading %s ...\n', a.tpFile);
tp = read_pattern_file(fullfile(a.folder, a.tpFile));

% -------------------------------------------------------------------------
% 1.2 Theoretical reference values (Balanis ch. 15 - reflector antennas)
% -------------------------------------------------------------------------
th.lambda  = cfg.c0 / a.freqHz;                    % wavelength [m]
% maximum directivity of a uniformly illuminated circular aperture of
% diameter D: D0 = (pi*D/lambda)^2 (aperture efficiency = 100 %)
th.Dmax_dBi = 10*log10((pi*a.D/th.lambda)^2);
% HPBW of a uniform circular aperture: 1.02*lambda/D rad = 58.4*lambda/D deg.
% A tapered illumination (here ~ -11 dB at the edge) widens the beam, a
% common rule of thumb for reflectors being ~70*lambda/D deg.
th.hpbwUniform_deg = 58.4 * th.lambda / a.D;
th.hpbwRule70_deg  = 70   * th.lambda / a.D;
% Beam squint of an offset paraboloid fed by a circularly polarised feed
% (Duan & Rahmat-Samii, IEEE TAP 1991): the beam moves off boresight, in
% the plane orthogonal to the offset plane (here the y-z plane, phi=+-90)
%     sin(psi_s) = lambda * sin(theta_offset) / (4*pi*F)
th.squint_deg = asind(th.lambda * sind(a.offsetDeg) / (4*pi*a.F));

fprintf('lambda = %.4f m, D/lambda = %.1f, Dmax(uniform) = %.2f dBi\n', ...
        th.lambda, a.D/th.lambda, th.Dmax_dBi);

% -------------------------------------------------------------------------
% 1.3 Theta-phi file: build the full plane cuts and compute the metrics
% -------------------------------------------------------------------------
% The file gives "half cuts" theta = 0..90 for phi = -180:45:180. Two half
% cuts at phi0 and phi0-180 form a full plane cut, which we describe with a
% signed theta in [-90, 90] deg (theta > 0 on the phi0 side).
cuts = build_plane_cuts(tp, cfg);
dTheta = tp.c1(2) - tp.c1(1);               % angular step of the file [deg]

for p = 1:numel(cuts)
    cuts(p).m = cut_metrics(cuts(p).theta, cuts(p).co, cfg.floorDb);
end

% best (interpolated) peak among the planes
pkVals = arrayfun(@(c) c.m.peakInterp, cuts);
[~, pBest] = max(pkVals);
tpm.peakGain_dBi   = cuts(pBest).m.peakInterp;
tpm.peakThetaSgn   = cuts(pBest).m.peakPosInterp;
[tpm.peakTheta_deg, tpm.peakPhi_deg] = signed_to_theta_phi( ...
                                tpm.peakThetaSgn, cuts(pBest).phi0);
% raw (not interpolated) maximum over every sample of the file
CO = tp.amp.(cfg.coComp);  XP = tp.amp.(cfg.xpComp);
[tpm.peakGainRaw_dBi, iMax] = max(CO(:));
[iPhi, iTh] = ind2sub(size(CO), iMax);
tpm.peakThetaRaw_deg = tp.c1(iTh);
tpm.peakPhiRaw_deg   = tp.c2(iPhi);

% HPBW, first null and side lobes per plane
tpm.hpbw_deg   = arrayfun(@(c) c.m.hpbw, cuts);
tpm.null_deg   = arrayfun(@(c) c.m.firstNull, cuts);
tpm.sll_dB     = arrayfun(@(c) c.m.sllRel, cuts);
tpm.samplesPerHpbw = mean(tpm.hpbw_deg) / dTheta;

% Cross-polar component.
% (a) global maximum over the whole file (all theta in [0,90] deg)
XPv = XP;  XPv(XPv <= cfg.floorDb) = -Inf;
[tpm.xpGlobal_dBi, iX] = max(XPv(:));
[iPhi, iTh] = ind2sub(size(XP), iX);
tpm.xpGlobalTheta_deg = tp.c1(iTh);
tpm.xpGlobalPhi_deg   = tp.c2(iPhi);
% (b) maximum inside the main lobe (theta <= mean first-null angle): this
%     is the cross-polar peak relevant for the link (it sits close to the
%     beam, where the co-polar gain is high)
tpm.mainLobeTheta_deg = mean(tpm.null_deg(~isnan(tpm.null_deg)));
mask = repmat(tp.c1(:).' <= tpm.mainLobeTheta_deg, numel(tp.c2), 1);
XPm = XPv;  XPm(~mask) = -Inf;
[tpm.xpMain_dBi, iX] = max(XPm(:));
[iPhi, iTh] = ind2sub(size(XP), iX);
tpm.xpMainTheta_deg = tp.c1(iTh);
tpm.xpMainPhi_deg   = tp.c2(iPhi);
% (c) cross-polar discrimination
tpm.XPD_main_dB   = tpm.peakGain_dBi - tpm.xpMain_dBi;
tpm.XPD_global_dB = tpm.peakGain_dBi - tpm.xpGlobal_dBi;
% (d) worst XPD (co - cross, point by point) inside the -3 dB beam
in3 = CO >= tpm.peakGain_dBi - 3;
tpm.XPD_worst3dB_dB = min(CO(in3) - XP(in3));
% (e) axial ratio at the raw peak sample, from the two CP components:
%     AR = (|EL| + |ER|) / (|EL| - |ER|)   (Balanis eq. 2-66 / 2-67)
tpm.AR_peak_dB = axial_ratio_dB(tpm.peakGainRaw_dBi, XP(iMax));

% -------------------------------------------------------------------------
% 1.4 uv file: metrics on the 2-D grid
% -------------------------------------------------------------------------
[U, V] = meshgrid(uv.c1, uv.c2);       % rows <-> v, columns <-> u
COuv = uv.amp.(cfg.coComp);
XPuv = uv.amp.(cfg.xpComp);
du   = uv.c1(2) - uv.c1(1);

% raw maximum
[uvm.peakGainRaw_dBi, iMax] = max(COuv(:));
[iv, iu] = ind2sub(size(COuv), iMax);
% sub-grid refinement: a parabola (in dB) through the 3 samples around the
% maximum along u and along v. For a separable quadratic surface the peak
% value is  y0 + (yu - y0) + (yv - y0).
[uPk, yu] = parabolic_peak(uv.c1, COuv(iv, :), iu);
[vPk, yv] = parabolic_peak(uv.c2, COuv(:, iu).', iv);
uvm.peakGain_dBi = yu + yv - uvm.peakGainRaw_dBi;
uvm.peakU = uPk;   uvm.peakV = vPk;
[uvm.peakTheta_deg, uvm.peakPhi_deg] = uv_to_theta_phi(uPk, vPk);
uvm.peakURaw = uv.c1(iu);  uvm.peakVRaw = uv.c2(iv);

% cross-polar: global over the uv window and inside the main lobe
XPuvV = XPuv;  XPuvV(XPuvV <= cfg.floorDb) = -Inf;
[uvm.xpGlobal_dBi, iX] = max(XPuvV(:));
uvm.xpGlobalU = U(iX);  uvm.xpGlobalV = V(iX);
[uvm.xpGlobalTheta_deg, uvm.xpGlobalPhi_deg] = uv_to_theta_phi(U(iX), V(iX));
inMain = sqrt(U.^2 + V.^2) <= sind(tpm.mainLobeTheta_deg);
XPm = XPuvV;  XPm(~inMain) = -Inf;
[uvm.xpMain_dBi, iX] = max(XPm(:));
uvm.xpMainU = U(iX);  uvm.xpMainV = V(iX);
[uvm.xpMainTheta_deg, uvm.xpMainPhi_deg] = uv_to_theta_phi(U(iX), V(iX));
uvm.XPD_main_dB = uvm.peakGain_dBi - uvm.xpMain_dBi;
in3 = COuv >= uvm.peakGain_dBi - 3;
uvm.nSamplesIn3dB = nnz(in3);
uvm.XPD_worst3dB_dB = min(COuv(in3) - XPuv(in3));

% HPBW from the uv grid: the 2-D pattern is extracted along the same phi
% planes used for the theta-phi file, so that the two files can be
% compared directly. Two steps:
%   1) the central part of the uv grid (|u|,|v| <= W, W covering the zoom
%      plots) is up-sampled 10 times with a 2-D cubic spline;
%   2) the plane cuts are taken from this fine grid by (bi)linear
%      interpolation, which is now accurate because the grid is fine.
% (A direct interp2(...,'spline') on scattered points is not supported by
% GNU Octave and 'cubic' on scattered points is extremely slow there.)
W   = max(a.uvZoom, 1.05*sind(a.thetaZoomDeg));     % half-width of window
iw  = find(abs(uv.c1) <= W + 3*du);                 % original samples used
uw  = uv.c1(iw);   vw = uv.c2(iw);
uf  = linspace(min(uw), max(uw), 10*(numel(uw)-1) + 1);   % fine axes
vf  = linspace(min(vw), max(vw), 10*(numel(vw)-1) + 1);
[Uf, Vf] = meshgrid(uf, vf);
COf = interp2(uw(:).', vw(:), COuv(iw, iw), Uf, Vf, 'spline');
XPf = interp2(uw(:).', vw(:), XPuv(iw, iw), Uf, Vf, 'spline');

thFine  = linspace(-asind(W), asind(W), 4001);
uvCuts  = struct([]);
for p = 1:numel(cfg.planesDeg)
    phi0 = cfg.planesDeg(p);
    uu = sind(thFine) * cosd(phi0);
    vv = sind(thFine) * sind(phi0);
    c.phi0  = phi0;
    c.theta = thFine;
    c.co    = interp2(Uf, Vf, COf, uu, vv, 'linear');
    c.xp    = interp2(Uf, Vf, XPf, uu, vv, 'linear');
    c.m     = cut_metrics(c.theta, c.co, cfg.floorDb);
    if isempty(uvCuts); uvCuts = c; else; uvCuts(p) = c; end
end
uvm.hpbw_deg = arrayfun(@(c) c.m.hpbw, uvCuts);
uvm.samplesPerHpbw = mean(uvm.hpbw_deg) / asind(du);
uvm.theta_step_deg = asind(du);

% -------------------------------------------------------------------------
% 1.5 Figures
% -------------------------------------------------------------------------
plot_uv_contours(a, uv, U, V, COuv, XPuv, uvm, tpm, cfg, outDir);
plot_tp_cuts(a, cuts, tpm, cfg, outDir);
plot_uv_vs_tp(a, cuts, uvCuts, uv, COuv, cfg, outDir);

% -------------------------------------------------------------------------
% 1.6 Collect results
% -------------------------------------------------------------------------
r.name   = a.name;
r.ant    = a;
r.theory = th;
r.tp     = tpm;
r.uv     = uvm;
r.cutsTp = cuts;
r.cutsUv = uvCuts;
r.planesDeg = cfg.planesDeg;
r.tpStep_deg = dTheta;
end


% =========================================================================
function P = read_pattern_file(fname)
%READ_PATTERN_FILE  Generic reader of the near/far-field output files.
%
%  P.title        first line of the file
%  P.c1  (N1x1)   values of the 1st coordinate (u or theta)
%  P.c2  (N2x1)   values of the 2nd coordinate (v or phi)
%  P.name2        name of the 2nd coordinate ('v' or 'phi')
%  P.amp.El, P.amp.Er    (N2 x N1) amplitudes in dB
%  P.ph.El,  P.ph.Er     (N2 x N1) phases in deg
%  -> row i corresponds to c2(i), column j to c1(j)

if ~exist(fname, 'file')
    error('File not found: %s', fname);
end
txt   = fileread(fname);
lines = regexp(txt, '\r?\n', 'split');
lines = lines(~cellfun(@isempty, strtrim(lines)));    % drop empty lines

% Header lines contain letters (other than the 'E' of the exponent)
isHdr  = ~cellfun(@isempty, regexp(lines, '[a-df-zA-DF-Z]', 'once'));
hdrIdx = find(isHdr);

P.title = strtrim(lines{1});
N = sscanf(lines{2}, '%d');                           % [N1; N2]
N1 = N(1);  N2 = N(2);

% list of the first coordinate: from line 3 up to the first block header
P.c1 = parse_numbers(lines(3:hdrIdx(2)-1));
if numel(P.c1) ~= N1
    error('%s: expected %d values of the 1st coordinate, found %d', ...
          fname, N1, numel(P.c1));
end

% data blocks
blk = hdrIdx(2:end);
blkEnd = [blk(2:end) - 1, numel(lines)];
P.c2 = nan(N2, 1);
P.amp.El = nan(N2, N1);  P.amp.Er = nan(N2, N1);
P.ph.El  = nan(N2, N1);  P.ph.Er  = nan(N2, N1);
i2 = 0;  lastVal = NaN;
for b = 1:numel(blk)
    h = lines{blk(b)};
    % "phi =  -180.0000  degrees, E-l amplitude" / "v = -0.4  E-r phase"
    tok = regexp(h, '^\s*(\w+)\s*=\s*([-+0-9.EeDd]+)', 'tokens', 'once');
    P.name2 = tok{1};
    val = str2double(tok{2});
    if b == 1 || val ~= lastVal          % new value of the 2nd coordinate
        i2 = i2 + 1;
        P.c2(i2) = val;
        lastVal = val;
    end
    if ~isempty(strfind(h, 'E-l')); comp = 'El'; else; comp = 'Er'; end
    data = parse_numbers(lines(blk(b)+1:blkEnd(b)));
    if numel(data) ~= N1
        error('%s: block "%s" has %d values instead of %d', ...
              fname, strtrim(h), numel(data), N1);
    end
    if ~isempty(strfind(h, 'amplitude'))
        P.amp.(comp)(i2, :) = data;
    else
        P.ph.(comp)(i2, :) = data;
    end
end
if i2 ~= N2
    error('%s: expected %d values of the 2nd coordinate, found %d', ...
          fname, N2, i2);
end
end


% =========================================================================
function x = parse_numbers(lineCell)
%PARSE_NUMBERS  Convert a group of text lines into a column of numbers.
% The files are written with fixed-width Fortran fields, so two negative
% numbers can be glued together, e.g. "0.8900023E-01-0.5656819".
% A blank is therefore inserted before every '-' that follows a digit
% (the '-' of an exponent follows an 'E' and is left untouched).
s = sprintf('%s ', lineCell{:});
s = regexprep(s, '(?<=[0-9.])-', ' -');
x = sscanf(s, '%f');
end


% =========================================================================
function cuts = build_plane_cuts(tp, cfg)
%BUILD_PLANE_CUTS  Join the half cuts phi0 and phi0-180 into full plane
%cuts with a signed theta in [-90, 90] deg.
th = tp.c1(:).';
cuts = struct([]);
for p = 1:numel(cfg.planesDeg)
    phi0 = cfg.planesDeg(p);
    iF = find_phi(tp.c2, phi0);           % forward half  (theta > 0)
    iB = find_phi(tp.c2, phi0 - 180);     % backward half (theta < 0)
    c.phi0  = phi0;
    c.theta = [-fliplr(th(2:end)), th];
    c.co = [fliplr(tp.amp.(cfg.coComp)(iB, 2:end)), tp.amp.(cfg.coComp)(iF, :)];
    c.xp = [fliplr(tp.amp.(cfg.xpComp)(iB, 2:end)), tp.amp.(cfg.xpComp)(iF, :)];
    if isempty(cuts); cuts = c; else; cuts(p) = c; end
end
end

function i = find_phi(phiList, phi)
% index of the angle phi in phiList, accounting for the 360 deg periodicity
d = mod(phiList(:) - phi + 180, 360) - 180;
[dmin, i] = min(abs(d));
if dmin > 1e-3
    error('phi = %g deg not present in the theta-phi file', phi);
end
end


% =========================================================================
function m = cut_metrics(x, y, floorDb)
%CUT_METRICS  Peak, HPBW, first nulls and first side lobes of a 1-D cut.
%  x : signed theta [deg] (monotonic),  y : co-polar gain [dBi]
y = y(:).';  x = x(:).';
y(y <= floorDb) = NaN;

[m.peakRaw, i0] = max(y);
m.peakPosRaw = x(i0);
[m.peakPosInterp, m.peakInterp] = parabolic_peak(x, y, i0);

% --- half power beamwidth: -3 dB crossings on both sides of the peak ---
lev = m.peakInterp - 3;
iR = find(y(i0:end) < lev, 1, 'first') + i0 - 1;   % first sample below
iL = find(y(1:i0)   < lev, 1, 'last');             % on each side
if isempty(iR) || isempty(iL)
    m.hpbw = NaN;  m.x3L = NaN;  m.x3R = NaN;
else
    m.x3R = level_crossing(x, y, iR-1, iR, lev);
    m.x3L = level_crossing(x, y, iL, iL+1, lev);
    m.hpbw = m.x3R - m.x3L;
end

% --- first nulls: walk away from the peak while the pattern decreases ---
iNR = i0;
while iNR < numel(y) && y(iNR+1) < y(iNR); iNR = iNR + 1; end
iNL = i0;
while iNL > 1 && y(iNL-1) < y(iNL); iNL = iNL - 1; end
m.nullR = x(iNR);  m.nullL = x(iNL);
m.firstNull = 0.5 * (abs(m.nullR - m.peakPosInterp) + ...
                     abs(m.nullL - m.peakPosInterp));   % mean null angle

% --- first side lobes: walk on from the nulls while the pattern rises ---
iSR = iNR;
while iSR < numel(y) && y(iSR+1) > y(iSR); iSR = iSR + 1; end
iSL = iNL;
while iSL > 1 && y(iSL-1) > y(iSL); iSL = iSL - 1; end
m.slR = y(iSR);  m.slL = y(iSL);
m.slPosR = x(iSR);  m.slPosL = x(iSL);
m.sllRel = max(m.slR, m.slL) - m.peakInterp;     % worst first SLL [dB]
end


% =========================================================================
function [xp, yp] = parabolic_peak(x, y, i0)
%PARABOLIC_PEAK  Refine a maximum with a parabola through 3 samples.
% Close to its maximum the main beam (in dB) is well approximated by a
% parabola: y = a*(x-xp)^2 + yp. With 3 samples the vertex is found in
% closed form. On the edges of the grid the raw sample is returned.
x = x(:).';  y = y(:).';
if i0 <= 1 || i0 >= numel(y) || any(isnan(y(i0-1:i0+1)))
    xp = x(i0);  yp = y(i0);  return
end
c = polyfit(x(i0-1:i0+1) - x(i0), y(i0-1:i0+1), 2);
if c(1) >= 0                       % not a maximum: keep the sample
    xp = x(i0);  yp = y(i0);  return
end
dx = -c(2) / (2*c(1));
xp = x(i0) + dx;
yp = polyval(c, dx);
end


% =========================================================================
function xc = level_crossing(x, y, i1, i2, lev)
%LEVEL_CROSSING  Abscissa where the cut crosses the level "lev" between
%samples i1 and i2. A cubic spline through the neighbouring samples is
%evaluated on a fine grid (more accurate than a straight line when the
%beam is sampled with few points), then a linear interpolation is used.
n  = numel(y);
w  = max(1, i1-3):min(n, i2+3);
w  = w(~isnan(y(w)));
xf = linspace(x(i1), x(i2), 2001);
yf = interp1(x(w), y(w), xf, 'spline');
k  = find(diff(sign(yf - lev)) ~= 0, 1, 'first');
if isempty(k)                       % fall back to linear interpolation
    xc = interp1([y(i1) y(i2)], [x(i1) x(i2)], lev);
else
    xc = xf(k) + (lev - yf(k)) * (xf(k+1) - xf(k)) / (yf(k+1) - yf(k));
end
end


% =========================================================================
function [theta, phi] = uv_to_theta_phi(u, v)
% u = sin(theta)cos(phi), v = sin(theta)sin(phi)
theta = asind(sqrt(u.^2 + v.^2));
phi   = atan2d(v, u);
end

function [theta, phi] = signed_to_theta_phi(thSigned, phi0)
% signed theta of a plane cut -> (theta >= 0, phi in (-180,180])
theta = abs(thSigned);
if thSigned >= 0; phi = phi0; else; phi = phi0 - 180; end
phi = mod(phi + 180, 360) - 180;
end

function ar = axial_ratio_dB(coDb, xpDb)
% Axial ratio from the two circular components (amplitudes in dB):
% AR = (|E_co| + |E_x|) / (|E_co| - |E_x|)
eco = 10.^(coDb/20);  ex = 10.^(xpDb/20);
ar = 20*log10((eco + ex) ./ abs(eco - ex));
end


% =========================================================================
%                               PLOTS
% =========================================================================
function plot_uv_contours(a, uv, U, V, CO, XP, uvm, tpm, cfg, outDir)
%PLOT_UV_CONTOURS  Contour plots of co- and cross-polar gain in the uv plane
%(full grid and zoom on the main beam).

pk  = uvm.peakGain_dBi;
xpk = uvm.xpGlobal_dBi;
coLev = floor(pk) - (60:-3:0);          % iso-gain levels every 3 dB
xpLev = floor(xpk) - (60:-3:0);

% Values below the lowest level are clipped to it: this avoids empty
% (white) regions and thousands of tiny contour polygons in the noisy far
% side-lobe region, which make contourf very slow.
COc = max(CO, min(coLev));
XPc = max(XP, min(xpLev));

for zoom = [false true]
    fig = figure('Color', 'w', 'Position', [60 60 1250 520]);
    if zoom
        lim = a.uvZoom;  sfx = 'zoom';
    else
        lim = max(abs(uv.c1));  sfx = 'full';
    end
    % only the part of the grid inside the plotted window is contoured
    ii = find(abs(uv.c1) <= lim + 2*(uv.c1(2) - uv.c1(1)));
    Uz = U(ii, ii);  Vz = V(ii, ii);
    % ---- co-polar ----
    subplot(1, 2, 1);
    contourf(Uz, Vz, COc(ii, ii), coLev, 'LineColor', 'none');  hold on;
    contour(Uz, Vz, CO(ii, ii), [pk-3 pk-3], 'w-', 'LineWidth', 1.5); % -3 dB
    plot(uvm.peakU, uvm.peakV, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    style_uv_axes(lim, coLev);
    title(sprintf(['%s  Co-polar (LHCP) [dBi]\npeak (+) %.2f dBi at (u,v) = ' ...
          '(%.4f, %.4f), white: -3 dB'], a.name, pk, nz(uvm.peakU), nz(uvm.peakV)));
    % ---- cross-polar ----
    subplot(1, 2, 2);
    contourf(Uz, Vz, XPc(ii, ii), xpLev, 'LineColor', 'none');  hold on;
    contour(Uz, Vz, CO(ii, ii), [pk-3 pk-3], 'w-', 'LineWidth', 1.5);
    plot(uvm.xpMainU, uvm.xpMainV, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
    style_uv_axes(lim, xpLev);
    title(sprintf(['%s  Cross-polar (RHCP) [dBi]\nmain-lobe peak (x) %.2f dBi at ' ...
          '(u,v) = (%.4f, %.4f)'], a.name, uvm.xpMain_dBi, nz(uvm.xpMainU), nz(uvm.xpMainV)));
    save_fig(fig, fullfile(outDir, sprintf('%s_uv_contour_%s', a.tag, sfx)), cfg);
end
end

function x = nz(x)
% avoid printing "-0.0000" for tiny negative numbers
x(abs(x) < 5e-5) = 0;
end

function style_uv_axes(lim, lev)
axis equal;  axis([-lim lim -lim lim]);
caxis([min(lev) max(lev)]);
colormap(jet(numel(lev) - 1));
colorbar;  grid on;  box on;
xlabel('u = sin\theta cos\phi');  ylabel('v = sin\theta sin\phi');
set(gca, 'Layer', 'top');
end


% -------------------------------------------------------------------------
function plot_tp_cuts(a, cuts, tpm, cfg, outDir)
%PLOT_TP_CUTS  Co- and cross-polar gain vs theta in the 4 phi planes:
%full range [-90, 90] deg and zoom on the main beam.
pk = tpm.peakGain_dBi;
for zoom = [false true]
    fig = figure('Color', 'w', 'Position', [60 60 1250 820]);
    for p = 1:numel(cuts)
        c = cuts(p);
        subplot(2, 2, p);
        plot(c.theta, c.co, 'b-', 'LineWidth', 1.2);  hold on;
        plot(c.theta, c.xp, 'r-', 'LineWidth', 1.0);
        if zoom
            lim = a.thetaZoomDeg;
            plot([-lim lim], (pk-3)*[1 1], 'k--');                  % -3 dB
            plot([c.m.x3L c.m.x3R], (c.m.peakInterp-3)*[1 1], 'ko', ...
                 'MarkerFaceColor', 'y');
            plot(c.m.peakPosInterp, c.m.peakInterp, 'k^', 'MarkerFaceColor', 'g');
            % lower limit: 70 dB below the peak, or lower if needed to
            % show the cross-polar peak inside the zoom window
            inZ = abs(c.theta) <= lim;
            yl = [min(pk-70, max(c.xp(inZ)) - 15), pk+3];
            leg = {'Co-pol (LHCP)', 'X-pol (RHCP)', '-3 dB level', ...
                   '-3 dB points', 'peak'};
            ttl = sprintf('%s  -  plane \\phi = %d deg\nHPBW = %.4f deg, peak %.2f dBi at \\theta = %.4f deg', ...
                   a.name, c.phi0, c.m.hpbw, c.m.peakInterp, nz(c.m.peakPosInterp));
            sfx = 'zoom';
        else
            lim = 90;  yl = [pk-100, pk+5];
            leg = {'Co-pol (LHCP)', 'X-pol (RHCP)'};
            ttl = sprintf('%s  -  plane \\phi = %d deg', a.name, c.phi0);
            sfx = 'full';
        end
        xlim([-lim lim]);  ylim(yl);  grid on;
        xlabel(sprintf('\\theta [deg]   (\\theta<0 \\leftrightarrow \\phi = %d deg)', ...
                       mod(c.phi0 - 180 + 180, 360) - 180));
        ylabel('Gain [dBi]');
        title(ttl);
        legend(leg, 'Location', 'northwest', 'FontSize', 7);
    end
    save_fig(fig, fullfile(outDir, sprintf('%s_theta_cuts_%s', a.tag, sfx)), cfg);
end

% overlay of the co-polar cuts of all planes (beam symmetry check)
fig = figure('Color', 'w', 'Position', [60 60 900 560]);
col = lines(numel(cuts));
for p = 1:numel(cuts)
    plot(cuts(p).theta, cuts(p).co, '-', 'Color', col(p,:), 'LineWidth', 1.3);
    hold on;
end
for p = 1:numel(cuts)
    plot(cuts(p).theta, cuts(p).xp, '--', 'Color', col(p,:), 'LineWidth', 1.0);
end
plot(a.thetaZoomDeg*[-1 1], (pk-3)*[1 1], 'k:');
xlim(a.thetaZoomDeg*[-1 1]);  ylim([pk-80, pk+3]);  grid on;
xlabel('\theta [deg] (signed)');  ylabel('Gain [dBi]');
title(sprintf('%s  -  co-polar (solid) and cross-polar (dashed), all planes', a.name));
legend([arrayfun(@(c) sprintf('co  \\phi=%d', c.phi0), cuts, 'UniformOutput', false), ...
        arrayfun(@(c) sprintf('x   \\phi=%d', c.phi0), cuts, 'UniformOutput', false)], ...
       'Location', 'eastoutside');
save_fig(fig, fullfile(outDir, sprintf('%s_theta_cuts_overlay', a.tag)), cfg);
end


% -------------------------------------------------------------------------
function plot_uv_vs_tp(a, cuts, uvCuts, uv, COuv, cfg, outDir)
%PLOT_UV_VS_TP  Consistency check between the two files: the theta-phi cut
%is compared with the same plane extracted from the uv grid (raw samples on
%the axes/diagonals of the grid, and spline interpolation between them).
fig = figure('Color', 'w', 'Position', [60 60 1250 820]);
n = numel(uv.c1);  ic = (n + 1) / 2;       % index of u = 0 / v = 0
for p = 1:numel(cuts)
    subplot(2, 2, p);
    plot(cuts(p).theta, cuts(p).co, 'b-', 'LineWidth', 1.3);  hold on;
    plot(uvCuts(p).theta, uvCuts(p).co, 'r--', 'LineWidth', 1.1);
    % raw uv samples lying exactly on this plane
    switch cuts(p).phi0
        case 0,   uu = uv.c1;  vv = zeros(size(uu));  g = COuv(ic, :);
        case 90,  vv = uv.c2;  uu = zeros(size(vv));  g = COuv(:, ic).';
        case 45,  uu = uv.c1;  vv = uv.c2;            g = diag(COuv).';
        case 135, uu = -uv.c1; vv = uv.c2;            g = diag(fliplr(COuv)).';
        otherwise, uu = [];  vv = [];  g = [];
    end
    if ~isempty(g)
        % signed theta of each sample (sign = side of the plane)
        s  = sign(uu(:).'*cosd(cuts(p).phi0) + vv(:).'*sind(cuts(p).phi0));
        ts = s .* asind(sqrt(uu(:).'.^2 + vv(:).'.^2));
        plot(ts, g, 'ko', 'MarkerSize', 4);
    end
    xlim(a.thetaZoomDeg*[-1 1]);
    ylim([cuts(p).m.peakInterp - 50, cuts(p).m.peakInterp + 3]);  grid on;
    xlabel('\theta [deg] (signed)');  ylabel('Co-polar gain [dBi]');
    title(sprintf('%s  \\phi = %d deg:  HPBW tp %.4f / uv %.4f deg', a.name, ...
          cuts(p).phi0, cuts(p).m.hpbw, uvCuts(p).m.hpbw));
    legend({'\theta-\phi file', 'uv file (spline interp.)', 'uv file samples'}, ...
           'Location', 'northwest', 'FontSize', 7);
end
save_fig(fig, fullfile(outDir, sprintf('%s_uv_vs_tp_comparison', a.tag)), cfg);
end


% -------------------------------------------------------------------------
function plot_normalized_comparison(results, cfg)
%PLOT_NORMALIZED_COMPARISON  Both reflectors on a common normalized axis.
% For a circular aperture the far field depends on theta only through
% X = (pi*D/lambda) * sin(theta) (Balanis ch. 12, eq. for circular
% apertures, k*a*sin(theta)). Two reflectors with the same f/D and edge
% taper should therefore give almost the same normalized pattern.
fig = figure('Color', 'w', 'Position', [60 60 1250 520]);
col = [0 0.35 0.8; 0.85 0.2 0.1];
for pl = 1:2
    subplot(1, 2, pl);
    p = find(results(1).planesDeg == 90*(pl-1));    % phi = 0 and phi = 90
    for k = 1:numel(results)
        c = results(k).cutsTp(p);
        X = pi * results(k).ant.D / results(k).theory.lambda * sind(c.theta);
        plot(X, c.co - c.m.peakInterp, '-', 'Color', col(k,:), 'LineWidth', 1.3);
        hold on;
    end
    for k = 1:numel(results)
        c = results(k).cutsTp(p);
        X = pi * results(k).ant.D / results(k).theory.lambda * sind(c.theta);
        plot(X, c.xp - c.m.peakInterp, ':', 'Color', col(k,:), 'LineWidth', 1.0);
    end
    xlim([-20 20]);  ylim([-100 3]);  grid on;
    xlabel('X = (\piD/\lambda) sin\theta');  ylabel('Normalized gain [dB]');
    title(sprintf('Normalized patterns, plane \\phi = %d deg', results(1).planesDeg(p)));
    legend({[results(1).name ' co'], [results(2).name ' co'], ...
            [results(1).name ' X'],  [results(2).name ' X']}, 'Location', 'south');
end
save_fig(fig, fullfile(cfg.outDir, 'comparison_S_vs_Ka_normalized'), cfg);
end


% -------------------------------------------------------------------------
function save_fig(fig, base, cfg)
drawnow;
if cfg.saveFigures
    print(fig, [base '.png'], '-dpng', sprintf('-r%d', cfg.dpi));
end
if cfg.closeFigs; close(fig); end
end


% =========================================================================
%                               REPORT
% =========================================================================
function write_report(R, fid)
%WRITE_REPORT  Print all the performance metrics (fid = 1 -> screen).
pr = @(varargin) fprintf(fid, varargin{:});
pr('\n');
pr('=====================================================================\n');
pr('              REFLECTOR RADIATION PATTERN - PERFORMANCE\n');
pr('=====================================================================\n');
pr(' Co-polar = LHCP (E-l), cross-polar = RHCP (E-r). Gains in dBi.\n');
pr(' "tp" = theta-phi file (9 phi cuts), "uv" = uv-grid file.\n');
pr(' Interpolated values: parabolic fit of the 3 samples around the max.\n');
for k = 1:numel(R)
    r = R(k);  a = r.ant;  th = r.theory;  t = r.tp;  u = r.uv;
    pr('\n---------------------------------------------------------------------\n');
    pr(' %s : D = %g m, F = %g m, f = %.3g GHz, lambda = %.4f m, D/lambda = %.1f\n', ...
       a.name, a.D, a.F, a.freqHz/1e9, th.lambda, a.D/th.lambda);
    pr('---------------------------------------------------------------------\n');
    pr(' THEORETICAL REFERENCES\n');
    pr('   Max directivity, uniform aperture (pi*D/lambda)^2 : %8.2f dBi\n', th.Dmax_dBi);
    pr('   HPBW uniform aperture 58.4*lambda/D              : %8.4f deg\n', th.hpbwUniform_deg);
    pr('   HPBW rule of thumb (tapered) 70*lambda/D         : %8.4f deg\n', th.hpbwRule70_deg);
    pr('   CP beam squint asin(lambda sin(aoff)/(4 pi F))  : %8.4f deg (first-order estimate)\n', th.squint_deg);
    pr('\n PEAK GAIN (co-polar)                         tp file        uv file\n');
    pr('   raw sample                [dBi]       %10.3f     %10.3f\n', t.peakGainRaw_dBi, u.peakGainRaw_dBi);
    pr('   interpolated              [dBi]       %10.3f     %10.3f\n', t.peakGain_dBi, u.peakGain_dBi);
    pr('   aperture efficiency        [%%]        %10.1f     %10.1f\n', ...
       100*10^((t.peakGain_dBi - th.Dmax_dBi)/10), 100*10^((u.peakGain_dBi - th.Dmax_dBi)/10));
    pr('\n PEAK LOCATION\n');
    pr('   tp raw sample   : theta = %.4f deg, phi = %.1f deg\n', t.peakThetaRaw_deg, t.peakPhiRaw_deg);
    pr('   tp interpolated : theta = %.4f deg, phi = %.1f deg\n', t.peakTheta_deg, t.peakPhi_deg);
    pr('   uv raw sample   : (u,v) = (%.5f, %.5f)\n', u.peakURaw, u.peakVRaw);
    pr('   uv interpolated : (u,v) = (%.5f, %.5f) -> theta = %.4f deg, phi = %.1f deg\n', ...
       u.peakU, u.peakV, u.peakTheta_deg, u.peakPhi_deg);
    pr('\n HALF POWER BEAMWIDTH [deg]\n');
    pr('   plane phi [deg] :'); pr('%10d', r.planesDeg); pr('      mean\n');
    pr('   tp file         :'); pr('%10.4f', t.hpbw_deg); pr('%10.4f\n', mean(t.hpbw_deg));
    pr('   uv file         :'); pr('%10.4f', u.hpbw_deg); pr('%10.4f\n', mean(u.hpbw_deg));
    pr('   samples per HPBW: tp %.1f (step %.4f deg), uv %.1f (step %.4f deg)\n', ...
       t.samplesPerHpbw, r.tpStep_deg, u.samplesPerHpbw, u.theta_step_deg);
    if u.samplesPerHpbw < 2
        pr('   WARNING: the uv grid has less than 2 samples per HPBW: the uv\n');
        pr('            peak gain, HPBW and X-pol/XPD values are NOT reliable\n');
        pr('            (the main beam falls between the grid points): use the tp file.\n');
    end
    pr('\n FIRST NULL / FIRST SIDE LOBE (tp file)\n');
    pr('   plane phi [deg] :'); pr('%10d', r.planesDeg); pr('\n');
    pr('   first null [deg]:'); pr('%10.4f', t.null_deg); pr('\n');
    pr('   SLL [dB]        :'); pr('%10.2f', t.sll_dB); pr('\n');
    pr('\n CROSS-POLAR (RHCP)\n');
    pr('   main lobe region: theta <= %.4f deg (mean first null)\n', t.mainLobeTheta_deg);
    pr('   tp main-lobe X-pol peak : %8.2f dBi at theta = %.4f deg, phi = %.1f deg\n', ...
       t.xpMain_dBi, t.xpMainTheta_deg, t.xpMainPhi_deg);
    pr('   uv main-lobe X-pol peak : %8.2f dBi at (u,v) = (%.4f, %.4f) -> theta = %.4f, phi = %.1f deg\n', ...
       u.xpMain_dBi, u.xpMainU, u.xpMainV, u.xpMainTheta_deg, u.xpMainPhi_deg);
    pr('   tp global X-pol peak    : %8.2f dBi at theta = %.4f deg, phi = %.1f deg\n', ...
       t.xpGlobal_dBi, t.xpGlobalTheta_deg, t.xpGlobalPhi_deg);
    pr('   uv global X-pol peak    : %8.2f dBi at (u,v) = (%.4f, %.4f) -> theta = %.4f, phi = %.1f deg\n', ...
       u.xpGlobal_dBi, u.xpGlobalU, u.xpGlobalV, u.xpGlobalTheta_deg, u.xpGlobalPhi_deg);
    pr('   XPD = co peak - X-pol main-lobe peak : tp %.2f dB, uv %.2f dB\n', t.XPD_main_dB, u.XPD_main_dB);
    pr('   worst XPD inside the -3 dB beam      : tp %.2f dB, uv %.2f dB\n', t.XPD_worst3dB_dB, u.XPD_worst3dB_dB);
    pr('   axial ratio at the beam peak (tp)    : %.2e dB\n', t.AR_peak_dB);
end
pr('\n');
end
