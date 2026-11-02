%% REFLECTOR RADIATION PATTERN POST-PROCESSING
%
%  Reads the far-field radiation pattern files of the two ideal offset
%  parabolic reflectors (S-band, D = 30 m and Ka-band, D = 10 m), plots the
%  radiation patterns and extracts the antenna performance metrics.
%
%  HOW TO RUN
%    Press "Run" (F5) to execute everything, or "Run Section"
%    (Ctrl+Enter) to execute one section at a time. All the intermediate
%    data stay in the workspace, so you can explore them.
%
%  STRUCTURE OF THE SCRIPT
%    SECTION 0  User configuration (files, frequencies, geometry, options)
%    SECTION 1  Reading of the pattern files
%    SECTION 2  [A] Theoretical expected values
%    SECTION 3  [B] Performance from the RAW file data (samples only)
%    SECTION 4  [C] Performance from the INTERPOLATED data (cubic spline of the complex field)
%    SECTION 5  Figures
%    SECTION 6  Report (Command Window + text file + .mat file)
%    LOCAL FUNCTIONS (at the end of the file)
%
%  ------------------------------------------------------------------------
%  INPUT FILES (two per reflector, in the folders S_band/ and Ka_band/)
%  ------------------------------------------------------------------------
%  *_uv : pattern on a regular (u,v) grid, u = v = [-0.4, 0.4], 401 points
%  *_tp : pattern on a (theta,phi) grid: theta = [0, 90] deg, 3601 points,
%         phi = [-180, 180] deg, 9 points (one half-cut every 45 deg)
%
%  Both files are in the antenna boresight frame (z axis = boresight), so
%  theta = 0 <=> u = v = 0, with the direction cosines (Balanis ch. 12):
%        u = sin(theta) cos(phi)          v = sin(theta) sin(phi)
%
%  File format (slide 2 of the ppt note):
%    line 1     : title (field components and coordinates)
%    line 2     : N1  N2  = number of values of the 1st and 2nd coordinate
%    next lines : the N1 values of the 1st coordinate (u or theta)
%    then, for every value of the 2nd coordinate (v or phi), 4 blocks:
%       "<2nd coord> = value  E-l amplitude"  -> N1 values [dB]
%       "<2nd coord> = value  E-l phase"      -> N1 values [deg]
%       "<2nd coord> = value  E-r amplitude"  -> N1 values [dB]
%       "<2nd coord> = value  E-r phase"      -> N1 values [deg]
%    E-l = left-hand circular (LHCP), E-r = right-hand circular (RHCP).
%    The feed is RHCP and the reflection inverts the sense of rotation,
%    so the reflector radiates LHCP:
%       CO-POLAR    component = E-l (LHCP)  -> gain G_co (or simply G)
%       CROSS-POLAR component = E-r (RHCP)  -> gain G_x
%    The amplitudes are absolute gains in dBi.
%
%  ------------------------------------------------------------------------
%  SYMBOLS USED IN THE RESULTS AND IN THE FIGURES
%  ------------------------------------------------------------------------
%  G        co-polar (LHCP) gain [dBi]; G_max = its peak value
%  G_x      cross-polar (RHCP) gain [dBi]
%  G0       maximum gain of a uniformly illuminated aperture (pi D/lambda)^2
%  eta      aperture efficiency  eta = G_max / G0
%  HPBW     half power beamwidth (-3 dB from G_max)
%  FML      first minor lobe: the first lobe after the first null of the
%           main beam (Balanis sec. 2.2: minor lobes = all lobes except the
%           main one)
%  FMLL     first minor lobe level = G(FML) - G_max [dB]
%  XPD      cross-polar discrimination = G_max - G_x,max [dB]
%  AR       axial ratio of the polarization ellipse, from the two circular
%           components:  AR = (|E_co| + |E_x|) / (|E_co| - |E_x|)
%           (0 dB = perfect circular polarization)
%
%  The code uses only base MATLAB functions (no toolbox needed).

clear; close all; clc;

%% ========================================================================
%  SECTION 0 - USER CONFIGURATION
%  ========================================================================
% Folder of this script: the data folders and the results are relative to it
rootDir = fileparts(mfilename('fullpath'));
if isempty(rootDir); rootDir = pwd; end   % (when running single sections)

cfg.outDir      = fullfile(rootDir, 'results'); % output folder
cfg.saveFigures = true;      % save every figure as PNG
cfg.dpi         = 150;       % resolution of the PNG files
cfg.plotFullUv  = true;      % contour of the whole uv grid (slow in Octave)
cfg.coComp      = 'El';      % co-polar component    (E-l = LHCP)
cfg.xpComp      = 'Er';      % cross-polar component (E-r = RHCP)
cfg.planesDeg   = [0 45 90 135];  % phi planes of the 1-D cuts [deg]
cfg.c0          = 299792458; % speed of light [m/s]
cfg.upTp        = 20;        % interpolation: theta step of the cuts / upTp
cfg.upUv        = 20;        % interpolation: uv grid step / upUv
cfg.minLobeRise = 1;         % [dB] a minimum followed by a rise smaller
                             % than this is a ripple, not a null
cfg.etaTaper    = 0.82;      % expected aperture efficiency for a ~-11 dB
                             % edge taper (optimum, Balanis sec. 15.4)

% Description of the two antennas (ppt note + supervisor e-mail)
%   freqHz       : operating frequency [Hz]
%   D, F         : projected aperture diameter and focal length [m]
%   offsetDeg    : offset angle of the feed (aoff in the ppt) [deg]
%   thetaZoomDeg : half-width of the zoomed theta cuts and of the
%                  interpolation window of the cuts [deg]
%   uvZoom       : half-width of the zoomed uv contour plots [-]
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

nAnt = numel(ant);
if ~exist(cfg.outDir, 'dir'); mkdir(cfg.outDir); end

%% ========================================================================
%  SECTION 1 - READING OF THE PATTERN FILES
%  ========================================================================
% For each antenna the struct dat(k) contains:
%   .uv, .tp       the raw content of the two files (see read_pattern_file)
%   .U, .V         uv grid (meshgrid), rows <-> v, columns <-> u
%   .Guv, .Gxuv    co- and cross-polar gain on the uv grid [dBi]
%   .Euv, .Exuv    co- and cross-polar complex field on the uv grid,
%                  E = 10^(G/20) exp(j phase)  (used for the interpolation)
%   .du            uv grid step,  .dth  theta step of the tp file [deg]
%   .cutsTp        full plane cuts (signed theta in [-90,90]) from tp file
%   .cutsUv        the uv samples lying on the same planes (axes and
%                  diagonals of the uv grid)
dat = struct([]);
for k = 1:nAnt
    fprintf('Reading %s files ...\n', ant(k).name);
    d = struct();
    d.uv   = read_pattern_file(fullfile(ant(k).folder, ant(k).uvFile));
    d.tp   = read_pattern_file(fullfile(ant(k).folder, ant(k).tpFile));
    [d.U, d.V] = meshgrid(d.uv.c1, d.uv.c2);
    d.Guv  = d.uv.amp.(cfg.coComp);
    d.Gxuv = d.uv.amp.(cfg.xpComp);
    d.Euv  = complex_field(d.Guv,  d.uv.ph.(cfg.coComp));
    d.Exuv = complex_field(d.Gxuv, d.uv.ph.(cfg.xpComp));
    d.du   = d.uv.c1(2) - d.uv.c1(1);
    d.dth  = d.tp.c1(2) - d.tp.c1(1);
    d.cutsTp = build_tp_plane_cuts(d.tp, cfg);
    d.cutsUv = build_uv_plane_cuts(d.uv, cfg);
    dat = add_to_array(dat, d, k);
end

%% ========================================================================
%  SECTION 2 - [A] THEORETICAL EXPECTED VALUES
%  ========================================================================
% Reference: uniformly illuminated circular aperture of diameter D
% (Balanis ch. 12), whose far field is the Airy pattern
%       E(theta) ~ 2 J1(X) / X ,   X = (pi D / lambda) sin(theta)
%   -3 dB points        : X = 1.6163
%   first null          : X = 3.8317  (first zero of J1)
%   first minor lobe    : X = 5.1356, level -17.57 dB
% A reflector has a tapered illumination (here ~ -11 dB at the edge): the
% beam is wider (rule of thumb HPBW ~ 70 lambda/D deg), the minor lobes are
% lower and the efficiency is below 100 % (about 82 % at the optimum edge
% taper of ~ -11 dB, Balanis sec. 15.4).
th = struct([]);
for k = 1:nAnt
    th = add_to_array(th, theory_values(ant(k), cfg), k);
end

%% ========================================================================
%  SECTION 3 - [B] PERFORMANCE FROM THE RAW FILE DATA
%  ========================================================================
% Only the samples written in the files are used: peaks, nulls and lobes
% are file samples. The only exception is the -3 dB point, which in general
% lies between two samples: it is found by a straight line between the two
% samples around -3 dB (the minimum needed to define a beamwidth).
raw = struct([]);
for k = 1:nAnt
    r = struct();
    % theta-phi file: all its samples are used for the 2-D searches
    r.tp = pattern_metrics(dat(k).cutsTp, tp_points(dat(k).tp, cfg), ...
                           th(k).G0, cfg);
    % uv file: all the grid samples are used for the 2-D searches
    r.uv = pattern_metrics(dat(k).cutsUv, ...
                           grid_points(dat(k).U, dat(k).V, dat(k).Guv, dat(k).Gxuv), ...
                           th(k).G0, cfg);
    r.stepTp = dat(k).dth;               % angular sampling steps [deg]
    r.stepUv = asind(dat(k).du);
    raw = add_to_array(raw, r, k);
end

%% ========================================================================
%  SECTION 4 - [C] PERFORMANCE FROM THE INTERPOLATED DATA
%  ========================================================================
% WHAT IS INTERPOLATED: the complex far field E = 10^(G/20) exp(j phase),
% built from the amplitude AND the phase written in the files; the gain is
% then G = 20 log10|E|. The far field of a finite aperture is a smooth
% (band-limited) function of the direction, so its real and imaginary
% parts are very well reproduced by cubic splines. Interpolating G in dB
% instead is not recommended: at the nulls G goes to -infinity and the
% spline oscillates (e.g. a cross-pol null of -92 dB between samples at
% -10 dB produces false cross-pol peaks).
% HOW:
%   * theta-phi file: every plane cut is resampled with interp1(...,
%     'spline') on a theta step cfg.upTp times finer, inside |theta| <=
%     thetaZoomDeg (main beam + first minor lobes). In a plane cut the two
%     half cuts phi0 and phi0-180 join continuously at theta = 0 (checked
%     on the phases of the files);
%   * uv file: the central part of the grid is resampled with
%     interp2(..., 'spline') on a grid cfg.upUv times finer; the plane
%     cuts are then taken from this fine grid.
% The SAME metric functions of section 3 are then applied to the fine data.
ip = struct([]);
for k = 1:nAnt
    a = ant(k);
    % --- theta-phi file: spline resampling of every plane cut ---
    stepTp = dat(k).dth / cfg.upTp;
    fineTp = struct([]);
    for p = 1:numel(dat(k).cutsTp)
        fineTp = add_to_array(fineTp, ...
            resample_cut(dat(k).cutsTp(p), a.thetaZoomDeg, stepTp), p);
    end
    % --- uv file: spline resampling of the central part of the grid ---
    W = max(a.uvZoom, sind(a.thetaZoomDeg));       % half-width of window
    fg = resample_uv_grid(dat(k).uv, dat(k).Euv, dat(k).Exuv, W, cfg.upUv);
    fineUv = struct([]);
    thMax  = asind(W);
    thf    = linspace(-thMax, thMax, round(2*thMax/stepTp) + 1);
    for p = 1:numel(cfg.planesDeg)
        c = struct();
        c.phi0  = cfg.planesDeg(p);
        c.theta = thf;
        uu = sind(thf) * cosd(c.phi0);
        vv = sind(thf) * sind(c.phi0);
        % linear interpolation of the complex field is enough here: the
        % grid is already fine
        c.co = to_dB(interp_c2(fg.U, fg.V, fg.E,  uu, vv, 'linear'));
        c.xp = to_dB(interp_c2(fg.U, fg.V, fg.Ex, uu, vv, 'linear'));
        fineUv = add_to_array(fineUv, c, p);
    end
    dat(k).fineTp = fineTp;   % kept for the figures
    dat(k).fineUv = fineUv;
    dat(k).fineGrid = fg;

    r = struct();
    r.tp = pattern_metrics(fineTp, cut_points(fineTp), th(k).G0, cfg);
    r.uv = pattern_metrics(fineUv, grid_points(fg.U, fg.V, fg.G, fg.Gx), ...
                           th(k).G0, cfg);
    % sampling check: number of file samples inside one HPBW
    r.samplesPerHpbwTp = mean(r.tp.hpbw) / dat(k).dth;
    r.samplesPerHpbwUv = mean(r.tp.hpbw) / asind(dat(k).du);
    ip = add_to_array(ip, r, k);
end

%% ========================================================================
%  SECTION 5 - FIGURES
%  ========================================================================
for k = 1:nAnt
    outDir = fullfile(cfg.outDir, ant(k).tag);
    if ~exist(outDir, 'dir'); mkdir(outDir); end
    if cfg.plotFullUv
        plot_uv_full(ant(k), dat(k), raw(k), cfg, outDir);
    end
    plot_uv_zoom(ant(k), dat(k), raw(k), ip(k), cfg, outDir);
    plot_cuts_full(ant(k), dat(k), raw(k), cfg, outDir);
    plot_cuts_zoom(ant(k), dat(k), ip(k), cfg, outDir);
    plot_cuts_overlay(ant(k), dat(k), ip(k), cfg, outDir);
    plot_axial_ratio(ant(k), dat(k), ip(k), cfg, outDir);
    plot_uv_vs_tp(ant(k), dat(k), ip(k), cfg, outDir);
end
plot_normalized_comparison(ant, dat, th, ip, cfg);

%% ========================================================================
%  SECTION 6 - REPORT
%  ========================================================================
write_report(1, ant, th, raw, ip, cfg);                        % screen
fid = fopen(fullfile(cfg.outDir, 'performance_summary.txt'), 'w');
write_report(fid, ant, th, raw, ip, cfg);                      % text file
fclose(fid);
results.ant = ant;  results.theory = th;  results.raw = raw;  results.interp = ip;
save(fullfile(cfg.outDir, 'performance_results.mat'), 'results', 'cfg');
fprintf('Figures and report saved in: %s\n', cfg.outDir);


%% ########################################################################
%                            LOCAL FUNCTIONS
%  ########################################################################

% =========================================================================
%  F1 - READING
% =========================================================================
function P = read_pattern_file(fname)
%READ_PATTERN_FILE  Generic reader of the near/far-field output files.
%  P.title          first line of the file
%  P.c1  (N1x1)     values of the 1st coordinate (u or theta)
%  P.c2  (N2x1)     values of the 2nd coordinate (v or phi)
%  P.amp.El, .Er    (N2 x N1) amplitudes [dB]
%  P.ph.El,  .Er    (N2 x N1) phases [deg]
%  Row i <-> c2(i), column j <-> c1(j).
if ~exist(fname, 'file')
    error('File not found: %s', fname);
end
txt = fileread(fname);
L   = regexp(txt, '\r?\n', 'split');
L   = L(~cellfun(@isempty, strtrim(L)));              % drop empty lines

% header lines contain letters (other than the 'E' of the exponent)
isHdr  = ~cellfun(@isempty, regexp(L, '[a-df-zA-DF-Z]', 'once'));
hdrIdx = find(isHdr);

P.title = strtrim(L{1});
N  = sscanf(L{2}, '%d');                              % [N1; N2]
N1 = N(1);  N2 = N(2);

% list of the 1st coordinate: from line 3 up to the first block header
P.c1 = parse_numbers(L(3:hdrIdx(2)-1));
if numel(P.c1) ~= N1
    error('%s: expected %d values of the 1st coordinate, found %d', ...
          fname, N1, numel(P.c1));
end

% data blocks
blk    = hdrIdx(2:end);
blkEnd = [blk(2:end) - 1, numel(L)];
P.c2 = nan(N2, 1);
P.amp.El = nan(N2, N1);  P.amp.Er = nan(N2, N1);
P.ph.El  = nan(N2, N1);  P.ph.Er  = nan(N2, N1);
i2 = 0;  lastVal = NaN;
for b = 1:numel(blk)
    h = L{blk(b)};
    % "phi =  -180.0000  degrees, E-l amplitude"  or  "v = -0.4  E-r phase"
    tok = regexp(h, '^\s*(\w+)\s*=\s*([-+0-9.EeDd]+)', 'tokens', 'once');
    val = str2double(tok{2});
    if b == 1 || val ~= lastVal          % new value of the 2nd coordinate
        i2 = i2 + 1;
        P.c2(i2) = val;
        lastVal  = val;
    end
    if ~isempty(strfind(h, 'E-l')); comp = 'El'; else; comp = 'Er'; end
    data = parse_numbers(L(blk(b)+1:blkEnd(b)));
    if numel(data) ~= N1
        error('%s: block "%s" has %d values instead of %d', ...
              fname, strtrim(h), numel(data), N1);
    end
    if ~isempty(strfind(h, 'amplitude'))
        P.amp.(comp)(i2, :) = data;
    else
        P.ph.(comp)(i2, :)  = data;
    end
end
if i2 ~= N2
    error('%s: expected %d values of the 2nd coordinate, found %d', ...
          fname, N2, i2);
end
end

function x = parse_numbers(lineCell)
%PARSE_NUMBERS  Text lines -> column of numbers.
% The files use fixed-width Fortran fields, so two numbers can be glued
% together, e.g. "0.8900023E-01-0.5656819". A blank is inserted before
% every '-' that follows a digit or a dot (the '-' of an exponent follows
% an 'E' and is left untouched).
s = sprintf('%s ', lineCell{:});
s = regexprep(s, '(?<=[0-9.])-', ' -');
x = sscanf(s, '%f');
end

% =========================================================================
%  F2 - DATA ORGANISATION (cuts and lists of points)
% =========================================================================
function cuts = build_tp_plane_cuts(tp, cfg)
%BUILD_TP_PLANE_CUTS  The tp file has half cuts theta = 0..90 at
%phi = -180:45:180. The half cuts phi0 and phi0-180 are joined into a full
%plane cut with a signed theta in [-90, 90] (theta > 0 on the phi0 side).
t = tp.c1(:).';
cuts = struct([]);
for p = 1:numel(cfg.planesDeg)
    phi0 = cfg.planesDeg(p);
    iF = find_phi(tp.c2, phi0);           % forward half  (theta > 0)
    iB = find_phi(tp.c2, phi0 - 180);     % backward half (theta < 0)
    c = struct();
    c.phi0  = phi0;
    c.theta = [-fliplr(t(2:end)), t];
    c.co = [fliplr(tp.amp.(cfg.coComp)(iB, 2:end)), tp.amp.(cfg.coComp)(iF, :)];
    c.xp = [fliplr(tp.amp.(cfg.xpComp)(iB, 2:end)), tp.amp.(cfg.xpComp)(iF, :)];
    cph  = [fliplr(tp.ph.(cfg.coComp)(iB, 2:end)), tp.ph.(cfg.coComp)(iF, :)];
    xph  = [fliplr(tp.ph.(cfg.xpComp)(iB, 2:end)), tp.ph.(cfg.xpComp)(iF, :)];
    c.E  = complex_field(c.co, cph);      % complex fields along the cut
    c.Ex = complex_field(c.xp, xph);
    cuts = add_to_array(cuts, c, p);
end
end

function i = find_phi(phiList, phi)
% index of the angle phi in phiList (360 deg periodicity)
dd = mod(phiList(:) - phi + 180, 360) - 180;
[dmin, i] = min(abs(dd));
if dmin > 1e-3
    error('phi = %g deg not present in the theta-phi file', phi);
end
end

function cuts = build_uv_plane_cuts(uv, cfg)
%BUILD_UV_PLANE_CUTS  uv samples lying exactly on the phi planes:
%u axis (phi = 0), v axis (phi = 90) and the two diagonals (45, 135).
%The grid must be square and symmetric (true for these files).
n  = numel(uv.c1);  ic = (n + 1) / 2;
G  = uv.amp.(cfg.coComp);  Gx = uv.amp.(cfg.xpComp);
sym = numel(uv.c2) == n && max(abs(uv.c1(:) + flipud(uv.c1(:)))) < 1e-6 ...
      && max(abs(uv.c1(:) - uv.c2(:))) < 1e-6;
if ~sym; error('build_uv_plane_cuts: the uv grid is not square/symmetric'); end
t = uv.c1(:).';                               % abscissa along the line
cuts = struct([]);
for p = 1:numel(cfg.planesDeg)
    c = struct();
    c.phi0 = cfg.planesDeg(p);
    switch c.phi0
        case 0,   r = abs(t);         co = G(ic, :);      xp = Gx(ic, :);
        case 90,  r = abs(t);         co = G(:, ic).';    xp = Gx(:, ic).';
        case 45,  r = sqrt(2)*abs(t); co = diag(G).';     xp = diag(Gx).';
        case 135, r = sqrt(2)*abs(t); co = diag(fliplr(G)).';
                  xp = diag(fliplr(Gx)).';
        otherwise
            error('build_uv_plane_cuts: plane %g deg not on the grid', c.phi0);
    end
    % on the 135 deg diagonal the point (u,v) = (-t, t) is at signed
    % distance +t along the direction (cos135, sin135): same sign as t
    c.theta = sign(t) .* asind(min(r, 1));
    c.co = co;  c.xp = xp;
    cuts = add_to_array(cuts, c, p);
end
end

function P = tp_points(tp, cfg)
%TP_POINTS  All the samples of the tp file as a list of points.
[PH, TH] = ndgrid(tp.c2, tp.c1);
P.th = TH(:);  P.ph = PH(:);
P.co = reshape(tp.amp.(cfg.coComp), [], 1);
P.xp = reshape(tp.amp.(cfg.xpComp), [], 1);
end

function P = grid_points(U, V, G, Gx)
%GRID_POINTS  All the points of a uv grid as a list of points.
P.th = asind(min(sqrt(U(:).^2 + V(:).^2), 1));
P.ph = atan2d(V(:), U(:));
P.co = G(:);  P.xp = Gx(:);
end

function P = cut_points(cuts)
%CUT_POINTS  All the points of a set of plane cuts as a list of points.
P.th = []; P.ph = []; P.co = []; P.xp = [];
for p = 1:numel(cuts)
    t  = cuts(p).theta(:);
    ph = cuts(p).phi0 * ones(size(t));
    ph(t < 0) = cuts(p).phi0 - 180;
    P.th = [P.th; abs(t)];
    P.ph = [P.ph; mod(ph + 180, 360) - 180];
    P.co = [P.co; cuts(p).co(:)];
    P.xp = [P.xp; cuts(p).xp(:)];
end
end

% =========================================================================
%  F3 - THEORY
% =========================================================================
function t = theory_values(a, cfg)
%THEORY_VALUES  Expected values for an aperture of diameter D.
t.lambda   = cfg.c0 / a.freqHz;               % wavelength [m]
t.DoL      = a.D / t.lambda;                  % electrical size
t.G0       = 10*log10((pi*a.D/t.lambda)^2);   % uniform aperture, eta = 1
t.etaTaper = cfg.etaTaper;                    % optimum-taper reflector
t.Gtaper   = t.G0 + 10*log10(cfg.etaTaper);
k = pi * a.D / t.lambda;                      % X = k sin(theta)
t.hpbwUniform  = 2 * asind(1.6163 / k);       % Airy -3 dB points
t.hpbwTaper    = 70 * t.lambda / a.D;         % rule of thumb [deg]
t.nullUniform  = asind(3.8317 / k);           % first zero of J1
t.fmlUniform   = asind(5.1356 / k);           % first minor lobe position
t.fmllUniform  = -17.57;                      % first minor lobe level [dB]
% Beam squint of an offset paraboloid with a circularly polarised feed
% (Duan & Rahmat-Samii, IEEE TAP 1991), in the plane orthogonal to the
% offset plane (phi = +-90 deg):  sin(psi) = lambda sin(aoff) / (4 pi F)
t.squint = asind(t.lambda * sind(a.offsetDeg) / (4*pi*a.F));
end

% =========================================================================
%  F4 - INTERPOLATION
% =========================================================================
function f = resample_cut(c, W, step)
%RESAMPLE_CUT  Cubic-spline resampling of a plane cut in |theta| <= W.
dth = abs(c.theta(2) - c.theta(1));
sel = abs(c.theta) <= W + 4*dth;              % samples used by the spline
f = struct();
f.phi0  = c.phi0;
f.theta = linspace(-W, W, round(2*W/step) + 1);
E  = interp_c1(c.theta(sel), c.E(sel),  f.theta);   % complex spline
Ex = interp_c1(c.theta(sel), c.Ex(sel), f.theta);
f.co = to_dB(E);
f.xp = to_dB(Ex);
end

function fg = resample_uv_grid(uv, E, Ex, W, up)
%RESAMPLE_UV_GRID  2-D cubic-spline resampling of the complex uv fields
%inside |u|,|v| <= W, on a grid "up" times finer.
du = uv.c1(2) - uv.c1(1);
iu = find(abs(uv.c1) <= W + 3*du);
iv = find(abs(uv.c2) <= W + 3*du);
uf = linspace(uv.c1(iu(1)), uv.c1(iu(end)), up*(numel(iu)-1) + 1);
vf = linspace(uv.c2(iv(1)), uv.c2(iv(end)), up*(numel(iv)-1) + 1);
[fg.U, fg.V] = meshgrid(uf, vf);
fg.E  = interp_c2(uv.c1(iu).', uv.c2(iv), E(iv, iu),  fg.U, fg.V, 'spline');
fg.Ex = interp_c2(uv.c1(iu).', uv.c2(iv), Ex(iv, iu), fg.U, fg.V, 'spline');
fg.G  = to_dB(fg.E);                          % gains [dBi]
fg.Gx = to_dB(fg.Ex);
end

function E = complex_field(GdB, phDeg)
% complex field from gain [dB] and phase [deg]: |E|^2 = gain (linear)
E = 10.^(GdB/20) .* exp(1i * phDeg * pi/180);
end

function G = to_dB(E)
% gain [dBi] from the complex field (floor at -200 dB to avoid -Inf)
G = 20*log10(max(abs(E), 1e-10));
end

function Ei = interp_c1(x, E, xi)
% 1-D cubic spline of a complex function (real and imaginary parts)
Ei = interp1(x, real(E), xi, 'spline') + 1i * interp1(x, imag(E), xi, 'spline');
end

function Ei = interp_c2(X, Y, E, Xi, Yi, method)
% 2-D interpolation of a complex function (real and imaginary parts)
Ei = interp2(X, Y, real(E), Xi, Yi, method) + 1i * interp2(X, Y, imag(E), Xi, Yi, method);
end

% =========================================================================
%  F5 - PERFORMANCE METRICS (same functions for raw and interpolated data)
% =========================================================================
function M = pattern_metrics(cuts, P, G0, cfg)
%PATTERN_METRICS  Performance of one pattern.
%  cuts : plane cuts (signed theta, co, xp)  -> HPBW, nulls, FML per plane
%  P    : list of all the available points (th, ph, co, xp) -> peak,
%         cross-polar peak, XPD and AR searched in 2-D
%  G0   : maximum gain of the uniform aperture -> eta

% ---- per-plane metrics ----
pl = struct([]);
for p = 1:numel(cuts)
    m = cut_metrics(cuts(p).theta, cuts(p).co, cuts(p).xp, cfg.minLobeRise);
    m.phi0 = cuts(p).phi0;
    pl = add_to_array(pl, m, p);
end
M.planes = pl;
M.hpbw   = [pl.hpbw];
M.null   = [pl.firstNull];
M.fmll   = [pl.fmll];

% ---- co-polar peak (2-D search over all points) ----
[M.G, i0]  = max(P.co);
M.th  = P.th(i0);  M.ph = P.ph(i0);
M.u   = sind(M.th) * cosd(M.ph);
M.v   = sind(M.th) * sind(M.ph);
M.eta = 10^((M.G - G0)/10);

% ---- worst first minor lobe among the planes ----
[~, q] = max(M.fmll);
M.fmlWorst.plane = pl(q).phi0;
M.fmlWorst.fmll  = pl(q).fmll;
M.fmlWorst.G     = pl(q).fmlG;
M.fmlWorst.theta = pl(q).fmlPos;

% ---- cross-polarization ----
% main-beam region: theta <= thNull (mean first-null angle of the planes)
M.thNull = mean(M.null);
inMain = P.th <= M.thNull;
xpm = P.xp;  xpm(~inMain) = -Inf;
[M.Gx, j] = max(xpm);                     % max G_x inside the main beam
M.GxTh = P.th(j);  M.GxPh = P.ph(j);
M.XPD  = M.G - M.Gx;
[M.GxAll, j] = max(P.xp);                 % max G_x over all the points
M.GxAllTh = P.th(j);  M.GxAllPh = P.ph(j);
% inside the -3 dB beam
in3 = P.co >= M.G - 3;
M.n3dB = nnz(in3);                        % number of points in the beam
M.XPDworst3 = min(P.co(in3) - P.xp(in3));
AR = axial_ratio_dB(P.co, P.xp);
M.ARpeak = AR(i0);
M.ARmax3 = max(AR(in3));
end

function m = cut_metrics(x, co, xp, minRise)
%CUT_METRICS  Peak, HPBW, first nulls, first minor lobes and cross-polar
%peak of one plane cut (x = signed theta [deg], co/xp in dBi). It works on
%the given samples; with interpolated (fine) data it gives the
%interpolated values.
x = x(:).';  y = co(:).';  yx = xp(:).';
n = numel(y);
[m.G, i0] = max(y);                       % peak of the cut
m.pkPos = x(i0);

% --- HPBW: -3 dB crossings, straight line between the 2 samples around ---
lev = m.G - 3;
iR = find(y(i0:end) < lev, 1, 'first') + i0 - 1;
iL = find(y(1:i0)   < lev, 1, 'last');
if isempty(iR) || isempty(iL)
    m.x3L = NaN;  m.x3R = NaN;
else
    m.x3R = x(iR-1) + (lev - y(iR-1)) * (x(iR) - x(iR-1)) / (y(iR) - y(iR-1));
    m.x3L = x(iL)   + (lev - y(iL))   * (x(iL+1) - x(iL)) / (y(iL+1) - y(iL));
end
m.hpbw = m.x3R - m.x3L;

% --- first nulls and first minor lobes, on both sides of the peak ---
[iNR, iSR] = null_and_lobe(y, i0, +1, minRise);
[iNL, iSL] = null_and_lobe(y, i0, -1, minRise);
m.nullL = x(iNL);  m.nullR = x(iNR);
m.firstNull = 0.5 * (m.nullR - m.nullL);  % half of null-to-null width
m.fmlPosL = x(iSL);  m.fmlGL = y(iSL);
m.fmlPosR = x(iSR);  m.fmlGR = y(iSR);
if m.fmlGR >= m.fmlGL                     % the highest of the two FMLs
    m.fmlPos = m.fmlPosR;  m.fmlG = m.fmlGR;
else
    m.fmlPos = m.fmlPosL;  m.fmlG = m.fmlGL;
end
m.fmll = m.fmlG - m.G;                    % first minor lobe level [dB]

% --- cross-polar peak inside the main lobe of this cut ---
[m.Gx, j] = max(yx(iNL:iNR));
m.GxPos = x(iNL + j - 1);
m.ARpeak = axial_ratio_dB(y(i0), yx(i0));
m.n = n;
end

function [iN, iS] = null_and_lobe(y, i0, dirn, minRise)
%NULL_AND_LOBE  Walk from the peak index i0 in direction dirn (+1/-1):
%down to the first null (local minimum), then up to the first minor lobe
%(local maximum). A minimum followed by a rise smaller than minRise dB is
%a ripple of the main beam and is skipped.
n = numel(y);
ok = @(i) i >= 1 && i <= n;
i = i0;
while true
    while ok(i+dirn) && y(i+dirn) <= y(i); i = i + dirn; end   % descend
    iN = i;
    j = i;
    while ok(j+dirn) && y(j+dirn) >= y(j); j = j + dirn; end   % ascend
    iS = j;
    if y(iS) - y(iN) >= minRise || ~ok(j+dirn); break; end
    i = j;                                % ripple: keep walking
end
end

function ar = axial_ratio_dB(coDb, xpDb)
%AXIAL_RATIO_DB  AR from the two circular components (Balanis sec. 2.12):
%  AR = (|E_co| + |E_x|) / (|E_co| - |E_x|),  in dB: 20 log10(AR)
eco = 10.^(coDb/20);  ex = 10.^(xpDb/20);
ar = 20*log10((eco + ex) ./ abs(eco - ex));
end

% =========================================================================
%  F6 - FIGURES
% =========================================================================
function plot_uv_full(a, d, r, cfg, outDir)
%PLOT_UV_FULL  Contours of G and G_x over the whole uv grid (raw data).
[GL, GxL] = contour_levels(r.uv.G, r.uv.GxAll);
fig = figure('Color', 'w', 'Position', [60 60 1300 560]);
lim = max(abs(d.uv.c1));
subplot(1, 2, 1);
uv_contour(d.U, d.V, d.Guv, GL, lim);  hold on;
plot(r.uv.u, r.uv.v, 'k+', 'MarkerSize', 10, 'LineWidth', 2);
title(sprintf('G (co-pol LHCP) [dBi], raw data\nG_{max} = %.2f dBi', r.uv.G));
subplot(1, 2, 2);
uv_contour(d.U, d.V, d.Gxuv, GxL, lim);  hold on;
title(sprintf('G_x (cross-pol RHCP) [dBi], raw data\nmax over the grid G_x = %.2f dBi', ...
      r.uv.GxAll));
fig_title(sprintf('%s reflector - uv plane, whole grid', a.name));
save_fig(fig, fullfile(outDir, [a.tag '_uv_contour_full']), cfg);
end

function plot_uv_zoom(a, d, r, ip, cfg, outDir)
%PLOT_UV_ZOOM  Zoom on the main beam: raw data (left) and spline
%interpolation (right); co-polar (top) and cross-polar (bottom).
% same colour scale for raw and interpolated panels
[GL, GxL] = contour_levels(max(r.uv.G, ip.uv.G), max(r.uv.Gx, ip.uv.Gx));
lim = a.uvZoom;  fg = d.fineGrid;
fig = figure('Color', 'w', 'Position', [40 40 1250 1050]);
src = {'raw samples', 'spline interpolation'};
for col = 1:2
    if col == 1
        U = d.U;  V = d.V;  G = d.Guv;  Gx = d.Gxuv;  M = r.uv;
    else
        U = fg.U; V = fg.V; G = fg.G;   Gx = fg.Gx;   M = ip.uv;
    end
    % ---- co-polar ----
    subplot(2, 2, col);
    uv_contour(U, V, G, GL, lim);  hold on;
    contour(U, V, G, (M.G-3)*[1 1], 'w-', 'LineWidth', 1.5);    % -3 dB
    draw_circle(sind(M.thNull), 'k--');                          % 1st null
    if col == 1; plot(d.U(:), d.V(:), '.', 'Color', [.3 .3 .3], 'MarkerSize', 4); end
    plot(M.u, M.v, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    title(sprintf('G (co-pol LHCP) [dBi] - %s\nG_{max} = %.2f dBi at \\theta = %.4f\\circ, \\phi = %.1f\\circ', ...
          src{col}, M.G, M.th, M.ph));
    % ---- cross-polar ----
    subplot(2, 2, 2 + col);
    uv_contour(U, V, Gx, GxL, lim);  hold on;
    contour(U, V, G, (M.G-3)*[1 1], 'w-', 'LineWidth', 1.5);
    draw_circle(sind(M.thNull), 'k--');
    if col == 1; plot(d.U(:), d.V(:), '.', 'Color', [.3 .3 .3], 'MarkerSize', 4); end
    xu = sind(M.GxTh) * cosd(M.GxPh);  xv = sind(M.GxTh) * sind(M.GxPh);
    plot(xu, xv, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
    title(sprintf('G_x (cross-pol RHCP) [dBi] - %s\nG_{x,max} = %.2f dBi at \\theta = %.4f\\circ, \\phi = %.1f\\circ', ...
          src{col}, M.Gx, M.GxTh, M.GxPh));
end
fig_title(sprintf('%s reflector - uv plane, zoom on the main beam', a.name));
fig_note(['+ : G_{max}     x : G_{x,max} (max cross-pol inside the main beam)     ' ...
          'white line : -3 dB contour of G     dashed circle : first null of G ' ...
          '(main-beam region)     grey dots : file samples']);
save_fig(fig, fullfile(outDir, [a.tag '_uv_contour_zoom']), cfg);
end

function [GL, GxL] = contour_levels(Gpk, Gxpk)
% iso-gain levels every 3 dB: 60 dB below the co-pol peak, 45 dB below
% the cross-pol peak
GL  = floor(Gpk)  - (60:-3:0);
GxL = floor(Gxpk) - (45:-3:0);
end

function uv_contour(U, V, G, lev, lim)
% filled contour of G on the uv plane; values below the lowest level are
% clipped to it (no white holes, faster contourf)
ii = find(abs(U(1, :)) <= lim * 1.05);
jj = find(abs(V(:, 1)) <= lim * 1.05);
contourf(U(jj, ii), V(jj, ii), max(G(jj, ii), min(lev)), lev, 'LineColor', 'none');
axis equal;  axis([-lim lim -lim lim]);
caxis([min(lev) max(lev)]);
colormap(jet(numel(lev) - 1));
colorbar;  grid on;  box on;  set(gca, 'Layer', 'top');
xlabel('u = sin\theta cos\phi');  ylabel('v = sin\theta sin\phi');
end

function draw_circle(rad, style)
t = linspace(0, 2*pi, 361);
plot(rad*cos(t), rad*sin(t), style, 'LineWidth', 1);
end

% -------------------------------------------------------------------------
function plot_cuts_full(a, d, r, cfg, outDir)
%PLOT_CUTS_FULL  G and G_x in the 4 planes, theta from -90 to 90 deg.
fig = figure('Color', 'w', 'Position', [60 60 1250 820]);
for p = 1:numel(d.cutsTp)
    c = d.cutsTp(p);
    subplot(2, 2, p);
    plot(c.theta, c.co, 'b-', 'LineWidth', 1.1);  hold on;
    plot(c.theta, c.xp, 'r-', 'LineWidth', 0.9);
    xlim([-90 90]);  ylim([r.tp.G - 100, r.tp.G + 5]);  grid on;
    xlabel(theta_label(c.phi0));  ylabel('Gain [dBi]');
    title(sprintf('plane \\phi = %d\\circ', c.phi0));
    legend({'G (co-pol LHCP)', 'G_x (cross-pol RHCP)'}, 'Location', 'northeast', 'FontSize', 7);
end
fig_title(sprintf('%s reflector - \\theta cuts from the \\theta-\\phi file (raw data)', a.name));
save_fig(fig, fullfile(outDir, [a.tag '_theta_cuts_full']), cfg);
end

function s = theta_label(phi0)
s = sprintf('\\theta [deg]   (\\theta < 0 : \\phi = %d\\circ)', mod(phi0, 360) - 180);
end

% -------------------------------------------------------------------------
function plot_cuts_zoom(a, d, ip, cfg, outDir)
%PLOT_CUTS_ZOOM  Zoom on the main beam of every plane: raw samples (dots),
%spline interpolation (lines) and all the metrics marked on the curves.
fig = figure('Color', 'w', 'Position', [40 40 1300 900]);
lim = a.thetaZoomDeg;
for p = 1:numel(d.cutsTp)
    c = d.cutsTp(p);  f = d.fineTp(p);  m = ip.tp.planes(p);
    subplot(2, 2, p);
    h1 = plot(f.theta, f.co, 'b-', 'LineWidth', 1.2);  hold on;
    h2 = plot(c.theta, c.co, 'b.', 'MarkerSize', 7);
    h3 = plot(f.theta, f.xp, 'r-', 'LineWidth', 1.0);
    h4 = plot(c.theta, c.xp, 'r.', 'MarkerSize', 7);
    plot([-lim lim], (m.G - 3)*[1 1], 'k:');                       % -3 dB
    h5 = plot(m.pkPos, m.G, 'k^', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    h6 = plot([m.x3L m.x3R], (m.G - 3)*[1 1], 'ko', 'MarkerFaceColor', 'y');
    h7 = plot([m.nullL m.nullR], interp1(f.theta, f.co, [m.nullL m.nullR]), ...
              'kv', 'MarkerFaceColor', 'c');
    h8 = plot([m.fmlPosL m.fmlPosR], [m.fmlGL m.fmlGR], 'kd', ...
              'MarkerFaceColor', 'm', 'MarkerSize', 7);
    h9 = plot(m.GxPos, m.Gx, 'ks', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    text(m.fmlPos, m.fmlG + 4, sprintf('FML: %.1f dB', m.fmll), ...
         'Color', [0.6 0 0.6], 'FontSize', 8, 'HorizontalAlignment', 'center', ...
         'FontWeight', 'bold');
    text(m.GxPos, m.Gx + 4, sprintf('G_{x,max} = %.1f dBi', m.Gx), ...
         'Color', [0.8 0 0], 'FontSize', 8, 'HorizontalAlignment', 'center', ...
         'FontWeight', 'bold');
    inZ = abs(f.theta) <= lim;
    ylim([min(m.G - 70, max(f.xp(inZ)) - 15), m.G + 6]);
    xlim([-lim lim]);  grid on;
    xlabel(theta_label(c.phi0));  ylabel('Gain [dBi]');
    title(sprintf(['plane \\phi = %d\\circ:  G_{max} = %.2f dBi,  HPBW = %.4f\\circ\n' ...
                   'FMLL = %.1f dB,  XPD = %.1f dB'], ...
          c.phi0, m.G, m.hpbw, m.fmll, m.G - m.Gx));
    if p == 1                % one legend for all, in a row at the bottom
        lg = legend([h1 h2 h3 h4 h5 h6 h7 h8 h9], ...
               {'G spline', 'G samples', 'G_x spline', 'G_x samples', 'G_{max}', ...
                '-3 dB points', 'first nulls', 'first minor lobes (FML)', ...
                'G_{x,max} in main lobe'}, 'Orientation', 'horizontal', 'FontSize', 8);
    end
end
set(lg, 'Units', 'normalized', 'Position', [0.03 0.005 0.94 0.03]);
fig_title(sprintf('%s reflector - \\theta cuts, zoom on the main beam', a.name));
save_fig(fig, fullfile(outDir, [a.tag '_theta_cuts_zoom']), cfg);
end

% -------------------------------------------------------------------------
function plot_cuts_overlay(a, d, ip, cfg, outDir)
%PLOT_CUTS_OVERLAY  All the planes on the same axes (beam symmetry).
fig = figure('Color', 'w', 'Position', [60 60 1000 600]);
col = lines(numel(d.fineTp));
lab = {};
for p = 1:numel(d.fineTp)
    f = d.fineTp(p);
    plot(f.theta, f.co, '-', 'Color', col(p,:), 'LineWidth', 1.3);  hold on;
    lab{end+1} = sprintf('G,   \\phi = %d\\circ', f.phi0);  %#ok<AGROW>
end
for p = 1:numel(d.fineTp)
    f = d.fineTp(p);
    plot(f.theta, f.xp, '--', 'Color', col(p,:), 'LineWidth', 1.0);
    lab{end+1} = sprintf('G_x, \\phi = %d\\circ', f.phi0);  %#ok<AGROW>
end
lim = a.thetaZoomDeg;
plot([-lim lim], (ip.tp.G - 3)*[1 1], 'k:');  lab{end+1} = '-3 dB';
xlim([-lim lim]);  ylim([ip.tp.G - 90, ip.tp.G + 3]);  grid on;
xlabel('\theta [deg] (signed)');  ylabel('Gain [dBi]');
title(sprintf('%s reflector - G (solid) and G_x (dashed) in all the planes (spline)', a.name));
legend(lab, 'Location', 'eastoutside');
save_fig(fig, fullfile(outDir, [a.tag '_theta_cuts_overlay']), cfg);
end

% -------------------------------------------------------------------------
function plot_axial_ratio(a, d, ip, cfg, outDir)
%PLOT_AXIAL_RATIO  AR versus theta in the main beam (spline data).
fig = figure('Color', 'w', 'Position', [60 60 1000 600]);
col = lines(numel(d.fineTp));
lab = {};
for p = 1:numel(d.fineTp)
    f = d.fineTp(p);  m = ip.tp.planes(p);
    in = f.theta >= m.nullL & f.theta <= m.nullR;      % main lobe only
    semilogy(f.theta(in), axial_ratio_dB(f.co(in), f.xp(in)), '-', ...
             'Color', col(p,:), 'LineWidth', 1.3);  hold on;
    lab{end+1} = sprintf('\\phi = %d\\circ', f.phi0);   %#ok<AGROW>
end
hb = mean(ip.tp.hpbw) / 2;
yl = ylim;  yl(1) = max(yl(1), 1e-4);      % the null at theta = 0 is not
ylim(yl);                                  % interesting: cut the axis
plot(-hb*[1 1], yl, 'k--');  plot(hb*[1 1], yl, 'k--');
lab{end+1} = '\pm HPBW/2';
xlim(1.1 * ip.tp.thNull * [-1 1]);  grid on;
xlabel('\theta [deg] (signed)');  ylabel('AR [dB]');
title(sprintf(['%s reflector - axial ratio AR inside the main lobe\n' ...
               'AR at G_{max} = %.2e dB,  max AR in the -3 dB beam = %.2e dB'], ...
      a.name, ip.tp.ARpeak, ip.tp.ARmax3));
legend(lab, 'Location', 'north');
save_fig(fig, fullfile(outDir, [a.tag '_axial_ratio']), cfg);
end

% -------------------------------------------------------------------------
function plot_uv_vs_tp(a, d, ip, cfg, outDir)
%PLOT_UV_VS_TP  The same plane from the two files: consistency check.
fig = figure('Color', 'w', 'Position', [60 60 1250 820]);
lim = a.thetaZoomDeg;
for p = 1:numel(d.fineTp)
    subplot(2, 2, p);
    plot(d.fineTp(p).theta, d.fineTp(p).co, 'b-', 'LineWidth', 1.3);  hold on;
    plot(d.fineUv(p).theta, d.fineUv(p).co, 'r--', 'LineWidth', 1.1);
    plot(d.cutsUv(p).theta, d.cutsUv(p).co, 'ko', 'MarkerSize', 4);
    xlim([-lim lim]);
    ylim([ip.tp.planes(p).G - 50, ip.tp.planes(p).G + 3]);  grid on;
    xlabel('\theta [deg] (signed)');  ylabel('G [dBi]');
    title(sprintf('plane \\phi = %d\\circ:  HPBW \\theta-\\phi file %.4f\\circ / uv file %.4f\\circ', ...
          d.fineTp(p).phi0, ip.tp.planes(p).hpbw, ip.uv.planes(p).hpbw));
    legend({'\theta-\phi file (spline)', 'uv file (spline)', 'uv file samples'}, ...
           'Location', 'northwest', 'FontSize', 7);
end
fig_title(sprintf('%s reflector - comparison of the two files', a.name));
save_fig(fig, fullfile(outDir, [a.tag '_uv_vs_tp_comparison']), cfg);
end

% -------------------------------------------------------------------------
function plot_normalized_comparison(ant, dat, th, ip, cfg)
%PLOT_NORMALIZED_COMPARISON  Both reflectors versus X = (pi D/lambda)
%sin(theta): for a circular aperture the pattern depends on theta only
%through X, so reflectors with the same f/D and edge taper overlap.
fig = figure('Color', 'w', 'Position', [60 60 1300 560]);
col = [0 0.35 0.8; 0.85 0.2 0.1];
for pl = 1:2
    subplot(1, 2, pl);
    p = find(cfg.planesDeg == 90*(pl-1));        % planes phi = 0 and 90
    lab = {};
    for k = 1:numel(ant)
        c = dat(k).cutsTp(p);
        X = pi * ant(k).D / th(k).lambda * sind(c.theta);
        plot(X, c.co - ip(k).tp.planes(p).G, '-', 'Color', col(k,:), 'LineWidth', 1.3);
        hold on;  lab{end+1} = [ant(k).name ' G'];   %#ok<AGROW>
    end
    for k = 1:numel(ant)
        c = dat(k).cutsTp(p);
        X = pi * ant(k).D / th(k).lambda * sind(c.theta);
        plot(X, c.xp - ip(k).tp.planes(p).G, ':', 'Color', col(k,:), 'LineWidth', 1.0);
        lab{end+1} = [ant(k).name ' G_x'];           %#ok<AGROW>
    end
    xlim([-20 20]);  ylim([-100 3]);  grid on;
    xlabel('X = (\piD/\lambda) sin\theta');  ylabel('G - G_{max} [dB]');
    title(sprintf('plane \\phi = %d\\circ', cfg.planesDeg(p)));
    legend(lab, 'Location', 'south');
end
fig_title('Normalized patterns of the two reflectors (raw \theta-\phi data)');
save_fig(fig, fullfile(cfg.outDir, 'comparison_S_vs_Ka_normalized'), cfg);
end

% -------------------------------------------------------------------------
function fig_title(str)
% figure title (sgtitle in MATLAB >= R2018b, a text box otherwise)
if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
    sgtitle(str, 'FontWeight', 'bold');
else
    annotation('textbox', [0 0.96 1 0.04], 'String', str, 'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
               'FontSize', 12, 'FitBoxToText', 'off');
end
end

function fig_note(str)
% explanation line at the bottom of a figure
annotation('textbox', [0 0 1 0.03], 'String', str, 'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', 'FontSize', 9, 'FitBoxToText', 'off');
end

function save_fig(fig, base, cfg)
drawnow;
if cfg.saveFigures
    print(fig, [base '.png'], '-dpng', sprintf('-r%d', cfg.dpi));
end
end

% =========================================================================
%  F7 - REPORT
% =========================================================================
function write_report(fid, ant, th, raw, ip, cfg)
%WRITE_REPORT  All the results, in the order: [A] theory, [B] raw data,
%[C] interpolated data, then a final summary. fid = 1 -> screen.
pr  = @(varargin) fprintf(fid, varargin{:});
bar = repmat('=', 1, 78);
sep = repmat('-', 1, 78);

pr('\n%s\n', bar);
pr('  REFLECTOR RADIATION PATTERNS - PERFORMANCE REPORT\n');
pr('%s\n', bar);
pr('  G      co-polar gain (LHCP, E-l)              G_max  its peak value\n');
pr('  G_x    cross-polar gain (RHCP, E-r)            G0     (pi D/lambda)^2\n');
pr('  eta    aperture efficiency = G_max / G0        HPBW   -3 dB beamwidth\n');
pr('  FML    first minor lobe (first lobe after the first null of the main beam)\n');
pr('  FMLL   first minor lobe level = G(FML) - G_max  [dB]\n');
pr('  XPD    cross-polar discrimination = G_max - G_x,max  [dB]\n');
pr('  AR     axial ratio (0 dB = perfect circular polarization)\n');
pr('  theta-phi file : 4 full plane cuts (phi = 0/45/90/135 deg)\n');
pr('  u-v file       : 401 x 401 grid, u,v in [-0.4, 0.4]\n');

for k = 1:numel(ant)
    a = ant(k);  t = th(k);
    pr('\n\n%s\n', bar);
    pr('  %s REFLECTOR   D = %g m, F = %g m, f = %g GHz\n', ...
       upper(a.name), a.D, a.F, a.freqHz/1e9);
    pr('%s\n', bar);

    % ------------------------------------------------------------- [A]
    pr('\n[A] THEORETICAL EXPECTED VALUES\n%s\n', sep);
    pr('  Uniformly illuminated circular aperture of diameter D (Airy pattern)\n');
    pr('    lambda                                    %10.4f m\n', t.lambda);
    pr('    D / lambda                                %10.1f\n', t.DoL);
    pr('    G0 = (pi D / lambda)^2     (eta = 100 %%)  %10.2f dBi\n', t.G0);
    pr('    HPBW                                      %10.4f deg\n', t.hpbwUniform);
    pr('    first null                                %10.4f deg\n', t.nullUniform);
    pr('    first minor lobe (FML) position           %10.4f deg\n', t.fmlUniform);
    pr('    first minor lobe level FMLL               %10.2f dB\n', t.fmllUniform);
    pr('  Reflector with ~ -11 dB edge taper (Balanis sec. 15.4)\n');
    pr('    eta (optimum edge taper)                  %10.0f %%\n', 100*t.etaTaper);
    pr('    expected G_max = G0 + 10 log10(eta)       %10.2f dBi\n', t.Gtaper);
    pr('    HPBW (rule of thumb 70 lambda / D)        %10.4f deg\n', t.hpbwTaper);
    pr('    FMLL                                      below -17.6 dB (lower than uniform)\n');
    pr('  Circular polarization on an offset reflector\n');
    pr('    beam squint asin(lambda sin(aoff)/(4 pi F))%9.4f deg  (towards phi = +-90 deg)\n', t.squint);
    pr('    ideal cross-polarization: very low, null on the beam axis\n');

    % ------------------------------------------------------------- [B]
    r = raw(k);
    pr('\n[B] VALUES FROM THE RAW FILE DATA (file samples only, no interpolation)\n%s\n', sep);
    pr('    (-3 dB points: straight line between the 2 samples around -3 dB)\n');
    table_header(pr);
    row(pr, 'sampling step (on the axes)', 'deg', '%.4f', r.stepTp, r.stepUv);
    write_metrics(pr, r.tp, r.uv, cfg, true);

    % ------------------------------------------------------------- [C]
    q = ip(k);
    pr('\n[C] VALUES FROM THE INTERPOLATED DATA (cubic spline of the complex field)\n%s\n', sep);
    pr('    E = 10^(G/20) exp(j phase) from amplitude and phase of the files,\n');
    pr('    real and imaginary parts interpolated with cubic splines, G = 20 log10|E|\n');
    pr('    theta-phi cuts resampled with a step %d times finer (|theta| <= %.2f deg)\n', ...
       cfg.upTp, a.thetaZoomDeg);
    pr('    uv grid resampled %d times finer around the beam\n', cfg.upUv);
    table_header(pr);
    write_metrics(pr, q.tp, q.uv, cfg, false);

    % ------------------------------------------------------- warnings
    pr('\n  Sampling of the main beam: %.1f samples per HPBW in the theta-phi file,\n', ...
       q.samplesPerHpbwTp);
    pr('  %.1f in the u-v file.\n', q.samplesPerHpbwUv);
    if q.samplesPerHpbwUv < 2
        pr('  WARNING: the u-v grid has less than 2 samples per HPBW: the main beam\n');
        pr('  falls between its points, so the u-v values (raw AND interpolated) of\n');
        pr('  G_max, HPBW, FML and cross-pol are NOT reliable. Use the theta-phi file.\n');
    end
end

% ----------------------------------------------------------------- summary
pr('\n\n%s\n', bar);
pr('  SUMMARY - best estimates (theta-phi file, spline interpolation)\n');
pr('%s\n', bar);
pr('  %-36s', '');
for k = 1:numel(ant); pr('%16s', ant(k).name); end;  pr('\n');
sumrow(pr, 'G_max                         [dBi]', '%16.2f', arrayfun(@(x) x.tp.G, ip));
sumrow(pr, 'expected G_max (eta = 82 %)   [dBi]', '%16.2f', [th.Gtaper]);
sumrow(pr, 'eta                             [%]', '%16.1f', 100*arrayfun(@(x) x.tp.eta, ip));
sumrow(pr, 'peak direction theta          [deg]', '%16.4f', arrayfun(@(x) x.tp.th, ip));
sumrow(pr, 'peak direction phi            [deg]', '%16.1f', arrayfun(@(x) x.tp.ph, ip));
sumrow(pr, 'HPBW (mean of the planes)     [deg]', '%16.4f', arrayfun(@(x) mean(x.tp.hpbw), ip));
sumrow(pr, 'FMLL (highest of the planes)   [dB]', '%16.2f', arrayfun(@(x) x.tp.fmlWorst.fmll, ip));
sumrow(pr, 'G_x,max in the main beam      [dBi]', '%16.2f', arrayfun(@(x) x.tp.Gx, ip));
sumrow(pr, 'XPD = G_max - G_x,max          [dB]', '%16.2f', arrayfun(@(x) x.tp.XPD, ip));
sumrow(pr, 'max AR in the -3 dB beam       [dB]', '%16.2e', arrayfun(@(x) x.tp.ARmax3, ip));
pr('\n');
end

function table_header(pr)
pr('  %-38s %6s %14s %14s\n', '', '', 'theta-phi file', 'u-v file');
end

function row(pr, label, unit, fmt, a, b)
% one line of the two-column table (theta-phi file | u-v file)
pr('  %-38s %6s %14s %14s\n', label, ['[' unit ']'], sprintf(fmt, a), sprintf(fmt, b));
end

function sumrow(pr, label, fmt, vals)
pr('  %-36s', label);  pr(fmt, vals);  pr('\n');
end

function write_metrics(pr, T, U, cfg, isRaw)
%WRITE_METRICS  Common block of results for [B] and [C].
pr('  CO-POLAR (LHCP)\n');
row(pr, 'G_max', 'dBi', '%.3f', T.G, U.G);
row(pr, 'eta = G_max / G0', '%', '%.1f', 100*T.eta, 100*U.eta);
row(pr, 'peak direction  theta', 'deg', '%.4f', T.th, U.th);
row(pr, '                phi', 'deg', '%.1f', T.ph, U.ph);
row(pr, '                u = sin(theta) cos(phi)', '-', '%.5f', T.u, U.u);
row(pr, '                v = sin(theta) sin(phi)', '-', '%.5f', T.v, U.v);
for p = 1:numel(cfg.planesDeg)
    row(pr, sprintf('HPBW, plane phi = %3d deg', cfg.planesDeg(p)), 'deg', '%.4f', ...
        T.hpbw(p), U.hpbw(p));
end
row(pr, 'HPBW, mean', 'deg', '%.4f', mean(T.hpbw), mean(U.hpbw));
for p = 1:numel(cfg.planesDeg)
    row(pr, sprintf('first null, plane phi = %3d deg', cfg.planesDeg(p)), 'deg', '%.4f', ...
        T.null(p), U.null(p));
end
pr('  FIRST MINOR LOBE (FML)\n');
for p = 1:numel(cfg.planesDeg)
    row(pr, sprintf('FML position, plane phi = %3d deg', cfg.planesDeg(p)), 'deg', '%.4f', ...
        T.planes(p).fmlPos, U.planes(p).fmlPos);
    row(pr, sprintf('FMLL,         plane phi = %3d deg', cfg.planesDeg(p)), 'dB', '%.2f', ...
        T.planes(p).fmll, U.planes(p).fmll);
end
row(pr, 'highest FMLL', 'dB', '%.2f', T.fmlWorst.fmll, U.fmlWorst.fmll);
row(pr, '   in plane phi', 'deg', '%d', T.fmlWorst.plane, U.fmlWorst.plane);
row(pr, '   at theta (signed)', 'deg', '%.4f', T.fmlWorst.theta, U.fmlWorst.theta);
row(pr, '   G(FML)', 'dBi', '%.2f', T.fmlWorst.G, U.fmlWorst.G);
pr('  CROSS-POLAR (RHCP)\n');
pr('    main-beam region = inside the first null of G (theta <= theta_null)\n');
row(pr, 'theta_null (mean of the planes)', 'deg', '%.4f', T.thNull, U.thNull);
row(pr, 'G_x,max  (max of G_x in main beam)', 'dBi', '%.2f', T.Gx, U.Gx);
row(pr, '   at theta', 'deg', '%.4f', T.GxTh, U.GxTh);
row(pr, '   at phi', 'deg', '%.1f', T.GxPh, U.GxPh);
row(pr, 'XPD = G_max - G_x,max', 'dB', '%.2f', T.XPD, U.XPD);
row(pr, 'worst XPD in -3 dB beam: min(G-G_x)', 'dB', '%.2f', T.XPDworst3, U.XPDworst3);
row(pr, 'AR at the peak of G', 'dB', '%.2e', T.ARpeak, U.ARpeak);
row(pr, 'max AR in the -3 dB beam', 'dB', '%.2e', T.ARmax3, U.ARmax3);
row(pr, 'points inside the -3 dB beam', '-', '%d', T.n3dB, U.n3dB);
if isRaw
    pr('    max of G_x over the WHOLE file (far from the beam, for information):\n');
    row(pr, 'G_x max, whole file', 'dBi', '%.2f', T.GxAll, U.GxAll);
    row(pr, '   at theta', 'deg', '%.4f', T.GxAllTh, U.GxAllTh);
    row(pr, '   at phi', 'deg', '%.1f', T.GxAllPh, U.GxAllPh);
end
end

% =========================================================================
%  F8 - UTILITIES
% =========================================================================
function S = add_to_array(S, s, k)
% S(k) = s for a struct array that may still be empty
if isempty(S); S = s; else; S(k) = s; end
end
