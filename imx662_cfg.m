function c = imx662_cfg()
%IMX662_CFG  Shared settings for all IMX662 analyses (lab conventions).
%   Units follow the lab standard used in Lupo (2025, MIT SM thesis):
%   ADU for raw signal, eV for X-ray energy, e- for charge, counts/s/mm^2
%   for count rates, conversion gain in ADU/eV.
    c.D   = '/Users/carissmamcgee/Downloads';
    c.OUT = fullfile(c.D, 'IMX662_results');
    c.zipDark = fullfile(c.D, 'IMX662_DarkCurrent_20250709.zip');
    c.zipLin  = fullfile(c.D, 'IMX662_Linearity_20250624.zip');
    c.zipFe   = fullfile(c.D, 'IMX662_Fe55_SpectrumData_20250625.zip');
    c.zipPo   = fullfile(c.D, 'IMX662_Po210_SpectrumData_20250626.zip');

    c.gainRef   = 0.2655;    % ADU/eV, lab value for IMX662 (Mn K-alpha), used until Fe-55 is re-fit
    c.epsSi     = 3.65;      % eV per electron-hole pair in silicon
    c.fano      = 0.115;     % Fano factor of silicon
    c.pixMM     = 2.9e-3;    % IMX662 pixel pitch [mm]
    c.nx = 1920; c.ny = 1080;
    c.areaMM2   = c.nx*c.ny*c.pixMM^2;   % optical area [mm^2] (~17.44)
    c.frameTime = 0.6;       % s, X-ray event-mode frame exposure
    c.adcMax    = 4095;      % 12-bit ADC

    % Reference line energies [eV] (lab values)
    c.E.SiKa = 1740; c.E.SiEsc = 4160; c.E.MnKa = 5900; c.E.MnKb = 6490;
    c.E.C = 277; c.E.O = 526; c.E.F = 677; c.E.Ti = 4510; c.E.MnKbEsc = 4750;

    % Colour scheme = the lab's thesis colours (Lupo 2025)
    c.col.all    = [0.00 0.00 0.00];   % all events: black
    c.col.n1     = [0.16 0.38 0.86];   % single-pixel events / measured points: blue
    c.col.n1234  = [0.18 0.62 0.24];   % multi-pixel (1-4) events: green
    c.col.hot    = [0.49 0.18 0.56];   % hot-pixel spectrum: purple
    c.col.fit    = [0.85 0.10 0.10];   % fits: red
    c.col.raw    = [0.90 0.12 0.12];   % raw light curve: red
    c.col.reb    = [0.10 0.25 0.80];   % re-binned light curve / Gaussian center: blue
    c.col.hist   = [1.00 0.70 0.30];   % light-curve histogram: orange
    c.col.bar    = [0.24 0.70 0.44];   % pixel-count / gain histograms: green bars
    c.col.orange = [1.00 0.60 0.10];   % saturation line, X-ray linearity fit: orange
    c.col.spat   = [0.20 0.65 0.25];   % spatial error, 5%/95% lines: green
    c.col.temp   = [0.55 0.25 0.75];   % temporal error: purple
    c.fs = 11;                        % base font size
    if ~exist(c.OUT, 'dir'), mkdir(c.OUT); end
end
