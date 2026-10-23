% WhaleCallLibrary - synthetic marine-mammal call generator for WAH-i Ocean.
%
% Presets are deliberately simple parametric approximations (tonal, FM sweep,
% or band-limited pulse) based on commonly reported frequency ranges and
% durations. They are meant to put the simulated signal in the right
% bandwidth/duration regime for TDOA performance, not to be faithful replicas.
% Adjust the parameters (or import a template) for your population.
%
% Usage
%   P   = WhaleCallLibrary.presets();               % struct array of presets
%   p   = WhaleCallLibrary.getPreset('Fin whale - 20 Hz pulse');
%   sig = WhaleCallLibrary.synthesize(spec, fs);    % column vector, unit peak
%   tpl = WhaleCallLibrary.loadTemplate(file);      % read a vector / wav
%   [fLo,fHi] = WhaleCallLibrary.estimateBand(sig, fs);
%
% spec fields: kind ('tonal'|'fm'|'pulse'|'template_wave'|'template_contour'),
%   f0, f1 (Hz; band edges for 'pulse'), dur (s), sweep ('linear'|'exponential'),
%   taper (fraction, Tukey), template (vector), templateFs (Hz, waveforms)

classdef WhaleCallLibrary
    methods (Static)

        function P = presets()
            % name, kind, f0, f1, dur(s), fs suggestion(Hz), sweep, SL (dB re 1uPa@1m), note
            D = {
                'Fin whale - 20 Hz pulse',            'fm',    23,    17,    1.0,    500,   'linear',      189, 'Downswept ~23->17 Hz, ~1 s (Watkins et al. 1987; Sirovic et al. 2007).'
                'Blue whale - NE Pacific B call',     'fm',    16.5,  15.0,  10.0,   250,   'linear',      186, 'Tonal ~16 Hz fundamental, ~10 s (McDonald et al. 2006).'
                'Blue whale - D call',                'fm',    90,    30,    2.0,    500,   'linear',      179, 'Variable downsweep ~90->30 Hz, 1-4 s.'
                'Antarctic blue whale - Z call (A)',  'fm',    28,    26,    8.0,    250,   'linear',      189, 'Approximates the ~27 Hz tonal first unit of the Z call.'
                'Humpback whale - song unit',         'fm',    350,   150,   1.5,    2000,  'exponential', 170, 'Generic FM song unit; real units span ~30 Hz-4 kHz.'
                'N. Atlantic right whale - upcall',   'fm',    90,    180,   1.0,    2000,  'linear',      172, 'Upsweep roughly 50-250 Hz, ~1 s.'
                'Minke whale - pulse',                'pulse', 100,   300,   0.25,   2000,  'linear',      165, 'Band-limited pulse approximation of N. Atlantic minke pulse trains.'
                'Sperm whale - click',                'pulse', 5000,  20000, 0.002,  96000, 'linear',      200, 'Broadband click, energy ~5-20 kHz (off-axis SL used).'
                'Delphinid - whistle',                'fm',    8000,  16000, 0.6,    48000, 'linear',      150, 'Upswept tonal whistle, generic delphinid.'
                'Beaked whale - FM pulse (Cuvier''s)', 'fm',   30000, 50000, 0.0003, 192000,'linear',      214, 'FM upswept pulse ~30-50 kHz, ~200-300 us (Zimmer et al. 2005).'
                'Custom - FM sweep',                  'fm',    100,   50,    1.0,    2000,  'linear',      180, 'User-defined sweep.'
                'Custom - tonal',                     'tonal', 50,    50,    1.0,    1000,  'linear',      180, 'Constant-frequency tone (uses Start frequency).'
                'Custom - pulse',                     'pulse', 1000,  5000,  0.005,  24000, 'linear',      180, 'Band-limited Gaussian pulse between Start/End frequency.'
                'Imported template - waveform',       'template_wave',    NaN, NaN, NaN, NaN, 'linear',  180, 'Waveform samples imported from file (csv/txt/mat/wav).'
                'Imported template - freq. contour',  'template_contour', NaN, NaN, NaN, NaN, 'linear',  180, 'Frequency contour (Hz) imported from file, stretched over Call duration.'
                };
            P = struct('name',D(:,1),'kind',D(:,2),'f0',D(:,3),'f1',D(:,4),'dur',D(:,5), ...
                'fs',D(:,6),'sweep',D(:,7),'SL',D(:,8),'note',D(:,9));
        end

        function names = presetNames()
            P = WhaleCallLibrary.presets();
            names = {P.name};
        end

        function p = getPreset(name)
            P = WhaleCallLibrary.presets();
            idx = find(strcmp({P.name}, name), 1);
            if isempty(idx), error('Unknown preset: %s', name); end
            p = P(idx);
        end

        %% ------------------------------------------------------------------
        function sig = synthesize(spec, fs)
            if ~isfield(spec,'taper') || isempty(spec.taper), spec.taper = 0.1; end
            switch spec.kind
                case 'tonal'
                    n = max(8, round(spec.dur*fs));
                    t = (0:n-1)'/fs;
                    sig = sin(2*pi*spec.f0*t);
                    sig = sig .* WhaleCallLibrary.tukey(n, spec.taper);
                case 'fm'
                    n = max(8, round(spec.dur*fs));
                    u = (0:n-1)'/max(n-1,1);
                    if strcmpi(spec.sweep,'exponential') && spec.f0 > 0 && spec.f1 > 0
                        f = spec.f0 * (spec.f1/spec.f0).^u;
                    else
                        f = spec.f0 + (spec.f1 - spec.f0)*u;
                    end
                    sig = WhaleCallLibrary.fromContour(f, fs);
                    sig = sig .* WhaleCallLibrary.tukey(n, spec.taper);
                case 'pulse'
                    fc = 0.5*(spec.f0 + spec.f1);
                    bw = max(abs(spec.f1 - spec.f0), 1);
                    sigma_f = bw/5;                     % band edges ~ +/-2.5 sigma (about -27 dB)
                    sigma_t = 1/(2*pi*sigma_f);
                    T = max(spec.dur, 8*sigma_t);
                    n = max(8, round(T*fs));
                    t = ((0:n-1)' - (n-1)/2)/fs;
                    sig = exp(-0.5*(t/sigma_t).^2) .* cos(2*pi*fc*t);
                case 'template_wave'
                    w = spec.template(:);
                    if isempty(w), error('No waveform template loaded.'); end
                    tfs = spec.templateFs;
                    if ~isfinite(tfs) || tfs <= 0, tfs = fs; end
                    if abs(tfs - fs) > 1e-6
                        tin = (0:numel(w)-1)'/tfs;
                        tout = (0:1/fs:tin(end))';
                        w = interp1(tin, w, tout, 'pchip', 0);   % band-limited enough for templates
                    end
                    w = w - mean(w);
                    sig = w .* WhaleCallLibrary.tukey(numel(w), min(spec.taper, 0.05));
                case 'template_contour'
                    fc = spec.template(:);
                    if isempty(fc), error('No frequency contour template loaded.'); end
                    n = max(8, round(spec.dur*fs));
                    f = interp1(linspace(0,1,numel(fc)), fc, linspace(0,1,n)', 'pchip');
                    f = max(f, 0);
                    sig = WhaleCallLibrary.fromContour(f, fs);
                    if isfield(spec,'templateAmp') && ~isempty(spec.templateAmp)
                        a = interp1(linspace(0,1,numel(spec.templateAmp)), spec.templateAmp(:), linspace(0,1,n)', 'linear');
                        sig = sig .* max(a, 0);
                    end
                    sig = sig .* WhaleCallLibrary.tukey(n, spec.taper);
                otherwise
                    error('Unknown call kind: %s', spec.kind);
            end
            sig = sig(:);
            pk = max(abs(sig));
            if pk > 0, sig = sig / pk; end
        end

        function sig = fromContour(f, fs)
            ph = 2*pi*cumsum(f(:))/fs;
            sig = sin(ph);
        end

        function w = tukey(n, r)
            if n < 2, w = ones(n,1); return; end
            r = min(max(r,0),1);
            w = ones(n,1);
            if r == 0, return; end
            x = (0:n-1)'/(n-1);
            a = x < r/2;
            w(a) = 0.5*(1 + cos(2*pi/r*(x(a) - r/2)));
            b = x >= 1 - r/2;
            w(b) = 0.5*(1 + cos(2*pi/r*(x(b) - 1 + r/2)));
        end

        %% ------------------------------------------------------------------
        function tpl = loadTemplate(file)
            % Returns struct(values, fs, amp, guess, file). fs is set for .wav.
            % Accepts a single vector (any length, e.g. 200 points) in csv/txt/mat,
            % a 2-column [freq amp] contour, or a .wav clip.
            [~,~,ext] = fileparts(file);
            tpl = struct('values',[],'fs',NaN,'amp',[],'guess','','file',file);
            switch lower(ext)
                case {'.wav','.flac','.ogg'}
                    [y, fs] = audioread(file);
                    tpl.values = y(:,1);
                    tpl.fs = fs;
                    tpl.guess = 'waveform';
                    return;
                case '.mat'
                    Sm = load(file);
                    fn = fieldnames(Sm);
                    A = [];
                    for k = 1:numel(fn)
                        v = Sm.(fn{k});
                        if isnumeric(v) && numel(v) >= 4, A = double(v); break; end
                    end
                    if any(strcmpi(fn,'fs')), tpl.fs = double(Sm.(fn{strcmpi(fn,'fs')})); end
                otherwise
                    A = OceanEnvironment.readNumericTable(file);
            end
            if isempty(A), error('No numeric data found in %s', file); end
            if isvector(A)
                v = A(:);
            elseif size(A,2) == 2 && size(A,1) > 2
                v = A(:,1); tpl.amp = A(:,2);
            elseif size(A,1) == 2 && size(A,2) > 2
                v = A(1,:)'; tpl.amp = A(2,:)';
            else
                v = A(:);
            end
            v = v(isfinite(v));
            if numel(v) < 4, error('Template needs at least 4 values.'); end
            tpl.values = v;
            % Heuristic guess: contour = all positive, smooth, typical Hz range
            dv = diff(v);
            smooth = median(abs(dv)) < 0.05*(max(v)-min(v)+eps) || all(dv>=0) || all(dv<=0);
            if all(v > 0) && min(v) >= 1 && smooth
                tpl.guess = 'contour';
            else
                tpl.guess = 'waveform';
            end
        end

        %% ------------------------------------------------------------------
        function [fLo, fHi] = estimateBand(sig, fs)
            % Band containing ~99 % of the call energy
            n = 2^nextpow2(max(numel(sig), 256)*2);
            X = abs(fft(sig(:), n)).^2;
            X = X(1:n/2+1);
            f = (0:n/2)'*fs/n;
            c = cumsum(X); c = c / max(c(end), eps);
            iLo = find(c >= 0.005, 1, 'first');
            iHi = find(c >= 0.995, 1, 'first');
            fLo = f(max(iLo,1)); fHi = f(max(iHi,1));
            bw = max(fHi - fLo, 2*fs/n);
            fLo = max(0, fLo - 0.1*bw);
            fHi = min(fs/2, fHi + 0.1*bw);
        end

        function [S, f, t] = spectrogramSimple(sig, fs)
            % Lightweight STFT (dB) without Signal Processing Toolbox
            sig = sig(:);
            nwin = max(16, 2^floor(log2(max(numel(sig)/16, 16))));
            nwin = min(nwin, 4096);
            hop = max(1, round(nwin/4));
            w = 0.5 - 0.5*cos(2*pi*(0:nwin-1)'/(nwin-1));
            if numel(sig) < nwin, sig = [sig; zeros(nwin-numel(sig),1)]; end
            nfr = floor((numel(sig)-nwin)/hop) + 1;
            nfft = 2^nextpow2(nwin*2);
            S = zeros(nfft/2+1, nfr);
            for k = 1:nfr
                seg = sig((k-1)*hop + (1:nwin)) .* w;
                Xk = fft(seg, nfft);
                S(:,k) = abs(Xk(1:nfft/2+1)).^2;
            end
            S = 10*log10(S + eps*max(S(:)) + eps);
            f = (0:nfft/2)'*fs/nfft;
            t = ((0:nfr-1)*hop + nwin/2)/fs;
        end
    end
end
