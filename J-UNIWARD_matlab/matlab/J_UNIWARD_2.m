function S_STRUCT = J_UNIWARD_2(coverPath, payload)

% -------------------------------------------------------------------------
% Copyright (c) 2013 DDE Lab, Binghamton University, NY.
% All Rights Reserved.
% -------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software for
% educational, research and non-profit purposes, without fee, and without a
% written agreement is hereby granted, provided that this copyright notice
% appears in all copies. The program is supplied "as is," without any
% accompanying services from DDE Lab. DDE Lab does not warrant the
% operation of the program will be uninterrupted or error-free. The
% end-user understands that the program was developed for research purposes
% and is advised not to rely exclusively on the program for any reason. In
% no event shall Binghamton University or DDE Lab be liable to any party
% for direct, indirect, special, incidental, or consequential damages,
% including lost profits, arising out of the use of this software. DDE Lab
% disclaims any warranties, and has no obligations to provide maintenance,
% support, updates, enhancements or modifications.
% -------------------------------------------------------------------------
% Contact: vojtech_holub@yahoo.com | fridrich@binghamton.edu | February
% 2013
%          http://dde.binghamton.edu/download/stego_algorithms/
% -------------------------------------------------------------------------
% This function simulates embedding using J-UNIWARD steganographic 
% algorithm.
% -------------------------------------------------------------------------
% Input:  coverPath ... path to the image
%         payload ..... payload in bits per non zero DCT coefficient
% Output: stego ....... resulting JPEG structure with embedded payload
% -------------------------------------------------------------------------

C_SPATIAL = double(imread(coverPath));
C_STRUCT = jpeg_read(coverPath);
C_COEFFS = C_STRUCT.coef_arrays{1};
C_QUANT = C_STRUCT.quant_tables{1};

wetConst = 10^13;
sgm = 2^(-6);
tic
%% Get 2D wavelet filters - Daubechies 8
% 1D high pass decomposition filter
hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, ...
        -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768];
% 1D low pass decomposition filter
lpdf = (-1).^(0:numel(hpdf)-1).*fliplr(hpdf);

F{1} = lpdf'*hpdf;
F{2} = hpdf'*lpdf;
F{3} = hpdf'*hpdf;
tic
%% Pre-compute impact in spatial domain when a jpeg coefficient is changed by 1
spatialImpact = cell(8, 8);
for bcoord_i=1:8
    for bcoord_j=1:8
        testCoeffs = zeros(8, 8);
        testCoeffs(bcoord_i, bcoord_j) = 1;
        spatialImpact{bcoord_i, bcoord_j} = idct2(testCoeffs)*C_QUANT(bcoord_i, bcoord_j);  %修改每一个DCT系数之后，计算其对应的8*8空域受到的影响
    end
end

%% Get 2D wavelet filters - Daubechies 8
% 1D high pass decomposition filter
hpdf = [-0.0544158422, 0.3128715909, -0.6756307363, 0.5853546837, 0.0158291053, -0.2840155430, -0.0004724846, 0.1287474266, 0.0173693010, -0.0440882539, ...
        -0.0139810279, 0.0087460940, 0.0048703530, -0.0003917404, -0.0006754494, -0.0001174768];
% 1D low pass decomposition filter
lpdf = (-1).^(0:numel(hpdf)-1).*fliplr(hpdf);

F{1} = lpdf'*hpdf;
F{2} = hpdf'*lpdf;
F{3} = hpdf'*hpdf;

%% Pre compute impact on wavelet coefficients when a jpeg coefficient is changed by 1
waveletImpact = cell(numel(F), 8, 8); %每个元素对应23*23的矩阵
for Findex = 1:numel(F)
    for bcoord_i=1:8
        for bcoord_j=1:8
            waveletImpact_temp = imfilter(spatialImpact{bcoord_i, bcoord_j}, F{Findex}, 'full');
            waveletImpact{Findex, bcoord_i, bcoord_j} = sum( abs(waveletImpact_temp(:)))/(23*23);
        end
    end   
end

avar_wavelet_Impact = zeros(8, 8);
for bcoord_i=1:8
    for bcoord_j=1:8
      avar_wavelet_Impact(bcoord_i, bcoord_j) = ( waveletImpact{1, bcoord_i, bcoord_j} + waveletImpact{2, bcoord_i, bcoord_j} + waveletImpact{3, bcoord_i, bcoord_j})/3; 
    end
end
avar_wavelet_Impact_matrix = repmat(avar_wavelet_Impact,[64 64]);                   
q_matrix = avar_wavelet_Impact_matrix;
%% Create reference cover wavelet coefficients (LH, HL, HH)
% Embedding should minimize their relative change. Computation uses mirror-padding
padSize = max([size(F{1})'; size(F{2})']);
C_SPATIAL_PADDED = padarray(C_SPATIAL, [padSize padSize], 'symmetric'); % pad image
nzAC = nnz(C_COEFFS)-nnz(C_COEFFS(1:8:end,1:8:end));

tic
avr = ones(23,23);
RC_2 = 0;
for i=1:numel(F)
    RC_temp = abs(imfilter(C_SPATIAL_PADDED, F{i}));
    RC_2 = RC_2 + imfilter(RC_temp, avr);    %%  (21,21), (29 29)
end
complex_block = zeros(64,64);
for block_R_M = 1:64
    for block_R_N = 1:64
        complex_block(block_R_M,block_R_N)  = RC_2((block_R_M-1)*8+21,(block_R_N-1)*8+21);
    end
end
J2 = complex_block (:);
J = ones(64,1)*J2';
J = col2im(J,[8 8], [512 512], 'distinct'); 
save('J','J')

rho = q_matrix./ ((sgm+J)/(23^2));
toc
rhoM1 = rho;
rhoP1 = rho;

rhoP1(rhoP1 > wetConst) = wetConst;
rhoP1(isnan(rhoP1)) = wetConst;    
rhoP1(C_COEFFS > 1023) = wetConst;
    
rhoM1(rhoM1 > wetConst) = wetConst;
rhoM1(isnan(rhoM1)) = wetConst;
rhoM1(C_COEFFS < -1023) = wetConst;
        
%% Embedding simulation
S_COEFFS = EmbeddingSimulator(C_COEFFS, rhoP1, rhoM1, round(payload * nzAC));

S_STRUCT = C_STRUCT;
S_STRUCT.coef_arrays{1} = S_COEFFS;

function [y, pChangeP1, pChangeM1] = EmbeddingSimulator(x, rhoP1, rhoM1, m)

    x = double(x);
    n = numel(x);
    
    lambda = calc_lambda(rhoP1, rhoM1, m, n);
    pChangeP1 = (exp(-lambda .* rhoP1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
    pChangeM1 = (exp(-lambda .* rhoM1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));

    randChange = rand(size(x));
    y = x;
    y(randChange < pChangeP1) = y(randChange < pChangeP1) + 1;
    y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
    
    function lambda = calc_lambda(rhoP1, rhoM1, message_length, n)

        l3 = 1e+3;
        m3 = double(message_length + 1);
        iterations = 0;
        while m3 > message_length
            l3 = l3 * 2;
            pP1 = (exp(-l3 .* rhoP1))./(1 + exp(-l3 .* rhoP1) + exp(-l3 .* rhoM1));
            pM1 = (exp(-l3 .* rhoM1))./(1 + exp(-l3 .* rhoP1) + exp(-l3 .* rhoM1));
            m3 = ternary_entropyf(pP1, pM1);
            iterations = iterations + 1;
            if (iterations > 10)
                lambda = l3;
                return;
            end
        end        
        
        l1 = 0; 
        m1 = double(n);        
        lambda = 0;
        
        alpha = double(message_length)/n;
        % limit search to 30 iterations
        % and require that relative payload embedded is roughly within 1/1000 of the required relative payload        
        while  (double(m1-m3)/n > alpha/1000.0 ) && (iterations<30)
            lambda = l1+(l3-l1)/2; 
            pP1 = (exp(-lambda .* rhoP1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
            pM1 = (exp(-lambda .* rhoM1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
            m2 = ternary_entropyf(pP1, pM1);
    		if m2 < message_length
    			l3 = lambda;
    			m3 = m2;
            else
    			l1 = lambda;
    			m1 = m2;
            end
    		iterations = iterations + 1;
        end
    end
    
    function Ht = ternary_entropyf(pP1, pM1)
        pP1 = pP1(:);
        pM1 = pM1(:);
        Ht = -(pP1.*log2(pP1))-(pM1.*log2(pM1))-((1-pP1-pM1).*log2(1-pP1-pM1));
        Ht(isnan(Ht)) = 0;
        Ht = sum(Ht);
    end

end

end
