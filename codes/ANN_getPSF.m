%% This is the core code for generating the theoretical PSF of dipoles
% This program is based on following papers:
% (1) Mortensen, K. I., Churchman, L. S., Spudich, J. A. & Flyvbjerg, H. Nat. Methods 7, 377-381 (2010).
% (2) Stallinga, S. & Rieger, B. Opt. Express 18, 24461-24476 (2010).
% (3) Patra, D., Gregor, I., Enderlein, J. & Sauer, M. Applied Physics Letters 87,(2005).
% (4) Irving, M. Biophys. J. 70, 1830-1835 (1996).
% The PSF of restricted dipole is implemented by Yongdeng Zhang & Lusheng GU
% 6th of May, 2012
%%*************************************************************************

function data = ANN_getPSF(info,self,disp_flag)
    if nargin<3
        disp_flag = 1;
    end
    data = ANN_getPSFp(info,self,disp_flag);
end