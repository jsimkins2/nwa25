# Core Masking



IO Layout 2,2


### Combining Output
mppnccombine from FRE-nc tools to combine the files together

worked for ocean_daily, not for ocean_5day

Error: cannot write variable "average_DT"'s values!

aka files too large

used this to split files in half by time variable

% ncks -d time,0,729 hgt.2006.nc hgt.2006-01.nc
% ncks -d time,730,1459 hgt.2006.nc hgt.2006-02.nc

still to large ( on triton at least)