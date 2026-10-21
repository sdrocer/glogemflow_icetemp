; Check of firnice_orient = 'table' with GloGEM's own code (orient_check.py writes the table, runs this and
; judges the output). Synthetic glacier of 30 bands, synthetic climate 1961-2020 at the cell.
;
; Run from the GloGEM tree under test:
;   GLOGEM_CONFIG=<a training config> HARNESS_DIR=<folder with orient_99999.dat> HARNESS_OUT=<file> \
;     idl -e ".r <path>/orient_harness.pro"
; Output: line 1 'tref' then nb x 12 reference temperatures; 'F' then nb x nd factors; 'months' then per month
; and band: month, band, air temperature of the thermal melt, melt without and with the factor, and whether the
; mass-balance state changed.

@procedures/initialise/settings.pro

openw, 2, getenv('HARNESS_OUT')
noval = -9999 & nb = 30
elev = 2600d + 70d * dindgen(nb) & elev0 = elev & gl = elev & hclim = 2500d
step = 10 & dpdz = 1.5d & mb_book = 0b & ar_gl = 0d
sur = dblarr(nb) & for i = 0, nb-1 do sur[i] = i mod 3
DDFsnow = 4.5d & DDFice = 9.0d
variab = dblarr(12, 31)
for mm = 0, 11 do for k = 0, 30 do variab[mm, k] = 3d * sin(2d * !dpi * (k + 0.3d * mm) / 31d)
submonth_variability = 'y'
firnice_forcing_toff = 'n' & firnice_thermal_lapse = 'reanalysis'
; climate series at the cell, 1961-2020
new = 'y'
nt = 60l * 12l
ryear = 1961l + lindgen(nt) / 12l & rmon = (lindgen(nt) mod 12l) + 1l
tempre = -2d - 9d * cos(2d * !dpi * (rmon - 1d) / 12d) + 0.03d * (ryear - 1990d) + 0.7d * sin(1.7d * dindgen(nt))
dtdz = -0.0058d - 0.0004d * cos(2d * !dpi * dindgen(12) / 12d)
id = ['99999'] & gg = [0l] & g = 0l
firnice_orient = 'table' & firnice_orient_dir = getenv('HARNESS_DIR')
@procedures/initialise/setup_firnice_orient.pro

printf, 2, 'tref'
for j = 0, nb-1 do printf, 2, forient_tref[j, *], format='(12(f22.15))'
printf, 2, 'F'
for j = 0, nb-1 do printf, 2, forient[j, *], format='(' + strtrim(n_elements(forient_dt), 2) + '(f22.15))'
printf, 2, 'months'
t_offset = 1.3d & pc = 0.08d
for m = 1, 12 do begin
  ; air temperature of each band from 5.5 K below to 7.5 K above its 1981-2020 mean of this month
  tg = forient_tref[*, m - 1] + (-5.5d + 13d * dindgen(nb) / (nb - 1)) + t_offset
  psg = dblarr(nb) & mel = psg & snowmel = mel & icemel = mel
  @procedures/processing/accumulation.pro
  @procedures/processing/meltmodel.pro
  keep = [psg, plg, mel, tg, tgs]
  firnice_orient = 'n'
  @procedures/processing/firnice_thermal_forcing.pro
  mel_n = mel_th & tl = ftw_tl
  firnice_orient = 'table'
  @procedures/processing/firnice_thermal_forcing.pro
  changed = array_equal(keep, [psg, plg, mel, tg, tgs]) ? 0 : 1
  for j = 0, nb-1 do printf, 2, m, j, tl[j], mel_n[j], mel_th[j], changed, format='(2(i4), 3(e26.17), i3)'
endfor
close, 2
end
