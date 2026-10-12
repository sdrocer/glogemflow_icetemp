; Checks of the thermal forcing switches with GloGEM's own code (thermal_forcing_check.py runs
; this and judges the output).
;
; Run from the GloGEM root:
;   GLOGEM_CONFIG=<a training config> HARNESS_PART=water|hl|latent HARNESS_OUT=<file> \
;     [HARNESS_W=<m w.e.> HARNESS_F=<f_ref>] idl -e ".r <path>/thermal_forcing_harness.pro"
;
; water:  firnice_thermal_forcing.pro against accumulation.pro and meltmodel.pro: identical with
;         t_offset = 0; with an offset, equal to the mass balance run without it; mass-balance state
;         untouched.
; hl:     firnice_perm_depth.pro with firnice_hl1980 = 'y' on a grid of accumulation and temperature
;         (constants taken from initialise_firnicetemp_spinup.pro).
; latent: one firn band, one month, conduction switched off (rf_dt tiny): the temperature change
;         from HARNESS_W m w.e. of water times HARNESS_F.

@procedures/initialise/settings.pro

part = getenv('HARNESS_PART')
openw, 2, getenv('HARNESS_OUT')

if part eq 'water' then begin
  noval = -9999 & nb = 30
  elev = 2000d + 100d * dindgen(nb) & gl = elev & gl[nb-1] = noval & hclim = 2000d
  step = 10 & dpdz = 1.5d & mb_book = 0b & ar_gl = 0d
  sur = dblarr(nb) & for i = 0, nb-1 do sur[i] = i mod 3 & sur[nb-1] = noval
  DDFsnow = 4.5d & DDFice = 9.0d
  variab = dblarr(12, 31)
  for mm = 0, 11 do for k = 0, 30 do variab[mm, k] = 3d * sin(2d * !dpi * (k + 0.3d * mm) / 31d)
  n_exact = 0l & n_fail = 0l & n_state = 0l & max_shift = 0d
  for pass = 0, 1 do begin
    submonth_variability = (pass eq 0) ? 'y' : 'n'
    for m = 1, 12 do begin
      for pat = 0, 2 do begin
        base = -15d + 27d * dindgen(nb) / (nb-1) + pat * 0.37d
        base[3] = 0d & base[4] = T_thres - 1 & base[5] = T_thres + 1 & base[6] = 0.5d & base[7] = 2.5d
        ; (1) t_offset = 0: the thermal forcing is the mass balance's, bit for bit
        t_offset = 0d & tg = base & pc = 0.08d
        psg = dblarr(nb) & mel = psg & snowmel = mel & icemel = mel
        @procedures/processing/accumulation.pro
        @procedures/processing/meltmodel.pro
        keep = [psg, plg, mel, tg, tgs]
        firnice_forcing_toff = 'n'
        @procedures/processing/firnice_thermal_forcing.pro
        if not array_equal(keep, [psg, plg, mel, tg, tgs]) then n_state++
        if array_equal(mel_th, mel) and array_equal(plg_th, plg) and array_equal(tgs_th, tgs) then n_exact++ else n_fail++
        firnice_forcing_toff = 'y'
        @procedures/processing/firnice_thermal_forcing.pro
        if not (array_equal(mel_th, mel) and array_equal(plg_th, plg) and array_equal(tgs_th, tgs)) then n_fail++
        ; (2) an offset of 2.67 K: equal to the mass balance without it (no value on a threshold)
        base2 = base + 0.0123d
        t_offset = 0d & tg = base2
        psg = dblarr(nb) & mel = psg & snowmel = mel & icemel = mel
        @procedures/processing/accumulation.pro
        @procedures/processing/meltmodel.pro
        ref = [mel, plg, tgs]
        t_offset = 2.67d & tg = base2 + t_offset
        psg = dblarr(nb) & mel = psg & snowmel = mel & icemel = mel
        @procedures/processing/accumulation.pro
        @procedures/processing/meltmodel.pro
        keep = [psg, plg, mel, tg, tgs]
        firnice_forcing_toff = 'n'
        @procedures/processing/firnice_thermal_forcing.pro
        if not array_equal(keep, [psg, plg, mel, tg, tgs]) then n_state++
        max_shift = max_shift > max(abs([mel_th, plg_th, tgs_th] - ref))
      endfor
    endfor
  endfor
  printf, 2, n_exact, n_fail, n_state, max_shift, fo='(3i8, e14.4)'
endif

if part eq 'hl' then begin
  ; the HL80_ constants exactly as initialise_firnicetemp_spinup.pro defines them
  src = 'procedures/initialise/initialise_firnicetemp_spinup.pro'
  lines_hl = strarr(file_lines(src))
  openr, 3, src & readf, 3, lines_hl & close, 3
  for k = 0, n_elements(lines_hl)-1 do $
    if strmid(strtrim(lines_hl[k], 2), 0, 5) eq 'HL80_' then void = execute(lines_hl[k])
  acc_g = [0.15d, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0]
  tk_g  = [248d, 253, 258, 263, 268, 272]
  nb = n_elements(acc_g) * n_elements(tk_g)
  hl_acc = dblarr(nb) & hl_tk = dblarr(nb) & hl_mask = bytarr(nb) + 1b
  for a = 0, n_elements(acc_g)-1 do for t = 0, n_elements(tk_g)-1 do begin
    hl_acc[a * n_elements(tk_g) + t] = acc_g[a] & hl_tk[a * n_elements(tk_g) + t] = tk_g[t]
  endfor
  firnice_perm_depth = dblarr(nb) & firnice_hl1980 = 'y'
  @procedures/initialise/firnice_perm_depth.pro
  for i = 0, nb-1 do printf, 2, hl_acc[i], hl_tk[i], firnice_perm_depth[i], fo='(f8.3, f9.2, f12.6)'
endif

if part eq 'latent' then begin
  k_ice_b = 0d & k_ice_a = 2.1d
  noval = -9999 & nb = 1
  gl = [3500d] & elev = [3500d] & elev0 = elev & hclim = 3500d
  thick = [200d] & slope = [5d]
  firn_th = [1d] & firn = [1d]
  sno = [0d] & mel = [double(getenv('HARNESS_W'))] & plg = [0d] & tgs = [-10d]
  tgs_th = tgs & mel_th = mel & plg_th = plg
  dtdz = dblarr(12) - 0.0065d & dtdz_icetemp = dtdz & m = 7 & ye = 0
  enable_advection = 'n' & advection_write = 'n' & use_flow_model_gl = 'n'
  enable_strain_heating = 'n' & firnice_melt_surface = 'n' & firnice_batch = 'n'
  firnice_write = ['n', 'n', 'n'] & fit_lacc_on = 0b & firnice_implicit = 'y'
  geothermal_flux = 0d
  firnice_adv_scale_b = [1d] & firnice_insul_scale_b = [1d]
  firnice_refreeze_frac_b = [double(getenv('HARNESS_F'))]
  firnice_perm_depth = [30d] & firnice_perm_frac_b = [1d]
  tl_fit = dblarr(1, total(fit_layers) + 1) - 10d & te_fit = tl_fit
  t_before = reform(tl_fit[0, *])
  rf_dt = 1d-12
  @procedures/processing/firnice_temperature_model.pro
  printf, 2, perm_limit, Lh_fit, fo='(i6, e24.16)'
  for j = 0, perm_limit + 1 do printf, 2, j, fit_dz[0, j], fit_dz[1, j], cap_fit[j], t_before[j], tl_fit[0, j], $
    fo='(i6, 2f10.3, e24.16, 2f22.15)'
endif

close, 2
exit
end
