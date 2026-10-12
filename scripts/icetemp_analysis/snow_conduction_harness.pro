; One ice band run to steady state by GloGEM's own firnice_temperature_model.pro, for the
; snow conduction test (snow_conduction_check.py computes the exact answer).
;
; Run from the GloGEM root:
;   GLOGEM_CONFIG=<any training config> HARNESS_SNO=0.3 HARNESS_IMPLICIT=y HARNESS_OUT=<file> \
;     idl -e ".r <path>/snow_conduction_harness.pro"
;
; Surface fixed at -30 C, geothermal flux 0.3 W/m2, constant ice conductivity (k_ice_b = 0),
; a fixed snowpack of HARNESS_SNO m w.e., no melt, no advection. Writes node depth and
; temperature of the resolved column.

@procedures/initialise/settings.pro

k_ice_b = 0d & k_ice_a = 2.1d          ; constant ice conductivity: exact steady answer exists
noval = -9999 & nb = 1
gl = [2615d] & elev = [2615d] & elev0 = elev & hclim = 2615d
thick = [29d] & slope = [5d]
firn_th = [0d] & firn = [0d]
sno = [double(getenv('HARNESS_SNO'))] & mel = [0d] & plg = [0d]
tgs = [-30d] & dtdz = dblarr(12) - 0.0065d & dtdz_icetemp = dtdz & m = 1 & ye = 0
tgs_th = tgs & mel_th = mel & plg_th = plg   ; the thermal model's energy inputs
enable_advection = 'n' & advection_write = 'n' & use_flow_model_gl = 'n'
enable_strain_heating = 'n' & firnice_melt_surface = 'n' & firnice_batch = 'n'
firnice_write = ['n', 'n', 'n'] & fit_lacc_on = 0b
firnice_implicit = getenv('HARNESS_IMPLICIT')
geothermal_flux = 0.3d
firnice_adv_scale_b = [1d] & firnice_insul_scale_b = [1d] & firnice_refreeze_frac_b = [1d]
firnice_perm_depth = [30d] & firnice_perm_frac_b = [1d]
tl_fit = dblarr(1, total(fit_layers) + 1) - 20d
te_fit = tl_fit

for it_h = 0l, 11999l do begin
@procedures/processing/firnice_temperature_model.pro
endfor

openw, 1, getenv('HARNESS_OUT')
for j = 0, tt-1 do printf, 1, fit_dz[1,j], tl_fit[0,j], fo='(f8.2, f14.8)'
close, 1
exit
end
