"""Size of the emulator's own uncertainty E across the parameter box (alps_s2 fixture): inflation
curve, trust radius, and E at a few parameter sets with their distance to the nearest run."""
import dill, numpy as np
from pathlib import Path
OUT = Path('/scratch_net/vierzack04_fourth/jabeer/tmp_icetemp_trees/step5b')
fx = dill.load(open(OUT / 'fixture.dill', 'rb'))
emu = fx['emulator']
print(f'p={emu.p}, LOO z-std {emu.loo_z_std_:.2f} (inflation {emu.variance_inflation_:.1f}x), '
      f'design spacing {emu.design_spacing_:.3f}, trusted distance {emu.max_trusted_distance_:.3f}, '
      f'degenerate curve {emu.inflation_curve_degenerate_}')
xs = emu._variance_inflation_curve.X_thresholds_; ys = emu._variance_inflation_curve.y_thresholds_
print('inflation curve (distance: factor):', ', '.join(f'{x:.2f}:{y:.1f}' for x, y in zip(xs, ys)))
for th in ([0.5245, 1.0893, 1.0], [0.33, 0.75, 1.0], [0.55, 1.0, 1.0], [0.3, 0.55, 1.0], [0.9, 1.8, 1.0]):
    mean, L = emu.predict_factor(th)
    ev = np.linalg.eigvalsh(L.T @ L)[::-1]
    nn = emu._nn_distance(np.atleast_2d(th))[0]
    print(f'theta {th}: nn distance {nn:.3f}, E trace/row {np.sum(L**2)/L.shape[0]:.3f} C^2, '
          f'top eigenvalues {np.round(ev[:3], 2)}')
