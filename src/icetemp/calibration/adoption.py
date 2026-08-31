"""
Adoption decision for the POLYTHERMAL-CLASSIFICATION deliverable: a two-stage
gate-then-bar rule.

This is deliberately SEPARATE from validation.Validator.summary's existing RMSE rule ("adopt KO
only if it beats both baselines on leave-one-out temperature RMSE"). The two answer different
questions about different deliverables, and conflating them would hide a real result:

  RMSE rule       -> "is KO good enough to predict ice TEMPERATURE at a point?"
  this rule       -> "is KO good enough to say WHERE polythermal glaciers are?" (the hazard
                     deliverable -- a coarser, categorical target)

Confirmed empirically across four real campaigns that these genuinely disagree: the 245-point
"Option C" campaign is simultaneously the WORST method on LOO RMSE (4.687 degC, losing to both
baselines) and the BEST on classification (62.5% accuracy, 40% polythermal recall, beating both
baselines on both). Reporting only the RMSE verdict would have thrown away a fit that is the
best available tool for the deliverable the project actually cares about.

WHY A GATE AT ALL, rather than just comparing classification scores: because a classification
score alone can be passed by a fit that is known-broken for unrelated reasons. Confirmed
against real campaign history -- the 145-point theta-expansion campaign produced a posterior at
dT_scale~3.93 with ZERO real training support within any reasonable radius (nearest real design
point 1.17 standardized units away; a trust-gate fallback bug, since fixed, let it through) and
its LOO RMSE nearly doubled -- yet it would have PASSED a classification-only bar (40% recall,
beating both baselines), because its systematic cold bias happened to align with this sample's
dominant classes. The gate exists specifically to stop that: a fit must first be shown
trustworthy on its own terms before its classification score is allowed to mean anything.

THRESHOLD PROVENANCE (important, and a real limitation): GATE_R_HAT_MAX, GATE_DENSITY_RADIUS
and GATE_MIN_NEARBY below are a REASONED STARTING POINT chosen against only three real historical
campaigns, not statistically fitted thresholds -- there is no way around that with this little
history. They separate the three known cases cleanly (see the table in check_gate's docstring),
which is evidence they are not absurd, not evidence they are optimal. Revisit once more
campaigns accumulate rather than treating them as settled.
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Split-chain R-hat ceiling. 1.05 is the conventional "chains have mixed" bar (stricter than
# the older 1.1 rule of thumb, looser than the 1.01 some modern references prefer) -- chosen so
# that a genuinely converged campaign passes comfortably while the two historically
# non-converged ones (R-hat 1.16-1.24) do not.
GATE_R_HAT_MAX = 1.05

# Local design-density check: how many REAL training points lie within GATE_DENSITY_RADIUS
# (standardized theta-space distance) of the posterior mean. This is a direct, self-contained
# measurement -- deliberately NOT a reuse of emulator.max_trusted_distance_, whose value depends
# on calibrate_variance's own fallback path and was itself the source of a real historical bug
# (see that method's "FALLBACK WHEN NO JUMP IS FOUND" docstring). A count of real neighbours
# cannot be inflated by a degenerate LOO-CV curve the way a fitted radius can.
GATE_DENSITY_RADIUS = 0.5
GATE_MIN_NEARBY = 5

# Likelihood-calibration bounds on Mahalanobis/n (see calibrator.mahalanobis_per_dof). ~1 means
# the residuals are the size the uncertainty budget claims. Campaign 7 scored 13.0 at its own
# posterior mode -- over-sharp by an order of magnitude, which is why its credible intervals came
# out 6-28x too narrow and why it excluded the whole high-recall region of theta-space with false
# confidence. [0.5, 2.0] admits an honest fit while refusing both over-sharp (>2) and
# absurdly-slack (<0.5, i.e. the uncertainty is padded until nothing can fail) covariance models.
MAHALANOBIS_BOUNDS = (0.5, 2.0)

# How far the MCMC's own mode may sit below an independent grid/multi-start search, in log-units.
# Campaign 7 reported r_hat 1.00 and ESS 15533 while sitting 223 log-units BELOW a gate-accepted
# grid point -- R-hat certifies that walkers agree with each other, NOT that they found the
# global mode, and all 32 started in a 1e-2*prior_std ball around one MAP. 5.0 is roughly the
# point beyond which a likelihood ratio stops being explicable as sampling noise.
GATE_MAX_MODE_GAP = 5.0


def local_design_density(emulator, theta, radius=GATE_DENSITY_RADIUS):
    """Number of real emulator training points within `radius` standardized units of `theta`.

    Standardization matches Emulator._nn_distance exactly (per-dimension, by theta_train_'s own
    std, with zero-variance dimensions treated as 1.0) so distances here mean the same thing
    they do everywhere else in this package.
    """
    theta_train = np.asarray(emulator.theta_train_, dtype=float)
    # theta_train_ holds only the emulator's ACTIVE columns (Emulator.active_params), which is
    # fewer than 3 whenever a parameter is held fixed -- callers still pass a full 3-vector, so
    # slice it the same way the emulator does rather than subtracting mismatched shapes.
    query = np.asarray(theta, dtype=float).reshape(1, -1)
    if hasattr(emulator, '_slice') and query.shape[1] != theta_train.shape[1]:
        query = emulator._slice(query)
    theta_std = theta_train.std(axis=0)
    theta_std[theta_std == 0] = 1.0
    dist = np.sqrt((((theta_train - query[0]) / theta_std) ** 2).sum(axis=1))
    return int((dist <= radius).sum())


@dataclass
class GateResult:
    """Stage 1: is this fit trustworthy ON ITS OWN TERMS, before any score is consulted?"""
    passed: bool
    r_hat_max: Optional[float]      # None = not available (see r_hat_ok)
    r_hat_ok: Optional[bool]         # None = COULD NOT BE EVALUATED, which never passes
    n_nearby: int
    density_ok: bool
    radius: float
    mahalanobis: Optional[float] = None       # None = not supplied -> NEVER passes
    mahalanobis_ok: Optional[bool] = None
    mode_gap: Optional[float] = None           # log-units below an independent grid best
    mode_gap_ok: Optional[bool] = None
    reasons: list = field(default_factory=list)

    def describe(self):
        lines = [f"  gate: {'PASS' if self.passed else 'FAIL'}"]
        if self.r_hat_ok is None:
            lines.append('    R-hat:   UNAVAILABLE (not persisted by step 4 -- cannot evaluate)')
        else:
            lines.append(f"    R-hat:   {self.r_hat_max:.3f} (max over free params) "
                          f"{'<=' if self.r_hat_ok else '>'} {GATE_R_HAT_MAX} -- "
                          f"{'ok' if self.r_hat_ok else 'NOT CONVERGED'}")
        lines.append(f"    support: {self.n_nearby} real training point(s) within {self.radius} "
                      f"standardized units of the posterior mean "
                      f"{'>=' if self.density_ok else '<'} {GATE_MIN_NEARBY} -- "
                      f"{'ok' if self.density_ok else 'POSTERIOR IS EXTRAPOLATING'}")
        if self.mahalanobis_ok is None:
            lines.append('    calib:   UNAVAILABLE (Mahalanobis/n not supplied -- cannot evaluate)')
        else:
            lo, hi = MAHALANOBIS_BOUNDS
            lines.append(f"    calib:   Mahalanobis/n = {self.mahalanobis:.3f} "
                          f"({'inside' if self.mahalanobis_ok else 'OUTSIDE'} [{lo}, {hi}]) -- "
                          f"{'ok' if self.mahalanobis_ok else 'LIKELIHOOD IS MIS-CALIBRATED'}")
        if self.mode_gap_ok is None:
            lines.append('    mode:    UNAVAILABLE (no independent grid search -- cannot evaluate)')
        else:
            lines.append(f"    mode:    posterior mode sits {self.mode_gap:.1f} log-units below "
                          f"the grid best (<= {GATE_MAX_MODE_GAP}) -- "
                          f"{'ok' if self.mode_gap_ok else 'MCMC MISSED THE GLOBAL MODE'}")
        return '\n'.join(lines)


def check_gate(emulator, theta_hat, r_hat_max=None, radius=GATE_DENSITY_RADIUS,
                min_nearby=GATE_MIN_NEARBY, r_hat_ceiling=GATE_R_HAT_MAX,
                mahalanobis=None, mode_gap=None,
                mahalanobis_bounds=MAHALANOBIS_BOUNDS, max_mode_gap=GATE_MAX_MODE_GAP):
    """Stage 1 of the adoption rule. Both conditions must hold:

      1. R-hat <= r_hat_ceiling for every free parameter (did the sampling converge?)
      2. >= min_nearby real training points within `radius` of the posterior mean
         (does the location it converged to actually have real evidence near it?)

    Both are kept even though, on the three campaigns available so far, the density check alone
    separates them -- they detect genuinely different failures. R-hat catches "the search never
    settled"; density catches "the search settled confidently somewhere we have no data". A
    future campaign could plausibly fail one while passing the other.

    r_hat_max=None means the value was not available (e.g. a campaign predating diagnostics
    persistence). This NEVER passes: an unverifiable gate is not a passed gate. Being unable to
    check convergence is a reason to withhold adoption, not to assume the best.

    Verified against real campaign history:
      campaign 2 (25-glacier, 100pt)   R-hat 1.16-1.18, 0 nearby -> FAIL (both conditions)
      campaign 3 (theta-exp, 145pt)    R-hat 1.00,      0 nearby -> FAIL (density only)
      campaign 4 (Option C, 245pt)     R-hat 1.00-1.01, 10 nearby -> PASS
    Campaign 3 is the case that justifies the whole gate: it would have passed a
    classification-only bar (see module docstring).
    """
    n_nearby = local_design_density(emulator, theta_hat, radius=radius)
    density_ok = n_nearby >= min_nearby
    r_hat_ok = None if r_hat_max is None else bool(r_hat_max <= r_hat_ceiling)
    lo, hi = mahalanobis_bounds
    mahalanobis_ok = (None if mahalanobis is None or not np.isfinite(mahalanobis)
                      else bool(lo <= mahalanobis <= hi))
    mode_gap_ok = None if mode_gap is None else bool(mode_gap <= max_mode_gap)

    reasons = []
    if r_hat_ok is None:
        reasons.append('R-hat unavailable (cannot verify convergence)')
    elif not r_hat_ok:
        reasons.append(f'R-hat {r_hat_max:.3f} > {r_hat_ceiling} (chains did not converge)')
    if not density_ok:
        reasons.append(f'only {n_nearby} real training point(s) within {radius} of the '
                        f'posterior mean (< {min_nearby}); the fit is extrapolating')
    if mahalanobis_ok is None:
        reasons.append('Mahalanobis/n unavailable (cannot verify the likelihood is calibrated)')
    elif not mahalanobis_ok:
        reasons.append(
            f'Mahalanobis/n {mahalanobis:.2f} outside [{lo}, {hi}] -- the likelihood is '
            f"{'OVER-SHARP' if mahalanobis > hi else 'over-slack'} and its credible intervals "
            f'are not trustworthy')
    if mode_gap_ok is None:
        reasons.append('no independent grid search (cannot verify the MCMC found the global mode)')
    elif not mode_gap_ok:
        reasons.append(f'the posterior mode sits {mode_gap:.1f} log-units below an independent '
                        f'grid best (> {max_mode_gap}); the MCMC did not find the global mode')

    return GateResult(
        passed=bool(r_hat_ok) and density_ok and bool(mahalanobis_ok) and bool(mode_gap_ok),
        r_hat_max=r_hat_max, r_hat_ok=r_hat_ok, n_nearby=n_nearby, density_ok=density_ok,
        radius=radius, mahalanobis=mahalanobis, mahalanobis_ok=mahalanobis_ok,
        mode_gap=mode_gap, mode_gap_ok=mode_gap_ok, reasons=reasons,
    )


# ── the classification bar ────────────────────────────────────────────────────────────
# WARM = the ice column contains non-cold ice somewhere. Both 'polythermal' and 'temperate'
# qualify: for the hazard deliverable ("where can we find polythermal-type glaciers") what
# matters is the PRESENCE of non-cold ice, and the temperate/polythermal split is a finer
# distinction the flowline model is not well placed to adjudicate at a single centreline column.
WARM_CLASSES = ('polythermal', 'temperate')

# ABSOLUTE FLOOR, both halves. Set 2026-08-20 by Janosch: "it would be good if both no. are well
# above 50 otherwise it is just like flipping a coin". Beating the baselines is NOT sufficient on
# its own -- the baselines are themselves weak, and a method can beat them while still being
# useless in absolute terms.
BAR_MIN_RECALL = 0.5        # of truly-warm glaciers, how many did we find
BAR_MIN_SPECIFICITY = 0.5   # of truly-cold glaciers, how many did we correctly leave alone


@dataclass
class BarResult:
    """Stage 2: does KO actually beat both baselines at the classification job, AND clear an
    absolute floor on both halves of the problem?"""
    passed: bool
    recall: dict            # method -> warm glaciers correctly found (sensitivity)
    specificity: dict       # method -> cold glaciers correctly left alone
    balanced: dict          # method -> (recall + specificity) / 2  <- the headline score
    accuracy: dict          # method -> plain accuracy. REPORTED ONLY, never part of the rule.
    polythermal_recall: dict  # method -> recall on strictly-'polythermal' truth, for continuity
    n_glaciers: int
    n_warm: int
    n_cold: int
    n_polythermal: int
    trivial_accuracy: float   # what "call everything cold" scores -- the honest baseline
    suffix: str
    floor_ok: dict

    def describe(self):
        lines = [f"  bar:  {'PASS' if self.passed else 'FAIL'}   [graded on '{self.suffix}']"]
        lines.append(f"    n={self.n_glaciers} ({self.n_warm} warm / {self.n_cold} cold; "
                      f"{self.n_polythermal} strictly polythermal)")
        lines.append(f"    BALANCED score (the rule):    " + ', '.join(
            f'{m}={self.balanced[m]:.1%}' for m in ('tier2', 'knn', 'ko')))
        lines.append(f"      warm found (recall):        " + ', '.join(
            f'{m}={self.recall[m]:.1%}' for m in ('tier2', 'knn', 'ko')))
        lines.append(f"      cold kept (specificity):    " + ', '.join(
            f'{m}={self.specificity[m]:.1%}' for m in ('tier2', 'knn', 'ko')))
        lines.append(f"    floor (both >{BAR_MIN_RECALL:.0%}): ko "
                      f"{'CLEARS' if self.floor_ok.get('ko') else 'FAILS'}")
        lines.append(f"    plain accuracy (context only):" + ', '.join(
            f' {m}={self.accuracy[m]:.1%}' for m in ('tier2', 'knn', 'ko'))
            + f'  [call-everything-cold scores {self.trivial_accuracy:.1%}]')
        return '\n'.join(lines)


def check_bar(results, methods=('tier2', 'knn', 'ko'), suffix='_class_real'):
    """Stage 2 of the adoption rule. KO must BOTH:

      1. clear an ABSOLUTE FLOOR -- recall > BAR_MIN_RECALL and specificity > BAR_MIN_SPECIFICITY;
      2. beat BOTH baselines on the BALANCED score, (recall + specificity) / 2.

    WHY BALANCED AND NOT PLAIN ACCURACY. Measured on the campaign-7 entity set: 83 of 116
    labelled entities are COLD, so a model that simply calls everything cold scores 71.6%
    accuracy while finding zero warm glaciers. Plain accuracy is therefore already above 50% for
    free and cannot express "better than a coin flip" at all; the best theta found in the
    2026-08-20 grid sweep scores 73.8% accuracy against that 71.6% trivial baseline -- a
    meaningless-looking 2-point gain -- while its balanced score is 74.5% against 50%. Accuracy
    is still REPORTED, next to the trivial baseline, so the imbalance is visible; it is not part
    of the rule. (The previous version of this docstring asserted the mix was "roughly half
    cold", which was simply wrong and is what made plain accuracy look defensible.)

    Recall alone remains trivially gamed by calling everything warm -- specificity is what stops
    that, and pairing them is why the balanced score is the right single number here.

    WHY suffix='_class_real'. `{m}_class` is written from physics.cp_model_single, the ANALYTICAL
    SURROGATE, while `{m}_class_real` is written from the EMULATOR -- the forward model KO's theta
    is actually fit against. Grading on the surrogate measures something the likelihood does not
    control, and the two genuinely disagree: on campaign 7 the LOO RMSE ordering INVERTS between
    them (ko worst at 3.682 in surrogate space, best at 2.033 in real space). Until 2026-08-20
    this function read `{m}_class` and `{m}_class_real` was written by validation.py and read
    NOWHERE. Pass suffix='_class' only to reproduce a pre-2026-08-20 campaign's verdict.
    """
    col = {m: f'{m}{suffix}' for m in methods}
    missing = [c for c in col.values() if c not in results.columns]
    if missing:
        raise KeyError(
            f'check_bar: missing column(s) {missing}. Grading silently on a different column is '
            f'exactly the defect this argument exists to prevent -- pass suffix explicitly.')

    truth = results['obs_class']
    is_warm = truth.isin(WARM_CLASSES)
    warm, cold = results[is_warm], results[~is_warm]
    poly = results[truth == 'polythermal']
    n = len(results)

    def _frac(sub, c, want_warm):
        if not len(sub):
            return float('nan')
        pred_warm = sub[c].isin(WARM_CLASSES)
        return float(pred_warm.mean() if want_warm else (~pred_warm).mean())

    recall = {m: _frac(warm, col[m], True) for m in methods}
    specificity = {m: _frac(cold, col[m], False) for m in methods}
    balanced = {m: (recall[m] + specificity[m]) / 2.0 for m in methods}
    accuracy = {m: float((results[col[m]] == truth).mean()) if n else float('nan')
                for m in methods}
    poly_recall = {m: float((poly[col[m]] == 'polythermal').mean()) if len(poly) else float('nan')
                   for m in methods}
    floor_ok = {m: bool(recall[m] > BAR_MIN_RECALL and specificity[m] > BAR_MIN_SPECIFICITY)
                for m in methods}

    passed = (
        floor_ok['ko']
        and balanced['ko'] > balanced['tier2']
        and balanced['ko'] > balanced['knn']
    )
    return BarResult(
        passed=bool(passed), recall=recall, specificity=specificity, balanced=balanced,
        accuracy=accuracy, polythermal_recall=poly_recall, n_glaciers=n,
        n_warm=int(is_warm.sum()), n_cold=int((~is_warm).sum()), n_polythermal=len(poly),
        trivial_accuracy=float((~is_warm).mean()) if n else float('nan'),
        suffix=suffix, floor_ok=floor_ok,
    )


@dataclass
class AdoptionDecision:
    adopt: bool
    gate: GateResult
    bar: Optional[BarResult]   # None when the gate failed -- the bar is deliberately not
                                # evaluated, so a failed fit's score never gets quoted as if
                                # it meant something.

    def describe(self):
        head = ('CLASSIFICATION ADOPTION DECISION (polythermal deliverable -- separate from the '
                'RMSE rule above):')
        body = [self.gate.describe()]
        if self.bar is None:
            body.append('  bar:  not evaluated (gate failed -- see reasons below)')
        else:
            body.append(self.bar.describe())
        verdict = f"  => {'ADOPT' if self.adopt else 'DO NOT ADOPT'} KO for polythermal classification"
        if not self.adopt and self.gate.reasons:
            verdict += '\n     because: ' + '; '.join(self.gate.reasons)
        elif not self.adopt and self.bar is not None and not self.bar.passed:
            why = []
            if not self.bar.floor_ok.get('ko'):
                why.append(f"KO is below the absolute floor (needs recall AND specificity "
                            f">{BAR_MIN_RECALL:.0%}; got {self.bar.recall['ko']:.1%} / "
                            f"{self.bar.specificity['ko']:.1%})")
            if not (self.bar.balanced['ko'] > self.bar.balanced['tier2']
                    and self.bar.balanced['ko'] > self.bar.balanced['knn']):
                why.append('KO did not beat both baselines on the balanced score')
            verdict += '\n     because: ' + '; '.join(why)
        return '\n'.join([head] + body + [verdict])


def decide(emulator, theta_hat, results, r_hat_max=None, suffix='_class_real', **gate_kwargs):
    """Full gate-then-bar decision. The bar is evaluated ONLY if the gate passes.

    `suffix` selects which classification column the bar grades -- see check_bar. The default
    grades the EMULATOR column, the forward model theta is actually fit against.
    """
    gate = check_gate(emulator, theta_hat, r_hat_max=r_hat_max, **gate_kwargs)
    if not gate.passed:
        return AdoptionDecision(adopt=False, gate=gate, bar=None)
    bar = check_bar(results, suffix=suffix)
    return AdoptionDecision(adopt=bar.passed, gate=gate, bar=bar)
