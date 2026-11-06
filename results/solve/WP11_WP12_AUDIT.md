# Auditing the (s, φ) discriminator: assumptions, recovery, and where it breaks

Response to the four requests: (1) audit the derivation, (2) test parameter recovery,
(3) propagate uncertainty, (4) establish what one independent rate measurement buys.

Scripts: `wp11_audit_identity.py`, `wp12_recovery.py`. Outputs: `wp11_audit.csv`, `wp12_recovery.csv`.

---

## 1. The audit — the invariance survives, one estimator does not

Four assumptions were at issue.

| assumption | relaxed how | does the identity survive? |
|---|---|---|
| A1 six polymerization steps | n_eff unknown | invariance **yes**; point estimate carries a 2× systematic |
| A2 equal k_pol at every step | six different k_pol | invariance **yes**; point estimate biased up to 2× |
| A3 single k_off | k_off differs between extension and translocation | invariance **yes**, and the result gets *stronger* |
| A4 dissociation at every step | — | invariance **yes** |

**The load-bearing claim is `d/dlog(k_tr)[log(s) − log(φ)] = 0`, and it survives every relaxation
tested.** That is what makes the discriminator usable: it does not depend on any rate being
constant across steps.

**What does not survive** is the literal inversion `k_off/k_pol = (s/φ)^(−1/6) − 1`.

- With heterogeneous k_pol the identity becomes `log(s/φ) = −Σᵢ log(1 + k_off/k_pol,ᵢ)`, so the
  inversion returns an **aggregate**, not a rate ratio. Measured bias when a single step is 10×
  slower: the naive estimate is off by a factor of ~2.
- With unknown n_eff, the estimate scales as `(s/φ)^(−1/n_eff)`: at s/φ = 0.99 the implied ratio
  runs from 0.00252 (n=4) to 0.00126 (n=8).

This is the same distinction you drew — identifies a ratio, not a mechanism — but it is now
quantified: the classification is reliable, the number is not.

### A3 is the useful surprise

Real heterogeneity in k_off is not a hypothetical here: the entire empirical claim under test is
that TPP1 leaves the *resting* k_off at 0.012 min⁻¹ while making dissociation during translocation
≥4× slower. Those are two different k_off values by construction. When the scheme is written that
way:

```
log(s/φ)          = Σᵢ log[k_pol,ᵢ/(k_pol,ᵢ + k_off_ext)]     depends only on EXTENSION
    d/dlog(k_tr)  = 0
    d/dlog(k_off_tr) = 0

φ/(1−φ)           = k_tr / k_off_tr                          depends only on TRANSLOCATION
```

The identity **splits** into two phase-specific ratios:

| channel | definition | reports on |
|---|---|---|
| **A = s/φ** | `Π k_pol,ᵢ/(k_pol,ᵢ + k_off_ext)` | retention during **extension** |
| **B = φ/(1−φ)** | `k_tr / k_off_tr` | retention during **translocation** |

That is strictly more informative than wp10 stated, and it is exactly the decomposition the
question needs — "does TPP1 hold DNA better" is two different claims depending on which phase.

---

## 2. Recovery — and a failure that is more useful than a pass

I first ran the recovery test with `k_pol = 2.0 min⁻¹`. **The classifier scored near zero.** That
was not a bug; it exposed the design's real weakness:

```
k_pol    k_off/k_pol    s/φ      after 4× k_off drop    ratio     detectable?
0.020       0.6000    0.0596          0.4323           7.253        yes
0.050       0.2400    0.2751          0.7050           2.563        yes
0.100       0.1200    0.5066          0.8375           1.653        yes
0.200       0.0600    0.7050          0.9145           1.297        yes
0.500       0.0240    0.8674          0.9647           1.112        yes
1.000       0.0120    0.9309          0.9822           1.055      marginal
2.000       0.0060    0.9647          0.9910           1.027         no
```

**The s/φ channel only has power when `k_off` is an appreciable fraction of `k_pol` — and
`k_off/k_pol` is precisely the unknown the design is trying to measure.** If k_pol is fast (≈2
min⁻¹), a 4× reduction in k_off moves s/φ by 2.7%, which no realistic assay resolves.

Two consequences, and both are design-relevant:

**(a) The design is conditional.** At k_pol = 0.10 min⁻¹ every signature is detectable at
optimistic-but-plausible precision (CV_s 2%, 6 replicates, 2000 molecules):

| mechanism | signature | effect | SE needed | SE achievable | detectable |
|---|---|---|---|---|---|
| retention | s/φ changes | 65.3% | 16.5% | 0.9% | yes |
| translocation | φ changes | 6.6% | 1.7% | 0.9% | yes |
| catalytic | s/φ changes | 47.3% | 11.9% | 0.9% | yes |
| combined | s/φ changes | 24.4% | 6.2% | 0.9% | yes |

At k_pol = 2.0 min⁻¹ none are detectable. The published record does not contain k_pol for human
telomerase in this framework, so **the design's statistical power is unknown until it is measured —
as the first experiment, not the last.**

**(b) The design is self-diagnosing, which rescues it.** The *baseline* value of s/φ is itself
informative:

| observed s/φ | implies k_off/k_pol | reading |
|---|---|---|
| ≈ 0.99 | < 0.002 | extension-phase retention negligible — TPP1 cannot be acting there; channel B is the only kinetic route left |
| ≈ 0.95 | ≈ 0.009 | channel has ~1.3–1.6× power per 4× k_off |
| ≤ 0.90 | > 0.018 | strong power; extension-phase retention is a live hypothesis |

So an s/φ near 1 is not a null result — it is a **positive exclusion** of extension-phase
retention, and it redirects the question to translocation. That is a better design than one that
can only return "significant / not significant".

---

## 3. Uncertainty — your point about φ → 1 is worse than stated, and it matters less

`B = φ/(1−φ)` has `dB/B = dφ/[φ(1−φ)]`, which diverges:

| φ | B | dφ = 0.01 | dφ = 0.005 | dφ = 0.002 |
|---|---|---|---|---|
| 0.50 | 1.0 | 4.0% | 2.0% | 0.8% |
| 0.80 | 4.0 | 6.3% | 3.1% | 1.3% |
| 0.90 | 9.0 | 11.1% | 5.6% | 2.2% |
| 0.95 | 19.0 | 21.1% | 10.5% | 4.2% |
| 0.99 | 99.0 | **101%** | 50.5% | 20.2% |

At the published φ = 0.99, a one-point error in φ doubles the implied ratio. **But the sampling is
not the problem** — binomial noise is cheap:

| φ | molecules for SE(dφ) = 0.005 |
|---|---|
| 0.90 | ~3,600 |
| 0.95 | ~1,900 |
| 0.99 | ~400 |

The binding constraint is instrumental, not statistical: whether the assay can resolve 0.99 from
0.995 at all. If it cannot, channel **B must be reported as a bound, not an estimate**, and the
analysis should lean on channel A — which does not saturate.

---

## 4. What one independent rate buys — asymmetry worth knowing

With (A, B) there are two ratios among four rates. Adding one absolute rate does **not** resolve
both phases symmetrically:

- **Extension phase: fully resolved.** With k_off_ext measured directly, k_pol follows exactly from
  A. Verified numerically at both scales, zero error.
- **Translocation phase: not resolved.** B = k_tr/k_off_tr still needs one of the two, and neither
  is easy to measure separately. φ is therefore almost always reported alone — which is why the
  translocation claim is **the least constrained part of this literature, not the best**.

Recommendation: spend the one extra measurement on **k_off_ext** (direct primer binding, resting
state), not on translocation. It converts channel A from a ratio into three absolute rates and
makes the extension-phase question closeable.

---

## 5. Revised decision table

| observation with TPP1 | supports | does not support |
|---|---|---|
| φ flat, s/φ changes | catalytic (k_pol up) or extension k_off down | which of the two |
| φ changes, s/φ flat | translocation only (k_tr/k_off_tr up) | that k_tr rose rather than k_off_tr falling |
| both change | retention and/or combined | a unique mechanism |
| neither changes | no change in these ratios | no TPP1 effect elsewhere |
| s/φ ≈ 0.99 at baseline | extension-phase retention is negligible | — (this is a positive exclusion) |

---

## 6. Status of the whole claim

| claim | status after this audit |
|---|---|
| `d/dlog(k_tr)[log(s/φ)] = 0` | **holds exactly**, robust to A1–A5 |
| split into extension / translocation channels | **holds**, and is an improvement on wp10 |
| `k_off/k_pol` from the inversion is a rate ratio | **fails** — aggregate, ±2× |
| the design has power at published effect sizes | **unknown** — depends on unmeasured k_pol |
| the design is self-diagnosing | **holds** — baseline s/φ tells you which regime you are in |

The honest summary: **the discriminator classifies mechanism class reliably and yields no
trustworthy rate.** That is what you said, and the audit now supports it quantitatively rather than
by assertion.
