# Computational & Biophysical Methodology (`telomerase-allosteric-switch`)

## 1. All-Atom Molecular Dynamics & Replica-Level Statistical Inference

### Simulation Setup (OpenMM 8.6.1)
All-atom molecular dynamics simulations were conducted using the OpenMM 8.6.1 biophysical simulation engine on the core TEN–TRAP–TPP1 junction:
- **Bound System (`junction_bound_fixed.pdb`)**: TERT TEN domain (residues 165–178), TERT TRAP domain (residues 750–795), and TPP1 TEL-patch segment (residues 165–185 and 210–225). Total: 1,547 protonated atoms.
- **Apo System (`junction_apo_fixed.pdb`)**: TERT TEN domain (residues 165–178) and TRAP domain (residues 750–795). Total: 944 protonated atoms.
- **Force Field**: AMBER `ff14SB` with Onufriev-Bashford-Case (OBC2) generalized Born implicit solvent model (`implicit/obc2.xml`).
- **Integration**: Langevin Middle Integrator ($T = 300\text{ K}$, friction coefficient $\gamma = 1.0\text{ ps}^{-1}$, integration time step $\Delta t = 2.0\text{ fs}$).
- **Equilibration & Production**: 150 steps L-BFGS energy minimization, followed by 1,500 integration steps ($3.0\text{ ps}$). A 300-step burn-in exclusion window ($0.6\text{ ps}$) was applied prior to production analysis.

### Avoiding Pseudoreplication: Replica-Level Statistics
Treating consecutive trajectory frames as independent data points constitutes pseudoreplication because adjacent timeframes exhibit high autocorrelation. Therefore, **independent simulation replicas ($N = 3$) initiated with distinct pseudorandom velocity seeds serve as the true statistical units**:

1. **Replicate Means**:
   $$\bar{x}_k = \frac{1}{M} \sum_{m=1}^M x_{k,m} \quad (\text{for replicate } k \in \{1, 2, 3\}, M = 60 \text{ production frames})$$
2. **Ensemble Statistics**:
   $$\mu = \frac{1}{N} \sum_{k=1}^N \bar{x}_k, \quad s^2 = \frac{1}{N - 1} \sum_{k=1}^N (\bar{x}_k - \mu)^2$$
3. **Two-Sample Independent t-Test ($df = 4$)**:
   $$t = \frac{\mu_{\text{bound}} - \mu_{\text{apo}}}{s_{\text{pooled}} \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}}, \quad s_{\text{pooled}} = \sqrt{\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}$$
4. **Hedges' $g$ Effect Size (Bias-Corrected for Small $N$)**:
   $$g = d \times \left(1 - \frac{3}{4(n_1 + n_2) - 9}\right), \quad d = \frac{\mu_{\text{bound}} - \mu_{\text{apo}}}{s_{\text{pooled}}}$$

---

## 2. Surface Plasmon Resonance (SPR) 1:1 Langmuir Kinetics

SPR binding sensorgrams measure real-time macromolecular interaction between immobilized TERT holoenzyme and injected TPP1 OB-fold analyte:

1. **Association Phase ($0 \le t \le t_{\text{inj}}$)**:
   $$R(t) = R_{\text{eq}} \left[1 - \exp\left(-(k_{\text{on}} C + k_{\text{off}}) t\right)\right], \quad R_{\text{eq}} = \frac{R_{\text{max}} C}{C + K_d}$$
   where $C$ is analyte concentration, $k_{\text{on}}$ is association rate ($\text{M}^{-1}\text{s}^{-1}$), $k_{\text{off}}$ is dissociation rate ($\text{s}^{-1}$), and $K_d = k_{\text{off}} / k_{\text{on}}$.
2. **Dissociation Phase ($t > t_{\text{inj}}$)**:
   $$R(t) = R_0 \exp\left(-k_{\text{off}} (t - t_0)\right)$$
3. **Mass Transport Limitation (MTL) & Quality Control**:
   - Flow rate: $30\text{ }\mu\text{L/min}$ to minimize diffusion boundary layer resistance.
   - Residual error criterion: $\chi^2 \le 0.10 \times R_{\text{max}}$.
   - Acceptance threshold: $K_d \le 300\text{ nM}$ with nanomolar recruitment preserved.

---

## 3. Quantitative RAP Gel Densitometry & Processivity Decay

Human telomerase extends telomeric DNA primers by adding tandem hexameric repeats ($5'\text{-GGTTAG-}3'$). Direct extension assays on 12% Urea-PAGE gels quantify this activity:

1. **Rolling-Ball Background Subtraction**:
   $$\text{Net}_i = \max(0.0, I_{\text{raw},i} - I_{\text{background}})$$
2. **$3i$ Radioactive Guanine Content Normalization**:
   Because each added repeat incorporates three $[\alpha\text{-}^{32}\text{P}]\text{-dGTP}$ molecules, longer products incorporate proportionally more radioactivity. The signal must be divided by $3i$ to reflect true molar stoichiometry:
   $$I_{\text{norm},i} = \frac{\text{Net}_i}{3i \times \text{Control}_{\text{recovery}}}$$
3. **Processivity per Repeat ($R_{p,i}$)**:
   According to Latrick & Cech (2010), the probability of extending beyond repeat $i$ without dissociating is:
   $$R_{p,i} = \frac{\sum_{j=i+1}^N I_{\text{norm},j}}{\sum_{j=i}^N I_{\text{norm},j}}$$
4. **Overall Processivity Factor ($P_{\text{overall}}$) & Fold-Stimulation ($\chi_{\text{RAP}}$)**:
   Log-linear fitting of normalized band intensities yields the decay rate per repeat ($k$):
   $$\ln(I_{\text{norm},i}) = -k \cdot i + C \implies P_{\text{overall}} = \exp(-k)$$
   Processivity stimulation fold is:
   $$\chi_{\text{RAP}} = \frac{P_{\text{overall}}(+\text{TPP1})}{P_{\text{overall}}(-\text{TPP1})}$$
   - Uncoupling Hit Criterion: $\chi_{\text{RAP}}(\text{K757A}) \le 1.8$ with basal catalysis $\ge 85\%$ of WT.
   - Rescue Criterion: $\chi_{\text{RAP}}(\text{K757R}) \ge 2.0$ with basal catalysis $\ge 85\%$ of WT.

---

## 4. Gaussian Network Model (GNM) Coarse-Grained Normal Modes

To compute equilibrium thermal fluctuations across the entire macromolecular complex without atomistic friction:

1. **Kirchhoff Adjacency Matrix ($\Gamma$)**:
   $$\Gamma_{ij} = \begin{cases} -1 & \text{if } i \neq j \text{ and } d_{ij} \le R_c \\ 0 & \text{if } i \neq j \text{ and } d_{ij} > R_c \\ -\sum_{k \neq i} \Gamma_{ik} & \text{if } i = j \end{cases}$$
   where $R_c = 10.0\text{ \AA}$ is the canonical inter-Cα cutoff.
2. **Covariance Matrix & Cross-Correlation**:
   $$C_{ij} = [\Gamma^{-1}]_{ij} = \sum_{k=1}^{N-1} \frac{1}{\lambda_k} \mathbf{v}_k \mathbf{v}_k^T$$
   $$\text{CC}(i,j) = \frac{C_{ij}}{\sqrt{C_{ii} C_{jj}}}$$
   Positive cross-correlation ($\text{CC} > 0.40$) between the TEN domain apex (`Tyr176`) and TRAP pore rim (`Lys757`–`Gln794`) demonstrates collective mechanical coupling upon TPP1 docking.

---

## 5. Identifiability Analysis of Processivity Observables

Added 2026-11-01. This section documents the methods behind
[`paper/MANUSCRIPT.md`](../../paper/MANUSCRIPT.md) and
[`results/causality/MECHANISM_CLASSIFICATION.md`](../results/causality/MECHANISM_CLASSIFICATION.md).
Everything here is verified by `scripts/mechanism_classification_checks.py`, which writes the
machine-readable values that the manuscript quotes, and by `tests/test_mechanism_classification.py`
and `tests/test_paper_manuscript.py`, which fail if the text drifts.

### 5.1 The model class

Per repeat the enzyme either extends or is lost, with primer binding in equilibrium:

$$s = \frac{k_{\text{ext}}}{k_{\text{ext}} + k_{\text{loss}} + R}, \qquad
R = \frac{k_{\text{off}} K_d}{K_d + [P]} = \frac{k_{\text{off}}^2}{k_{\text{off}} + [P]\,k_{\text{on}}}$$

and the observable is $P(\geq n) = s^{\,n-1}$ for a single population. $k_{\text{loss}}$ is an
**effective** rate — the loss that more primer cannot rescue — not a molecular state. The second
form of $R$ is an exact identity, not an approximation, and it is what makes §5.2 work.

### 5.2 Symbolic verification (SymPy)

Each of these is computed, not asserted:

| quantity | how |
|---|---|
| scale invariance | `simplify(s(λθ) − s(θ))` must be exactly `0` |
| B–C2 crossing | solve `s_B = s_C2` for `P`; the result is **degree 1 in P**, so at most one crossing |
| degenerate case | `α = β` numerator must simplify to `α·k_off²(α−1)` |
| loss-term rewrite | `simplify(R − k_off²/(k_off + P·k_on))` must be `0` |
| matched-*K*_d scaling | `simplify(R(λθ) − λR(θ))` must be `0` at fixed *K*_d |
| mixture monotonicity | `dr/du = s(ρ−1)/(1+u)²`, so *r* has no interior extremum |

The displayed formula for $P^*$ is checked twice: against `sympy.solve` and against a 20 001-point
numeric scan that counts sign changes of $s_B - s_{C2}$ in $[P]$. The scan must find **one** crossing
iff $P^* > 0$ and **zero** otherwise. This is not decoration — an earlier draft of the manuscript
displayed the numerator with the opposite sign, which turns every real crossing into a non-crossing.

### 5.3 Identifiability

**Structural.** The Jacobian $J_{ij} = \partial s(P_i)/\partial\theta_j$ is evaluated with
**complex-step** derivatives, $\operatorname{Im}[s(\theta_j + ih)]/h$ with $h = 10^{-20}$, because
finite differences produce unreliable rank here. Rank is declared from the relative singular values
against a $10^{-8}$ cut-off, and is checked to be grid-independent across three primer grids.

**Practical.** A log-space Fisher matrix $F = J^\top W J$ with $W = 1/(0.15\,s)^2$, condition number
$\sqrt{\lambda_{\max}/\lambda_{\min}}$, and relative standard errors from $\mathrm{pinv}(F)$. Two
parametrisations are reported — ratio form (3 parameters) and absolute four-rate form (4) — because
condition numbers from different parametrisations are not comparable, and the manuscript states
which is which in every row.

### 5.4 Mixture ratio and the tail-slope audit

For a two-species mixture the adjacent-band ratio is

$$r(n) = s\,\frac{1+\kappa\rho^{n}}{1+\kappa\rho^{n-1}}, \qquad
\kappa = \frac{\phi(1-s_2)}{(1-\phi)(1-s)}, \quad \rho = \frac{s_2}{s}$$

evaluated in log space ($\log(1+e^x)$ via `logaddexp`) because $\rho^n$ overflows for $n$ in the
hundreds. Three quantities are extracted: the tail slope over a given band window, the implied
$k_{\text{loss}}/k_{\text{ext}} = (1-s)/s$, and the onset length
$n^* = 1 + \ln(1/\kappa)/\ln\rho$ at which the minority component overtakes the majority.

**Standing rule.** The tail slope is **forbidden** as the basis for any mechanistic claim. It is a
function of the read window, not of the distribution: for a true $s = 0.60$ it returns 0.610 over
$n = 3\text{–}10$ and 0.950 over $n = 40\text{–}80$, a 12-fold change in the implied loss ratio from
identical data. Mechanistic conclusions must come from a full-distribution fit with the mixture
weight $\phi$ free.

### 5.5 What these methods cannot do

They cannot establish a biological mechanism, and no result from this section may be cited as
having done so. The analysis is exact *within the model class*; whether the model describes a real
gel is an empirical question, and the test for it is per-band densitometry of an existing gel, not
another simulation.
