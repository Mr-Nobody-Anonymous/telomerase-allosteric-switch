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
