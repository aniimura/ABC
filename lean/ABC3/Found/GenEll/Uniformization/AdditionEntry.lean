/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import ABC3.Found.GenEll.LatticeCurve
import ABC3.Found.GenEll.WeierstrassODE
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Uniformization.Surjective

/-!
# 一様化 —— Laurent の入口・加法定理の欠損関数・`z ≡ −w` の側・零点勘定の第一の煉瓦

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★Laurent の入口——加法定理へ -/

/-- ★★★★★**`℘(z) − 1/z²` は mathlib の `℘[Λ − 0]` そのもの**（原点で解析的）。 -/
theorem weierstrassP_sub_invSq (P : PeriodPair) (z : ℂ) :
    P.weierstrassP z - 1 / z ^ 2 = P.weierstrassPExcept 0 z := by
  have h := P.weierstrassPExcept_add ⟨0, P.lattice.zero_mem⟩ z
  simp only [sub_zero] at h
  rw [← h]
  simp

/-- ★★★★★**`℘′(z) + 2/z³` は `℘′[Λ − 0]`**。 -/
theorem derivWeierstrassP_add_invCube (P : PeriodPair) (z : ℂ) :
    P.derivWeierstrassP z + 2 / z ^ 3 = P.derivWeierstrassPExcept 0 z := by
  have h := P.derivWeierstrassPExcept_sub ⟨0, P.lattice.zero_mem⟩ z
  simp only [sub_zero] at h
  rw [← h]
  ring

/-- ★★★★★★**`z²·℘(z)` の解析接続**——原点で `1`。

★これが `℘(z) = z⁻² + O(z²)` の Lean 上の姿である
（`weierstrassPExcept` は原点で解析的で値 `0`）。 -/
noncomputable def laurentB (P : PeriodPair) (z : ℂ) : ℂ :=
  1 + z ^ 2 * P.weierstrassPExcept 0 z

/-- ★★★★★★**`z³·℘′(z)` の解析接続**——原点で `−2`。 -/
noncomputable def laurentA (P : PeriodPair) (z : ℂ) : ℂ :=
  -2 + z ^ 3 * P.derivWeierstrassPExcept 0 z

@[simp] theorem laurentB_zero (P : PeriodPair) : laurentB P 0 = 1 := by simp [laurentB]

@[simp] theorem laurentA_zero (P : PeriodPair) : laurentA P 0 = -2 := by simp [laurentA]

/-- ★★★★★★`z ≠ 0` では `laurentB P z = z²·℘(z)`。 -/
theorem laurentB_eq (P : PeriodPair) (z : ℂ) (hz : z ≠ 0) :
    laurentB P z = z ^ 2 * P.weierstrassP z := by
  have h := P.weierstrassPExcept_add ⟨0, P.lattice.zero_mem⟩ z
  simp only [sub_zero] at h
  simp only [laurentB, ← h]
  have h0 : (1 : ℂ) / 0 ^ 2 = 0 := by norm_num
  rw [h0]
  field_simp
  ring

/-- ★★★★★★`z ≠ 0` では `laurentA P z = z³·℘′(z)`。 -/
theorem laurentA_eq (P : PeriodPair) (z : ℂ) (hz : z ≠ 0) :
    laurentA P z = z ^ 3 * P.derivWeierstrassP z := by
  simp only [laurentA, ← derivWeierstrassP_add_invCube]
  field_simp
  ring

/-- ★★★★★★★★`laurentA`・`laurentB` は原点で解析的。 -/
theorem analyticAt_laurentB (P : PeriodPair) : AnalyticAt ℂ (laurentB P) 0 := by
  refine analyticAt_const.add (((analyticAt_id).pow 2).mul ?_)
  exact ((P.differentiableOn_weierstrassPExcept 0).analyticOnNhd
    P.isOpen_compl_lattice_sdiff) 0 (by simp)

theorem analyticAt_laurentA (P : PeriodPair) : AnalyticAt ℂ (laurentA P) 0 := by
  refine analyticAt_const.add (((analyticAt_id).pow 3).mul ?_)
  exact ((P.differentiableOn_derivWeierstrassPExcept 0).analyticOnNhd
    P.isOpen_compl_lattice_sdiff) 0 (by simp)

def laurentB.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(z²℘(z) の解析接続——加法定理の Laurent の入口。★無条件)",
    sectionId := "genell-lemma-3-5" }

def laurentA.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(z³℘′(z) の解析接続——加法定理の Laurent の入口。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★`M(z)`——`laurent_cancel` の右辺の因子。 -/
noncomputable def laurentM (P : PeriodPair) (x y z : ℂ) : ℂ :=
  4 * P.weierstrassPExcept 0 z + 8 * (P.weierstrassPExcept 0 z - x)
  + 4 * z ^ 2 * (P.weierstrassPExcept 0 z - x) ^ 2
  + 8 * z ^ 2 * P.weierstrassPExcept 0 z * (P.weierstrassPExcept 0 z - x)
  + 4 * z ^ 4 * P.weierstrassPExcept 0 z * (P.weierstrassPExcept 0 z - x) ^ 2
  + 4 * z * (P.derivWeierstrassPExcept 0 z - y)
  - z ^ 4 * (P.derivWeierstrassPExcept 0 z - y) ^ 2

/-- ★★★★★★★★★★★★★★★★★★★★★★**加法定理の極の打ち消し（原点）**——★純粋な恒等式。

    `4·B·(B − z²·x) − (A − z³·y)² = z² · M(z)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `Skeleton/GenEll/AdditionTheorem.lean`（第 606）の証明の核である。
`x = ℘(w)`・`y = ℘′(w)` と置くと、`z ≠ 0` では `B = z²℘(z)`・`A = z³℘′(z)` なので

    `℘(z) + ℘(w) + ℘(z+w) − (1/4)·((℘′(z)−℘′(w))/(℘(z)−℘(w)))²`
      `= ℘(z+w) + ℘(w) + [4·B·(B − z²℘(w))² − (A − z³℘′(w))²] / (4·z²·(B − z²℘(w))²)`

であり、★**分子が `z²` でくくれる**（本定理）ので `z²` が約されて原点で有界になる。

☆`B(0) = 1`・`A(0) = −2` なので `4·1·1 − 4 = 0`——**それが打ち消しの正体**である。

☆残るのは (1) この式から「`F` は原点の除いた近傍で解析関数と一致する」を出すこと、
(2) `z ≡ −w` での極の打ち消し（そちらは `℘` の 2 階の Taylor が要る）。 -/

theorem laurent_cancel (P : PeriodPair) (z x y : ℂ) :
    4 * laurentB P z * (laurentB P z - z ^ 2 * x) ^ 2
        - (laurentA P z - z ^ 3 * y) ^ 2
      = z ^ 2 * laurentM P x y z := by
  simp only [laurentA, laurentB, laurentM]
  ring

/-- ★★★★★`M(0) = −8·x`——これが原点での値を `0` にする。 -/
@[simp] theorem laurentM_zero (P : PeriodPair) (x y : ℂ) : laurentM P x y 0 = -8 * x := by
  simp [laurentM]

/-- ★★★★★★★★★★★★**打ち消しの正体**——`z = 0` での値。

`4·B(0)·B(0)² − A(0)² = 4·1·1 − (−2)² = 0`。 -/
theorem laurent_cancel_zero (P : PeriodPair) (x y : ℂ) :
    4 * laurentB P 0 * (laurentB P 0 - (0:ℂ) ^ 2 * x) ^ 2
      - (laurentA P 0 - (0:ℂ) ^ 3 * y) ^ 2 = 0 := by
  simp
  norm_num

def laurent_cancel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(加法定理の極の打ち消し——分子が z² でくくれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★加法定理の欠損関数 -/

/-- ★★★★★★★★**加法定理の欠損** `F_w(z)`。

    `F_w(z) ≔ ℘(z+w) + ℘(z) + ℘(w) − (1/4)·((℘′(z)−℘′(w))/(℘(z)−℘(w)))²`

★`Skeleton/GenEll/AdditionTheorem.lean`（第 606）の `weierstrassP_add` は
**`F_w ≡ 0`** と同値である。 -/
noncomputable def addDefect (P : PeriodPair) (w z : ℂ) : ℂ :=
  P.weierstrassP (z + w) + P.weierstrassP z + P.weierstrassP w
    - ((P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4

/-- ★★★★★★★★**原点の近くでの解析的な姿**——`z²` が約された形。 -/
noncomputable def addDefectNear (P : PeriodPair) (w z : ℂ) : ℂ :=
  P.weierstrassP (z + w) + P.weierstrassP w
    + laurentM P (P.weierstrassP w) (P.derivWeierstrassP w) z
      / (4 * (laurentB P z - z ^ 2 * P.weierstrassP w) ^ 2)

/-- ★★★★★★★★★★★★★★★★★★**`z ≠ 0` では両者は一致する**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`℘(z) = B/z²`・`℘′(z) = A/z³` を入れると分母の `z²` が `laurent_cancel` の
`z²` と約されるからである。☆**これが「原点の極が消える」ことの中身**である。 -/
theorem addDefect_eq_near (P : PeriodPair) (w z : ℂ) (hz : z ≠ 0)
    (hpz : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    addDefect P w z = addDefectNear P w z := by
  have hB := laurentB_eq P z hz
  have hA := laurentA_eq P z hz
  have hcan := laurent_cancel P z (P.weierstrassP w) (P.derivWeierstrassP w)
  rw [hB, hA] at hcan
  have hz2 : z ^ 2 ≠ 0 := pow_ne_zero _ hz
  have hkey : z ^ 4 * (4 * P.weierstrassP z
        * (P.weierstrassP z - P.weierstrassP w) ^ 2
      - (P.derivWeierstrassP z - P.derivWeierstrassP w) ^ 2)
      = laurentM P (P.weierstrassP w) (P.derivWeierstrassP w) z := by
    refine mul_left_cancel₀ hz2 ?_
    linear_combination hcan
  simp only [addDefect, addDefectNear, hB]
  field_simp
  linear_combination hkey

/-- ★★★★★★★★★★★★★★★★**原点での値は `0`**——`M(0) = −8℘(w)`・`B(0) = 1` から

    `℘(w) + ℘(w) + (−8℘(w))/4 = 2℘(w) − 2℘(w) = 0`。

★★☆**これが「加法定理が原点で成り立つ」ことである**——`F_w` の解析接続は原点で消える。 -/
@[simp] theorem addDefectNear_zero (P : PeriodPair) (w : ℂ) :
    addDefectNear P w 0 = 0 := by
  simp only [addDefectNear, laurentM_zero, laurentB_zero, zero_add]
  norm_num
  ring

/-- ★★★★★★★★★★★★**`addDefectNear` は原点で解析的**（`w ∉ Λ` のとき）。 -/
theorem analyticAt_addDefectNear (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    AnalyticAt ℂ (addDefectNear P w) 0 := by
  have hpw : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w)) 0 := by
    have hf : AnalyticAt ℂ (fun z : ℂ => z + w) 0 := analyticAt_id.add analyticAt_const
    have hg : AnalyticAt ℂ P.weierstrassP ((fun z : ℂ => z + w) 0) := by
      simpa using P.analyticOnNhd_weierstrassP w hw
    exact AnalyticAt.comp (f := fun z : ℂ => z + w) (x := 0) hg hf
  have he : AnalyticAt ℂ (P.weierstrassPExcept 0) 0 :=
    ((P.differentiableOn_weierstrassPExcept 0).analyticOnNhd
      P.isOpen_compl_lattice_sdiff) 0 (by simp)
  have hf' : AnalyticAt ℂ (P.derivWeierstrassPExcept 0) 0 :=
    ((P.differentiableOn_derivWeierstrassPExcept 0).analyticOnNhd
      P.isOpen_compl_lattice_sdiff) 0 (by simp)
  have hM : AnalyticAt ℂ (laurentM P (P.weierstrassP w) (P.derivWeierstrassP w)) 0 := by
    unfold laurentM
    fun_prop (disch := assumption)
  have hD : AnalyticAt ℂ (fun z : ℂ => 4 * (laurentB P z - z ^ 2 * P.weierstrassP w) ^ 2) 0 := by
    unfold laurentB
    fun_prop (disch := assumption)
  have hDne : (fun z : ℂ => 4 * (laurentB P z - z ^ 2 * P.weierstrassP w) ^ 2) 0 ≠ 0 := by
    simp
  exact (hpw.add analyticAt_const).add (hM.div hD hDne)

def addDefect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(加法定理の欠損 F_w——F_w ≡ 0 が加法定理そのもの)",
    sectionId := "genell-lemma-3-5" }

def addDefect_eq_near.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(原点の極が消える——z² が約される。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★`z ≡ −w` の側——`q` の零点の位数 -/

/-- ★★★★★★★★**`z ≡ −w` の側の鍵となる関数**

    `q(t) ≔ 2·(℘(t−w) − ℘(w)) − t·(℘′(t−w) − ℘′(w))`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`t ≔ z + w` と置くと `F_w(t−w)` の極は `q` が `t = 0` で 3 位の零点を持つことで消える
（`u ≔ ℘(t−w) − ℘(w)`・`v ≔ ℘′(t−w) − ℘′(w)`・`û ≔ u/t` として `2û − v = q/t`、
`4û² − v² = (2û−v)(2û+v)`）。 -/
noncomputable def addQ (P : PeriodPair) (w t : ℂ) : ℂ :=
  2 * (P.weierstrassP (t - w) - P.weierstrassP w)
    - t * (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)

/-- ★★★★★★**`q(0) = 0`**——`℘` は偶だから `℘(−w) = ℘(w)`。 -/
@[simp] theorem addQ_zero (P : PeriodPair) (w : ℂ) : addQ P w 0 = 0 := by
  simp only [addQ, zero_sub, zero_mul, sub_zero, P.weierstrassP_neg]
  ring

/-- ★★★★★★★★**`q` の 1 階導関数** `q′(t) = ℘′(t−w) + ℘′(w) − t·℘″(t−w)`。 -/
theorem hasDerivAt_addQ (P : PeriodPair) (w t : ℂ) (ht : t - w ∉ P.lattice) :
    HasDerivAt (addQ P w)
      (P.derivWeierstrassP (t - w) + P.derivWeierstrassP w
        - t * deriv P.derivWeierstrassP (t - w)) t := by
  have h1 : HasDerivAt (fun s : ℂ => P.weierstrassP (s - w))
      (P.derivWeierstrassP (t - w)) t :=
    HasDerivAt.comp_sub_const t w (hasDerivAt_weierstrassP P ht)
  have h2 : HasDerivAt (fun s : ℂ => P.derivWeierstrassP (s - w))
      (deriv P.derivWeierstrassP (t - w)) t :=
    HasDerivAt.comp_sub_const t w (hasDerivAt_derivWeierstrassP P ht)
  have h3 : HasDerivAt (fun s : ℂ => s * (P.derivWeierstrassP (s - w) - P.derivWeierstrassP w))
      (1 * (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)
        + t * deriv P.derivWeierstrassP (t - w)) t :=
    (hasDerivAt_id t).mul (h2.sub_const _)
  have h4 := ((h1.sub_const (P.weierstrassP w)).const_mul (2:ℂ)).sub h3
  have hval : (2:ℂ) * P.derivWeierstrassP (t - w)
      - (1 * (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)
        + t * deriv P.derivWeierstrassP (t - w))
      = P.derivWeierstrassP (t - w) + P.derivWeierstrassP w
        - t * deriv P.derivWeierstrassP (t - w) := by ring
  rw [← hval]
  exact h4

/-- ★★★★★★★★★★**`q′(0) = 0`**——`℘′` が奇であることがちょうど効く。

    `q′(0) = ℘′(−w) + ℘′(w) − 0 = −℘′(w) + ℘′(w) = 0` -/
theorem deriv_addQ_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    deriv (addQ P w) 0 = 0 := by
  have hnw : (0 : ℂ) - w ∉ P.lattice := by
    intro hc
    exact hw (by simpa using neg_mem hc)
  have h := hasDerivAt_addQ P w 0 hnw
  rw [h.deriv]
  simp only [zero_sub, zero_mul, sub_zero, P.derivWeierstrassP_neg]
  ring

def addQ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(z ≡ −w の側の鍵——q は t = 0 で 3 位の零点)",
    sectionId := "genell-lemma-3-5" }

def deriv_addQ_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q′(0) = 0——℘′ が奇であることが効く。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★**`q` の 1 階導関数（閉じた形）**——`Found/GenEll/WeierstrassODE.lean` の
`deriv_derivWeierstrassP`（`deriv ℘′ = 6℘² − g₂/2`）を入れた形。 -/
noncomputable def addQ' (P : PeriodPair) (w t : ℂ) : ℂ :=
  P.derivWeierstrassP (t - w) + P.derivWeierstrassP w
    - t * (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)

theorem hasDerivAt_addQ_closed (P : PeriodPair) (w t : ℂ) (ht : t - w ∉ P.lattice) :
    HasDerivAt (addQ P w) (addQ' P w t) t := by
  have h := hasDerivAt_addQ P w t ht
  rw [deriv_derivWeierstrassP P ht] at h
  exact h

/-- ★★★★★★★★★★★★**`q″(t) = −12·t·℘(t−w)·℘′(t−w)`**——★**`t` の因子が残る**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`deriv ℘′ = 6℘² − g₂/2` の項がちょうど打ち消し合うので `t` の因子だけが残る。
☆これが第 611 で「自動的に消える」と書いた段である。 -/
theorem hasDerivAt_addQ' (P : PeriodPair) (w t : ℂ) (ht : t - w ∉ P.lattice) :
    HasDerivAt (addQ' P w)
      (-(12 * t * P.weierstrassP (t - w) * P.derivWeierstrassP (t - w))) t := by
  have hp : HasDerivAt (fun s : ℂ => P.weierstrassP (s - w))
      (P.derivWeierstrassP (t - w)) t :=
    HasDerivAt.comp_sub_const t w (hasDerivAt_weierstrassP P ht)
  have h2 : HasDerivAt (fun s : ℂ => P.derivWeierstrassP (s - w))
      (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2) t := by
    have h := HasDerivAt.comp_sub_const t w (hasDerivAt_derivWeierstrassP P ht)
    rwa [deriv_derivWeierstrassP P ht] at h
  have hsq : HasDerivAt (fun s : ℂ => 6 * P.weierstrassP (s - w) ^ 2 - P.g₂ / 2)
      (6 * (2 * P.weierstrassP (t - w) ^ 1 * P.derivWeierstrassP (t - w))) t :=
    ((hp.pow 2).const_mul 6).sub_const _
  have hmul : HasDerivAt (fun s : ℂ => s * (6 * P.weierstrassP (s - w) ^ 2 - P.g₂ / 2))
      (1 * (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)
        + t * (6 * (2 * P.weierstrassP (t - w) ^ 1 * P.derivWeierstrassP (t - w)))) t :=
    (hasDerivAt_id t).mul hsq
  have h4 := (h2.add_const (P.derivWeierstrassP w)).sub hmul
  have hval : (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)
      - (1 * (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)
        + t * (6 * (2 * P.weierstrassP (t - w) ^ 1 * P.derivWeierstrassP (t - w))))
      = -(12 * t * P.weierstrassP (t - w) * P.derivWeierstrassP (t - w)) := by
    ring
  rw [← hval]
  exact h4

/-- ★★★★★★★★★★★★**`q″(0) = 0`**。 -/
theorem iteratedDeriv_two_addQ_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    iteratedDeriv 2 (addQ P w) 0 = 0 := by
  have hnw : (0 : ℂ) - w ∉ P.lattice := fun hc => hw (by simpa using neg_mem hc)
  have hopen : IsOpen {t : ℂ | t - w ∉ P.lattice} := by
    have : {t : ℂ | t - w ∉ P.lattice} = (fun t : ℂ => t - w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
    rw [this]
    exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
  have hnhds : {t : ℂ | t - w ∉ P.lattice} ∈ nhds (0 : ℂ) := hopen.mem_nhds hnw
  have heq : deriv (addQ P w) =ᶠ[nhds (0:ℂ)] addQ' P w := by
    filter_upwards [hnhds] with t ht
    exact (hasDerivAt_addQ_closed P w t ht).deriv
  rw [iteratedDeriv_succ, iteratedDeriv_one, heq.deriv_eq,
    (hasDerivAt_addQ' P w 0 hnw).deriv]
  ring

/-- ★★★★★★★★★★★★★★★★★★**`q` は `t = 0` で 3 位の零点**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`q(0) = q′(0) = q″(0) = 0`（第 612・第 614）から、mathlib の
`natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero` で位数が `≥ 3` と分かる。

☆これで `2û − v = q/t` が `2` 位の零点になり、
`4û² − v² = (2û−v)(2û+v)` が `2` 位で消えて **`z ≡ −w` の極が打ち消される**。 -/
theorem three_le_analyticOrderAt_addQ (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (hana : AnalyticAt ℂ (addQ P w) 0) :
    (3 : ℕ) ≤ analyticOrderAt (addQ P w) 0 := by
  rw [natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero hana]
  intro i hi
  interval_cases i
  · simpa using addQ_zero P w
  · rw [iteratedDeriv_one]; exact deriv_addQ_zero P w hw
  · exact iteratedDeriv_two_addQ_zero P w hw

def hasDerivAt_addQ'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q″(t) = −12t·℘(t−w)·℘′(t−w)——t の因子が残る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def three_le_analyticOrderAt_addQ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q は t = 0 で 3 位の零点——z ≡ −w の極が打ち消される)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★零点勘定への第一の煉瓦 -/

/-- ★★★★★★★★★★★★★★★★★★★★**楕円関数の周期平行四辺形の境界積分は消える**。

    `∮_{∂D} f dz = 0`（`D` は `a` を頂点とする周期平行四辺形）

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★機構は**周期性だけ**——Cauchy の定理は要らない:

* `f(a + ω₂ + tω₁) = f(a + tω₁)`（`ω₂ ∈ Λ`）なので上辺と下辺が打ち消し合う
* `f(a + ω₁ + tω₂) = f(a + tω₂)`（`ω₁ ∈ Λ`）なので左辺と右辺が打ち消し合う

★★★★☆**これが楕円関数の零点勘定の第一の煉瓦である**。
残る半分は**留数定理**（`∮ f = 2πi·Σ res`）であり、それと合わせて

    `Σ_{D} res(f) = 0`

が出る。`f = g′/g` に当てると **`#零点 = #極`**（偏角の原理）になり、
`g = ℘ − c` なら「`℘` は各値をちょうど 2 回取る」が従う
——それが `Skeleton/GenEll/AdditionTheorem.lean`（第 615）で同定した唯一の入口である。

☆mathlib は軸平行な長方形での Cauchy（`integral_boundary_rect_eq_zero_of_...`）を持つが、
一般の格子の平行四辺形には当たらない（2026-08-29 に測定）。 -/
theorem elliptic_boundary_integral_zero (P : PeriodPair) (f : ℂ → ℂ)
    (hper : ∀ (z : ℂ), ∀ l ∈ P.lattice, f (z + l) = f z) (a : ℂ) :
    (∫ t in (0:ℝ)..1, f (a + (t : ℂ) * P.ω₁) * P.ω₁)
      + (∫ t in (0:ℝ)..1, f (a + P.ω₁ + (t : ℂ) * P.ω₂) * P.ω₂)
      - (∫ t in (0:ℝ)..1, f (a + P.ω₂ + (t : ℂ) * P.ω₁) * P.ω₁)
      - (∫ t in (0:ℝ)..1, f (a + (t : ℂ) * P.ω₂) * P.ω₂) = 0 := by
  have hω₁ : P.ω₁ ∈ P.lattice := Submodule.subset_span (by simp)
  have hω₂ : P.ω₂ ∈ P.lattice := Submodule.subset_span (by simp)
  have h1 : ∀ t : ℝ, f (a + P.ω₂ + (t : ℂ) * P.ω₁) = f (a + (t : ℂ) * P.ω₁) := by
    intro t
    have hz : a + P.ω₂ + (t : ℂ) * P.ω₁ = (a + (t : ℂ) * P.ω₁) + P.ω₂ := by ring
    rw [hz, hper _ _ hω₂]
  have h2 : ∀ t : ℝ, f (a + P.ω₁ + (t : ℂ) * P.ω₂) = f (a + (t : ℂ) * P.ω₂) := by
    intro t
    have hz : a + P.ω₁ + (t : ℂ) * P.ω₂ = (a + (t : ℂ) * P.ω₂) + P.ω₁ := by ring
    rw [hz, hper _ _ hω₁]
  simp only [h1, h2]
  ring

def elliptic_boundary_integral_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(楕円関数の境界積分は消える——周期性だけから。★無条件)",
    sectionId := "genell-lemma-3-5" }

def elliptic_boundary_integral_zero.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残る半分は留数定理(∮ f = 2πi·Σ res)である。" ++
       "★合わせて Σ_D res(f) = 0 が出て、f = g′/g に当てると #零点 = #極(偏角の原理)、" ++
       "g = ℘ − c なら「℘ は各値をちょうど 2 回取る」が従う" ++
       "——Skeleton/GenEll/AdditionTheorem.lean(第 615)で同定した唯一の入口である") 21,
    .implicitStep
      ("☆mathlib の在庫(2026-08-29): 軸平行な長方形での Cauchy" ++
       "(Complex.integral_boundary_rect_eq_zero_of_differentiable_on_off_countable)・" ++
       "Jensen の公式・Nevanlinna 理論(ValueDistribution/)・MeromorphicOn.divisor はある。" ++
       "☆偏角の原理と楕円関数の零点和 = 極和は無い。" ++
       "★一般の格子の平行四辺形には長方形版の Cauchy は当たらない" ++
       "(変数変換が ℝ-線型で正則性を壊す)") 13 ]

/-- ★★★★★★★★**`2f + a ∈ Λ` なら `℘(f+a) = ℘(f)`**——`℘` が偶であることから。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★`G(z) ≔ ℘(z+a) − ℘(z)` は `z ↦ −a−z` について奇であり、その**不動点**
（`2f ≡ −a (mod Λ)` を満たす `f`——`Λ/2Λ` の分だけ **4 個**ある）で消える。
☆これが「`℘(z+a) = ℘(z)` の解は `a ∉ Λ` でも 4 個ある」という事実であり、
★★★零点勘定（`#零点 = #極`）と組み合わせると
**`℘` が各値を 2 回しか取らない**ことの矛盾を作る材料になる
（`G` の極は `Λ` と `−a+Λ` に 2 位ずつ＝計 4）。 -/
theorem weierstrassP_shift_eq_of_two_add_mem (P : PeriodPair) (f a : ℂ)
    (h : 2 * f + a ∈ P.lattice) :
    P.weierstrassP (f + a) = P.weierstrassP f := by
  have h1 : P.weierstrassP (f + a) = P.weierstrassP (-(f + a)) := (P.weierstrassP_neg _).symm
  have h2 : -(f + a) = f + (-(2 * f + a)) := by ring
  rw [h1, h2, P.weierstrassP_add_coe f ⟨-(2 * f + a), neg_mem h⟩]

open Filter Topology Bornology Metric in
/-- ★★★★★★★★★★★★**`℘` の周期群はちょうど `Λ`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★機構: `℘(z+a) = ℘(z)` が恒等的に成り立つなら、`℘` は `−a` でも極を持つ
（`℘(z+a)` が `z = −a` で極だから）。★★しかし `−a ∉ Λ` なら `℘` は `−a` で解析的
——連続なので有界な極限を持ち、`cobounded` へ発散することと両立しない。

☆これで `G(z) ≔ ℘(z+a) − ℘(z)` が `a ∉ Λ` のとき恒等的に `0` でないことが分かる
——零点勘定の議論で「`G ≢ 0`」を言うのに要る。 -/
theorem mem_lattice_of_weierstrassP_periodic (P : PeriodPair) (a : ℂ)
    (h : ∀ z, P.weierstrassP (z + a) = P.weierstrassP z) : a ∈ P.lattice := by
  by_contra hc
  have hna : -a ∉ P.lattice := fun hm => hc (by simpa using neg_mem hm)
  have hcont : ContinuousAt P.weierstrassP (-a) :=
    (P.analyticOnNhd_weierstrassP (-a) hna).continuousAt
  have hord : meromorphicOrderAt P.weierstrassP 0 < 0 := by
    rw [P.order_weierstrassP 0 P.lattice.zero_mem]; decide
  have h1 : Tendsto P.weierstrassP (𝓝[≠] (0:ℂ)) (cobounded ℂ) :=
    tendsto_cobounded_of_meromorphicOrderAt_neg hord
  have hshift : Tendsto (fun z : ℂ => z + a) (𝓝[≠] (-a)) (𝓝[≠] (0:ℂ)) := by
    rw [tendsto_nhdsWithin_iff]
    refine ⟨?_, ?_⟩
    · have ht : Tendsto (fun z : ℂ => z + a) (𝓝 (-a)) (𝓝 ((-a) + a)) :=
        (continuous_id.add continuous_const).tendsto _
      simpa using ht.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with z hz
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff] at hz ⊢
      intro hcc
      exact hz (by linear_combination hcc)
  have h3 : Tendsto P.weierstrassP (𝓝[≠] (-a)) (cobounded ℂ) := by
    have h2 := h1.comp hshift
    have hfun : (P.weierstrassP ∘ fun z : ℂ => z + a) = P.weierstrassP := by
      funext z; exact h z
    rwa [hfun] at h2
  have h4 : Tendsto P.weierstrassP (𝓝[≠] (-a)) (𝓝 (P.weierstrassP (-a))) :=
    hcont.continuousWithinAt
  exact (h4.not_tendsto (disjoint_nhds_cobounded _)) h3

def weierstrassP_shift_eq_of_two_add_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2f + a ∈ Λ なら ℘(f+a) = ℘(f)——℘ が偶だから。★無条件)",
    sectionId := "genell-lemma-3-5" }

def mem_lattice_of_weierstrassP_periodic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘ の周期群はちょうど Λ。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
