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
import ABC3.Found.GenEll.Uniformization.G2G3

/-!
# 一様化 —— 代表系の像は `⟨Q⟩`・Vélu の `B` の和・組み立て・代表系の具体形・出典

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★代表系の像は `⟨Q⟩` -/

open scoped Classical in
/-- ★★★★★★★★`Φ(k·z) = k • Φ(z)`（`k : ℕ`）。 -/
theorem uniformMap_natCast_mul (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ)
    (k : ℕ) : uniformMap P hΔ ((k : ℂ) * z) = k • uniformMap P hΔ z := by
  rw [show ((k : ℂ) * z) = k • z by simp [nsmul_eq_mul], ← uniformHom_apply,
    map_nsmul, uniformHom_apply]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**代表系 `T = {0, z₀, …, (l−1)z₀}` の像は `⟨Q⟩ = {0, Q, …, (l−1)Q}`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★☆**これが Galois 降下の鍵である**——解析側で選んだ代表系 `T` の像は
`Q` が生成する巡回部分群そのものなので、`veluQuotientFull` に渡す点集合 `S` は
**`Q` だけで決まる**（`z₀` の選び方に依らない）。
したがって `Q` が `L`-有理なら `S` も `L`-有理である。 -/
theorem image_uniformMap_veluRep (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z₀ : ℂ)
    (l : ℕ) :
    ((Finset.range l).image (fun k : ℕ => (k : ℂ) * z₀)).image (uniformMap P hΔ)
      = (Finset.range l).image (fun k : ℕ => k • uniformMap P hΔ z₀) := by
  rw [Finset.image_image]
  refine Finset.image_congr ?_
  intro k _
  exact uniformMap_natCast_mul P hΔ z₀ k

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**`T∖{0}` の像は `⟨Q⟩∖{O}`**。

★`0 < k < l` で `k • Q ≠ 0`（`Q` の位数がちょうど `l` だから）。 -/
theorem uniformMap_ne_zero_of_mem_erase (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q) {k : ℕ} (hk0 : k ≠ 0) (hkl : k < l) :
    uniformMap P hΔ ((k : ℂ) * z₀) ≠ 0 := by
  rw [uniformMap_natCast_mul, hz₀]
  intro hc
  have hdvd : addOrderOf Q ∣ k := addOrderOf_dvd_of_nsmul_eq_zero hc
  rw [hQ] at hdvd
  have := Nat.le_of_dvd (Nat.pos_of_ne_zero hk0) hdvd
  omega

def image_uniformMap_veluRep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系 T の像は ⟨Q⟩——Galois 降下の鍵。★無条件)",
    sectionId := "genell-lemma-3-5" }

def uniformMap_ne_zero_of_mem_erase.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(0 < k < l なら k·z₀ の像は O でない。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**代表系の上で `f(℘)·℘′` の和は消える**

    `Σ_{w ∈ T∖{0}} f(℘_Λ(w))·℘′_Λ(w) = 0`   （任意の `f : ℂ → ℂ`）

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★第 670 の一般化。`ν w ≡ −w` は `T∖{0}` の対合で、
`℘` は偶（`℘(ν w) = ℘(w)`）・`℘′` は奇（`℘′(ν w) = −℘′(w)`）だから
`S = Σ f(℘(ν w))·℘′(ν w) = −S`。

★★☆`f = 1` なら第 670（`Σ ℘′ = 0`）、`f = id` なら `Σ ℘·℘′ = 0`。
後者は Vélu の `Σ B·x = 0`（第 688 の仮説）に要る。 -/
theorem sum_mul_derivWeierstrassP_rep_eq_zero (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) (f : ℂ → ℂ) :
    ∑ w ∈ T.erase 0, f (P.weierstrassP w) * P.derivWeierstrassP w = 0 := by
  classical
  have huniq : ∀ {w v : ℂ}, w ∈ T → v ∈ T → w - v ∈ P.lattice → w = v :=
    fun hw hv hd => rep_sub_mem_lattice_imp_eq P P' T hT hrep hw hv hd
  have hex : ∀ w ∈ T, ∃ v, v ∈ T ∧ w + v ∈ P.lattice := by
    intro w hw
    obtain ⟨v, hv, hv2, -⟩ := hrep w (hT w hw)
    exact ⟨v, hv, hv2⟩
  choose! ν hνT hνΛ using hex
  have hνe : ∀ w ∈ T.erase 0, ν w ∈ T.erase 0 := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have hw0 : w ≠ 0 := Finset.ne_of_mem_erase hw
    refine Finset.mem_erase.2 ⟨?_, hνT w hw'⟩
    intro hc
    refine hw0 (huniq hw' h0T ?_)
    have hz := hνΛ w hw'
    rw [hc, add_zero] at hz
    simpa using hz
  have hinvol : ∀ w ∈ T.erase 0, ν (ν w) = w := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have h1 := hνΛ w hw'
    have h2 := hνΛ (ν w) (hνT w hw')
    refine huniq (hνT (ν w) (hνT w hw')) hw' ?_
    have hd := P.lattice.sub_mem h2 h1
    rw [show ν w + ν (ν w) - (w + ν w) = ν (ν w) - w by ring] at hd
    exact hd
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0, ν w = ν v → w = v := by
    intro w hw v hv he
    rw [← hinvol w hw, ← hinvol v hv, he]
  have hodd : ∀ w ∈ T, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w)
      = -(f (P.weierstrassP w) * P.derivWeierstrassP w) := by
    intro w hw
    have hl : w + ν w ∈ P.lattice := hνΛ w hw
    have he : ν w = -w + (w + ν w) := by ring
    have hp : P.weierstrassP (ν w) = P.weierstrassP w := by
      rw [he, P.weierstrassP_add_coe (-w) ⟨w + ν w, hl⟩, P.weierstrassP_neg]
    have hpd : P.derivWeierstrassP (ν w) = -P.derivWeierstrassP w := by
      rw [he, P.derivWeierstrassP_add_coe (-w) ⟨w + ν w, hl⟩, P.derivWeierstrassP_neg]
    rw [hp, hpd]
    ring
  have hinjOn : Set.InjOn ν ↑(T.erase 0) := fun w hw v hv he =>
    hinj w (Finset.mem_coe.1 hw) v (Finset.mem_coe.1 hv) he
  have himg : (T.erase 0).image ν = T.erase 0 :=
    Finset.eq_of_subset_of_card_le
      (fun v hv => by
        obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 hv
        exact hνe w hw)
      (le_of_eq (Finset.card_image_of_injOn hinjOn).symm)
  have h1 : ∑ v ∈ T.erase 0, f (P.weierstrassP v) * P.derivWeierstrassP v
      = ∑ w ∈ T.erase 0, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w) := by
    conv_lhs => rw [← himg]
    exact Finset.sum_image (fun w hw v hv he => hinj w hw v hv he)
  have h2 : ∑ w ∈ T.erase 0, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w)
      = -∑ w ∈ T.erase 0, f (P.weierstrassP w) * P.derivWeierstrassP w := by
    have hc : ∑ w ∈ T.erase 0, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w)
        = ∑ w ∈ T.erase 0, (-(f (P.weierstrassP w) * P.derivWeierstrassP w)) :=
      Finset.sum_congr rfl (fun w hw => hodd w (Finset.mem_of_mem_erase hw))
    rw [hc, Finset.sum_neg_distrib]
  have h3 : (2 : ℂ) * ∑ w ∈ T.erase 0,
      f (P.weierstrassP w) * P.derivWeierstrassP w = 0 := by
    linear_combination h1.trans h2
  simpa using h3

def sum_mul_derivWeierstrassP_rep_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の上で f(℘)·℘′ の和は消える——第 670 の一般化。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★Vélu の `B` の和が消えること -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**`Σ_S B = 0`**（`S` は解析側の代表点集合）。

`latticeCurve` では `a₁ = a₃ = 0` なので `B(Q) = 2·y_Q = ℘′(w)`。
第 670（`Σ ℘′(w) = 0`）そのものである。 -/
theorem sum_veluB_image_eq_zero (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    ∑ Q ∈ (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)),
        (2 * Q.2 + (latticeCurve P).a₁ * Q.1 + (latticeCurve P).a₃) = 0 := by
  have hnot : ∀ w ∈ T.erase 0, w ∉ P.lattice :=
    fun w hw => rep_notMem_lattice P P' T h0T hT hrep hw
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0,
      (latticePointX P w, latticePointY P w)
        = (latticePointX P v, latticePointY P v) → w = v := by
    intro w hw v hv he
    refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw)
      (Finset.mem_of_mem_erase hv) ?_
    exact sub_mem_lattice_of_point_eq P (hnot w hw) (hnot v hv)
      (congrArg Prod.fst he) (congrArg Prod.snd he)
  rw [Finset.sum_image hinj]
  have hstep : ∀ w ∈ T.erase 0,
      (2 * latticePointY P w + (latticeCurve P).a₁ * latticePointX P w
        + (latticeCurve P).a₃) = P.derivWeierstrassP w := by
    intro w _
    simp only [latticePointY, latticePointX, latticeCurve]
    ring
  rw [Finset.sum_congr rfl hstep]
  exact sum_derivWeierstrassP_rep_eq_zero P P' T h0T hT hrep

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**`Σ_S B·x = 0`**——第 690 で `f = id`。 -/
theorem sum_veluBx_image_eq_zero (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    ∑ Q ∈ (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)),
        (2 * Q.2 + (latticeCurve P).a₁ * Q.1 + (latticeCurve P).a₃) * Q.1 = 0 := by
  have hnot : ∀ w ∈ T.erase 0, w ∉ P.lattice :=
    fun w hw => rep_notMem_lattice P P' T h0T hT hrep hw
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0,
      (latticePointX P w, latticePointY P w)
        = (latticePointX P v, latticePointY P v) → w = v := by
    intro w hw v hv he
    refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw)
      (Finset.mem_of_mem_erase hv) ?_
    exact sub_mem_lattice_of_point_eq P (hnot w hw) (hnot v hv)
      (congrArg Prod.fst he) (congrArg Prod.snd he)
  rw [Finset.sum_image hinj]
  have hstep : ∀ w ∈ T.erase 0,
      (2 * latticePointY P w + (latticeCurve P).a₁ * latticePointX P w
          + (latticeCurve P).a₃) * latticePointX P w
        = (fun z : ℂ => z) (P.weierstrassP w) * P.derivWeierstrassP w := by
    intro w _
    simp only [latticePointY, latticePointX, latticeCurve]
    ring
  rw [Finset.sum_congr rfl hstep]
  exact sum_mul_derivWeierstrassP_rep_eq_zero P P' T h0T hT hrep (fun z => z)

def sum_veluB_image_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Σ_S B = 0——解析側の代表点集合で。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_veluBx_image_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Σ_S B·x = 0——第 690 で f = id。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★組み立て -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**位数 `l` の点から、`Vélu` の商に要るデータをすべて出す**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

第 679 に `Σ_S B = 0`・`Σ_S B·x = 0`（第 691）を足した形。
★☆この 2 つが `veluQuotientFull` を変数変換で移すのに要る（第 688・第 693）。 -/
theorem exists_veluQuotientFull_data_of_torsion (P : PeriodPair)
    (hΔ : latticeDisc P ≠ 0) {Q : (latticeCurve P).toAffine.Point} {l : ℕ}
    (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (P' : PeriodPair) (A B C D : ℤ) (S : Finset (ℂ × ℂ)),
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      S.card = l - 1 ∧
      (∑ q ∈ S, (2 * q.2 + (latticeCurve P).a₁ * q.1 + (latticeCurve P).a₃) = 0) ∧
      (∑ q ∈ S, (2 * q.2 + (latticeCurve P).a₁ * q.1 + (latticeCurve P).a₃) * q.1 = 0) ∧
      latticeCurve P' = veluQuotientFull (latticeCurve P) S := by
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' (intCast_mul_mem_lattice_iff P hΔ hQ hz₀)
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  have hnot : ∀ w ∈ T.erase 0, w ∉ P.lattice :=
    fun w hw => rep_notMem_lattice P P' T h0T hT hrep hw
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0,
      (latticePointX P w, latticePointY P w)
        = (latticePointX P v, latticePointY P v) → w = v := by
    intro w hw v hv he
    refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw)
      (Finset.mem_of_mem_erase hv) ?_
    exact sub_mem_lattice_of_point_eq P (hnot w hw) (hnot v hv)
      (congrArg Prod.fst he) (congrArg Prod.snd he)
  refine ⟨P', A, B, C, D,
    (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)),
    h1, h2, hdet, ?_,
    sum_veluB_image_eq_zero P P' T h0T hT hrep,
    sum_veluBx_image_eq_zero P P' T h0T hT hrep, ?_⟩
  · rw [Finset.card_image_of_injOn (fun w hw v hv he =>
      hinj w (Finset.mem_coe.1 hw) v (Finset.mem_coe.1 hv) he),
      Finset.card_erase_of_mem h0T, hcard]
  · exact latticeCurve_eq_veluQuotientFull P P' T h0T hT hrep
      (g₂_isogeny P P' T h0T hT hrep hvelu) (g₃_isogeny P P' T h0T hT hrep hvelu)

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**一意化の変数変換を外した形——`Lemma 3.5` の `ℂ` 側の最終形**

`C • W = latticeCurve P`（`W` は `E.map σ`、`C` は一意化の変数変換）のとき、
位数 `l` の点 `Q` から

    latticeCurve P′ = C • veluQuotientFull W S_W

★★★★★★☆**これが `htFalt_isogeny_le_of_velu`（第 678）に渡す形そのものである**
——`E′ ≔ veluQuotientFull W S_W` と取れば `C` は両側で同じ、すなわち `α = 1`。 -/
theorem exists_vc_veluQuotientFull_of_torsion (P : PeriodPair)
    (hΔ : latticeDisc P ≠ 0) (W : WeierstrassCurve ℂ) (C : VariableChange ℂ)
    (hCW : C • W = latticeCurve P)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l)
    (hQ : addOrderOf Q = l) :
    ∃ (P' : PeriodPair) (A B Cc D : ℤ) (SW : Finset (ℂ × ℂ)),
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (Cc : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * Cc).natAbs = l ∧
      latticeCurve P' = C • veluQuotientFull W SW := by
  obtain ⟨P', A, B, Cc, D, S, h1, h2, hdet, hcard, hB, hBx, hEq⟩ :=
    exists_veluQuotientFull_data_of_torsion P hΔ hl hQ
  refine ⟨P', A, B, Cc, D, S.image (vcInvPair C), h1, h2, hdet, ?_⟩
  rw [hEq, ← hCW]
  refine veluQuotientFull_eq_vc_pullback C W S ?_ ?_
  · rw [hCW]; exact hB
  · rw [hCW]; exact hBx

def exists_veluQuotientFull_data_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(位数 l の点から Vélu の商に要るデータをすべて出す。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_vc_veluQuotientFull_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一意化の変数変換を外した形——ℂ 側の最終形。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★**`Φ(z)` の座標は `(℘(z), ℘′(z)/2)`**（`z ∉ Λ`）。

★★☆これで解析側の代表点集合が「`Φ` の像の座標」として書ける——
`Q` が生成する部分群（第 680）の座標そのものである。 -/
theorem uniformMap_coords (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {z : ℂ}
    (hz : z ∉ P.lattice) :
    pointCoords (uniformMap P hΔ z) = (latticePointX P z, latticePointY P z) := by
  rw [uniformMap_of_notMem P hΔ hz]
  rfl

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**解析側の代表点集合は `⟨Q⟩∖{O}` の座標**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これで `S` が `Q` だけで決まることが**座標の言葉で**書けた——
`z₀` の選び方にも代表系 `T` の取り方にも依らない。
☆したがって `Q` が `L`-有理なら `S ⊆ L × L` である。 -/
theorem pointCoords_uniformMap_natCast_mul (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {z₀ : ℂ} {k : ℕ} (hk : ((k : ℂ) * z₀) ∉ P.lattice) :
    pointCoords (uniformMap P hΔ ((k : ℂ) * z₀))
      = (latticePointX P ((k : ℂ) * z₀), latticePointY P ((k : ℂ) * z₀)) :=
  uniformMap_coords P hΔ hk

def uniformMap_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ(z) の座標は (℘(z), ℘′(z)/2)。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★代表系の具体形 -/

open scoped Classical in
/-- ★★★★★★★★**`T∖{0}` は `{k·z₀ : 0 < k < l}`**。

★`k·z₀ = 0` は `k = 0` のときだけ（`k·z₀ ∈ Λ ⟺ l ∣ k`）。 -/
theorem velu_rep_erase_eq (P : PeriodPair) (z₀ : ℂ) (l : ℕ) (hl : 0 < l)
    (hord : ∀ k : ℤ, (k : ℂ) * z₀ ∈ P.lattice ↔ (l : ℤ) ∣ k) :
    ((Finset.range l).image (fun k : ℕ => (k : ℂ) * z₀)).erase 0
      = ((Finset.range l).erase 0).image (fun k : ℕ => (k : ℂ) * z₀) := by
  have hne : ∀ k : ℕ, k ≠ 0 → k < l → ((k : ℂ) * z₀) ≠ 0 := by
    intro k hk0 hkl hc
    have hmem : ((k : ℤ) : ℂ) * z₀ ∈ P.lattice := by
      have : ((k : ℤ) : ℂ) * z₀ = (k : ℂ) * z₀ := by push_cast; ring
      rw [this, hc]
      exact P.lattice.zero_mem
    have hdvd := (hord _).1 hmem
    have hdvdN : l ∣ k := by exact_mod_cast hdvd
    have := Nat.le_of_dvd (Nat.pos_of_ne_zero hk0) hdvdN
    omega
  ext x
  simp only [Finset.mem_erase, Finset.mem_image, Finset.mem_range]
  constructor
  · rintro ⟨hx0, k, hk, rfl⟩
    refine ⟨k, ⟨?_, hk⟩, rfl⟩
    rintro rfl
    exact hx0 (by simp)
  · rintro ⟨k, ⟨hk0, hkl⟩, rfl⟩
    exact ⟨hne k hk0 hkl, k, hkl, rfl⟩

open scoped Classical in
/-- ★★★★★★★★★★★★**`k·Q` の座標は `(℘(k z₀), ℘′(k z₀)/2)`**（`0 < k < l`）。

★★☆これで解析側の代表点集合が `⟨Q⟩∖{O}` の座標そのものだと分かる。 -/
theorem pointCoords_nsmul_eq (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q) {k : ℕ} (hk0 : k ≠ 0) (hkl : k < l) :
    pointCoords (k • Q)
      = (latticePointX P ((k : ℂ) * z₀), latticePointY P ((k : ℂ) * z₀)) := by
  have hnot : ((k : ℂ) * z₀) ∉ P.lattice := by
    intro hc
    have hmem : ((k : ℤ) : ℂ) * z₀ ∈ P.lattice := by
      have he : ((k : ℤ) : ℂ) * z₀ = (k : ℂ) * z₀ := by push_cast; ring
      rw [he]; exact hc
    have hdvd := (intCast_mul_mem_lattice_iff P hΔ hQ hz₀ (k : ℤ)).1 hmem
    have hdvdN : l ∣ k := by exact_mod_cast hdvd
    have := Nat.le_of_dvd (Nat.pos_of_ne_zero hk0) hdvdN
    omega
  rw [← hz₀, ← uniformMap_natCast_mul]
  exact uniformMap_coords P hΔ hnot

def velu_rep_erase_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(T∖{0} は {k·z₀ : 0 < k < l}。★無条件)",
    sectionId := "genell-lemma-3-5" }

def pointCoords_nsmul_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(k·Q の座標は (℘(k z₀), ℘′(k z₀)/2)。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（`ℂ` 側・点集合を `Q` で決めた形）**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

位数 `l` の点 `Q` に対し、**点集合を `⟨Q⟩∖{O}` の座標に固定した形**で

    latticeCurve P′ = veluQuotientFull (latticeCurve P)
      (((range l).erase 0).image (fun k => coords(k·Q)))

★★★★★☆**これで `S` が `Q` だけで決まる**——`z₀` の選び方にも代表系 `T` の
取り方にも依らない。☆したがって `Q` が `L`-有理なら `S ⊆ L × L` であり、
第 697 の `embPoint` で各 `σ` へ送ったものと一致する。 -/
theorem exists_veluQuotientFull_zmultiples (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l)
    (hQ : addOrderOf Q = l) :
    ∃ (P' : PeriodPair) (A B C D : ℤ),
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      (∑ q ∈ ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)),
        (2 * q.2 + (latticeCurve P).a₁ * q.1 + (latticeCurve P).a₃) = 0) ∧
      (∑ q ∈ ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)),
        (2 * q.2 + (latticeCurve P).a₁ * q.1 + (latticeCurve P).a₃) * q.1 = 0) ∧
      latticeCurve P'
        = veluQuotientFull (latticeCurve P)
            (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) := by
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  have hord := intCast_mul_mem_lattice_iff P hΔ hQ hz₀
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' hord
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  -- ★点集合が `⟨Q⟩∖{O}` の座標に等しいこと
  have hS : (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w))
      = ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)) := by
    rw [hTdef, velu_rep_erase_eq P z₀ l hl hord, Finset.image_image]
    refine Finset.image_congr ?_
    intro k hk
    rw [Finset.mem_coe, Finset.mem_erase, Finset.mem_range] at hk
    exact (pointCoords_nsmul_eq P hΔ hQ hz₀ hk.1 hk.2).symm
  refine ⟨P', A, B, C, D, h1, h2, hdet, ?_, ?_, ?_⟩
  · rw [← hS]; exact sum_veluB_image_eq_zero P P' T h0T hT hrep
  · rw [← hS]; exact sum_veluBx_image_eq_zero P P' T h0T hT hrep
  · rw [← hS]
    exact latticeCurve_eq_veluQuotientFull P P' T h0T hT hrep
      (g₂_isogeny P P' T h0T hT hrep hvelu) (g₃_isogeny P P' T h0T hT hrep hvelu)

def exists_veluQuotientFull_zmultiples.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ℂ 側・点集合を ⟨Q⟩∖{O} の座標に固定した形。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（`ℂ` 側・一意化の変数変換を外した最終形）**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`C • W = latticeCurve P`（`W = E.map σ`、`C` は一意化の変数変換）で
`W` 上の位数 `l` の点 `Q` があれば

    latticeCurve P′ = C • veluQuotientFull W (⟨Q⟩∖{O} の座標)

★★★★★★☆**これが `htFalt_isogeny_le_of_velu`（第 678）に渡す形そのものである。**
`E′ ≔ veluQuotientFull E S` と取れば、`veluQuotientFull` は底変換と可換（第 679）
なので各 `σ` で `C σ • (E′.map σ) = latticeCurve (P′ σ)` になる——
すなわち `α = 1`。 -/
theorem exists_periodPair_veluQuotientFull (W : WeierstrassCurve ℂ) [W.IsElliptic]
    (C : VariableChange ℂ) [(C • W).IsElliptic]
    (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (hCW : C • W = latticeCurve P)
    {Q : W.toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (P' : PeriodPair) (A B Cc D : ℤ),
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (Cc : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * Cc).natAbs = l ∧
      latticeCurve P' = C • veluQuotientFull W
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) := by
  have hQ₂ : addOrderOf
      (hCW ▸ vcPoint C W Q : (latticeCurve P).toAffine.Point) = l := by
    rw [addOrderOf_congr_curve, addOrderOf_vcPoint]
    exact hQ
  obtain ⟨P', A, B, Cc, D, h1, h2, hdet, hB, hBx, hEq⟩ :=
    exists_veluQuotientFull_zmultiples P hΔ hl hQ₂
  have hS₂ : ((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords
          (k • (hCW ▸ vcPoint C W Q : (latticeCurve P).toAffine.Point)))
      = (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))).image
          (fun q : ℂ × ℂ => (vcX C q.1, vcY C q.1 q.2)) := by
    rw [Finset.image_image]
    refine Finset.image_congr ?_
    intro k hk
    rw [Finset.mem_coe, Finset.mem_erase, Finset.mem_range] at hk
    have hkne : k • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hk.1 hk.2
    simp only [Function.comp_apply]
    rw [← nsmul_congr_curve hCW, pointCoords_congr_curve, ← vcPoint_nsmul,
      pointCoords_vcPoint_of_ne C W hkne]
  rw [hS₂] at hB hBx hEq
  rw [← hCW] at hB hBx hEq
  refine ⟨P', A, B, Cc, D, h1, h2, hdet, ?_⟩
  rw [hEq, veluQuotientFull_eq_vc_pullback C W _ hB hBx, image_vcInvPair_vcPair]

def exists_periodPair_veluQuotientFull.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ℂ 側・一意化の変数変換を外した最終形。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
