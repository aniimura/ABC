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
import ABC3.Found.GenEll.Uniformization.VeluAnalytic

/-!
# Surjective —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★平行移動が誘導する置換 -/

open Classical in
/-- ★★★★★平行移動 `w ↦ w + w₀` が代表系 `T` に誘導する写像。

★`-(w + w₀)` の代表を選ぶ（`Λ` を法として `w + w₀` と合同な `T` の元）。 -/
noncomputable def shiftRep (P : PeriodPair) (T : Finset ℂ) (w₀ : ℂ) (w : ℂ) : ℂ :=
  if h : ∃ v ∈ T, -(w + w₀) + v ∈ P.lattice then h.choose else 0

/-- ★★★★★★**代表系の一意性**——`hrep` の言い換え。 -/
theorem rep_unique (P P' : PeriodPair) (T : Finset ℂ)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (p : ℂ) (hp : p ∈ P'.lattice) (v v' : ℂ) (hv : v ∈ T) (hv' : v' ∈ T)
    (h1 : p + v ∈ P.lattice) (h2 : p + v' ∈ P.lattice) : v = v' := by
  obtain ⟨u, hu, -, huniq⟩ := hrep p hp
  have e1 : v = u := by by_contra hc; exact huniq v hv hc h1
  have e2 : v' = u := by by_contra hc; exact huniq v' hv' hc h2
  rw [e1, e2]

/-- ★★★★★★★★★★★★★★★★**平行移動は代表系を置換する**——`veluAnalyticX_shift` の
仮定 `σ` が**代表系の定義から作れる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★単射性は代表系の一意性から、全射性は `T` が有限だから出る。 -/
theorem exists_veluShiftPerm (P P' : PeriodPair) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (w₀ : ℂ) (hw₀ : w₀ ∈ P'.lattice) :
    ∃ σ : ℂ → ℂ, (∀ w ∈ T, σ w ∈ T) ∧ (∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
      ∧ (∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v) ∧ (∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) := by
  classical
  have hex : ∀ w ∈ T, ∃ v ∈ T, -(w + w₀) + v ∈ P.lattice := by
    intro w hw
    have hmem : -(w + w₀) ∈ P'.lattice := neg_mem (P'.lattice.add_mem (hT w hw) hw₀)
    obtain ⟨v, hv, hv2, -⟩ := hrep _ hmem
    exact ⟨v, hv, hv2⟩
  have hmemT : ∀ w ∈ T, shiftRep P T w₀ w ∈ T := by
    intro w hw
    rw [shiftRep, dif_pos (hex w hw)]
    exact (hex w hw).choose_spec.1
  have hinj : ∀ w ∈ T, ∀ w' ∈ T, shiftRep P T w₀ w = shiftRep P T w₀ w' → w = w' := by
    intro w hw w' hw' h
    have s1 := (hex w hw).choose_spec.2
    have s2 := (hex w' hw').choose_spec.2
    rw [shiftRep, dif_pos (hex w hw)] at h
    rw [shiftRep, dif_pos (hex w' hw')] at h
    rw [h] at s1
    have hdiff : -w + w' ∈ P.lattice := by
      have hs := P.lattice.sub_mem s2 s1
      have he : (-(w' + w₀) + (hex w' hw').choose) - (-(w + w₀) + (hex w' hw').choose)
          = -w' + w := by ring
      rw [he] at hs
      have hn := neg_mem hs
      simpa using hn
    have hzero : -w + w ∈ P.lattice := by simpa using P.lattice.zero_mem
    exact rep_unique P P' T hrep (-w) (neg_mem (hT w hw)) w w' hw hw' hzero hdiff
  refine ⟨shiftRep P T w₀, hmemT, hinj, ?_, ?_⟩
  · intro v hv
    obtain ⟨a, ha, hav⟩ := Finset.surj_on_of_inj_on_of_card_le (s := T) (t := T)
      (fun a _ => shiftRep P T w₀ a) (fun a ha => hmemT a ha)
      (fun a₁ a₂ ha₁ ha₂ h => hinj a₁ ha₁ a₂ ha₂ h) le_rfl v hv
    exact ⟨a, ha, hav.symm⟩
  · intro w hw
    have s1 := (hex w hw).choose_spec.2
    rw [shiftRep, dif_pos (hex w hw)]
    have he : w + w₀ - (hex w hw).choose = -(-(w + w₀) + (hex w hw).choose) := by ring
    rw [he]
    exact neg_mem s1

/-- ★★★★★★★★★★★★★★★★★★★★**`veluAnalyticX` は `Λ′`-周期的**——★**無条件**
（代表系であることだけから）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 599 では `σ` を仮定として受けていたが、第 602 でそれが**代表系の定義から作れた**。 -/
theorem veluAnalyticX_periodic (P P' : PeriodPair) (T : Finset ℂ) (c : ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (z : ℂ) (l : ℂ) (hl : l ∈ P'.lattice) :
    veluAnalyticX P T c (z + l) = veluAnalyticX P T c z := by
  obtain ⟨σ, hmem, hinj, hsurj, hshift⟩ := exists_veluShiftPerm P P' T hT hrep l hl
  exact veluAnalyticX_shift P T c l σ hmem hinj hsurj hshift z

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側
——代表系であることだけから**。

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − Σ_{w ∈ T∖{0}} ℘_Λ(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**仮定は `T` が `Λ′/Λ` の代表系であること、それだけである**:

| 仮定 | 内容 |
|---|---|
| `hle` | `Λ ⊆ Λ′` |
| `h0T` | `0 ∈ T` |
| `hT` | `T ⊆ Λ′` |
| `hrep` | `p ∈ Λ′` ならちょうど 1 つの `w₀ ∈ T` が `p + w₀ ∈ Λ` |

★★★★★周期性・極の打ち消し・正規化はすべて塞がった
（第 598 Liouville、第 599 周期性、第 600 極、第 601 正規化、第 602 置換の構成）。

☆残るのは、この `℘` の等式を `Found/GenEll/Velu.lean` の**代数側 `veluXGen`**
（第 591）へ翻訳することだけである——`℘` の加法定理が要る（mathlib に無い）。 -/
theorem weierstrassP_eq_velu_of_rep (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
  weierstrassP_eq_velu P P' hle T h0T hT hrep
    (fun y l hl => veluAnalyticX_periodic P P' T _ hT hrep y l hl) z

def exists_veluShiftPerm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(平行移動は代表系を置換する——σ は代表系の定義から作れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu_of_rep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側——仮定は T が Λ′/Λ の代表系であることだけ)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu_of_rep.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆この ℘ の等式を Found/GenEll/Velu.lean の代数側 veluXGen(第 591)へ" ++
       "翻訳すること。★℘ の加法定理が要る(mathlib に無い、2026-08-29 に測定)") 8,
    .implicitStep
      ("★★★★到達点(2026-08-29、第 602): Vélu の公式の解析側が" ++
       "「T が Λ′/Λ の代表系である」だけから従う形になった。" ++
       "第 598 Liouville・第 599 周期性・第 600 極・第 601 正規化・第 602 置換の構成") 9 ]

/-! ## ★★★★★★★★★★★★★★★★2-捩れ点 -/

/-- ★★★★★★★★★★★★**2-捩れ点では `℘′` が消える**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★証明は 3 行: `℘′` は奇であり、`2z ∈ Λ` なら `−z ≡ z (mod Λ)` なので
`℘′(z) = ℘′(−z) = −℘′(z)`。

★★`Found/GenEll/Velu.lean` の `veluV2`（2-捩れの場合分け）と
`velu2_omega` の仮定 `2y₀ + a₁x₀ + a₃ = 0` は、`latticeCurve` では
`a₁ = a₃ = 0` なので `2·(℘′(z)/2) = ℘′(z) = 0`——**本定理そのもの**である。 -/
theorem derivWeierstrassP_eq_zero_of_two_mem (P : PeriodPair) (z : ℂ)
    (h2 : 2 * z ∈ P.lattice) : P.derivWeierstrassP z = 0 := by
  have hneg : P.derivWeierstrassP (-z) = -P.derivWeierstrassP z := P.derivWeierstrassP_neg z
  have hper : P.derivWeierstrassP (z + (-(2 * z))) = P.derivWeierstrassP z :=
    P.derivWeierstrassP_add_coe z ⟨-(2 * z), neg_mem h2⟩
  have hz : z + (-(2 * z)) = -z := by ring
  rw [hz, hneg] at hper
  linear_combination -hper / 2

/-- ★★★★★★★★★★**2-捩れ点の `y` 座標は `0`**——`Velu.lean` の 2-捩れの場合分けの中身。 -/
theorem latticePointY_eq_zero_of_two_mem (P : PeriodPair) (z : ℂ)
    (h2 : 2 * z ∈ P.lattice) : latticePointY P z = 0 := by
  simp [latticePointY, derivWeierstrassP_eq_zero_of_two_mem P z h2]

/-- ★★★★★★★★★★★★**2-捩れ点の `x` 座標は `4x³ − g₂x − g₃` の根**。

★`℘′² = 4℘³ − g₂℘ − g₃` の左辺が `0` になるから。
☆`latticeCurve P = ⟨0,0,0,−g₂/4,−g₃/4⟩` の 2-捩れ点はちょうどここである。 -/
theorem cubic_eq_zero_of_two_mem (P : PeriodPair) (z : ℂ) (hz : z ∉ P.lattice)
    (h2 : 2 * z ∈ P.lattice) :
    4 * (latticePointX P z) ^ 3 - P.g₂ * (latticePointX P z) - P.g₃ = 0 := by
  have hsq := P.derivWeierstrassP_sq z hz
  rw [derivWeierstrassP_eq_zero_of_two_mem P z h2] at hsq
  simp only [latticePointX]
  linear_combination -hsq

def derivWeierstrassP_eq_zero_of_two_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2-捩れ点では ℘′ が消える——℘′ は奇だから。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
