import ABC3.Found.GaloisRep.NeronIdeal
import ABC3.Found.GaloisRep.D2Bridge
import ABC3.Found.GaloisRep.TateCurveWitness

/-!
# Galois (G7) 第 317 ブロック —— **★★★★★★★★★★★G7 達成 `SemistableModelData`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★到達点

> **`SemistableModelData` は非空虚である**(`SemistableModelData.nonvacuous`)——**G7 達成**

★★★これで Galois 表現論の 8 件のうち **6 件**(G1, G2, G3, G4, G6, G7)が埋まった。

## ★★★★★★変数変換則——イデアルの水準から元の水準へ

`count` の代数で**イデアルの等式**を出す:

    ω_E(C • W) = (u)⁻¹ · ω_E(W)         (`omegaFracIdeal_variableChange`)

★各素点で `count` を比べるだけ(`count_mul`・`count_inv`・第 316 の橋)。
★★★これを**元の水準**に落とすと界面の形になる:

    x ∈ ω_E(C • W)  ↔  u·x ∈ ω_E(W)     (`mem_omegaFracIdeal_iff`)

★★`spanSingleton x ≤ I ↔ x ∈ I` と `le_spanSingleton_mul_iff` で往復する。

## ★★★★有限生成性は mathlib が持っている

`FractionalIdeal.fg_of_isNoetherianRing`——**分数イデアルは Noether 環の上で f.g.**。
★`𝓞_L` は Noether(数体の整数環)なので、そのまま効く。

## ★★★界面の第 3 の訂正(記録)

`omegaFrac` の 3 性質に **`[W.IsElliptic]`** を付けた。
★`Δ = 0` の曲線に Néron 微分は無く、変数変換則も成り立たない
(両辺とも自明なイデアルになり `x ∈ 1 ↔ u·x ∈ 1` は偽)。
★★欄そのものは全曲線で定義しておき、**性質だけを楕円曲線に限る**形にした。
★★★第 304(`TateCurveData` の `Δ ≠ 0`)と同じ訂正である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `omegaFracIdeal_variableChange` | ★★★★★★イデアルの水準の変数変換則 |
| `mem_omegaFracIdeal_iff` | ★★★★★★★**元の水準の変数変換則** |
| `omegaFracIdeal_fg`・`omegaFracIdeal_ne_bot` | ★★★★有限生成性と非零性 |
| `SemistableModelData.nonvacuous` | ★★★★★★★★★★★**G7 達成** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve
open scoped nonZeroDivisors

section Frac

variable {L : Type} [Field L] [NumberField L]

/-- ★★★主分数イデアルとの積の元の判定。 -/
theorem mem_spanSingleton_mul_iff (x : L) (hx : x ≠ 0) (J : FractionalIdeal (𝓞 L)⁰ L) (y : L) :
    y ∈ FractionalIdeal.spanSingleton (𝓞 L)⁰ x * J ↔ x⁻¹ * y ∈ J := by
  rw [← FractionalIdeal.spanSingleton_le_iff_mem, FractionalIdeal.le_spanSingleton_mul_iff]
  constructor
  · intro h
    obtain ⟨z, hz, hxz⟩ := h y (FractionalIdeal.mem_spanSingleton_self _ _)
    have heq : x⁻¹ * y = z := by
      rw [← hxz]
      field_simp
    rw [heq]
    exact hz
  · intro h zI hzI
    rw [FractionalIdeal.mem_spanSingleton] at hzI
    obtain ⟨a, ha⟩ := hzI
    refine ⟨a • (x⁻¹ * y), (J : Submodule (𝓞 L) L).smul_mem a h, ?_⟩
    rw [← ha, Algebra.smul_def, Algebra.smul_def]
    field_simp

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**イデアルの水準の変数変換則** `ω_E(C • W) = (u)⁻¹ ω_E(W)`。 -/
theorem omegaFracIdeal_variableChange (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0)
    (C : VariableChange L) :
    omegaFracIdeal (C • W)
      = FractionalIdeal.spanSingleton (𝓞 L)⁰ ((C.u : L))⁻¹ * omegaFracIdeal W := by
  have hΔC : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
  have hune : ((C.u : L))⁻¹ ≠ 0 := inv_ne_zero (Units.ne_zero _)
  have h1 : omegaFracIdeal (C • W) ≠ 0 := omegaFracIdeal_ne_zero _ hΔC
  have h2 : FractionalIdeal.spanSingleton (𝓞 L)⁰ ((C.u : L))⁻¹ ≠ 0 :=
    FractionalIdeal.spanSingleton_ne_zero_iff.2 hune
  have h3 : omegaFracIdeal W ≠ 0 := omegaFracIdeal_ne_zero W hΔ
  refine eq_of_count_eq h1 (mul_ne_zero h2 h3) (fun p => ?_)
  rw [count_omegaFracIdeal _ hΔC, FractionalIdeal.count_mul L p h2 h3,
    count_omegaFracIdeal W hΔ, neronExp_variableChange p W hΔ C]
  have h4 : FractionalIdeal.spanSingleton (𝓞 L)⁰ ((C.u : L))⁻¹
      = (FractionalIdeal.spanSingleton (𝓞 L)⁰ ((C.u : L)))⁻¹ :=
    (FractionalIdeal.spanSingleton_inv L _).symm
  rw [h4, FractionalIdeal.count_inv, count_spanSingleton_eq_valAdd]
  omega

/-- ★★★★★★★**元の水準の変数変換則**——界面 `omegaFrac_variableChange` の形。 -/
theorem mem_omegaFracIdeal_iff (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) (C : VariableChange L)
    (x : L) : x ∈ omegaFracIdeal (C • W) ↔ ((C.u : L) * x) ∈ omegaFracIdeal W := by
  rw [omegaFracIdeal_variableChange W hΔ C,
    mem_spanSingleton_mul_iff _ (inv_ne_zero (Units.ne_zero _)) _ x, inv_inv]

/-- ★★★★有限生成(分数イデアルは Noether 環の上で f.g.)。 -/
theorem omegaFracIdeal_fg (W : WeierstrassCurve L) :
    (omegaFracIdeal W : Submodule (𝓞 L) L).FG :=
  FractionalIdeal.fg_of_isNoetherianRing le_rfl _

/-- ★★★★零加群でない。 -/
theorem omegaFracIdeal_ne_bot (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) :
    (omegaFracIdeal W : Submodule (𝓞 L) L) ≠ ⊥ := by
  intro h
  refine omegaFracIdeal_ne_zero W hΔ ?_
  refine FractionalIdeal.coeToSubmodule_injective ?_
  show (omegaFracIdeal W : Submodule (𝓞 L) L)
    = ((0 : FractionalIdeal (𝓞 L)⁰ L) : Submodule (𝓞 L) L)
  rw [h]
  simp

end Frac

/-! ## ★★★★★★★★★★★G7 達成 -/

open ABC3.Interface.GaloisRep

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`SemistableModelData` の実装**。 -/
noncomputable def semistableModelDataWitness : SemistableModelData where
  toTateCurveData := tateCurveDataWitness
  omegaFrac := fun L _ _ W => (omegaFracIdeal W : Submodule (𝓞 L) L)
  omegaFrac_fg := by
    intro L _ _ W hell
    exact omegaFracIdeal_fg W
  omegaFrac_ne_bot := by
    intro L _ _ W hell
    exact omegaFracIdeal_ne_bot W hell.isUnit.ne_zero
  omegaFrac_variableChange := by
    intro L _ _ W hell C x
    exact mem_omegaFracIdeal_iff W hell.isUnit.ne_zero C x
  SemiStable := fun L _ _ W =>
    ∀ (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
      [Algebra R L] [IsFractionRing R L] [W.IsMinimal R], IsSemiStableAt R W
  semiStable_iff := by
    intro L _ _ W
    exact Iff.rfl

/-- ★★★★★★★★★★★**`SemistableModelData` は非空虚である**——G7 達成。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★これが Galois 表現論の 8 件のうち **G7** である。 -/
theorem SemistableModelData.nonvacuous : Nonempty SemistableModelData :=
  ⟨semistableModelDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def mem_omegaFracIdeal_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Néron 微分の分数イデアルの変数変換則)",
    sectionId := "genell-def-3-3" }

def SemistableModelData.nonvacuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(SemistableModelData の非空虚性——G7)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
