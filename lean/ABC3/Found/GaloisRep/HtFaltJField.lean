/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtFaltJ
import ABC3.Meta.Claim

/-!
# `ht^Falt` は `j` だけで決まる —— **体が違ってもよい形**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★これは何か

`HtFaltJ.lean`（`§9-1166`、第 739）は**同じ体の上で**「`j` が同じなら `ht^Falt` が同じ」
を証明した。本ファイルは**共通の体 `L₃` へ上げられるなら体が違ってもよい**形にする。

    L₁ ⊆ L₃ ⊇ L₂,  E₁/L₁ 半安定,  E₂/L₂ 半安定,  j(E₁) = j(E₂)（`L₃` の中で）
    ⟹ ht^Falt(E₁) = ht^Falt(E₂)

## ★★★機構——3 段

1. `ht^Falt(E₁) = ht^Falt(E₁×L₃)`（`htFaltOf_baseChange_of_semistable`、`§9-1165`、第 738、
   仮説は**下での半安定性だけ**に弱めた——`§9-1167`、第 740）
2. 同様に `E₂`
3. `L₃` の上で `j` が同じなので `ht^Falt(E₁×L₃) = ht^Falt(E₂×L₃)`
   （`htFaltOf_congr_j_of_maxJ`、`§9-1166`、第 739）

★★2 と 3 を繋ぐ鍵が `minDeltaExp_eq_maxJ_baseChange`（第 740）である
——上では `SemistableAt` そのものは持ち上げられないが、
**`v_P(Δ_min) = max(0, −v_P(j))` という等式は持ち上げられる**。

## ☆残っていること

`SSCurve`（`Found/GenEll/EllModuliObjects.lean`）2 つに対して `L₃` を実際に作る段
——`ℂ` の中で `emb₁(L₁)` と `emb₂(L₂)` が生成する部分体を取り、
それが数体であることと本ファイルの型クラスをすべて満たすことを確かめる。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve

section Common

variable (L₁ L₂ L₃ : Type) [Field L₁] [NumberField L₁] [Field L₂] [NumberField L₂]
  [Field L₃] [NumberField L₃]
  [Algebra L₁ L₃] [Algebra (𝓞 L₁) L₃] [IsScalarTower (𝓞 L₁) L₁ L₃]
  [IsScalarTower (𝓞 L₁) (𝓞 L₃) L₃] [Module.Finite (𝓞 L₁) (𝓞 L₃)] [IsScalarTower ℚ L₁ L₃]
  [Algebra L₂ L₃] [Algebra (𝓞 L₂) L₃] [IsScalarTower (𝓞 L₂) L₂ L₃]
  [IsScalarTower (𝓞 L₂) (𝓞 L₃) L₃] [Module.Finite (𝓞 L₂) (𝓞 L₃)] [IsScalarTower ℚ L₂ L₃]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**共通の体へ上げられるなら、半安定な曲線の `ht^Falt` は `j` だけで決まる**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★これが `EllModuliData` の `faltingsHeight : EllClass → ℝ` の
**well-defined 性そのもの**である（`EllClass := ℂ`、`cls E = j(E)`）。 -/
theorem htFaltOf_congr_j_of_common (E₁ : WeierstrassCurve L₁) [E₁.IsElliptic]
    (E₂ : WeierstrassCurve L₂) [E₂.IsElliptic]
    (h₁ : ∀ p : HeightOneSpectrum (𝓞 L₁), SemistableAt p E₁)
    (h₂ : ∀ p : HeightOneSpectrum (𝓞 L₂), SemistableAt p E₂)
    (hj : algebraMap L₁ L₃ E₁.j = algebraMap L₂ L₃ E₂.j) :
    htFaltOf L₁ E₁ = htFaltOf L₂ E₂ := by
  have e1 : htFaltOf L₃ (E₁.baseChange L₃) = htFaltOf L₁ E₁ :=
    htFaltOf_baseChange_of_semistable L₁ L₃ E₁ h₁
  have e2 : htFaltOf L₃ (E₂.baseChange L₃) = htFaltOf L₂ E₂ :=
    htFaltOf_baseChange_of_semistable L₂ L₃ E₂ h₂
  have hj1 : (E₁.baseChange L₃).j = algebraMap L₁ L₃ E₁.j :=
    WeierstrassCurve.map_j (W := E₁) (f := algebraMap L₁ L₃)
  have hj2 : (E₂.baseChange L₃).j = algebraMap L₂ L₃ E₂.j :=
    WeierstrassCurve.map_j (W := E₂) (f := algebraMap L₂ L₃)
  have hjj : (E₁.baseChange L₃).j = (E₂.baseChange L₃).j := by rw [hj1, hj2, hj]
  have key := htFaltOf_congr_j_of_maxJ (E₁.baseChange L₃) (E₂.baseChange L₃)
    (fun P => minDeltaExp_eq_maxJ_baseChange L₁ L₃
      (HeightOneSpectrumUnder (A := 𝓞 L₁) P) P E₁ (h₁ _))
    (fun P => minDeltaExp_eq_maxJ_baseChange L₂ L₃
      (HeightOneSpectrumUnder (A := 𝓞 L₂) P) P E₂ (h₂ _)) hjj
  rw [← e1, ← e2, key]

/-- ★★★★★★★★★★★★★★**`deg∞` も同じ**——共通の体へ上げられるなら `j` だけで決まる。 -/
theorem degInfOf_congr_j_of_common (E₁ : WeierstrassCurve L₁) [E₁.IsElliptic]
    (E₂ : WeierstrassCurve L₂) [E₂.IsElliptic]
    (h₁ : ∀ p : HeightOneSpectrum (𝓞 L₁), SemistableAt p E₁)
    (h₂ : ∀ p : HeightOneSpectrum (𝓞 L₂), SemistableAt p E₂)
    (hj : algebraMap L₁ L₃ E₁.j = algebraMap L₂ L₃ E₂.j) :
    degInfOf L₁ E₁ = degInfOf L₂ E₂ := by
  have e1 : degInfOf L₃ (E₁.baseChange L₃) = degInfOf L₁ E₁ :=
    degInfOf_baseChange_of_semistable L₁ L₃ E₁ h₁
  have e2 : degInfOf L₃ (E₂.baseChange L₃) = degInfOf L₂ E₂ :=
    degInfOf_baseChange_of_semistable L₂ L₃ E₂ h₂
  have hj1 : (E₁.baseChange L₃).j = algebraMap L₁ L₃ E₁.j :=
    WeierstrassCurve.map_j (W := E₁) (f := algebraMap L₁ L₃)
  have hj2 : (E₂.baseChange L₃).j = algebraMap L₂ L₃ E₂.j :=
    WeierstrassCurve.map_j (W := E₂) (f := algebraMap L₂ L₃)
  have hjj : (E₁.baseChange L₃).j = (E₂.baseChange L₃).j := by rw [hj1, hj2, hj]
  have key := degInfOf_congr_j_of_maxJ (E₁.baseChange L₃) (E₂.baseChange L₃)
    (fun P => minDeltaExp_eq_maxJ_baseChange L₁ L₃
      (HeightOneSpectrumUnder (A := 𝓞 L₁) P) P E₁ (h₁ _))
    (fun P => minDeltaExp_eq_maxJ_baseChange L₂ L₃
      (HeightOneSpectrumUnder (A := 𝓞 L₂) P) P E₂ (h₂ _)) hjj
  rw [← e1, ← e2, key]

end Common

/-! ## ★出典の紐付け(`.src`) -/

def htFaltOf_congr_j_of_common.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(共通の体へ上げられるなら半安定な曲線の ht^Falt は j だけで決まる。★無条件)",
    sectionId := "genell-prop-3-4" }

def htFaltOf_congr_j_of_common.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは SSCurve 2 つに対して L₃ を実際に作る段——ℂ の中で " ++
       "emb₁(L₁) と emb₂(L₂) が生成する部分体を取り、それが数体であることと " ++
       "本ファイルの型クラスをすべて満たすことを確かめる") 6 ]

end ABC3.Found.GaloisRep
