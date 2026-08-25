/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithFundamental
import ABC3.Found.Divisor.ArithFrobenioid
import ABC3.Found.FrdI.Def24ScTransport

/-!
# 節点 `rlf-flat`（`S = ℝ≥0`）—— 算術ではトレース写像で通る

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

## ★★何のために要るか

`Theorem 6.4, (i)` の 5 圏のうち **`𝒞^rlf`** を算術で出すには
`Φ^rlf : MonoidOn 𝒟` の条件 (a)

    f : M ↪ N  ⟹  id ⊗ f : ℝ≥0 ⊗_ℕ M ↪ ℝ≥0 ⊗_ℕ N

が要る（節点 `rlf-flat`）。★一般には「`ℝ≥0` が ℕ 上平坦」という
凸幾何の議論（濾過余極限・単体錐）になるが、

★★**`f` に `k` 倍の引き戻しがあれば平坦性は要らない**
（`isCharacteristicallyInjective_scMap_of_nsmul_retraction`、`Def24ScTransport.lean`）。

## ★★★算術での引き戻し —— トレース写像

`Φ.map α` は素点の引き戻し `arithExtend`（局所次数倍つき）なので、
**上にある素点についての和**

    tr := Finsupp.mapDomain resPlace

を取ると、相対基本等式 `Σ_{V | v} localDeg V = [M:L]`（`ArithFundamental.lean`）から

    tr ∘ arithExtend = [M : L] · id

が出る（`mapDomain_arithExtend`）。★`arithDivGroup` を保つことは
相対ノルムの冪表示 `N(𝔓) = N(𝔭)^f` による（`mapDomain_mem_arithDivGroup`）。
-/

namespace ABC3.Found.Divisor

open CategoryTheory ABC3.Found.FrdI
open scoped NNReal

universe v u w

/-- ★★★★★★**`ArithDatum` に `k` 倍の引き戻しがあれば条件 (a) が出る**。

★これが節点 `rlf-flat`（`S = ℝ≥0`）の迂回路の、`ArithDatum` の水準での形である。 -/
theorem arithDatum_isCharInj_scMap {D : Type u} [Category.{v} D]
    (Δ : ArithDatum.{v, u, w} D) (hD : IsOfFSMType D)
    {A B : D} (α : B ⟶ A)
    (hsharpN : IsSharp (ScT ℝ≥0 (effR (Δ.grp B))))
    (tr : (Δ.primes B →₀ ℝ) →+ (Δ.primes A →₀ ℝ))
    (htr_mem : ∀ {y : Δ.primes B →₀ ℝ}, y ∈ Δ.grp B → tr y ∈ Δ.grp A)
    (htr_nonneg : ∀ {y : Δ.primes B →₀ ℝ}, (∀ s, 0 ≤ y s) → ∀ s, 0 ≤ (tr y) s)
    (k : ℕ) (hk : 0 < k)
    (htr : ∀ x : Δ.primes A →₀ ℝ, tr (Δ.pull α x) = k • x) :
    IsCharacteristicallyInjective (scMap (S := ℝ≥0) ((Δ.phi hD).map α)) := by
  classical
  refine isCharacteristicallyInjective_scMap_of_nsmul_retraction hsharpN hk
    (g := AddMonoidHom.codRestrict (tr.comp (effR (Δ.grp B)).subtype) _
      (fun y => mem_effR.mpr ⟨htr_mem (mem_effR.mp y.2).1,
        htr_nonneg (mem_effR.mp y.2).2⟩)) ?_
  intro x
  exact Subtype.ext (htr _)

/-! ### ★出典の紐付け -/

def arithDatum_isCharInj_scMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — ArithDatum に k 倍の引き戻しがあれば Φ^rlf の条件 (a) が出る",
    sectionId := "frdi-thm-6-4" }

def arithDatum_isCharInj_scMap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isCharacteristicallyInjective_scMap_of_nsmul_retraction(迂回路の一般形)"
      (.inProject "ABC3" "ABC3.Found.FrdI.isCharacteristicallyInjective_scMap_of_nsmul_retraction") 51,
    .citation "[ABC3]" "mapDomain_arithExtend(トレース写像が [M:L] 倍の引き戻し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.mapDomain_arithExtend") 114,
    .citation "[ABC3]" "mapDomain_mem_arithDivGroup(arithDivGroup を保つ)"
      (.inProject "ABC3" "ABC3.Found.Divisor.mapDomain_mem_arithDivGroup") 114,
    .implicitStep
      ("★節点 rlf-flat(S = ℝ≥0)の凸幾何の議論(ℝ≥0 が ℕ 上の有限生成自由半加群の" ++
       "濾過余極限であること)は、§6 の算術の実例では要らない") 114 ]

end ABC3.Found.Divisor
