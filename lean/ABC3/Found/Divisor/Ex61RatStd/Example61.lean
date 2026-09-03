/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex61Units
import ABC3.Found.FrdI.Thm64RatStd
import ABC3.Found.Divisor.Ex61RatStd.Proposition55

/-!
# Ex61RatStd —— `[FrdI] Example 6.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI ABC3.Meta
universe u
variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]
  (DK : Set (PrimeDivisorPt V))
  [∀ L : FinSub V.functionField Kbar, IsLocallyNoetherian (normObj V L)]
  [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)]
  [IsGalois V.functionField Kbar]
  (hkq : IsKQCartier V DK
    (fun (L : FinSub V.functionField Kbar) _ => normObj_isNormalScheme V L))

/-! ### ★★★★★★★★項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
(`tools/frdi-progress.mjs` の規則)。下の 2 つは、条がすべて閉じたので置く。 -/

/-- ★★★★★★★★**[FrdI] Example 6.1** —— 条がすべて実装された。

| 原文の主張 | 宣言 |
|---|---|
| `Φ(L)^gp ⊆ ℤ[D_L]` は `V[L]` 上の Cartier 因子の群 | `cartierOnDL`(`NormKQC.lean`) |
| `Φ(L)^pf = ℚ≥0[D_L] ⊆ ℚ[D_L] = (Φ(L)^pf)^gp` | `CartierPf.lean` |
| `Φ(L)` は perf-factorial | `isPerfFactorial_effSub` |
| `Prime(Φ(L)) ≃ D_L` | `effSubPrimeEquiv` |
| 台は `D_L` の有限部分集合ちょうど | `exists_effSub_support_eq` |
| `B(L) → Φ(L)^gp` と `L` についての関手性 | `divBHom` / `divBHom_bHom` / `bMonoidFunctor` |
| `Theorem 5.2, (ii)` で model Frobenioid `C_{V,K̄,D_K}` | `ex61Frobenioid` |
| isotropic 型 | `ex61Frobenioid_isotropicType` |
| birationally Frobenius-normalized 型 | `ex61Frobenioid_biratFrobNormalizedType` |
| 射は `(Spec L → Spec M, d, L^⊗d → M|_{V[L]})` の 3 つ組 | `Ex61Morph.lean` |
| `𝒪^×(A) = 𝒪^▷(A) = k_L^×` | `ex61_otri_le_otimes` / `SchemeHartogs.lean` | -/
def ex_6_1.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
