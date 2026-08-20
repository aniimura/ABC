/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Main
import ABC3.Found.SixExp.Setup
import ABC3.Meta.Claim

/-!
# six exponentials theorem

★★★★★★★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーンの頂点。
[FrdI] `Lemma 6.5, (ii)` が原文で Lang の定理に送っているものである。

> `x₁, x₂` が ℚ 上一次独立、`y₁, y₂, y₃` が ℚ 上一次独立な複素数なら、
> `exp(xⱼ·yᵢ)`(6 個)の少なくとも 1 つは超越数である。

## ★証明の骨(背理法)

6 個がすべて代数的だと仮定する。

1. **数体を作る**(`Setup.lean`)—— `K = ℚ(6 個)` は有限次、したがって数体。
   共通分母 `b ∈ 𝓞_K` を取り、`Ê j i = b·E j i ∈ 𝓞_K` とする。
2. **矛盾**(`Main.lean` の `sixExp_contradiction`)——
   Siegel で補助関数の係数を取り、外挿で**すべての格子点で消える**ことを示し、
   指標の一次独立性で**係数がすべて 0**になって Siegel の「0 でない」に反する。

★数え上げ(`ParamsGap.lean`)がちょうど合うのは **`L' ≈ N^{3/2}`** のときで、
そこが `y` を 3 個取る理由である(零点は `N³` 個、増大は `N` の 1 乗)。
-/

namespace ABC3.Found.SixExp

open Complex NumberField

/-- ★★★★★★★**six exponentials theorem**(Lang / Ramachandra)。

`x₁, x₂` が ℚ 上一次独立、`y₁, y₂, y₃` が ℚ 上一次独立な複素数なら、
`exp(xⱼ·yᵢ)`(6 個)の少なくとも 1 つは超越数である。 -/
theorem six_exponentials {x : Fin 2 → ℂ} {y : Fin 3 → ℂ}
    (hx : LinearIndependent ℚ x) (hy : LinearIndependent ℚ y) :
    ∃ (j : Fin 2) (i : Fin 3), Transcendental ℚ (Complex.exp (x j * y i)) := by
  classical
  by_contra hcon
  push_neg at hcon
  set v : Fin 2 → Fin 3 → ℂ := sixExpVals x y with hvdef
  have halg : ∀ j i, IsAlgebraic ℚ (v j i) := by
    intro j i
    have h := hcon j i
    rw [Transcendental, not_not] at h
    exact h
  haveI hfd : FiniteDimensional ℚ ↥(sixField v) := sixField_finiteDimensional v halg
  haveI hnf : NumberField ↥(sixField v) := sixField_numberField v halg
  haveI : DecidableEq (↥(sixField v) →+* ℂ) := Classical.decEq _
  set σ : ↥(sixField v) →+* ℂ := sixFieldHom v with hσdef
  obtain ⟨b, hb, hbspec⟩ :=
    exists_common_denom ↥(sixField v) (fun j i => sixFieldElt v j i)
  choose Eh hEh using hbspec
  refine sixExp_contradiction σ x y hx hy b hb Eh ?_
  intro j i
  rw [hEh j i, map_mul, hσdef, sixFieldHom_elt]

/-! ### ★出典の紐付け -/

/-- ★locator —— [FrdI] `Lemma 6.5, (ii)` が引く Lang の定理。 -/
def six_exponentials.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116, item := "Lemma 6.5, (ii) — six exponentials theorem",
    sectionId := "frdi-lemma-6-5" }

end ABC3.Found.SixExp
