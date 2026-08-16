import ABC3.Meta.Claim
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic

/-!
# [NCBelyi] Lemma 2.2 —— Belyi 型多項式は合成で閉じる(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.3。
**260 dpi 目視確認 2026-08-17**。

原文 (NCBelyi p.3):
> Lemma 2.2. (Belyi Maps Noncritical at Prescribed Rational Points)

## ★★原文の証明は 6 行の粗描である

> Indeed, we induct on the cardinality |S| of S and apply Lemma 2.1 [with,
> say, C = 2] to the set λ · S ⊆ P1Q, for some appropriate positive rational number
> λ. ... Thus, by applying the induction hypothesis and **composing** the resulting
> morphisms P1Q → PQ1 , we conclude the existence of an “f”, “φ”

★**「composing」が構造の要**である——帰納法の各段で得た射を合成して
最終的な `φ` を作る。合成しても
(a) `φ(S) ⊆ {0,1,∞}` と (c) `ℙ¹\{0,1,∞}` 上不分岐 が保たれねばならない。

## ★本ファイルが取るもの

> **`IsBelyiPoly` は合成で閉じる**

`IsBelyiPoly f` を「(1) 非定数、(2) `f(0), f(1) ∈ {0,1}`、
(3) `ℂ` 上の臨界**値**が `{0,1}` に入る」と定めると、
`g ∘ f` もまた `IsBelyiPoly` である。

★★**証明は連鎖律 `(g∘f)′ = g′(f)·f′` だけ**である。
積が `0` になるのは
- `f′(x) = 0` —— このとき `f(x) ∈ {0,1}` なので `g(f(x)) ∈ {0,1}`
  (**ここで条件 (2) が効く**)
- `g′(f(x)) = 0` —— このとき `g(f(x)) ∈ {0,1}`

★★★**条件 (2)(`f(0), f(1) ∈ {0,1}`)は (3) だけでは出ない**——
合成の第 1 の場合でちょうど要る。
★原文が (a) と (c) を**両方**要求する理由がここにある。

## ★∞ の扱い

原文は `ℙ¹` の射として述べるが、多項式写像は `∞ ↦ ∞` で全分岐する。
`∞ ∈ {0,1,∞}` なので (a)(c) はどちらも満たされ、
**有限部分だけを見れば足りる**。本ファイルはその有限部分を扱う。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-- **Belyi 型多項式** —— 原文 (a)(c) の有限部分。

- `deg_pos` —— 非定数(原文「nonconstant polynomial」)
- `maps_zero` / `maps_one` —— `f({0,1}) ⊆ {0,1}`(原文 (a) の一部)
- `crit` —— 臨界**値**が `{0,1}` に入る(原文 (c)) -/
structure IsBelyiPoly (f : ℚ[X]) : Prop where
  deg_pos : 0 < f.natDegree
  maps_zero : f.eval 0 = 0 ∨ f.eval 0 = 1
  maps_one : f.eval 1 = 0 ∨ f.eval 1 = 1
  crit : ∀ x : ℂ, (derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
    (f.map (algebraMap ℚ ℂ)).eval x = 0 ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1

/-- ★`ℚ` での評価を `ℂ` へ運ぶ。 -/
theorem eval_map_ratCast (f : ℚ[X]) (a : ℚ) :
    (f.map (algebraMap ℚ ℂ)).eval (algebraMap ℚ ℂ a) = algebraMap ℚ ℂ (f.eval a) := by
  rw [eval_map, eval₂_at_apply]

/-- ★★★**Belyi 型多項式は合成で閉じる**。

原文 (NCBelyi p.4):
> existence of an “f”, “φ” as in the statement of the present Lemma 2.2.

★★証明は**連鎖律 1 回**である。
`(g∘f)′ = g′(f)·f′` が `0` になるのは `f′(x) = 0` か `g′(f(x)) = 0` のときで、
前者では `f(x) ∈ {0,1}` から `g(f(x)) ∈ {0,1}`(**条件 (2) が効く**)、
後者では `g` の条件 (3) から直ちに従う。 -/
theorem IsBelyiPoly.comp {f g : ℚ[X]} (hf : IsBelyiPoly f) (hg : IsBelyiPoly g) :
    IsBelyiPoly (g.comp f) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [natDegree_comp]
    exact Nat.mul_pos hg.deg_pos hf.deg_pos
  · rw [eval_comp]
    rcases hf.maps_zero with h | h <;> rw [h]
    · exact hg.maps_zero
    · exact hg.maps_one
  · rw [eval_comp]
    rcases hf.maps_one with h | h <;> rw [h]
    · exact hg.maps_zero
    · exact hg.maps_one
  · intro x hx
    rw [Polynomial.map_comp] at hx ⊢
    rw [derivative_comp, eval_mul, eval_comp] at hx
    rw [eval_comp]
    rcases mul_eq_zero.1 hx with h1 | h2
    · -- `f′(x) = 0` なので `f(x) ∈ {0,1}`
      rcases hf.crit x h1 with h | h <;> rw [h]
      · have hz := eval_map_ratCast g 0
        rcases hg.maps_zero with hg0 | hg0 <;> rw [hg0] at hz <;>
          simp only [map_zero, map_one] at hz
        · exact Or.inl hz
        · exact Or.inr hz
      · have ho := eval_map_ratCast g 1
        rcases hg.maps_one with hg1 | hg1 <;> rw [hg1] at ho <;>
          simp only [map_zero, map_one] at ho
        · exact Or.inl ho
        · exact Or.inr ho
    · -- `g′(f(x)) = 0`
      exact hg.crit _ h2

/-! ## ★出典の紐付け(`.src`) -/

def IsBelyiPoly.comp.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 3,
    item := "Lemma 2.2(合成で閉じることのみ)",
    sectionId := "ncbelyi-lemma-2-2" }

end ABC3.Found.NCBelyi
