/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.MobiusRedDeg
import ABC3.Found.NCBelyi.ConjStable
import ABC3.Found.NCBelyi.CoeffBound
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★`Lemma 2.4` の (b) —— 1 段ぶんの分離（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> bounded by some fixed expression in d0. Thus, for a suitable choice of C, it follows

## ★★★★★★★★★★★★★★★★★★★★★★これは何か

`Found/NCBelyi/Lemma24Package.lean` の測定で、`Lemma 2.4` に残るのは
**(b)（分離）** だけになった:

    `∀ x ∈ S,  h(x) ≠ h(β)`   かつ   `∀ w, h′(w) = 0 → h(w) ≠ h(β)`

★★★**本ファイルはその「1 段ぶん」を取る**——原文の

> Thus, for a suitable choice of C, it follows that f0(β) /∈S′

そのものである。

## ★機構 —— 部品を 3 つ繋ぐだけ

| 部品 | 場所 |
|---|---|
| `S` が `Gal`-安定かつ `‖α‖ ≤ 1` なら `f₀ ≔ minpoly x₀` の根も単位円板の中 | `ConjStable.lean` の `norm_roots_minpoly_le_one` |
| `f₀′` の根も単位円板の中（Gauss–Lucas） | `CoeffBound.lean` の `norm_root_derivative_le_one` |
| `‖β‖ > 3` なら `f₀(β) ≠ f₀(x)`（`‖x‖ ≤ 1` または `f₀′(x) = 0`） | `CoeffBound.lean` の `eval_ne_of_root_or_root_derivative` |

★★`C = 3` で足りることは第 405-406 で測ってある（原文の「`d₀` に相対的に十分大きい `C`」は
`f₀` がモニックで根が単位円板の中である以上、**`d₀` に依らない**）。

## ★★★★★これで `Lemma 2.4` に残るのは配管 2 つ

1. **`ℚ`-Möbius は `Gal`-安定性を保つ**——`σ` が `ℚ` 係数なら
   `conjSet (σ x) = σ (conjSet x)`。★これがあれば
   `exists_rat_normalization`（`|α| ≤ 1`、`|β| ≥ 4`）で正規化したあとも
   `separation_step` の仮定がそろう。
   ★★測度の側（`redDeg`／`maxRedDeg`）は `MobiusRedDeg.lean` で済んでいる。
2. **合成を `ℙ¹` の有理写像として組み立てる**——原文の `f(x) ∈ ℚ(x)` は
   Möbius と多項式の交互合成である。★`Separation.lean` の `P1C` がその受け皿。

★★★★★★**数値の核（本ファイル）と測度（`MobiusRedDeg.lean`）は取れた。
残りは配管である。**
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★★★★★★★★★★★★★★★★★★★★★1 段ぶんの分離 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**原文 p.5 の `f0(β) ∉ S′`**。

原文 (NCBelyi p.5):
> bounded by some fixed expression in d0. Thus, for a suitable choice of C, it follows

`S` が `Gal(ℚ̄/ℚ)`-安定で `‖α‖ ≤ 1`（`∀ α ∈ S`）、`x₀ ∈ S`、`‖β‖ > 3` のとき、
`f₀ ≔ minpoly ℚ x₀` は

    `f₀(β) ≠ f₀(x)`   （`x ∈ S` または `f₀′(x) = 0` なるすべての `x`）

★すなわち **`f₀(β) ∉ f₀(S) ∪ f₀(S₀)` = `S′`** である。

★★`Gal`-安定性が効くのは**`f₀` の根がすべて `S` に入る**ため
（原文の『without loss of generality, that S is Gal-stable』はここへ効く）。
★★★`f₀′` の根が単位円板の中であることは Gauss–Lucas による。 -/
theorem separation_step (S : Finset ℂ) (hstab : IsConjStable S)
    (hnorm : ∀ x ∈ S, ‖x‖ ≤ 1)
    {x₀ : ℂ} (hx₀ : x₀ ∈ S) (hint : IsIntegral ℚ x₀)
    (hdeg : 0 < (minpoly ℚ x₀).natDegree)
    {β : ℂ} (hβ : 3 < ‖β‖) :
    ∀ x : ℂ, (x ∈ S ∨ aeval x (derivative (minpoly ℚ x₀)) = 0) →
      aeval β (minpoly ℚ x₀) ≠ aeval x (minpoly ℚ x₀) := by
  intro x hx
  set f : ℂ[X] := (minpoly ℚ x₀).map (algebraMap ℚ ℂ) with hf
  have hmonic : f.Monic := (minpoly.monic hint).map _
  have hfdeg : 0 < f.natDegree := by
    rw [hf, Polynomial.natDegree_map_eq_of_injective (algebraMap ℚ ℂ).injective]
    exact hdeg
  have hroots : ∀ z ∈ f.roots, ‖z‖ ≤ 1 := norm_roots_minpoly_le_one hstab hnorm hx₀ hint
  have hx' : ‖x‖ ≤ 1 ∨ (derivative f).eval x = 0 := by
    rcases hx with h | h
    · exact Or.inl (hnorm x h)
    · refine Or.inr ?_
      rw [hf, Polynomial.derivative_map, Polynomial.eval_map, ← Polynomial.aeval_def]
      exact h
  have h2 := eval_ne_of_root_or_root_derivative f hmonic hfdeg hroots hβ hx'
  rw [hf, Polynomial.eval_map, Polynomial.eval_map, ← Polynomial.aeval_def,
    ← Polynomial.aeval_def] at h2
  exact h2

/-- ★★★★★★**`S′` の言葉で書いた形** —— `f₀(β) ∉ f₀(S) ∪ E`。

★`E` は `f₀′` の根の像（原文の `f₀(S₀)`）であり、
`exists_descend_data` が返す `E` と同じものである。 -/
theorem separation_step_not_mem (S : Finset ℂ) (hstab : IsConjStable S)
    (hnorm : ∀ x ∈ S, ‖x‖ ≤ 1)
    {x₀ : ℂ} (hx₀ : x₀ ∈ S) (hint : IsIntegral ℚ x₀)
    (hdeg : 0 < (minpoly ℚ x₀).natDegree)
    {β : ℂ} (hβ : 3 < ‖β‖) :
    aeval β (minpoly ℚ x₀)
      ∉ (S.image (fun x => aeval x (minpoly ℚ x₀)))
        ∪ ((((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)).roots.toFinset).image
            (fun w => aeval w (minpoly ℚ x₀))) := by
  classical
  have hd0 : derivative (minpoly ℚ x₀) ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := minpoly ℚ x₀)).1 hc
    omega
  have hmapne : ((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)) ≠ 0 :=
    (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
  intro hmem
  rcases Finset.mem_union.1 hmem with hm | hm
  · obtain ⟨x, hxS, hx⟩ := Finset.mem_image.1 hm
    exact separation_step S hstab hnorm hx₀ hint hdeg hβ x (Or.inl hxS) hx.symm
  · obtain ⟨w, hw, hx⟩ := Finset.mem_image.1 hm
    rw [Multiset.mem_toFinset] at hw
    have hwr : aeval w (derivative (minpoly ℚ x₀)) = 0 := by
      have h := (Polynomial.mem_roots hmapne).1 hw
      simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using h
    exact separation_step S hstab hnorm hx₀ hint hdeg hβ w (Or.inr hwr) hx.symm

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def separation_step.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(1 段ぶんの分離——f₀(β) ≠ f₀(x))",
    sectionId := "ncbelyi-lemma-2-4" }

def separation_step_not_mem.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(f₀(β) ∉ S′——原文の for a suitable choice of C)",
    sectionId := "ncbelyi-lemma-2-4" }

def separation_step_not_mem.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "eval_ne_of_root_or_root_derivative(C = 3 で足りる、第 405-406)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.eval_ne_of_root_or_root_derivative") 3,
    .citation "[ABC3]" "norm_roots_minpoly_le_one(Gal-安定なら f₀ の根も単位円板の中、第 410)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.norm_roots_minpoly_le_one") 3,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 原文 p.5 の『Thus, for a suitable choice of C, " ++
       "it follows that f0(β) ∉ S′』は、部品 3 つ" ++
       "(Gal-安定性 ⟹ f₀ の根が単位円板／Gauss–Lucas ⟹ f₀′ の根も単位円板／" ++
       "‖β‖ > 3 ⟹ 値が分かれる)を繋ぐだけである。" ++
       "★C = 3 で足りることは第 405-406 で測ってある" ++
       "——原文の『d₀ に相対的に十分大きい C』は d₀ に依らない") 6,
    .implicitStep
      ("★★★★★これで Lemma 2.4 に残るのは**配管 2 つ**である: " ++
       "(1) ℚ-Möbius が Gal-安定性を保つこと(conjSet (σ x) = σ (conjSet x))" ++
       "——測度の側は MobiusRedDeg.lean で済んでいる。" ++
       "(2) 合成を ℙ¹ の有理写像として組み立てること" ++
       "(原文の f(x) ∈ ℚ(x) は Möbius と多項式の交互合成。Separation.lean の P1C が受け皿)。" ++
       "★★数値の核(本ファイル)と測度(MobiusRedDeg.lean)は取れた") 7 ]

end ABC3.Found.NCBelyi
