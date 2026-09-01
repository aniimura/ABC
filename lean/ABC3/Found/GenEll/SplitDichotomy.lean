/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.QuadTwist

/-!
# 第 963 ブロック —— **★★★★★★★★★★★★★★★★分裂するか、捧れば分裂するか**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——(D3) の (a)

mathlib の `HasSplitMultiplicativeReduction` は剰余体で 2 次式

    `c₄X² + a₁c₄X - (54b₆ - 3b₂b₄ + a₂c₄)`

が分裂することを要求する。★だが**分裂しなくても良い**——
非平方で捧れば分裂する（第 938）からである。

☆剰余体は有限体で、標数が `2` でなければ**非平方元が存在する**
（`FiniteField.exists_nonsquare`）。それを持ち上げればよい。

★これで `Lemma 3.5` の分裂/非分裂の場合分けが**型の上で閉じる**。
☆非分裂の側は `minDeltaExp_eq_mul_of_nonsplit`（第 929）が受ける。

| 定理 | 内容 |
|---|---|
| `splits_or_exists_twist_splits` | ★★★★★★★★★★★★★★★★**二者择一** |
-/

namespace ABC3.Found.GenEll

open Polynomial

/-- ★★★★★★★★★★★★★★★★**分裂するか、ある捧りで分裂するか**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 963）**——これが (D3) の (a) である。
☆場合分けは `∃ x, c₄x² = K` かどうかだけ:

* ある ⇒ `a₁ = 0`（`IsCharNeTwoNF`）なので 2 次式は `c₄X² - K`、
  根 `x` を持つので分裂（`splits_quadratic_of_root`）
* ない ⇒ 非平方元 `a` を取り（`FiniteField.exists_nonsquare`）、
  全射性で持ち上げて第 938 を当てる -/
theorem splits_or_exists_twist_splits {R k : Type} [CommRing R] [Field k] [Fintype k]
    [DecidableEq k]
    (φ : R →+* k) (hφ : Function.Surjective φ) (V : WeierstrassCurve R) [V.IsCharNeTwoNF]
    (hchar : ringChar k ≠ 2)
    (hc : φ V.c₄ ≠ 0)
    (hK : φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄) ≠ 0) :
    Splits (Polynomial.map φ
        (C V.c₄ * X ^ 2 + C (V.a₁ * V.c₄) * X
          - C (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)))
      ∨ ∃ d : R, φ d ≠ 0 ∧ Splits (Polynomial.map φ
        (C (quadTwist V d).c₄ * X ^ 2
          + C ((quadTwist V d).a₁ * (quadTwist V d).c₄) * X
          - C (54 * (quadTwist V d).b₆ - 3 * (quadTwist V d).b₂ * (quadTwist V d).b₄
              + (quadTwist V d).a₂ * (quadTwist V d).c₄))) := by
  by_cases hns : ∃ x : k, φ V.c₄ * x ^ 2 = φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)
  · left
    obtain ⟨x, hx⟩ := hns
    have ha1 : V.a₁ = 0 := ‹V.IsCharNeTwoNF›.a₁
    have hmap : Polynomial.map φ
        (C V.c₄ * X ^ 2 + C (V.a₁ * V.c₄) * X
          - C (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄))
        = C (φ V.c₄) * X ^ 2 + C (0 : k) * X
          - C (φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)) := by
      rw [ha1]
      simp
    rw [hmap]
    exact splits_quadratic_of_root _ _ _ hc x (by rw [← hx]; ring)
  · right
    obtain ⟨a, ha⟩ := FiniteField.exists_nonsquare (F := k) hchar
    obtain ⟨d, rfl⟩ := hφ a
    have hd0 : φ d ≠ 0 := by
      intro h
      exact ha (by rw [h]; exact ⟨0, by ring⟩)
    exact ⟨d, hd0, splits_quadTwist_of_not_isSquare φ V d hc hd0 hK hns ha⟩

/-! ## ★★★★★★★★★★★★★★★★第 979 —— 剰余標数 2 でも二者択一は成り立つ

★第 963 は `ringChar k ≠ 2` を受けていた——`FiniteField.exists_nonsquare` が
それを要求するからである。☆だが**標数 2 では非平方元が無い代わりに、
2 次式が必ず分裂する**:

有限体の標数 2 では `x ↦ x²` は全単射（`FiniteField.isSquare_of_char_two`）なので
`∃ x, c₄x² = K` が常に解ける。`IsCharNeTwoNF` により `a₁ = 0` だから
2 次式は `c₄X² − K` であり、根をもつので分裂する。

★★したがって**二者択一は標数の仮定なしに成り立つ**。 -/

/-- ★**標数 2 の有限体では `c·x² = K` が常に解ける**——`x ↦ x²` が全単射だから。 -/
theorem exists_sq_of_charTwo {k : Type} [Field k] [Finite k] (h2 : ringChar k = 2)
    (c K : k) (hc : c ≠ 0) : ∃ x : k, c * x ^ 2 = K := by
  obtain ⟨y, hy⟩ := FiniteField.isSquare_of_char_two h2 (K / c)
  refine ⟨y, ?_⟩
  have hy2 : y ^ 2 = K / c := by rw [hy]; ring
  rw [hy2]
  field_simp

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**標数の仮定なしの二者択一**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 979）**——第 963 から `ringChar k ≠ 2` が落ちた。
☆標数 2 では非平方元が無い代わりに 2 次式が必ず分裂するからである。 -/
theorem splits_or_exists_twist_splits' {R k : Type} [CommRing R] [Field k] [Fintype k]
    [DecidableEq k]
    (φ : R →+* k) (hφ : Function.Surjective φ) (V : WeierstrassCurve R) [V.IsCharNeTwoNF]
    (hc : φ V.c₄ ≠ 0)
    (hK : φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄) ≠ 0) :
    Splits (Polynomial.map φ
        (C V.c₄ * X ^ 2 + C (V.a₁ * V.c₄) * X
          - C (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)))
      ∨ ∃ d : R, φ d ≠ 0 ∧ Splits (Polynomial.map φ
        (C (quadTwist V d).c₄ * X ^ 2
          + C ((quadTwist V d).a₁ * (quadTwist V d).c₄) * X
          - C (54 * (quadTwist V d).b₆ - 3 * (quadTwist V d).b₂ * (quadTwist V d).b₄
              + (quadTwist V d).a₂ * (quadTwist V d).c₄))) := by
  by_cases hchar : ringChar k = 2
  · left
    obtain ⟨x, hx⟩ := exists_sq_of_charTwo hchar (φ V.c₄)
      (φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)) hc
    have ha1 : V.a₁ = 0 := ‹V.IsCharNeTwoNF›.a₁
    have hmap : Polynomial.map φ
        (C V.c₄ * X ^ 2 + C (V.a₁ * V.c₄) * X
          - C (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄))
        = C (φ V.c₄) * X ^ 2 + C (0 : k) * X
          - C (φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)) := by
      rw [ha1]; simp
    rw [hmap]
    exact splits_quadratic_of_root _ _ _ hc x (by rw [← hx]; ring)
  · exact splits_or_exists_twist_splits φ hφ V hchar hc hK

def exists_sq_of_charTwo.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(標数 2 の有限体では c·x² = K が常に解ける。★無条件)",
    sectionId := "genell-lemma-3-5" }

def splits_or_exists_twist_splits'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(標数の仮定なしの二者択一。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★出典の紐付け(`.src`) -/

def splits_or_exists_twist_splits.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(分裂するか、ある捧りで分裂するか。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
