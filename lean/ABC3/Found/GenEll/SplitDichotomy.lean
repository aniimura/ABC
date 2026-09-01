/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.QuadTwist
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction

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

/-! ## ★★★★★★★★★★★★第 980 —— 整モデルを `IsCharNeTwoNF` に正規化する

★第 979 の二者択一は `[V.IsCharNeTwoNF]`（`a₁ = a₃ = 0`）を受ける。
☆`2` が可逆な環では変数変換 1 回でその形にできる
（mathlib の `exists_variableChange_isCharNeTwoNF`、`[Invertible 2]` を要求）。

★完備化の整数環で `2` が可逆になるのは `v_p(2) = 0`、すなわち `p ∤ 2` のときである
（第 953 の `isUnit_natCast_of_valAdd_eq_zero` に `n = 2` を入れる）。
☆`p ∣ 2` では別の手当て（不分岐 2 次拡大）が要る——第 979 の測定を見よ。 -/

/-- ★★★★★★★★★★★★**`2` が単元なら `IsCharNeTwoNF` に正規化できる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆mathlib の `exists_variableChange_isCharNeTwoNF` は `[Invertible 2]` を受けるので、
`IsUnit (2 : R)` から `Invertible` を作って渡すだけである。 -/
theorem exists_variableChange_isCharNeTwoNF_of_isUnit_two {R : Type} [CommRing R]
    (W : WeierstrassCurve R) (h2 : IsUnit (2 : R)) :
    ∃ C : WeierstrassCurve.VariableChange R, (C • W).IsCharNeTwoNF := by
  haveI : Invertible (2 : R) := h2.invertible
  exact W.exists_variableChange_isCharNeTwoNF

def exists_variableChange_isCharNeTwoNF_of_isUnit_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2 が単元なら IsCharNeTwoNF に正規化できる。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★第 981 —— 整モデルの `IsCharNeTwoNF` は体の側で決まる

★第 979 の二者択一は `[V.IsCharNeTwoNF]` を **整モデル** `V = integralModel R W` について受ける。
☆第 980 は `R` の中で正規化する道（`2` が `R` で単元、すなわち `p ∤ 2`）だった。

★★だが `integralModel_a₁_eq`（mathlib）は
`algebraMap R K (integralModel R W).a₁ = W.a₁` を与える。
`algebraMap R K` は単射だから、**`W.a₁ = 0` なら整モデルの `a₁` も `0`** である。

☆つまり正規化は**体 `K`（あるいは `L`）の側でやればよい**——
標数 0 なら `2` は必ず可逆だから、剰余標数の条件は要らない。
★残る問題は「体の側で正規化した曲線が `p` で整のままか」であり、
`p ∣ 2` ではそこが崩れうる（第 979・980 の測定を見よ）。 -/

/-- ★★★★★★★★★★★★**体の側が `a₁ = a₃ = 0` なら整モデルもそうである**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`integralModel_a₁_eq`・`integralModel_a₃_eq`（mathlib）と
`algebraMap R K` の単射性だけである。★剰余標数の条件は要らない。 -/
theorem isCharNeTwoNF_integralModel {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [WeierstrassCurve.IsIntegral R W]
    (ha1 : W.a₁ = 0) (ha3 : W.a₃ = 0) :
    (WeierstrassCurve.integralModel R W).IsCharNeTwoNF := by
  refine ⟨?_, ?_⟩
  · refine IsFractionRing.injective R K ?_
    rw [WeierstrassCurve.integralModel_a₁_eq R W, ha1, map_zero]
  · refine IsFractionRing.injective R K ?_
    rw [WeierstrassCurve.integralModel_a₃_eq R W, ha3, map_zero]

def isCharNeTwoNF_integralModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(体の側が a₁ = a₃ = 0 なら整モデルもそう。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★第 982 —— 仮説は `c₄` が単元であることだけ

★第 979 はまだ `hK : φ(54b₆ − 3b₂b₄ + a₂c₄) ≠ 0` を受けていた。
☆だが **`φK = 0` なら 2 次式は `c₄X²`（`a₁ = 0` だから）で、`X = 0` を根にもつ**——
やはり分裂する。

★★したがって二者択一に要るのは **`φ c₄ ≠ 0`（＝乗法還元）だけ**である。
☆これは `HasMultiplicativeReduction` の `multiplicativeReduction` フィールドそのものなので、
第 976 が与える。 -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**二者択一に要るのは `c₄` が単元であることだけ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 982）**——第 963 から `ringChar k ≠ 2`（第 979）が落ち、
さらに `hK`（本ブロック）も落ちた。
☆残る仮説は `φ c₄ ≠ 0` と `IsCharNeTwoNF` の 2 つだけである。 -/
theorem splits_or_exists_twist_splits'' {R k : Type} [CommRing R] [Field k] [Fintype k]
    [DecidableEq k]
    (φ : R →+* k) (hφ : Function.Surjective φ) (V : WeierstrassCurve R) [V.IsCharNeTwoNF]
    (hc : φ V.c₄ ≠ 0) :
    Splits (Polynomial.map φ
        (C V.c₄ * X ^ 2 + C (V.a₁ * V.c₄) * X
          - C (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)))
      ∨ ∃ d : R, φ d ≠ 0 ∧ Splits (Polynomial.map φ
        (C (quadTwist V d).c₄ * X ^ 2
          + C ((quadTwist V d).a₁ * (quadTwist V d).c₄) * X
          - C (54 * (quadTwist V d).b₆ - 3 * (quadTwist V d).b₂ * (quadTwist V d).b₄
              + (quadTwist V d).a₂ * (quadTwist V d).c₄))) := by
  by_cases hK : φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄) = 0
  · left
    have ha1 : V.a₁ = 0 := ‹V.IsCharNeTwoNF›.a₁
    have hmap : Polynomial.map φ
        (C V.c₄ * X ^ 2 + C (V.a₁ * V.c₄) * X
          - C (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄))
        = C (φ V.c₄) * X ^ 2 + C (0 : k) * X
          - C (φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)) := by
      rw [ha1]; simp
    rw [hmap]
    exact splits_quadratic_of_root _ _ _ hc 0 (by rw [hK]; ring)
  · exact splits_or_exists_twist_splits' φ hφ V hc hK

def splits_or_exists_twist_splits''.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(二者択一に要るのは c₄ が単元であることだけ。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★出典の紐付け(`.src`) -/

def splits_or_exists_twist_splits.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(分裂するか、ある捧りで分裂するか。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
