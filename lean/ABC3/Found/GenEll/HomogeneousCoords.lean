/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottPacked
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★点の同次座標を**構成する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★これは何か —— 座標は仮定でなく構成物である

`§9-940` まで、`Proposition 1.4, (iv)` は点ごとの**同次座標** `x : Fin (N+1) → 𝓞_F` を
仮定として受けていた。★★本ファイルはそれを**構成する**。

## ★★★機構 —— 生成点で読んで分母を払う

1. `Spec F` は**1 点**（`F` は体）なので、`ℙᴺ` のどれかのチャート `D₊(X_i)` に必ず入る
   （`§9-C2b` の `exists_chart_range`）
2. そのチャートで正規化座標 `c_k ≔ (x_k/x_i)(点) ∈ F` が読める（`projPointCoord`）
   ——`c_i = 1` である
3. ★`𝓞_F` は `F` の中で分母を払える（`IsLocalization.exist_integer_multiples_of_finite`）
   ——`b ∈ 𝓞_F∖{0}` を取って `x_k ≔ b·c_k ∈ 𝓞_F` と置く
4. すると `x_i = b ≠ 0` かつ **`x_k = c_k·x_i`**

★★★★これが原文の「同次座標」であり、`§9-940` の `hw`・`hcv` の**生成点版**である。

## ★これで何が残ったか

★★★残るのは**各素点での読み替え**である:

* 有限素点 `Q`——`v_Q(x_j)` が最小の `j` を取れば `𝔞_Q = (x_j)` かつ `x_0/x_j ∈ (𝓞_F)_Q`
* 無限素点 `v`——生成点の関係を `σ_v` で送るだけ（★`x_i ≠ 0` なのでチャートは同じ `i` でよい）

★どちらも本ファイルの `x_k = c_k·x_i` から出る形である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★★★★★★★★★★★★★★同次座標の構成 -/

/-- ★★★★★★★★★★★★★★★**`ℙᴺ` の `F`-点は `𝓞_F` の同次座標を持つ**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は 3 段:
* `Spec F` は 1 点なのでどれかのチャートに入る（`exists_chart_range`）
* そこで正規化座標 `c_k` が読める（`c_i = 1`）
* ★`𝓞_F` は `F` の中で分母を払える——`x_k ≔ b·c_k` と置けば `x_i = b ≠ 0`

★★これが原文の「同次座標」であり、`§9-940` の `hw`・`hcv` の**生成点版**である。 -/
theorem exists_homogeneous_coords (N : ℕ) (F : Type) [Field F] [NumberField F]
    (q : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    ∃ (i : Fin (N + 1)) (hx : Set.range q.base ⊆ Set.range (chartA N ℤ i).base)
      (x : Fin (N + 1) → 𝓞 F), x i ≠ 0 ∧ x ≠ 0 ∧
        ∀ k, ((x k : F)) = projPointCoord N ℤ F q i hx k * ((x i : F)) := by
  obtain ⟨i, hxr⟩ := exists_chart_range N ℤ F q
  refine ⟨i, hxr, ?_⟩
  set c : Fin (N + 1) → F := fun k => projPointCoord N ℤ F q i hxr k with hc
  obtain ⟨b, hb⟩ := IsLocalization.exist_integer_multiples_of_finite
    (nonZeroDivisors (𝓞 F)) (S := F) c
  choose x hx using fun k => hb k
  have hci : c i = 1 := projPointCoord_self N ℤ F q i hxr
  have hxi : ((x i : F)) = ((b : 𝓞 F) : F) := by
    have h := hx i
    rw [hci, Algebra.smul_def, mul_one] at h
    simpa using h
  have hbne : ((b : 𝓞 F) : F) ≠ 0 := by
    simpa using (nonZeroDivisors.coe_ne_zero b)
  have hine : x i ≠ 0 := by
    intro h0
    apply hbne
    rw [← hxi, h0]
    simp
  refine ⟨x, hine, ?_, fun k => ?_⟩
  · intro h0
    exact hine (by rw [h0]; rfl)
  · have h := hx k
    rw [Algebra.smul_def] at h
    have h2 : ((x k : F)) = ((b : 𝓞 F) : F) * c k := by simpa using h
    rw [h2, hxi]
    ring

/-! ## ★★★★★★★★★★★★スキームの点に付ける -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★**`X` の `𝓞_F`-点は `ψ` を通して同次座標を持つ**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`exists_homogeneous_coords` を**生成点** `Spec F ⟶ Spec 𝓞_F ⟶ X ⟶ ℙᴺ` に当てるだけである。
★★`§9-940` の点ごとの座標 `x p` は、これで**構成物になった**。 -/
theorem exists_homogeneous_coords_of_point {X : Scheme.{0}}
    (N : ℕ) (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X)
    (ψ : X ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :
    ∃ (i : Fin (N + 1))
      (hx : Set.range (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ ψ).base
        ⊆ Set.range (chartA N ℤ i).base)
      (x : Fin (N + 1) → 𝓞 F), x i ≠ 0 ∧ x ≠ 0 ∧
        ∀ k, ((x k : F))
          = projPointCoord N ℤ F
              (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ ψ) i hx k
            * ((x i : F)) :=
  exists_homogeneous_coords N F _

/-! ## ★出典の紐付け(`.src`) -/

def exists_homogeneous_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ℙᴺ の F-点は 𝓞_F の同次座標を持つ)",
    sectionId := "genell-prop-1-4" }

def exists_homogeneous_coords_of_point.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X の 𝓞_F-点は ψ を通して同次座標を持つ)",
    sectionId := "genell-prop-1-4" }

def exists_homogeneous_coords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_chart_range(F-点はどれかのチャートに入る、§9-C2b)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_chart_range") 2,
    .citation "[mathlib]" "IsLocalization.exist_integer_multiples_of_finite(分母を払う)"
      (.inMathlib "IsLocalization.exist_integer_multiples_of_finite") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-940 まで**仮定**だった点ごとの同次座標 x は" ++
       "**構成物である**。Spec F は 1 点なのでどれかのチャートに入り、" ++
       "そこで正規化座標 c_k(c_i = 1)が読め、𝓞_F は F の中で分母を払えるから、" ++
       "x_k ≔ b·c_k と置けば x_i = b ≠ 0 かつ x_k = c_k·x_i である") 5,
    .implicitStep
      ("★残るのは各素点での読み替えである——" ++
       "有限素点 Q では v_Q(x_j) が最小の j を取れば 𝔞_Q = (x_j) かつ x_0/x_j ∈ (𝓞_F)_Q、" ++
       "無限素点 v では生成点の関係を σ_v で送るだけ" ++
       "(x_i ≠ 0 なのでチャートは同じ i でよい)") 4 ]

end ABC3.Found.GenEll
