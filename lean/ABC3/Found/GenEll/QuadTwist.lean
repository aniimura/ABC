/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.JScale
import Mathlib.AlgebraicGeometry.EllipticCurve.NormalForms
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction
import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic

/-!
# 第 920 ブロック —— **★★★★★★★★★★★★二次の捧り（quadratic twist）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★これは何か

`Lemma 3.5` に残っているのは「`j` が同じ**分裂**乗法還元の対を 1 つ作る」
ことである（第 919 `minDeltaExp_eq_mul_of_twist`）。
★それを与えるのが**二次の捧り**である。

`a₁ = a₃ = 0` の形（mathlib の `IsCharNeTwoNF`、`2` が可逆ならいつでも取れる）で

    `W^d : a₂ ↦ d·a₂,  a₄ ↦ d²·a₄,  a₆ ↦ d³·a₆`

とすると、`c₄ ↦ d²c₄`・`c₆ ↦ d³c₆`・`Δ ↦ d⁶Δ` となり、
★★**`j` は変わらない**。

☆`j = Δ⁻¹c₄³`（第 913）に直しておけば、`d⁻⁶·d⁶ = 1` の計算だけである。

| 定義・定理 | 内容 |
|---|---|
| `quadTwist` | ★二次の捧り |
| `quadTwist_c₄` / `_c₆` / `_Δ` | 重み `d²`・`d³`・`d⁶` |
| `quadTwist_j` | ★★★★**`j` は変わらない** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve

variable {F : Type} [Field F]

/-- ★**二次の捧り**（`a₁ = a₃ = 0` の形で）。 -/
def quadTwist {A : Type} [CommRing A] (W : WeierstrassCurve A) (d : A) :
    WeierstrassCurve A where
  a₁ := 0
  a₂ := d * W.a₂
  a₃ := 0
  a₄ := d ^ 2 * W.a₄
  a₆ := d ^ 3 * W.a₆

instance quadTwist_isCharNeTwoNF {A : Type} [CommRing A] (W : WeierstrassCurve A) (d : A) :
    (quadTwist W d).IsCharNeTwoNF := ⟨rfl, rfl⟩

@[simp] theorem quadTwist_a₂ {A : Type} [CommRing A] (W : WeierstrassCurve A) (d : A) :
    (quadTwist W d).a₂ = d * W.a₂ := rfl

@[simp] theorem quadTwist_a₄ {A : Type} [CommRing A] (W : WeierstrassCurve A) (d : A) :
    (quadTwist W d).a₄ = d ^ 2 * W.a₄ := rfl

@[simp] theorem quadTwist_a₆ {A : Type} [CommRing A] (W : WeierstrassCurve A) (d : A) :
    (quadTwist W d).a₆ = d ^ 3 * W.a₆ := rfl

/-- ★★`c₄` は `d²` 倍。 -/
theorem quadTwist_c₄ {A : Type} [CommRing A] (W : WeierstrassCurve A) [W.IsCharNeTwoNF]
    (d : A) : (quadTwist W d).c₄ = d ^ 2 * W.c₄ := by
  rw [WeierstrassCurve.c₄_of_isCharNeTwoNF, WeierstrassCurve.c₄_of_isCharNeTwoNF,
    quadTwist_a₂, quadTwist_a₄]
  ring

/-- ★★`c₆` は `d³` 倍。 -/
theorem quadTwist_c₆ {A : Type} [CommRing A] (W : WeierstrassCurve A) [W.IsCharNeTwoNF]
    (d : A) : (quadTwist W d).c₆ = d ^ 3 * W.c₆ := by
  rw [WeierstrassCurve.c₆_of_isCharNeTwoNF, WeierstrassCurve.c₆_of_isCharNeTwoNF,
    quadTwist_a₂, quadTwist_a₄, quadTwist_a₆]
  ring

/-- ★★`Δ` は `d⁶` 倍。 -/
theorem quadTwist_Δ {A : Type} [CommRing A] (W : WeierstrassCurve A) [W.IsCharNeTwoNF]
    (d : A) : (quadTwist W d).Δ = d ^ 6 * W.Δ := by
  rw [WeierstrassCurve.Δ_of_isCharNeTwoNF, WeierstrassCurve.Δ_of_isCharNeTwoNF,
    quadTwist_a₂, quadTwist_a₄, quadTwist_a₆]
  ring

/-- ★★★捧りも楼円曲線である（`d ≠ 0` なら）。 -/
theorem quadTwist_isElliptic (W : WeierstrassCurve F) [W.IsCharNeTwoNF] [W.IsElliptic]
    {d : F} (hd : d ≠ 0) : (quadTwist W d).IsElliptic := by
  refine ⟨?_⟩
  rw [quadTwist_Δ]
  refine (isUnit_iff_ne_zero).2 ?_
  exact mul_ne_zero (pow_ne_zero 6 hd) W.isUnit_Δ.ne_zero

/-- ★★★★**捧りで `j` は変わらない**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`j = Δ⁻¹c₄³`（第 913）なので `(d⁶Δ)⁻¹(d²c₄)³ = Δ⁻¹c₄³` である。 -/
theorem quadTwist_j (W : WeierstrassCurve F) [W.IsCharNeTwoNF] [W.IsElliptic]
    {d : F} (hd : d ≠ 0) [(quadTwist W d).IsElliptic] :
    (quadTwist W d).j = W.j := by
  rw [j_eq_inv_Delta_mul, j_eq_inv_Delta_mul, quadTwist_Δ, quadTwist_c₄]
  have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  field_simp

/-! ## ★★★★★★★★分裂性は「剰余体に根がある」だけ -/

open Polynomial in
/-- ★★**2 次式は根があれば分裂する**——mathlib の `Splits.of_degree_eq_two` を
`HasSplitMultiplicativeReduction` のフィールの形に合わせたもの。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これで「分裂乗法還元」の残る中身は
**剰余体の中で 2 次式の根を 1 つ見つける**ことだけになる。 -/
theorem splits_quadratic_of_root {k : Type} [Field k] (A B C₀ : k) (hA : A ≠ 0)
    (x : k) (hx : A * x ^ 2 + B * x - C₀ = 0) :
    Polynomial.Splits (Polynomial.C A * Polynomial.X ^ 2
      + Polynomial.C B * Polynomial.X - Polynomial.C C₀) := by
  refine Polynomial.Splits.of_degree_eq_two (x := x) ?_ ?_
  · have : (Polynomial.C A * Polynomial.X ^ 2 + Polynomial.C B * Polynomial.X
        - Polynomial.C C₀)
        = Polynomial.C A * Polynomial.X ^ 2 + Polynomial.C B * Polynomial.X
          + Polynomial.C (-C₀) := by
      rw [map_neg]
      ring
    rw [this]
    exact Polynomial.degree_quadratic hA
  · simp only [Polynomial.eval_sub, Polynomial.eval_add, Polynomial.eval_mul,
      Polynomial.eval_pow, Polynomial.eval_C, Polynomial.eval_X]
    exact hx

/-! ## ★★★★★★★★有限体では非平方 × 非平方 = 平方 -/

/-- ★★★★★★**有限体（奇標数）では非平方同士の積は平方**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが「`d` を非平方単数に取れば捧りは分裂になる」の理由である。
☆二次指標は乗法的なので `(−1)·(−1) = 1`。 -/
theorem isSquare_mul_of_not_isSquare {k : Type} [Field k] [Fintype k] [DecidableEq k]
    {a d : k} (ha : a ≠ 0) (hd : d ≠ 0)
    (hna : ¬ IsSquare a) (hnd : ¬ IsSquare d) : IsSquare (d * a) := by
  have h1 : quadraticChar k a = -1 := quadraticChar_neg_one_iff_not_isSquare.2 hna
  have h2 : quadraticChar k d = -1 := quadraticChar_neg_one_iff_not_isSquare.2 hnd
  have h3 : quadraticChar k (d * a) = 1 := by
    rw [map_mul, h1, h2]
    norm_num
  exact (quadraticChar_one_iff_isSquare (mul_ne_zero hd ha)).1 h3

/-! ## ★★★★★★★★★★Vélu の商は捧りと可換 -/

/-- ★★★★★★★★★★**Vélu の商は捧りと可換する**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `veluCurve (W^d) (d²v) (d³w) = (veluCurve W v w)^d`

☆`v` は重み 4、`w` は重み 6 なので、捧り（`c₄ ↦ d²c₄`）では
`v ↦ d²v`・`w ↦ d³w` となる。
★`b₂(W^d) = d·b₂(W)`（`a₁ = a₃ = 0` の形）が鍵である。

★これが「`E′^d = (E^d)/H^d`」の曲線の水準での中身であり、
非分裂の降下に要る最後の道具である。 -/
theorem veluCurve_quadTwist {A : Type} [CommRing A] (W : WeierstrassCurve A)
    [W.IsCharNeTwoNF] (v w d : A) :
    veluCurve (quadTwist W d) (d ^ 2 * v) (d ^ 3 * w)
      = quadTwist (veluCurve W v w) d := by
  have hb₂ : (quadTwist W d).b₂ = d * W.b₂ := by
    rw [WeierstrassCurve.b₂_of_isCharNeTwoNF, WeierstrassCurve.b₂_of_isCharNeTwoNF,
      quadTwist_a₂]
    ring
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · rfl
  · rfl
  · rfl
  · show (quadTwist W d).a₄ - 5 * (d ^ 2 * v) = d ^ 2 * (W.a₄ - 5 * v)
    rw [quadTwist_a₄]
    ring
  · show (quadTwist W d).a₆ - (quadTwist W d).b₂ * (d ^ 2 * v) - 7 * (d ^ 3 * w)
      = d ^ 3 * (W.a₆ - W.b₂ * v - 7 * w)
    rw [quadTwist_a₆, hb₂]
    ring

/-! ## ★★捧りは底変換と可換 -/

/-- ★★**捧りは底変換と可換する**。

☆定義が係数の多項式なので `ext` と `map_mul`・`map_pow` だけである。 -/
theorem quadTwist_map {A B : Type} [CommRing A] [CommRing B] (f : A →+* B)
    (W : WeierstrassCurve A) (d : A) :
    (quadTwist W d).map f = quadTwist (W.map f) (f d) := by
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · show f (0 : A) = 0
    exact map_zero f
  · show f (d * W.a₂) = f d * f W.a₂
    exact map_mul f d W.a₂
  · show f (0 : A) = 0
    exact map_zero f
  · show f (d ^ 2 * W.a₄) = f d ^ 2 * f W.a₄
    rw [map_mul, map_pow]
  · show f (d ^ 3 * W.a₆) = f d ^ 3 * f W.a₆
    rw [map_mul, map_pow]

/-! ## ★★★★★★★★★★★★捧った Vélu の商 -/

/-- ★★★★★★★★★★★★**`(E/H)^d = veluCurve (E^d) (d²v) (d³w)`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**これが非分裂の降下の鍵である**——
捧った商は、捧った曲線の `veluCurve` として**点集合を介さずに**書ける。
☆したがって `√d` を体に追加する必要がない——
`tateParam_quot_mu`（第 883）が消費するのは `veluCurve` の形だからである。 -/
theorem quadTwist_veluQuotientFull (W : WeierstrassCurve F) [W.IsCharNeTwoNF]
    (S : Finset (F × F)) (d : F) :
    quadTwist (veluQuotientFull W S) d
      = veluCurve (quadTwist W d) (d ^ 2 * veluVFull W S) (d ^ 3 * veluWFull W S) := by
  rw [veluQuotientFull, veluCurve_quadTwist]

/-- ★★★★捧った商を `veluCurve` の形で取る（仮説を `hE′` で受けた形）。 -/
theorem quadTwist_eq_veluCurve (W E' : WeierstrassCurve F) [W.IsCharNeTwoNF]
    (S : Finset (F × F)) (d : F) (hE' : E' = veluQuotientFull W S) :
    quadTwist E' d
      = veluCurve (quadTwist W d) (d ^ 2 * veluVFull W S) (d ^ 3 * veluWFull W S) := by
  rw [hE', quadTwist_veluQuotientFull]

/-! ## ★★★★★★★★整モデルは一意である -/

section IntegralModel

variable {R : Type} [CommRing R] {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★**整モデルは一意である**——`algebraMap R K` が単射なら。

☆mathlib の `integralModel` は `.choose` で定義されているが、
底変換が一致する `R` 上の曲線は係数ごとに一意に決まる。
★これで「`integralModel R (W^d)` は何か」を手で指定できる。 -/
theorem integralModel_eq_of_map_eq (hinj : Function.Injective (algebraMap R K))
    (W : WeierstrassCurve K) [WeierstrassCurve.IsIntegral R W] (V : WeierstrassCurve R)
    (h : V.map (algebraMap R K) = W) :
    WeierstrassCurve.integralModel R W = V := by
  have h1 : (WeierstrassCurve.integralModel R W).map (algebraMap R K) = W :=
    WeierstrassCurve.baseChange_integralModel_eq R W
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_ <;>
    refine hinj ?_
  · show ((WeierstrassCurve.integralModel R W).map (algebraMap R K)).a₁
      = (V.map (algebraMap R K)).a₁
    rw [h1, h]
  · show ((WeierstrassCurve.integralModel R W).map (algebraMap R K)).a₂
      = (V.map (algebraMap R K)).a₂
    rw [h1, h]
  · show ((WeierstrassCurve.integralModel R W).map (algebraMap R K)).a₃
      = (V.map (algebraMap R K)).a₃
    rw [h1, h]
  · show ((WeierstrassCurve.integralModel R W).map (algebraMap R K)).a₄
      = (V.map (algebraMap R K)).a₄
    rw [h1, h]
  · show ((WeierstrassCurve.integralModel R W).map (algebraMap R K)).a₆
      = (V.map (algebraMap R K)).a₆
    rw [h1, h]

/-- ★★★★★★★★**捧りの整モデルは整モデルの捧り**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが `HasSplitMultiplicativeReduction` の `Splits` 条件を
捧りの側へ送るための配管である。 -/
theorem integralModel_quadTwist (hinj : Function.Injective (algebraMap R K))
    (W : WeierstrassCurve K) [WeierstrassCurve.IsIntegral R W] (d : R)
    [WeierstrassCurve.IsIntegral R (quadTwist W (algebraMap R K d))] :
    WeierstrassCurve.integralModel R (quadTwist W (algebraMap R K d))
      = quadTwist (WeierstrassCurve.integralModel R W) d := by
  refine integralModel_eq_of_map_eq hinj _ _ ?_
  rw [quadTwist_map]
  congr 1
  exact WeierstrassCurve.baseChange_integralModel_eq R W

end IntegralModel

/-- ★★★★★★★★**非分裂の 2 次式は非平方で捧ると根をもつ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★捧りで 2 次式の係数は `c ↦ d²c`・`K ↦ d³K` となるので、
判定量 `K/c` は `d` 倍される。有限体では非平方×非平方 = 平方（第 922）なので、
`d` を非平方に取れば捧った側は根をもつ。 -/
theorem exists_sq_of_twist_of_not_isSquare {k : Type} [Field k] [Fintype k] [DecidableEq k]
    (c d K : k) (hc : c ≠ 0) (hd0 : d ≠ 0) (hK : K ≠ 0)
    (hnsK : ¬ ∃ x : k, c * x ^ 2 = K) (hnd : ¬ IsSquare d) :
    ∃ x : k, (d ^ 2 * c) * x ^ 2 = d ^ 3 * K := by
  have hane : K * c⁻¹ ≠ 0 := mul_ne_zero hK (inv_ne_zero hc)
  have hna : ¬ IsSquare (K * c⁻¹) := by
    rintro ⟨y, hy⟩
    refine hnsK ⟨y, ?_⟩
    have : K * c⁻¹ = y ^ 2 := by rw [hy]; ring
    field_simp at this
    linear_combination -this
  obtain ⟨x, hx⟩ := isSquare_mul_of_not_isSquare hane hd0 hna hnd
  refine ⟨x, ?_⟩
  have hx2 : d * (K * c⁻¹) = x ^ 2 := by rw [hx]; ring
  have : (d ^ 2 * c) * x ^ 2 = d ^ 2 * c * (d * (K * c⁻¹)) := by rw [hx2]
  rw [this]
  field_simp

/-! ## ★★★★★★`Splits` の 2 次式の係数は `d²`・`d³` 倍 -/

/-- ★★★★★★**分裂性の 2 次式の定数項は `d³` 倍**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆mathlib の `HasSplitMultiplicativeReduction` のフィールに現れる
`54b₆ − 3b₂b₄ + a₂c₄` が、捧りでちょうど `d³` 倍になる。
★`c₄` は `d²` 倍（第 920）なので、判定量は `d` 倍される——第 936 の形である。 -/
theorem quadTwist_splitConst {A : Type} [CommRing A] (V : WeierstrassCurve A)
    [V.IsCharNeTwoNF] (d : A) :
    54 * (quadTwist V d).b₆ - 3 * (quadTwist V d).b₂ * (quadTwist V d).b₄
        + (quadTwist V d).a₂ * (quadTwist V d).c₄
      = d ^ 3 * (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄) := by
  rw [WeierstrassCurve.b₆_of_isCharNeTwoNF, WeierstrassCurve.b₆_of_isCharNeTwoNF,
    WeierstrassCurve.b₂_of_isCharNeTwoNF, WeierstrassCurve.b₂_of_isCharNeTwoNF,
    WeierstrassCurve.b₄_of_isCharNeTwoNF, WeierstrassCurve.b₄_of_isCharNeTwoNF,
    WeierstrassCurve.c₄_of_isCharNeTwoNF, WeierstrassCurve.c₄_of_isCharNeTwoNF,
    quadTwist_a₂, quadTwist_a₄, quadTwist_a₆]
  ring

/-- ★★**分裂性の 2 次式の `X` の係数は消える**（`a₁ = 0` だから）。 -/
@[simp] theorem quadTwist_splitLin {A : Type} [CommRing A] (V : WeierstrassCurve A) (d : A) :
    (quadTwist V d).a₁ * (quadTwist V d).c₄ = 0 := by
  show (0 : A) * _ = 0
  ring

/-! ## ★★★★★★★★★★★★非分裂の曲線は非平方で捧ると分裂になる -/

section SplitTwist

open Polynomial

variable {R : Type} [CommRing R] {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★★★★★**非分裂の曲線を非平方単数で捧ると分裂になる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが `Lemma 3.5` の非分裂の降下の残る中身である。
`HasSplitMultiplicativeReduction` の `Splits` フィールそのものを、
捧った整モデルについて出す。

☆部品はすべて揃っている:
`quadTwist_c₄`（920）・`splits_quadratic_of_root`（921）・
`exists_sq_of_twist_of_not_isSquare`（936）・`quadTwist_splitConst`（937）。 -/
theorem splits_quadTwist_of_not_isSquare {k : Type} [Field k] [Fintype k] [DecidableEq k]
    (φ : R →+* k) (V : WeierstrassCurve R) [V.IsCharNeTwoNF] (d : R)
    (hc : φ V.c₄ ≠ 0) (hd0 : φ d ≠ 0)
    (hK : φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄) ≠ 0)
    (hns : ¬ ∃ x : k, φ V.c₄ * x ^ 2 = φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄))
    (hnd : ¬ IsSquare (φ d)) :
    Polynomial.Splits (Polynomial.map φ
      (Polynomial.C (quadTwist V d).c₄ * Polynomial.X ^ 2
        + Polynomial.C ((quadTwist V d).a₁ * (quadTwist V d).c₄) * Polynomial.X
        - Polynomial.C (54 * (quadTwist V d).b₆
            - 3 * (quadTwist V d).b₂ * (quadTwist V d).b₄
            + (quadTwist V d).a₂ * (quadTwist V d).c₄))) := by
  obtain ⟨x, hx⟩ := exists_sq_of_twist_of_not_isSquare (φ V.c₄) (φ d)
    (φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)) hc hd0 hK hns hnd
  have hmap : Polynomial.map φ
      (Polynomial.C (quadTwist V d).c₄ * Polynomial.X ^ 2
        + Polynomial.C ((quadTwist V d).a₁ * (quadTwist V d).c₄) * Polynomial.X
        - Polynomial.C (54 * (quadTwist V d).b₆
            - 3 * (quadTwist V d).b₂ * (quadTwist V d).b₄
            + (quadTwist V d).a₂ * (quadTwist V d).c₄))
      = Polynomial.C (φ d ^ 2 * φ V.c₄) * Polynomial.X ^ 2
        + Polynomial.C (0 : k) * Polynomial.X
        - Polynomial.C (φ d ^ 3
            * φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄)) := by
    rw [quadTwist_splitConst, quadTwist_splitLin, quadTwist_c₄]
    simp only [Polynomial.map_sub, Polynomial.map_add, Polynomial.map_mul,
      Polynomial.map_pow, Polynomial.map_C, Polynomial.map_X, Polynomial.map_zero,
      map_mul, map_pow, map_zero]
  rw [hmap]
  refine splits_quadratic_of_root _ _ _ ?_ x ?_
  · exact mul_ne_zero (pow_ne_zero 2 hd0) hc
  · have : φ d ^ 2 * φ V.c₄ * x ^ 2
        = φ d ^ 3 * φ (54 * V.b₆ - 3 * V.b₂ * V.b₄ + V.a₂ * V.c₄) := hx
    linear_combination this

end SplitTwist

/-! ## ★出典の紐付け(`.src`) -/

def integralModel_eq_of_map_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(整モデルは一意である。★無条件)",
    sectionId := "genell-lemma-3-5" }

def integralModel_quadTwist.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(捧りの整モデルは整モデルの捧り。★無条件)",
    sectionId := "genell-lemma-3-5" }

def quadTwist_veluQuotientFull.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5((E/H)^d = veluCurve (E^d) (d²v) (d³w)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def quadTwist_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(捧りは底変換と可換する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluCurve_quadTwist.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は捧りと可換する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def splits_quadTwist_of_not_isSquare.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(非分裂の曲線を非平方で捧ると分裂になる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def quadTwist_splitConst.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(分裂性の 2 次式の定数項は d³ 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_sq_of_twist_of_not_isSquare.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(非分裂の 2 次式は非平方で捧ると根をもつ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isSquare_mul_of_not_isSquare.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有限体では非平方同士の積は平方。★無条件)",
    sectionId := "genell-lemma-3-5" }

def splits_quadratic_of_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2 次式は根があれば分裂する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def quadTwist.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(二次の捧りの定義。★無条件)",
    sectionId := "genell-lemma-3-5" }

def quadTwist_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(捧りで j は変わらない。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
