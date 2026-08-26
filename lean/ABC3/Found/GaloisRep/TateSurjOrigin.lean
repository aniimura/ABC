import ABC3.Found.GaloisRep.TateFormalParam

/-!
# Galois (G6) 第 288 ブロック —— **★★★★★★★★★★原点近傍での全射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点——葉 (e) が閉じた

> `1/x ∈ 𝔪` で `(x,y)` が曲線上なら、母数 `u`(`1 − u ∈ 𝔪`)が在って
> `X(u) = x`、`Y(u) = y`(`tate_surjective_origin`)

★★★これで葉 (e)(全射性)が**3 領域すべてで**閉じた。

## ★★★★★★★`z` が `s` を決める——縮小は要らない

`z := x/y`、`s := −1/y` と置くと、曲線の式は

    s·(1 + z + a₄·z·s − a₆·s²) = −z³

★★★★同じ `z` に対して `s, s' ∈ 𝔪` が 2 つあれば、差を括って

    (s − s')·[1 + (𝔪 の元)] = 0

★★**`1 + 𝔪` は単元**なので `s = s'`。★★★★★不動点定理も逐次近似も要らなかった——
**`z³` の側が同じなら差が消える**という一行である。

★仮定に `s, s' ∈ 𝔪` すら要らない(`z ∈ 𝔪` と `a₆ ∈ 𝔪` だけで括弧が単元になる)。

## ★★★★★★分母を払った形が `s` を `R` に運ぶ

`s = −1/Y` は一見 `K` の元だが、`Y = YE/(1−u)³` で `YE` は単元(第 280)なので

    s = −(1−u)³ · YE⁻¹ ∈ 𝔪

★★これが第 275 の分母払いの 4 度目の働きである。

## ★★組み立て

| 段 | 内容 |
|---|---|
| 1 | `tate_formal_param`(第 287)で `z = x/y ∈ 𝔪` |
| 2 | `s := −(z·α) ∈ 𝔪`(`α = 1/x`)——**積だから `𝔪` に入る** |
| 3 | `exists_tateZ_eq`(第 276)で `tateZ(u) = −z` なる `u` |
| 4 | 両者の `s` が同じ 3 次関係式を満たす → `s` 一致 |
| 5 | `s·y = −1 = s·Y` → `y = Y`、`z·y = x`・`z·Y = X` → `x = X` |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_s_unique` | ★★★★★★★**`z` は `s` を決める** |
| `tate_s_relation` | ★★★★★★曲線上の点が満たす 3 次関係式 |
| `tateS`・`tateS_mul_tateYK` | ★★★★★原点近傍の点の `s` |
| `tate_surjective_origin` | ★★★★★★★★★★**原点近傍での全射性** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★★★★★★`z` は `s` を決める -/

section Unique

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★★★**形式母数 `z` は `s = −1/y` を決める**——`1 + 𝔪` が単元だから。 -/
theorem tate_s_unique {z s s' A B : R} (hz : z ∈ I) (hs : s ∈ I) (hs' : s' ∈ I)
    (hA : A ∈ I) (hB : B ∈ I)
    (h : s * (1 + z + A * z * s - B * s ^ 2) = -z ^ 3)
    (h' : s' * (1 + z + A * z * s' - B * s' ^ 2) = -z ^ 3) : s = s' := by
  have hfac : (s - s') * (1 + (z + A * z * (s + s') - B * (s ^ 2 + s * s' + s' ^ 2))) = 0 := by
    linear_combination h - h'
  have hmem : z + A * z * (s + s') - B * (s ^ 2 + s * s' + s' ^ 2) ∈ I := by
    refine Ideal.sub_mem _ (Ideal.add_mem _ hz ?_) (Ideal.mul_mem_right _ _ hB)
    exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ hz)
  have hu : IsUnit (1 + (z + A * z * (s + s') - B * (s ^ 2 + s * s' + s' ^ 2))) :=
    isUnit_add_mem isUnit_one hmem
  exact sub_eq_zero.1 ((hu.mul_left_eq_zero).1 hfac)

end Unique

/-! ## ★★★★★★曲線上の点が満たす 3 次関係式 -/

section Relation

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★曲線上の点が満たす `s` の 3 次関係式。 -/
theorem tate_s_relation (hinj : Function.Injective (algebraMap R K))
    (q : R) (hq : q ∈ I) (x y : K) (z₀ s₀ : R)
    (hzy : algebraMap R K z₀ * y = x) (hsy : algebraMap R K s₀ * y = -1)
    (he : y ^ 2 + x * y = x ^ 3 + algebraMap R K ((tateCurveAt q hq).a₄) * x
      + algebraMap R K ((tateCurveAt q hq).a₆)) :
    s₀ * (1 + z₀ + (tateCurveAt q hq).a₄ * z₀ * s₀ - (tateCurveAt q hq).a₆ * s₀ ^ 2)
      = -z₀ ^ 3 := by
  set A := (tateCurveAt q hq).a₄ with hA
  set B := (tateCurveAt q hq).a₆ with hB
  set S := algebraMap R K s₀ with hS
  set Z := algebraMap R K z₀ with hZ
  have hSx : S * x = -Z := by rw [← hzy]; linear_combination Z * hsy
  refine hinj ?_
  simp only [map_mul, map_sub, map_add, map_one, map_pow, map_neg, ← hS, ← hZ]
  have e1 : S ^ 3 * y ^ 2 = S := by linear_combination (S * (S * y) - S) * hsy
  have e2 : S ^ 3 * (x * y) = S * Z := by
    linear_combination (S ^ 2 * x) * hsy + (-S) * hSx
  have e3 : S ^ 3 * x ^ 3 = -Z ^ 3 := by
    linear_combination ((S * x) ^ 2 - (S * x) * Z + Z ^ 2) * hSx
  have e4 : S ^ 3 * x = -(S ^ 2 * Z) := by linear_combination S ^ 2 * hSx
  linear_combination S ^ 3 * he - e1 - e2 + e3 + algebraMap R K A * e4

/-! ## ★★★★★原点近傍の点の `s` -/

/-- ★★★★★原点近傍の点の `s = −1/Y`——分母を払えば `R` の元。 -/
noncomputable def tateS (u w q : R) (hq : q ∈ I) : R :=
  -((1 - u) ^ 3 * Ring.inverse (tateYpairE u w q hq))

theorem tateS_mem (u w q : R) (hq : q ∈ I) (hu : 1 - u ∈ I) : tateS u w q hq ∈ I := by
  rw [tateS]
  refine neg_mem (Ideal.mul_mem_right _ _ ?_)
  rw [pow_succ]
  exact Ideal.mul_mem_left _ _ hu

theorem tateS_mul_tateYK (u w q : R) (hq : q ∈ I) (hu : 1 - u ∈ I)
    (hne : algebraMap R K (1 - u) ≠ 0) :
    algebraMap R K (tateS u w q hq) * tateYK (K := K) u w q hq = -1 := by
  have hYE := isUnit_tateYpairE u w q hq hu
  have hYne : algebraMap R K (tateYpairE u w q hq) ≠ 0 := (hYE.map (algebraMap R K)).ne_zero
  have hinv : algebraMap R K (Ring.inverse (tateYpairE u w q hq))
      = (algebraMap R K (tateYpairE u w q hq))⁻¹ := by
    have h2 := congrArg (algebraMap R K) (Ring.mul_inverse_cancel _ hYE)
    rw [map_mul, map_one] at h2
    exact (inv_eq_of_mul_eq_one_right h2).symm
  rw [tateS, tateYK, map_neg, map_mul, map_pow, hinv]
  field_simp

/-! ## ★★★★★★★★★★原点近傍での全射性 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**原点近傍での全射性**——葉 (e) が 3 領域すべてで閉じた。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_surjective_origin (hinj : Function.Injective (algebraMap R K)) (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (q : R) (hq : q ∈ I) (x y : K) (hx0 : x ≠ 0) (α₀ : R) (hα : α₀ ∈ I)
    (hαx : algebraMap R K α₀ = x⁻¹)
    (he : y ^ 2 + x * y = x ^ 3 + algebraMap R K ((tateCurveAt q hq).a₄) * x
      + algebraMap R K ((tateCurveAt q hq).a₆)) :
    ∃ u : R, 1 - u ∈ I ∧ algebraMap R K (1 - u) ≠ 0 ∧
      tateXK u (wOf q u) q hq = x ∧ tateYK (K := K) u (wOf q u) q hq = y := by
  have hax : algebraMap R K α₀ * x = 1 := by rw [hαx, inv_mul_cancel₀ hx0]
  obtain ⟨z₀, hz₀, hzy⟩ := tate_formal_param hinj hI hquad q hq x y hx0 α₀ hα hαx he
  set s₀ := -(z₀ * α₀) with hs₀def
  have hs₀ : s₀ ∈ I := neg_mem (Ideal.mul_mem_right _ _ hz₀)
  have hsy : algebraMap R K s₀ * y = -1 := by
    rw [hs₀def, map_neg, map_mul]
    linear_combination (-(algebraMap R K α₀)) * hzy - hax
  have hrel := tate_s_relation hinj q hq x y z₀ s₀ hzy hsy he
  have hz0ne : z₀ ≠ 0 := by
    intro h
    rw [h, map_zero, zero_mul] at hzy
    exact hx0 hzy.symm
  obtain ⟨u, hu, hZu⟩ := exists_tateZ_eq (-z₀) q (neg_mem hz₀) hq
  have hu1 : (1 : R) - u ≠ 0 := by
    intro h
    have hzero : tateZ u (wOf q u) q hq = 0 := by
      rw [tateZ, show u - 1 = -(1 - u) by ring, h, neg_zero, zero_mul, zero_mul]
    rw [hZu] at hzero
    exact hz0ne (by linear_combination -hzero)
  have hne : algebraMap R K (1 - u) ≠ 0 := fun h => hu1 (hinj (by rw [h, map_zero]))
  have hSY := tateS_mul_tateYK (K := K) u (wOf q u) q hq hu hne
  have hZK : algebraMap R K (-(tateZ u (wOf q u) q hq)) * tateYK (K := K) u (wOf q u) q hq
      = tateXK u (wOf q u) q hq := by
    have hYK0 : tateYK (K := K) u (wOf q u) q hq ≠ 0 := by
      intro h
      rw [h, mul_zero] at hSY
      exact (by norm_num : (0 : K) ≠ -1) hSY
    rw [map_neg, tateZ_map u (wOf q u) q hq hu hne]
    field_simp
  have heK := tate_equationK (K := K) u (wOf q u) q hq
    (mul_wOf q u (isUnit_of_one_sub_mem hu)) (isUnit_one_sub (wOf_mem u hq)) hne
  have hrel' := tate_s_relation hinj q hq (tateXK u (wOf q u) q hq)
    (tateYK (K := K) u (wOf q u) q hq) (-(tateZ u (wOf q u) q hq))
    (tateS u (wOf q u) q hq) hZK hSY heK
  rw [hZu, neg_neg] at hrel'
  have hss : s₀ = tateS u (wOf q u) q hq :=
    tate_s_unique hz₀ hs₀ (tateS_mem u (wOf q u) q hq hu)
      (tateCurveAt_a4_mem (I := I) q hq) (tateCurveAt_a6_mem (I := I) q hq) hrel hrel'
  have hs0ne : algebraMap R K s₀ ≠ 0 := by
    intro h
    rw [h, zero_mul] at hsy
    exact (by norm_num : (0 : K) ≠ -1) hsy
  have hy : y = tateYK (K := K) u (wOf q u) q hq := by
    refine mul_left_cancel₀ hs0ne ?_
    rw [hsy, hss, hSY]
  refine ⟨u, hu, hne, ?_, hy.symm⟩
  rw [← hZK, hZu, neg_neg, ← hy, hzy]

end Relation

/-! ## ★出典の紐付け(`.src`) -/

def tate_s_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——z は s を決める)",
    sectionId := "genell-def-3-3" }

def tate_surjective_origin.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍での全射性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
