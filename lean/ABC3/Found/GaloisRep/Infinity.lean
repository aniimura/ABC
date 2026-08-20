import ABC3.Found.GaloisRep.RedHom

/-!
# Galois (G5) 第 161 ブロック —— **★★★★★★無限遠に落ちる点の付値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`w(x) > 1` なら `w(y)² = w(x)³`

第 160 で `redPoint` が加法を保つことを示したが、それは
**両方の点の座標が付値環に入る場合**である。★`n • S` との可換性には
「無限遠に落ちる点の全体が部分群である」ことが要り、その入口が本ブロックである。

## ★★★★機構——第 145 と同じ形

曲線上の点では**平方完成**

    (2y + a₁x + a₃)² = 4x³ + (a₁²+4a₂)x² + (4a₄+2a₁a₃)x + (4a₆+a₃²)

が成り立つ(`linear_combination 4·heq` の 1 行)。
★`w(x) > 1` なら右辺は `4x³` が単独で最大なので `w(RHS) = w(x)³`。
★★左辺は `w(2y+a₁x+a₃)²`。したがって `w(2y+a₁x+a₃)² = w(x)³ > w(x)²` から
`w(2y+a₁x+a₃) > w(x)` となり、`a₁x + a₃` の項は落ちて `w(y) = w(2y+a₁x+a₃)`。

★★★これは第 145(生成点の場合)と**同じ形**である——
そこでは `z² = Ψ₂Sq(x)` と `deg Ψ₂Sq = 3` を使った。

## ★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psi2_sq_eq` | ★★★★★★**曲線上の点での平方完成**(任意の可換環) |
| `val_psi_sq` | ★★★★★`w(2y+a₁x+a₃)² = w(x)³` |
| `val_y_of_psi` | ★★★★`w(y) = w(2y+a₁x+a₃)` |
| `val_y_sq_eq_val_x_cube` | ★★★★★★**`w(y)² = w(x)³`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-- ★★★★★★**曲線上の点での平方完成**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`linear_combination 4·heq` の 1 行で出る(任意の可換環)。 -/
theorem psi2_sq_eq {R : Type} [CommRing R] (W : WeierstrassCurve.Affine R) {x y : R}
    (heq : W.Equation x y) :
    (2 * y + W.a₁ * x + W.a₃) ^ 2
      = 4 * x ^ 3 + (W.a₁ ^ 2 + 4 * W.a₂) * x ^ 2
        + (4 * W.a₄ + 2 * W.a₁ * W.a₃) * x + (4 * W.a₆ + W.a₃ ^ 2) := by
  rw [equation_iff] at heq
  linear_combination 4 * heq

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)

/-- ★★★★★`w(2y+a₁x+a₃)² = w(x)³`——右辺は `4x³` が単独で最大。 -/
theorem val_psi_sq (h2 : IsUnit (2 : F))
    {x y : W.FunctionField} (heq : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hx : 1 < v.valuation W.FunctionField x) :
    (v.valuation W.FunctionField
        (2 * y + algebraMap F W.FunctionField W.a₁ * x
          + algebraMap F W.FunctionField W.a₃)) ^ 2
      = (v.valuation W.FunctionField x) ^ 3 := by
  set w := v.valuation W.FunctionField with hw
  set A := algebraMap F W.FunctionField with hA
  have h4 : (4 : F) ≠ 0 := by
    have h44 : (4 : F) = 2 * 2 := by norm_num
    rw [h44]; exact mul_ne_zero h2.ne_zero h2.ne_zero
  have hw4 : w (4 : W.FunctionField) = 1 := by
    have he : (4 : W.FunctionField) = A (4 : F) := by rw [hA, map_ofNat]
    rw [he, hA, IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v ((isUnit_iff_ne_zero.2 h4).map
      (algebraMap F W.CoordinateRing))
  have hw2 : w (2 : W.FunctionField) = 1 := by
    have he : (2 : W.FunctionField) = A (2 : F) := by rw [hA, map_ofNat]
    rw [he, hA, IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v (h2.map (algebraMap F W.CoordinateRing))
  have hid := psi2_sq_eq (W.map A) heq
  simp only [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₃,
    WeierstrassCurve.map_a₄, WeierstrassCurve.map_a₆] at hid
  have hgroup : 4 * x ^ 3 + (A W.a₁ ^ 2 + 4 * A W.a₂) * x ^ 2
      + (4 * A W.a₄ + 2 * A W.a₁ * A W.a₃) * x + (4 * A W.a₆ + A W.a₃ ^ 2)
      = 4 * x ^ 3 + ((A W.a₁ ^ 2 + 4 * A W.a₂) * x ^ 2
        + (4 * A W.a₄ + 2 * A W.a₁ * A W.a₃) * x + (4 * A W.a₆ + A W.a₃ ^ 2)) := by ring
  rw [hgroup] at hid
  have hx1 : w x < w x ^ 2 := by
    have hp := pow_lt_pow_right₀ hx (by norm_num : 1 < 2); rwa [pow_one] at hp
  have hx2 : w x ^ 2 < w x ^ 3 := pow_lt_pow_right₀ hx (by norm_num : 2 < 3)
  have hrest : w ((A W.a₁ ^ 2 + 4 * A W.a₂) * x ^ 2
      + (4 * A W.a₄ + 2 * A W.a₁ * A W.a₃) * x + (4 * A W.a₆ + A W.a₃ ^ 2)) < w x ^ 3 := by
    refine lt_of_le_of_lt (le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)) hx2
    · refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
      · rw [Valuation.map_mul, Valuation.map_pow]
        calc w (A W.a₁ ^ 2 + 4 * A W.a₂) * w x ^ 2 ≤ 1 * w x ^ 2 := by
              refine mul_le_mul' ?_ le_rfl
              refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
              · rw [Valuation.map_pow]
                exact pow_le_one' (by rw [hA]; exact valuation_algebraMap_field W v _) 2
              · rw [Valuation.map_mul, hw4, one_mul, hA]
                exact valuation_algebraMap_field W v _
          _ = w x ^ 2 := one_mul _
      · rw [Valuation.map_mul]
        calc w (4 * A W.a₄ + 2 * A W.a₁ * A W.a₃) * w x ≤ 1 * w x := by
              refine mul_le_mul' ?_ le_rfl
              refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
              · rw [Valuation.map_mul, hw4, one_mul, hA]
                exact valuation_algebraMap_field W v _
              · rw [Valuation.map_mul, Valuation.map_mul, hw2, one_mul, hA]
                calc w (algebraMap F W.FunctionField W.a₁)
                      * w (algebraMap F W.FunctionField W.a₃) ≤ 1 * 1 :=
                      mul_le_mul' (valuation_algebraMap_field W v _)
                        (valuation_algebraMap_field W v _)
                    _ = 1 := one_mul 1
          _ = w x := one_mul _
          _ ≤ w x ^ 2 := le_of_lt hx1
    · refine le_trans (le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)) (le_of_lt hx1)
      · rw [Valuation.map_mul, hw4, one_mul, hA]
        exact le_trans (valuation_algebraMap_field W v _) (le_of_lt hx)
      · rw [Valuation.map_pow]
        exact le_trans (pow_le_one' (by rw [hA]; exact valuation_algebraMap_field W v _) 2)
          (le_of_lt hx)
  rw [← Valuation.map_pow, hid, Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_mul, hw4, one_mul, Valuation.map_pow]
      exact ne_of_gt hrest), Valuation.map_mul, hw4, one_mul, Valuation.map_pow,
    max_eq_left (le_of_lt hrest)]

/-- ★★★★`w(y) = w(2y+a₁x+a₃)` から `w(y)² = w(x)³`。 -/
theorem val_y_of_psi (h2 : IsUnit (2 : F))
    {x y : W.FunctionField}
    (hx : 1 < v.valuation W.FunctionField x)
    (hpsi : (v.valuation W.FunctionField
        (2 * y + algebraMap F W.FunctionField W.a₁ * x
          + algebraMap F W.FunctionField W.a₃)) ^ 2
      = (v.valuation W.FunctionField x) ^ 3) :
    (v.valuation W.FunctionField y) ^ 2 = (v.valuation W.FunctionField x) ^ 3 := by
  set w := v.valuation W.FunctionField with hw
  set A := algebraMap F W.FunctionField with hA
  set Z := 2 * y + A W.a₁ * x + A W.a₃ with hZ
  have hx2 : w x ^ 2 < w x ^ 3 := pow_lt_pow_right₀ hx (by norm_num : 2 < 3)
  have hZgt : w x < w Z := by
    by_contra hcon
    rw [not_lt] at hcon
    have hle : w Z ^ 2 ≤ w x ^ 2 := pow_le_pow_left' hcon 2
    rw [hpsi] at hle
    exact absurd hle (not_le.2 hx2)
  have h2u : w (2 : W.FunctionField) = 1 := by
    have h2e : (2 : W.FunctionField) = A (2 : F) := by rw [hA, map_ofNat]
    rw [h2e, hA, IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v (h2.map (algebraMap F W.CoordinateRing))
  have hrest : w (A W.a₁ * x + A W.a₃) ≤ w x := by
    refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
    · rw [Valuation.map_mul, hA]
      calc w (algebraMap F W.FunctionField W.a₁) * w x ≤ 1 * w x :=
            mul_le_mul' (valuation_algebraMap_field W v _) le_rfl
        _ = w x := one_mul _
    · rw [hA]
      exact le_trans (valuation_algebraMap_field W v _) (le_of_lt hx)
  have hdec : (2 : W.FunctionField) * y = Z - (A W.a₁ * x + A W.a₃) := by rw [hZ]; ring
  have h2y : w ((2 : W.FunctionField) * y) = w Z := by
    rw [hdec, sub_eq_add_neg, Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_neg]
      exact ne_of_gt (lt_of_le_of_lt hrest hZgt)), Valuation.map_neg,
      max_eq_left (le_of_lt (lt_of_le_of_lt hrest hZgt))]
  rw [Valuation.map_mul, h2u, one_mul] at h2y
  rw [h2y, hpsi]

/-- ★★★★★★**無限遠に落ちる点では `w(y)² = w(x)³`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが「無限遠に落ちる点の全体が部分群である」ことの入口になる。 -/
theorem val_y_sq_eq_val_x_cube (h2 : IsUnit (2 : F))
    {x y : W.FunctionField} (heq : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hx : 1 < v.valuation W.FunctionField x) :
    (v.valuation W.FunctionField y) ^ 2 = (v.valuation W.FunctionField x) ^ 3 :=
  val_y_of_psi W v h2 hx (val_psi_sq W v h2 heq hx)

/-! ## ★出典の紐付け(`.src`) -/

def val_y_sq_eq_val_x_cube.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——無限遠に落ちる点での付値の関係)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
