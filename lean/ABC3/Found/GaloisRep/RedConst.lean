import ABC3.Found.GaloisRep.Reduction

/-!
# Galois (G5) 第 153 ブロック —— **★★★★★★還元は環準同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★還元 `red` の代数的性質

第 152 で `w t ≤ 1 ⟹ ∃! d ∈ F, w(t − d) < 1` が取れた。
★その `d` を `red(t)` と書くと、`red` は付値環の上で**環準同型**である:

| 性質 | 証明 |
|---|---|
| `red(s + t) = red s + red t` | `(s+t) − (rs+rt) = (s−rs) + (t−rt)`、超距離で max < 1 |
| `red(s·t) = red s · red t` | `st − rs·rt = s(t−rt) + rt(s−rs)`、`w s ≤ 1`・`w(rt) ≤ 1` |
| `red(−t) = −red t` | 付値は符号を見ない |
| `red(s/t) = red s / red t`(`w t = 1`) | `red(s/t)·red(t) = red(s)` |
| `red(d) = d`(`d ∈ F`) | 定数は動かない |

★★★**`w t = 1` なら `red t ≠ 0`** も出る——これが加法公式の場合分けで
「傾きの分母が還元後も 0 でない」ことを言うのに使う。

## ★★★★これが加法公式の場合分けの土台である

`slope x₁ x₂ y₁ y₂ = (y₁ − y₂)/(x₁ − x₂)`(`x₁ ≠ x₂` のとき)に `red` を当てると、
`w(x₁ − x₂) = 1` すなわち `red x₁ ≠ red x₂` の場合には

    red(slope) = (red y₁ − red y₂)/(red x₁ − red x₂)

となり、**還元先の傾きと一致する**。★残るのは `red x₁ = red x₂` の場合である。

## ★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `redConst` | ★★★★★還元 |
| `redConst_spec` / `redConst_eq` | ★★★★特徴づけ |
| `redConst_algebraMap` | ★定数は動かない |
| `redConst_add` / `redConst_mul` | ★★★★★★**環準同型** |
| `redConst_neg` / `redConst_sub` | ★★ |
| `redConst_ne_zero_of_eq_one` | ★★★**`w t = 1` なら `red t ≠ 0`** |
| `redConst_div` | ★★★★除算(傾きに要る) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)
  {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

open Classical in
/-- ★★★★★**還元**——付値環の元を定数へ。 -/
noncomputable def redConst (t : W.FunctionField) : F :=
  if ht : v.valuation W.FunctionField t ≤ 1 then
    Classical.choose (exists_const_of_valuation_le_one W v h hv ht) else 0

theorem redConst_spec {t : W.FunctionField} (ht : v.valuation W.FunctionField t ≤ 1) :
    v.valuation W.FunctionField
      (t - algebraMap F W.FunctionField (redConst W v h hv t)) < 1 := by
  rw [redConst, dif_pos ht]
  exact Classical.choose_spec (exists_const_of_valuation_le_one W v h hv ht)

theorem redConst_eq {t : W.FunctionField} (ht : v.valuation W.FunctionField t ≤ 1) {d : F}
    (hd : v.valuation W.FunctionField (t - algebraMap F W.FunctionField d) < 1) :
    redConst W v h hv t = d :=
  const_unique W v (redConst_spec W v h hv ht) hd

/-! ## ★★★★★★環準同型であること -/

theorem redConst_algebraMap (d : F) :
    redConst W v h hv (algebraMap F W.FunctionField d) = d := by
  refine redConst_eq W v h hv (valuation_algebraMap_field W v d) ?_
  rw [sub_self, Valuation.map_zero]
  exact zero_lt_one

/-- ★★★★★★**還元は加法を保つ**。 -/
theorem redConst_add {s t : W.FunctionField}
    (hs : v.valuation W.FunctionField s ≤ 1) (ht : v.valuation W.FunctionField t ≤ 1) :
    redConst W v h hv (s + t) = redConst W v h hv s + redConst W v h hv t := by
  refine redConst_eq W v h hv (le_trans (Valuation.map_add _ s t) (max_le hs ht)) ?_
  have heq : s + t - algebraMap F W.FunctionField (redConst W v h hv s + redConst W v h hv t)
      = (s - algebraMap F W.FunctionField (redConst W v h hv s))
        + (t - algebraMap F W.FunctionField (redConst W v h hv t)) := by
    rw [map_add]; ring
  rw [heq]
  exact lt_of_le_of_lt (Valuation.map_add _ _ _)
    (max_lt (redConst_spec W v h hv hs) (redConst_spec W v h hv ht))

/-- ★★★★★★**還元は乗法を保つ**。 -/
theorem redConst_mul {s t : W.FunctionField}
    (hs : v.valuation W.FunctionField s ≤ 1) (ht : v.valuation W.FunctionField t ≤ 1) :
    redConst W v h hv (s * t) = redConst W v h hv s * redConst W v h hv t := by
  set w := v.valuation W.FunctionField with hw
  have hst : w (s * t) ≤ 1 := by
    rw [Valuation.map_mul]
    calc w s * w t ≤ 1 * 1 := mul_le_mul' hs ht
      _ = 1 := one_mul 1
  refine redConst_eq W v h hv hst ?_
  have heq : s * t - algebraMap F W.FunctionField (redConst W v h hv s * redConst W v h hv t)
      = s * (t - algebraMap F W.FunctionField (redConst W v h hv t))
        + algebraMap F W.FunctionField (redConst W v h hv t)
          * (s - algebraMap F W.FunctionField (redConst W v h hv s)) := by
    rw [map_mul]; ring
  rw [heq]
  refine lt_of_le_of_lt (Valuation.map_add _ _ _) (max_lt ?_ ?_)
  · rw [Valuation.map_mul]
    calc w s * w (t - algebraMap F W.FunctionField (redConst W v h hv t))
        ≤ 1 * w (t - algebraMap F W.FunctionField (redConst W v h hv t)) := mul_le_mul' hs le_rfl
      _ = w (t - algebraMap F W.FunctionField (redConst W v h hv t)) := one_mul _
      _ < 1 := redConst_spec W v h hv ht
  · rw [Valuation.map_mul]
    calc w (algebraMap F W.FunctionField (redConst W v h hv t))
          * w (s - algebraMap F W.FunctionField (redConst W v h hv s))
        ≤ 1 * w (s - algebraMap F W.FunctionField (redConst W v h hv s)) :=
          mul_le_mul' (valuation_algebraMap_field W v _) le_rfl
      _ = w (s - algebraMap F W.FunctionField (redConst W v h hv s)) := one_mul _
      _ < 1 := redConst_spec W v h hv hs

theorem redConst_neg {t : W.FunctionField} (ht : v.valuation W.FunctionField t ≤ 1) :
    redConst W v h hv (-t) = -redConst W v h hv t := by
  refine redConst_eq W v h hv (by rwa [Valuation.map_neg]) ?_
  have heq : -t - algebraMap F W.FunctionField (-redConst W v h hv t)
      = -(t - algebraMap F W.FunctionField (redConst W v h hv t)) := by rw [map_neg]; ring
  rw [heq, Valuation.map_neg]
  exact redConst_spec W v h hv ht

theorem redConst_sub {s t : W.FunctionField}
    (hs : v.valuation W.FunctionField s ≤ 1) (ht : v.valuation W.FunctionField t ≤ 1) :
    redConst W v h hv (s - t) = redConst W v h hv s - redConst W v h hv t := by
  rw [sub_eq_add_neg, redConst_add W v h hv hs (by rwa [Valuation.map_neg]),
    redConst_neg W v h hv ht, sub_eq_add_neg]

/-- ★★★**`w t = 1` なら `red t ≠ 0`**——傾きの分母に効く。 -/
theorem redConst_ne_zero_of_eq_one {t : W.FunctionField}
    (ht : v.valuation W.FunctionField t = 1) : redConst W v h hv t ≠ 0 := by
  intro h0
  have hs := redConst_spec W v h hv (le_of_eq ht)
  rw [h0, map_zero, sub_zero, ht] at hs
  exact absurd hs (lt_irrefl 1)

/-- ★★★★**還元は除算を保つ**(分母の付値が 1 のとき)。 -/
theorem redConst_div {s t : W.FunctionField}
    (hs : v.valuation W.FunctionField s ≤ 1) (ht : v.valuation W.FunctionField t = 1) :
    redConst W v h hv (s / t) = redConst W v h hv s / redConst W v h hv t := by
  have htne : t ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at ht
    exact absurd ht (by simp)
  have hdiv : v.valuation W.FunctionField (s / t) ≤ 1 := by
    rw [div_eq_mul_inv, Valuation.map_mul, map_inv₀, ht, inv_one, mul_one]
    exact hs
  have hkey : redConst W v h hv (s / t) * redConst W v h hv t = redConst W v h hv s := by
    rw [← redConst_mul W v h hv hdiv (le_of_eq ht), div_mul_cancel₀ _ htne]
  field_simp [redConst_ne_zero_of_eq_one W v h hv ht] at hkey ⊢
  exact hkey

/-! ## ★出典の紐付け(`.src`) -/

def redConst_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——還元が環準同型であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
