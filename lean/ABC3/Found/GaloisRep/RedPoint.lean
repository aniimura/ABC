import ABC3.Found.GaloisRep.RedConst

/-!
# Galois (G5) 第 154 ブロック —— **★★★★★★点の還元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★点の還元が定義できた

第 153 で `red` が環準同型だと分かった。★したがって

    (W ⊗ F(W)).Equation x y  かつ  w x ≤ 1, w y ≤ 1
      ⟹  W.Equation (red x) (red y)

が出る(方程式の係数は `F` にあるので `red` で動かない)。
★★楕円曲線では方程式を満たす点は自動的に非特異(`equation_iff_nonsingular`)なので、

    redPoint : (W ⊗ F(W)).Point → W.Point

が定義できる。★★★座標が付値環に入らない点は `0` に送る(無限遠へ落ちる点)。

## ★★★★残るのは加法性だけ

    redPoint (S₁ + S₂) = redPoint S₁ + redPoint S₂

★`red x₁ ≠ red x₂` の場合は第 153 の `redConst_div` で傾きがそのまま還元される。
★★残るのは `red x₁ = red x₂` の場合(還元先で 2 倍公式に切り替わる)である。

## ★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `val_mul_le` / `val_pow_le` / `val_add_le` | ★付値が 1 以下であることの保存 |
| `redConst_pow` | ★★`red(t^k) = (red t)^k` |
| `equation_redConst` | ★★★★★**還元した座標は方程式を満たす** |
| `redPoint` | ★★★★★★**点の還元** |
| `redPoint_zero` / `redPoint_some` | ★計算則 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)

/-! ## ★付値が 1 以下であることの保存 -/

theorem val_mul_le {s t : W.FunctionField} (hs : v.valuation W.FunctionField s ≤ 1)
    (ht : v.valuation W.FunctionField t ≤ 1) :
    v.valuation W.FunctionField (s * t) ≤ 1 := by
  rw [Valuation.map_mul]
  calc v.valuation W.FunctionField s * v.valuation W.FunctionField t ≤ 1 * 1 :=
        mul_le_mul' hs ht
    _ = 1 := one_mul 1

theorem val_pow_le {t : W.FunctionField} (ht : v.valuation W.FunctionField t ≤ 1) (k : ℕ) :
    v.valuation W.FunctionField (t ^ k) ≤ 1 := by
  rw [Valuation.map_pow]; exact pow_le_one' ht k

theorem val_add_le {s t : W.FunctionField} (hs : v.valuation W.FunctionField s ≤ 1)
    (ht : v.valuation W.FunctionField t ≤ 1) :
    v.valuation W.FunctionField (s + t) ≤ 1 :=
  le_trans (Valuation.map_add _ _ _) (max_le hs ht)

variable {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

theorem redConst_pow {t : W.FunctionField} (ht : v.valuation W.FunctionField t ≤ 1) (k : ℕ) :
    redConst W v h hv (t ^ k) = (redConst W v h hv t) ^ k := by
  induction k with
  | zero =>
    rw [pow_zero, pow_zero]
    have h1 : (1 : W.FunctionField) = algebraMap F W.FunctionField 1 := (map_one _).symm
    rw [h1, redConst_algebraMap]
  | succ k ih =>
    rw [pow_succ, pow_succ, redConst_mul W v h hv (val_pow_le W v ht k) ht, ih]

/-! ## ★★★★★還元した座標は方程式を満たす -/

/-- ★★★★★**還元した座標は方程式を満たす**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★方程式の係数は `F` にあるので `red` で動かない。 -/
theorem equation_redConst {x y : W.FunctionField}
    (heq : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hx : v.valuation W.FunctionField x ≤ 1) (hy : v.valuation W.FunctionField y ≤ 1) :
    W.Equation (redConst W v h hv x) (redConst W v h hv y) := by
  set A := algebraMap F W.FunctionField with hA
  rw [equation_iff] at heq ⊢
  simp only [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₃,
    WeierstrassCurve.map_a₄, WeierstrassCurve.map_a₆] at heq
  have hA1 : v.valuation W.FunctionField (A W.a₁) ≤ 1 := valuation_algebraMap_field W v _
  have hA2 : v.valuation W.FunctionField (A W.a₂) ≤ 1 := valuation_algebraMap_field W v _
  have hA3 : v.valuation W.FunctionField (A W.a₃) ≤ 1 := valuation_algebraMap_field W v _
  have hA4 : v.valuation W.FunctionField (A W.a₄) ≤ 1 := valuation_algebraMap_field W v _
  have hA6 : v.valuation W.FunctionField (A W.a₆) ≤ 1 := valuation_algebraMap_field W v _
  have hq := congrArg (redConst W v h hv) heq
  rw [redConst_add W v h hv
        (val_add_le W v (val_pow_le W v hy 2) (val_mul_le W v (val_mul_le W v hA1 hx) hy))
        (val_mul_le W v hA3 hy),
      redConst_add W v h hv (val_pow_le W v hy 2) (val_mul_le W v (val_mul_le W v hA1 hx) hy),
      redConst_pow W v h hv hy 2,
      redConst_mul W v h hv (val_mul_le W v hA1 hx) hy,
      redConst_mul W v h hv hA1 hx,
      redConst_mul W v h hv hA3 hy,
      redConst_add W v h hv
        (val_add_le W v (val_add_le W v (val_pow_le W v hx 3) (val_mul_le W v hA2
          (val_pow_le W v hx 2))) (val_mul_le W v hA4 hx)) hA6,
      redConst_add W v h hv (val_add_le W v (val_pow_le W v hx 3) (val_mul_le W v hA2
          (val_pow_le W v hx 2))) (val_mul_le W v hA4 hx),
      redConst_add W v h hv (val_pow_le W v hx 3) (val_mul_le W v hA2 (val_pow_le W v hx 2)),
      redConst_pow W v h hv hx 3,
      redConst_mul W v h hv hA2 (val_pow_le W v hx 2),
      redConst_pow W v h hv hx 2,
      redConst_mul W v h hv hA4 hx,
      hA, redConst_algebraMap, redConst_algebraMap, redConst_algebraMap,
      redConst_algebraMap, redConst_algebraMap] at hq
  exact hq

/-! ## ★★★★★★点の還元 -/

open Classical in
/-- ★★★★★★**点の還元**——座標が付値環に入らない点は `0` に送る。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★楕円曲線では方程式を満たす点は自動的に非特異である。 -/
noncomputable def redPoint [W.IsElliptic] (S : (W.map (algebraMap F W.FunctionField)).Point) :
    W.Point :=
  match S with
  | 0 => 0
  | Point.some x y hns =>
      if hxy : v.valuation W.FunctionField x ≤ 1 ∧ v.valuation W.FunctionField y ≤ 1 then
        Point.some (redConst W v h hv x) (redConst W v h hv y)
          (equation_iff_nonsingular.mp (equation_redConst W v h hv hns.1 hxy.1 hxy.2))
      else 0

theorem redPoint_zero [W.IsElliptic] : redPoint W v h hv 0 = 0 := rfl

theorem redPoint_some [W.IsElliptic] {x y : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular x y)
    (hx : v.valuation W.FunctionField x ≤ 1) (hy : v.valuation W.FunctionField y ≤ 1) :
    redPoint W v h hv (Point.some x y hns)
      = Point.some (redConst W v h hv x) (redConst W v h hv y)
        (equation_iff_nonsingular.mp (equation_redConst W v h hv hns.1 hx hy)) := by
  rw [redPoint, dif_pos ⟨hx, hy⟩]

/-! ## ★出典の紐付け(`.src`) -/

def redPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——点の還元写像)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
