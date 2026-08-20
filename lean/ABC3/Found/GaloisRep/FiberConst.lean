import ABC3.Found.GaloisRep.TranslateCount

/-!
# Galois (G5) 第 171 ブロック —— **★★★★★★★★差が `n` 等分点なら指数は等しい**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★`e_v` の一定性が使える形になった

第 170 の `count_translate_eq` は平行移動 `τ` を仮定として要求する。★本ブロックで
**「2 つの点の差が `n` 等分点である」だけを仮定にした形**にする:

    n·(Q_{v'} − Q_v) = 0   ⟹   count_v(μ r) = count_{v'}(μ r)

★★差が 0 の場合は `v = v'` になる(`xyIdeal` は点を決める)。0 でなければ
第 124 の `exists_translateAut_all` で `τ` を取り、第 168・170 を繋ぐだけ。

## ★★★★★★場合 A の判定が点の言葉になった

第 167 で `hnn` が定理になったので、さらに

    n·Q_v ≠ 0   ⟹   ∀ r, w_v(μ r) ≤ 1

が言える(`valuation_mulByN_le_one`)。★`redPoint(n·生成点) = n·Q_v`(第 163・164)と
`redPoint_eq_zero_iff`(第 163)を繋ぐだけである。
★★**場合 A / 場合 B の区別が「`n·Q_v` が 0 かどうか」になった**——
これで台が `[n]⁻¹(P)` と `E[n]∖{O}` の 2 つに分かれることが見える。

## ★★★残るのは総和の組み立てだけ

    S := Σ_v e_v · Q_v
       = α·(Σ_{[n]⁻¹(P)} Q) − β·(Σ_{T ∈ E[n]∖{O}} T)
       = α·(n·P) − β·0 = 0

★**α と β が一致する必要はない**(分岐指数の大域的な一定性は不要)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_some_of_ne_zero` | ★0 でない点はアフィン点 |
| `hom_eq_curveHom` | ★★★★★生成元の像が一致すれば `curveHom` |
| `valuation_hom_le_one` | ★★★★★★座標が整なら像はすべて整 |
| `redPoint_mulByN` | ★★★★★**`red([n] の引き戻し) = n·Q_v`** |
| `valuation_mulByN_le_one` | ★★★★★★**`n·Q_v ≠ 0` なら `hnn`** |
| `count_eq_of_sub_torsion` | ★★★★★★★★**差が `n` 等分点なら指数は等しい** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-- ★0 でない点はアフィン点として書ける。 -/
theorem exists_some_of_ne_zero (S : W.Point) (hS : S ≠ 0) :
    ∃ (a b : F) (hb : W.Nonsingular a b), S = Point.some a b hb := by
  match S with
  | 0 => exact absurd rfl hS
  | Point.some a b hb => exact ⟨a, b, hb, rfl⟩

/-- ★★★★★**生成元の像が一致する環準同型は `curveHom` である**。 -/
theorem hom_eq_curveHom {S : Type} [CommRing S] [Algebra F S] (μ : W.CoordinateRing →+* S)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F S d)
    {x y : S} (heq : (W.map (algebraMap F S)).Equation x y)
    (hμx : μ (genX W) = x) (hμy : μ (genY W) = y) (r : W.CoordinateRing) :
    μ r = curveHom W heq r := by
  have hext : μ = curveHom W heq :=
    coordinateRing_hom_ext _ _ (fun d => by rw [hμF d, curveHom_algebraMap])
      (by rw [hμx, curveHom_genX]) (by rw [hμy, curveHom_genY])
  exact congrFun (congrArg (fun f => (f : W.CoordinateRing →+* S).toFun) hext) r

variable [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)

/-- ★★★★★★座標が付値環に入れば `hnn` は自動。 -/
theorem valuation_hom_le_one (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {x y : W.FunctionField}
    (heq : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hμx : μ (genX W) = x) (hμy : μ (genY W) = y)
    (hx : v.valuation W.FunctionField x ≤ 1) (hy : v.valuation W.FunctionField y ≤ 1)
    (r : W.CoordinateRing) :
    v.valuation W.FunctionField (μ r) ≤ 1 := by
  rw [hom_eq_curveHom W μ hμF heq hμx hμy r]
  exact valuation_curveHom_le_one W v heq hx hy r

variable [DecidableEq F] [W.IsElliptic]
  {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

include hv in
/-- ★★★★★**`[n]` の引き戻しの還元は `n·Q_v`**。 -/
theorem redPoint_mulByN (h2 : IsUnit (2 : F)) (n : ℕ)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns) :
    redPoint W v h hv (Point.some xn yn hns)
      = n • Point.some c y₀ (equation_iff_nonsingular.mp h) := by
  rw [← hμP, redPoint_nsmul W v h hv h2, redPoint_generic W v h hv]

include hv in
/-- ★★★★★★**`n·Q_v ≠ 0` なら `hnn` が成り立つ**——場合 A の判定が点の言葉になる。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 167 で `hnn` が定理になったので、`redPoint_eq_zero_iff` と繋ぐだけで出る。 -/
theorem valuation_mulByN_le_one (h2 : IsUnit (2 : F)) (n : ℕ)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hne : n • Point.some c y₀ (equation_iff_nonsingular.mp h) ≠ 0)
    (r : W.CoordinateRing) :
    v.valuation W.FunctionField (μ r) ≤ 1 := by
  have hred := redPoint_mulByN W v h hv h2 n hns hμP
  have hxle : v.valuation W.FunctionField xn ≤ 1 := by
    by_contra hcon
    rw [not_le] at hcon
    rw [(redPoint_eq_zero_iff W v h hv hns).2 hcon] at hred
    exact hne hred.symm
  exact valuation_hom_le_one W v μ hμF hns.1 hμx hμy hxle
    (val_y_le_of_val_x_le W v hns.1 hxle) r

end ABC3.Found.GaloisRep

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [IsAlgClosed F] [Infinite F] [DecidableEq F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]
  [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★★★★★**差が `n` 等分点なら指数は等しい**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★差が 0 なら `v = v'`。0 でなければ第 124 の `exists_translateAut_all` で `τ` を取り、
第 168(`τ ∘ μ = μ`)と第 170(`e_v` の一定性)を繋ぐ。 -/
theorem count_eq_of_sub_torsion (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (v v' : HeightOneSpectrum W.CoordinateRing)
    {c y₀ c' y₀' : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    (h' : W.Equation c' y₀')
    (hv' : v'.asIdeal = CoordinateRing.XYIdeal W c' (Polynomial.C y₀'))
    (n : ℕ) (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hdiff : n • (Point.some c' y₀' (equation_iff_nonsingular.mp h')
      - Point.some c y₀ (equation_iff_nonsingular.mp h)) = 0)
    (r : W.CoordinateRing) (hz : μ r ≠ 0) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ r))
      = FractionalIdeal.count W.FunctionField v'
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ r)) := by
  by_cases hDz : Point.some c' y₀' (equation_iff_nonsingular.mp h')
      - Point.some c y₀ (equation_iff_nonsingular.mp h) = 0
  · have hpt : Point.some c' y₀' (equation_iff_nonsingular.mp h')
        = Point.some c y₀ (equation_iff_nonsingular.mp h) := sub_eq_zero.1 hDz
    have hcc : c' = c := by injection hpt
    have hyy : y₀' = y₀ := by injection hpt
    subst hcc; subst hyy
    rw [HeightOneSpectrum.ext (hv.trans hv'.symm)]
  · obtain ⟨x₀, y₀T, hQD, hDm⟩ := exists_some_of_ne_zero W _ hDz
    obtain ⟨τ, hx, hy⟩ := exists_translateAut_all W h4 hQD
    have hsum : Point.some c y₀ (equation_iff_nonsingular.mp h) + Point.some x₀ y₀T hQD
        = Point.some c' y₀' (equation_iff_nonsingular.mp h') := by
      rw [← hDm]; abel
    have hT : n • Point.some x₀ y₀T hQD = 0 := by rw [← hDm]; exact hdiff
    have hcomp := aut_comp_mulByN W τ hQD hx hy n hT μ hμF hns hμx hμy hμP
    have hτz : τ (μ r) = μ r :=
      congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun) hcomp) r
    exact count_translate_eq W h2 v v' h hv h' hv' τ hQD hx hy hsum hz hτz

/-! ## ★出典の紐付け(`.src`) -/

def valuation_mulByN_le_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——場合 A の判定が n·Q_v ≠ 0 であること)",
    sectionId := "genell-thm-3-8" }

def count_eq_of_sub_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——差が n 等分点なら指数が等しいこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
