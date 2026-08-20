import ABC3.Found.GaloisRep.PlaceTransport

/-!
# Galois (G5) 第 170 ブロック —— **★★★★★★★★`e_v` の一定性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★D2 の最後の道具がそろった

第 168(`τ ∘ μ = μ`)と第 169(素点の輸送)を繋ぐと

    count_v(z) = count_{v'}(z)      (τ z = z、Q_{v'} = Q_v + T)

が出る。★`z := μ f_P` に当てれば **`e_v` はファイバー上で一定**である。

### ★★★★★輸送の仮説を第 167 でそろえる

第 169 の `count_comp_eq` は次の 4 つを要求する:

| 仮説 | 出どころ |
|---|---|
| `w_v(τ(a)) ≤ 1` (∀a ∈ F[W]) | 第 167 `valuation_curveHom_le_one` |
| `w_v(τ(a)) < 1 ⟺ a ∈ v'` | 第 167 `pullbackPrime_curveHom` |
| 同じものを `τ.symm`・`v'` について | 下の `autFF_symm_coord` |

★`τ ∘ algebraMap` が `curveHom` であることは第 119 の `coordinateRing_hom_ext` から
(`aut_algebraMap_eq_curveHom`)。
★★`τ(coordX)` の付値が 1 以下であることは、第 163 の `redPoint_eq_zero_iff` と
`redPoint(生成点 + T) = Q_v + T` から出る(`Q_v + T ≠ 0` が効く)。

### ★★★★`τ.symm` は `−T` の平行移動である

`τ.symm ∘ τ = id` を点の側で読むと `τ.symm(生成点) = 生成点 − T`。★したがって
`τ.symm(coordX) = translateX W x₀ (negY x₀ y₀T)`(`autFF_symm_coord`)。
★★`exists_translateAut_all` を `−T` に当て直す必要はない。

## ★★★残るのは総和の組み立てだけ

    S := Σ_v e_v · Q_v = 0

★台はファイバー `[n]⁻¹(P)` と `E[n] ∖ {O}` の 2 つ。**係数が一致する必要はない**——
`Σ_{ファイバー} Q = n·P = 0`(第 150)と `Σ_{T ∈ E[n]∖{O}} T = 0`(第 150)が
それぞれ 0 だからである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `autFF_symm_autFF` / `autFF_toFF` | ★★自己同型の基本性質 |
| `autFF_symm_coord` | ★★★★★★**`τ.symm` は `−T` の平行移動** |
| `redPoint_toFF` | ★★定数点の還元はそれ自身 |
| `redPoint_aut_generic` | ★★★★★★**`red(τ(生成点)) = Q_v + T`** |
| `aut_algebraMap_eq_curveHom` | ★★★★★`τ ∘ algebraMap = curveHom` |
| `aut_transport_hyps` | ★★★★★★★輸送の仮説がそろう |
| `count_translate_eq` | ★★★★★★★★**`e_v` の一定性** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-! ## ★★自己同型の基本性質 -/

/-- ★★逆向きの自己同型は元に戻す。 -/
theorem autFF_symm_autFF (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    (S : (W.map (algebraMap F W.FunctionField)).Point) :
    autFF W τ.symm (autFF W τ S) = S := by
  match S with
  | 0 => rw [map_zero, map_zero]
  | Point.some x y h =>
    obtain ⟨h1, e1⟩ := autFF_some W τ h
    obtain ⟨h2, e2⟩ := autFF_some W τ.symm h1
    rw [e1, e2]
    exact point_some_congr (τ.symm_apply_apply x) (τ.symm_apply_apply y) h2 h

/-- ★★定数点は自己同型で固定される。 -/
theorem autFF_toFF (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    autFF W τ (toFF W (Point.some x₀ y₀ hQ)) = toFF W (Point.some x₀ y₀ hQ) := by
  rw [toFF_some]
  obtain ⟨h1, e1⟩ := autFF_some W τ (mapNonsingular W hQ)
  rw [e1]
  exact point_some_congr (τ.commutes x₀) (τ.commutes y₀) h1 (mapNonsingular W hQ)

/-- ★★★★★★**`τ.symm` は `−T` の平行移動である**。 -/
theorem autFF_symm_coord (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀)
    (hx : τ (coordX W) = translateX W x₀ y₀) (hy : τ (coordY W) = translateY W x₀ y₀) :
    τ.symm (coordX W) = translateX W x₀ (W.negY x₀ y₀)
      ∧ τ.symm (coordY W) = translateY W x₀ (W.negY x₀ y₀) := by
  have hQ' : W.Nonsingular x₀ (W.negY x₀ y₀) := (nonsingular_neg x₀ y₀).mpr hQ
  have h1 : autFF W τ.symm (autFF W τ (ABC3.Found.GaloisRep.genericPoint W))
      = ABC3.Found.GaloisRep.genericPoint W := autFF_symm_autFF W τ _
  rw [autFF_generic W τ hQ hx hy, map_add, autFF_toFF W τ.symm hQ] at h1
  have h2 : autFF W τ.symm (ABC3.Found.GaloisRep.genericPoint W)
      = ABC3.Found.GaloisRep.genericPoint W + toFF W (Point.some x₀ (W.negY x₀ y₀) hQ') := by
    have hneg : toFF W (Point.some x₀ (W.negY x₀ y₀) hQ') = -toFF W (Point.some x₀ y₀ hQ) := by
      rw [← map_neg, Point.neg_some]
    rw [hneg, ← sub_eq_add_neg, eq_sub_iff_add_eq]
    exact h1
  rw [generic_add_toFF W hQ'] at h2
  obtain ⟨h3, e3⟩ := autFF_some W τ.symm (nonsingular_coord W)
  rw [ABC3.Found.GaloisRep.genericPoint, e3] at h2
  exact ⟨by injection h2, by injection h2⟩

/-- ★★★★★**自己同型と座標環の合成は `curveHom` である**。 -/
theorem aut_algebraMap_eq_curveHom (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    (heq : (W.map (algebraMap F W.FunctionField)).Equation (τ (coordX W)) (τ (coordY W)))
    (a : W.CoordinateRing) :
    τ (algebraMap W.CoordinateRing W.FunctionField a) = curveHom W heq a := by
  have hext : (τ.toAlgHom.toRingHom).comp (algebraMap W.CoordinateRing W.FunctionField)
      = curveHom W heq := by
    refine coordinateRing_hom_ext _ _ (fun d => ?_) ?_ ?_
    · rw [RingHom.comp_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe,
        ← IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField d,
        curveHom_algebraMap]
      exact τ.commutes d
    · rw [RingHom.comp_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, curveHom_genX]
      rfl
    · rw [RingHom.comp_apply, AlgHom.toRingHom_eq_coe, RingHom.coe_coe, curveHom_genY]
      rfl
  exact congrFun (congrArg (fun f => (f : W.CoordinateRing →+* W.FunctionField).toFun) hext) a

/-! ## ★★★★★★輸送の仮説 -/

variable [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★定数点の還元はそれ自身。 -/
theorem redPoint_toFF (v : HeightOneSpectrum W.CoordinateRing) {c y₀ : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    {x₀ y₀' : F} (hQ : W.Nonsingular x₀ y₀') :
    redPoint W v h hv (toFF W (Point.some x₀ y₀' hQ)) = Point.some x₀ y₀' hQ := by
  rw [toFF_some, redPoint_some W v h hv (mapNonsingular W hQ)
    (valuation_algebraMap_field W v _) (valuation_algebraMap_field W v _)]
  exact point_some_congr (redConst_algebraMap W v h hv x₀)
    (redConst_algebraMap W v h hv y₀') _ hQ

/-- ★★★★★★**平行移動した生成点の還元は `Q_v + T`**。 -/
theorem redPoint_aut_generic (v : HeightOneSpectrum W.CoordinateRing)
    {c y₀ : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    (h2 : IsUnit (2 : F)) (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀T : F} (hQ : W.Nonsingular x₀ y₀T)
    (hx : τ (coordX W) = translateX W x₀ y₀T) (hy : τ (coordY W) = translateY W x₀ y₀T) :
    redPoint W v h hv (autFF W τ (ABC3.Found.GaloisRep.genericPoint W))
      = Point.some c y₀ (equation_iff_nonsingular.mp h) + Point.some x₀ y₀T hQ := by
  rw [autFF_generic W τ hQ hx hy, redPoint_add W v h hv h2,
    redPoint_generic W v h hv, redPoint_toFF W v h hv hQ]

/-- ★★★★★★★**輸送の仮説がそろう**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Q_v + T ≠ 0`(アフィン点)であることが `w_v(τ coordX) ≤ 1` を保証する。 -/
theorem aut_transport_hyps (v : HeightOneSpectrum W.CoordinateRing)
    {c y₀ : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    (h2 : IsUnit (2 : F)) (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀T : F} (hQ : W.Nonsingular x₀ y₀T)
    (hx : τ (coordX W) = translateX W x₀ y₀T) (hy : τ (coordY W) = translateY W x₀ y₀T)
    {c' y₀' : F} (hQ' : W.Nonsingular c' y₀')
    (hsum : Point.some c y₀ (equation_iff_nonsingular.mp h) + Point.some x₀ y₀T hQ
      = Point.some c' y₀' hQ') :
    (∀ a : W.CoordinateRing,
        v.valuation W.FunctionField (τ (algebraMap W.CoordinateRing W.FunctionField a)) ≤ 1)
      ∧ (∀ a : W.CoordinateRing,
        v.valuation W.FunctionField (τ (algebraMap W.CoordinateRing W.FunctionField a)) < 1
          ↔ a ∈ CoordinateRing.XYIdeal W c' (Polynomial.C y₀')) := by
  obtain ⟨hns1, hpt⟩ := autFF_some W τ (nonsingular_coord W)
  have hred : redPoint W v h hv (Point.some (τ (coordX W)) (τ (coordY W)) hns1)
      = Point.some c' y₀' hQ' := by
    rw [← hpt, ← ABC3.Found.GaloisRep.genericPoint,
      redPoint_aut_generic W v h hv h2 τ hQ hx hy, hsum]
  have hxle : v.valuation W.FunctionField (τ (coordX W)) ≤ 1 := by
    by_contra hcon
    rw [not_le] at hcon
    rw [(redPoint_eq_zero_iff W v h hv hns1).2 hcon] at hred
    exact absurd hred.symm (by simp)
  have hyle : v.valuation W.FunctionField (τ (coordY W)) ≤ 1 :=
    val_y_le_of_val_x_le W v hns1.1 hxle
  rw [redPoint_some W v h hv hns1 hxle hyle] at hred
  have hcx : redConst W v h hv (τ (coordX W)) = c' := by injection hred
  have hcy : redConst W v h hv (τ (coordY W)) = y₀' := by injection hred
  have heq1 := aut_algebraMap_eq_curveHom W τ hns1.1
  constructor
  · intro a
    rw [heq1 a]
    exact valuation_curveHom_le_one W v hns1.1 hxle hyle a
  · intro a
    rw [heq1 a]
    have hpb := pullbackPrime_curveHom W v h hv hns1.1 hxle hyle
    rw [hcx, hcy] at hpb
    rw [← hpb]
    exact Iff.rfl

/-! ## ★★★★★★★★`e_v` の一定性 -/

/-- ★★★★★★★★**平行移動で不変な関数の指数は、平行移動した素点でも同じ**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 168 の `τ ∘ μ = μ` を `z := μ f_P` に当てれば **`e_v` はファイバー上で一定**。 -/
theorem count_translate_eq (h2 : IsUnit (2 : F))
    (v v' : HeightOneSpectrum W.CoordinateRing)
    {c y₀ c' y₀' : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    (h' : W.Equation c' y₀')
    (hv' : v'.asIdeal = CoordinateRing.XYIdeal W c' (Polynomial.C y₀'))
    (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀T : F} (hQ : W.Nonsingular x₀ y₀T)
    (hx : τ (coordX W) = translateX W x₀ y₀T) (hy : τ (coordY W) = translateY W x₀ y₀T)
    (hsum : Point.some c y₀ (equation_iff_nonsingular.mp h) + Point.some x₀ y₀T hQ
      = Point.some c' y₀' (equation_iff_nonsingular.mp h'))
    {z : W.FunctionField} (hz : z ≠ 0) (hτz : τ z = z) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ z)
      = FractionalIdeal.count W.FunctionField v'
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ z) := by
  obtain ⟨hle, hlt⟩ := aut_transport_hyps W v h hv h2 τ hQ hx hy
    (equation_iff_nonsingular.mp h') hsum
  have hQ' : W.Nonsingular x₀ (W.negY x₀ y₀T) := (nonsingular_neg x₀ y₀T).mpr hQ
  obtain ⟨hx', hy'⟩ := autFF_symm_coord W τ hQ hx hy
  have hsum' : Point.some c' y₀' (equation_iff_nonsingular.mp h')
      + Point.some x₀ (W.negY x₀ y₀T) hQ'
      = Point.some c y₀ (equation_iff_nonsingular.mp h) := by
    rw [← hsum, ← Point.neg_some hQ, add_assoc, add_neg_cancel, add_zero]
  obtain ⟨hle', hlt'⟩ := aut_transport_hyps W v' h' hv' h2 τ.symm hQ' hx' hy'
    (equation_iff_nonsingular.mp h) hsum'
  rw [← hv'] at hlt
  rw [← hv] at hlt'
  have hmain := count_comp_eq τ.toRingEquiv v v' hle hlt hle' hlt' hz
  rw [show τ.toRingEquiv z = z from hτz] at hmain
  exact hmain

/-! ## ★出典の紐付け(`.src`) -/

def aut_transport_hyps.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——素点の輸送の仮説がそろうこと)",
    sectionId := "genell-thm-3-8" }

def count_translate_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——分岐指数がファイバー上で一定であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
