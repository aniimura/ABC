import ABC3.Found.GaloisRep.FiberConst

/-!
# Galois (G5) 第 172 ブロック —— **★★★★★素点と点の 1 対 1 対応**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★総和を点の側で取るための足場

D2 の残りは

    Σ_v e_v · Q_v = 0

の組み立てだけである。★これを**点の側の有限和**に移すには、
素点 `v` と点 `Q_v` の間の 1 対 1 対応を関数として持つ必要がある:

* `pointOf v` —— 第 138 の `exists_point_of_heightOneSpectrum` からの選択
* `placeOf S hS` —— 第 138 の `pointSpectrum`(場合分けで定義)

★★往復は `xyIdeal_inj`(第 164)で閉じる:

    pointOf (placeOf S hS) = S,      placeOf (pointOf v) _ = v

## ★★★★類の言い換え

第 165 の `classGroup_primeUnit` を `pointOf` の形に直しておく:

    [v] = toClass (pointOf v)

★これで第 166 の `classGroup_rootUnit_eq_toClass` に `Q := pointOf` を渡せる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_point_asIdeal` | ★素点はつねに点の `XYIdeal` |
| `pointOf` / `pointOf_ne_zero` | ★★★素点に対応する点 |
| `placeOf` / `placeOf_some` | ★★★点に対応する素点 |
| `pointOf_placeOf` / `placeOf_pointOf` | ★★★★★**往復が恒等** |
| `classGroup_primeUnit_pointOf` | ★★★★★`[v] = toClass (pointOf v)` |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [IsAlgClosed F] [DecidableEq F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]
  [inst : IsDedekindDomain W.CoordinateRing]

/-! ## ★★★素点に対応する点 -/

/-- ★素点はつねに点の `XYIdeal` である(第 138 の言い換え)。 -/
theorem exists_point_asIdeal (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) :
    ∃ c y₀ : F, W.Equation c y₀
      ∧ v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀) := by
  obtain ⟨c, y₀, heq, hv⟩ := exists_point_of_heightOneSpectrum h2 W v
  exact ⟨c, y₀, heq, by rw [hv]; rfl⟩

/-- ★素点に対応する点の `x` 座標(選択)。 -/
noncomputable def pointC (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) : F :=
  Classical.choose (exists_point_asIdeal W h2 v)

/-- ★素点に対応する点の `y` 座標(選択)。 -/
noncomputable def pointY (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) : F :=
  Classical.choose (Classical.choose_spec (exists_point_asIdeal W h2 v))

theorem pointEq (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) :
    W.Equation (pointC W h2 v) (pointY W h2 v) :=
  (Classical.choose_spec (Classical.choose_spec (exists_point_asIdeal W h2 v))).1

theorem pointAsIdeal (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) :
    v.asIdeal = CoordinateRing.XYIdeal W (pointC W h2 v)
      (Polynomial.C (pointY W h2 v)) :=
  (Classical.choose_spec (Classical.choose_spec (exists_point_asIdeal W h2 v))).2

/-- ★★★**素点に対応する点**。 -/
noncomputable def pointOf (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) :
    W.Point :=
  Point.some (pointC W h2 v) (pointY W h2 v)
    (equation_iff_nonsingular.mp (pointEq W h2 v))

theorem pointOf_eq (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) :
    pointOf W h2 v = Point.some (pointC W h2 v) (pointY W h2 v)
      (equation_iff_nonsingular.mp (pointEq W h2 v)) := rfl

theorem pointOf_ne_zero (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) :
    pointOf W h2 v ≠ 0 := by
  rw [pointOf_eq]; simp

/-! ## ★★★点に対応する素点 -/

/-- ★★★**点に対応する素点**。 -/
noncomputable def placeOf (S : W.Point) (hS : S ≠ 0) : HeightOneSpectrum W.CoordinateRing :=
  match S, hS with
  | 0, hh => absurd rfl hh
  | Point.some _ _ hns, _ => pointSpectrum W hns.1

theorem placeOf_some {x y : F} (hns : W.Nonsingular x y) (hS : Point.some x y hns ≠ 0) :
    (placeOf W (Point.some x y hns) hS).asIdeal
      = CoordinateRing.XYIdeal W x (Polynomial.C y) := rfl

/-! ## ★★★★★往復が恒等 -/

/-- ★★★★★**`pointOf ∘ placeOf = id`**——`xyIdeal_inj`(第 164)から。 -/
theorem pointOf_placeOf (h2 : IsUnit (2 : F)) (S : W.Point) (hS : S ≠ 0) :
    pointOf W h2 (placeOf W S hS) = S := by
  match S, hS with
  | 0, hh => exact absurd rfl hh
  | Point.some x y hns, hh =>
    set v := placeOf W (Point.some x y hns) hh with hvdef
    have hid : v.asIdeal = CoordinateRing.XYIdeal W x (Polynomial.C y) := placeOf_some W hns hh
    have h1 : CoordinateRing.XYIdeal W (pointC W h2 v) (Polynomial.C (pointY W h2 v))
        = CoordinateRing.XYIdeal W x (Polynomial.C y) := by
      rw [← pointAsIdeal W h2 v, hid]
    obtain ⟨hx, hy⟩ := xyIdeal_inj W hns.1 h1
    rw [pointOf_eq]
    exact point_some_congr hx hy _ hns

/-- ★★★★★**`placeOf ∘ pointOf = id`**。 -/
theorem placeOf_pointOf (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing) :
    placeOf W (pointOf W h2 v) (pointOf_ne_zero W h2 v) = v :=
  HeightOneSpectrum.ext (pointAsIdeal W h2 v).symm

/-- ★★★★★**素点の類は `pointOf` の類である**。 -/
theorem classGroup_primeUnit_pointOf (h2 : IsUnit (2 : F))
    (v : HeightOneSpectrum W.CoordinateRing) :
    ClassGroup.mk W.FunctionField (primeUnit W.FunctionField v)
      = Additive.toMul (Point.toClass (pointOf W h2 v)) :=
  classGroup_primeUnit W v (equation_iff_nonsingular.mp (pointEq W h2 v)) (pointAsIdeal W h2 v)

/-! ## ★出典の紐付け(`.src`) -/

def pointOf_placeOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——素点と点の 1 対 1 対応)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
