import ABC3.Found.GaloisRep.TranslateAut

/-!
# Galois (G5) 第 118 ブロック —— **★★★★★★点は環準同型を与える(一般形)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★平行移動で使った道具は**点一般**に効く

第 115・117 ブロックの `translateHom` は、実は次の一般形の特別な場合であった:

    関数体上の曲線の点 `(x, y)`  ⟹  環準同型 `F[W] →+* F(W)`,  `x ↦ x`, `y ↦ y`

★これを `pointHom` として切り出す。★★すると **`[n]` の引き戻しも同じ道で出る**
——`n • (生成点)` の座標に `pointHom` を当てるだけである。

## ★★★★葉 2 が構成できた

`Skeleton/GaloisRep/WeilFunctionField.lean` の葉 2(乗法 `[n]` の引き戻し)について、
**環準同型そのものは本ブロックで `Found` に入る**(`exists_mulByNHom`)。

★残るのは 2 つ:
- **生成点が捩れ点でないこと**(`n • genericPoint ≠ 0`)——`.needs` に記録
- 単射性——`pointHom_injective_of_transcendental` に帰着済み

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pointHom` | ★★★★★★**点 ⟹ 環準同型**(一般形) |
| `pointHom_genX` / `_genY` | ★★★生成元の像は点の座標 |
| `pointHom_injective_of_transcendental` | ★★★★★単射性は `x` の超越性に帰着 |
| `exists_mulByNHom` | ★★★★★★**`[n]` の引き戻し** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★★★★★点から環準同型へ -/

/-- ★★★★★★**関数体上の曲線の点は座標環からの環準同型を与える**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★方程式を満たすことがそのまま `AdjoinRoot.lift` の仮説になる。 -/
noncomputable def pointHom (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y) :
    W.CoordinateRing →+* W.FunctionField :=
  AdjoinRoot.lift (Polynomial.eval₂RingHom (algebraMap F W.FunctionField) x) y (by
    rw [Polynomial.eval₂_eval₂RingHom_apply, ← WeierstrassCurve.Affine.map_polynomial]
    exact h)

theorem pointHom_genX (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y) :
    pointHom W h (genX W) = x := by
  rw [pointHom, genX, CoordinateRing.mk, AdjoinRoot.lift_mk]
  simp

theorem pointHom_genY (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y) :
    pointHom W h (genY W) = y := by
  rw [pointHom, genY, CoordinateRing.mk, AdjoinRoot.lift_mk]
  simp

theorem pointHom_polynomial (W : WeierstrassCurve.Affine F) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Equation x y) (p : Polynomial F) :
    pointHom W h (algebraMap (Polynomial F) W.CoordinateRing p)
      = Polynomial.eval₂ (algebraMap F W.FunctionField) x p := by
  rw [algebraMap_polynomial_coordinateRing, pointHom, CoordinateRing.mk, AdjoinRoot.lift_mk,
    Polynomial.eval₂_C, Polynomial.coe_eval₂RingHom]

/-- ★★★★★**単射性は `x` 座標の超越性に帰着する**(一般形)。

★`F[W]` が `F[X]` 上整であることから `Ideal.eq_bot_of_comap_eq_bot` が効く。 -/
theorem pointHom_injective_of_transcendental (W : WeierstrassCurve.Affine F)
    {x y : W.FunctionField} (h : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (htr : ∀ p : Polynomial F, p ≠ 0 →
      Polynomial.eval₂ (algebraMap F W.FunctionField) x p ≠ 0) :
    Function.Injective (pointHom W h) := by
  rw [RingHom.injective_iff_ker_eq_bot]
  refine Ideal.eq_bot_of_comap_eq_bot (R := Polynomial F) ?_
  refine le_antisymm (fun p hp => ?_) bot_le
  rw [Ideal.mem_comap, RingHom.mem_ker, pointHom_polynomial] at hp
  by_contra hne
  exact htr p (by simpa using hne) hp

/-! ## ★★★★★★葉 2 —— 乗法 `[n]` の引き戻し -/

theorem exists_coords_of_ne_zero {W : WeierstrassCurve.Affine F} (P : W.Point) (hP : P ≠ 0) :
    ∃ (x y : F) (h : W.Nonsingular x y), P = Point.some x y h := by
  cases P with
  | zero => exact absurd rfl hP
  | some x y h => exact ⟨x, y, h, rfl⟩

/-- ★★★★★★**乗法 `[n]` の座標環からの引き戻し**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`n • (生成点)` の座標に `pointHom` を当てるだけである
——**平行移動と同じ道が `[n]` にもそのまま効いた**。
★★仮定 `n • genericPoint ≠ 0`(生成点が捩れ点でないこと)は
`Skeleton/GaloisRep/WeilFunctionField.lean` に葉として記録する。 -/
theorem exists_mulByNHom (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ)
    (hn : n • genericPoint W ≠ 0) :
    ∃ (x y : W.FunctionField) (h : (W.map (algebraMap F W.FunctionField)).Nonsingular x y),
      n • genericPoint W = Point.some x y h ∧
      ∃ μ : W.CoordinateRing →+* W.FunctionField,
        μ (genX W) = x ∧ μ (genY W) = y := by
  obtain ⟨x, y, h, hP⟩ := exists_coords_of_ne_zero (n • genericPoint W) hn
  exact ⟨x, y, h, hP, pointHom W h.1, pointHom_genX W h.1, pointHom_genY W h.1⟩

/-! ## ★出典の紐付け(`.src`) -/

def pointHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——曲線の点から座標環の環準同型を作る)",
    sectionId := "genell-thm-3-8" }

def exists_mulByNHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——乗法 [n] の座標環からの引き戻し)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
