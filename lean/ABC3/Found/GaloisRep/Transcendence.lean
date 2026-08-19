import ABC3.Found.GaloisRep.Translate

/-!
# Galois (G5) 第 116 ブロック —— **★★★★★単射性を超越性に帰着する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★葉 1 に残っていたもの

第 115 ブロックで環準同型 `τ_Q : F[W] →+* F(W)` が構成できた。
★分数体へ延ばすには**単射性**が要る。それを次のように削る:

    `F[W]` は `F[X]` 上**整**(階数 2 の自由加群)
      ⟹ `Ideal.eq_bot_of_comap_eq_bot` により
        「核 ∩ `F[X]` = 0」から「核 = 0」が出る
      ⟹ 残るのは **`translateX` が `F` 上超越的**であることだけ

★★これで葉が一段細くなった。

## ★★★★足場(2026-08-20 実測)

| mathlib | 役割 |
|---|---|
| `Ideal.eq_bot_of_comap_eq_bot` | 整拡大で (0) の上のイデアルは (0) |
| `Module.Finite.of_basis` + `CoordinateRing.basis` | `F[W]` が `F[X]` 上階数 2 の自由加群 |
| `Algebra.IsIntegral.of_finite` | 有限 ⟹ 整 |
| `CoordinateRing.smul_basis_eq_zero` | 生成点の `x` の超越性 |
| `Polynomial.hom_eval₂` | 評価と環準同型の交換 |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `genX_transcendental` | ★★★★**生成点の `x` は `F` 上超越的** |
| `coordX_transcendental` | ★★★★関数体でも同じ |
| `coordinateRing_isIntegral` | ★★★`F[W]` は `F[X]` 上整 |
| `translateHom_injective_of_transcendental` | ★★★★★**単射性は超越性に帰着する** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★★★生成点の超越性 -/

/-- ★`F[X]` の元を生成点で評価すると `mk (C p)` になる。 -/
theorem eval₂_genX (W : WeierstrassCurve.Affine F) (p : Polynomial F) :
    Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W) p
      = CoordinateRing.mk W (Polynomial.C p) := by
  have hext : Polynomial.eval₂RingHom (algebraMap F W.CoordinateRing) (genX W)
      = (CoordinateRing.mk W).comp (Polynomial.C (R := Polynomial F)) := by
    refine Polynomial.ringHom_ext ?_ ?_
    · intro a
      simp only [Polynomial.coe_eval₂RingHom, Polynomial.eval₂_C, RingHom.coe_comp,
        Function.comp_apply]
      rfl
    · simp only [Polynomial.coe_eval₂RingHom, Polynomial.eval₂_X, RingHom.coe_comp,
        Function.comp_apply]
      rfl
  exact congrArg (fun h : Polynomial F →+* W.CoordinateRing => h p) hext

/-- ★★★★**生成点の `x` 座標は `F` 上超越的である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`F[W]` が `F[X]` 上階数 2 の自由加群であること(`smul_basis_eq_zero`)から出る。 -/
theorem genX_transcendental (W : WeierstrassCurve.Affine F) {p : Polynomial F} (hp : p ≠ 0) :
    Polynomial.eval₂ (algebraMap F W.CoordinateRing) (genX W) p ≠ 0 := by
  intro h
  rw [eval₂_genX] at h
  have h0 : p • (1 : W.CoordinateRing)
      + (0 : Polynomial F) • CoordinateRing.mk W Polynomial.X = 0 := by
    rw [zero_smul, add_zero, CoordinateRing.smul, mul_one]
    exact h
  exact hp (CoordinateRing.smul_basis_eq_zero h0).1

/-! ## ★★★整拡大 -/

/-- ★★`F[W]` は `F[X]` 上有限生成加群である(階数 2 の自由加群)。 -/
instance coordinateRing_finite (W : WeierstrassCurve.Affine F) :
    Module.Finite (Polynomial F) W.CoordinateRing :=
  Module.Finite.of_basis (CoordinateRing.basis W)

/-- ★★★`F[W]` は `F[X]` 上整である。 -/
instance coordinateRing_isIntegral (W : WeierstrassCurve.Affine F) :
    Algebra.IsIntegral (Polynomial F) W.CoordinateRing :=
  Algebra.IsIntegral.of_finite _ _

/-- ★★★★**関数体の座標関数 `x` も `F` 上超越的である**。 -/
theorem coordX_transcendental (W : WeierstrassCurve.Affine F) {p : Polynomial F} (hp : p ≠ 0) :
    Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) p ≠ 0 := by
  intro h
  refine genX_transcendental W hp ?_
  have hinj : Function.Injective (algebraMap W.CoordinateRing W.FunctionField) :=
    IsFractionRing.injective W.CoordinateRing W.FunctionField
  refine hinj ?_
  rw [map_zero, ← h, Polynomial.hom_eval₂, ← IsScalarTower.algebraMap_eq]
  rfl

/-! ## ★★★★★単射性の帰着 -/

theorem algebraMap_polynomial_coordinateRing (W : WeierstrassCurve.Affine F) (p : Polynomial F) :
    algebraMap (Polynomial F) W.CoordinateRing p = CoordinateRing.mk W (Polynomial.C p) := rfl

theorem translateHom_polynomial (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) (p : Polynomial F) :
    translateHom W hQ (algebraMap (Polynomial F) W.CoordinateRing p)
      = Polynomial.eval₂ (algebraMap F W.FunctionField) (translateX W x₀ y₀) p := by
  rw [algebraMap_polynomial_coordinateRing, translateHom, CoordinateRing.mk, AdjoinRoot.lift_mk,
    Polynomial.eval₂_C, Polynomial.coe_eval₂RingHom]

/-- ★★★★★**平行移動の環準同型の単射性は `translateX` の超越性に帰着する**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`F[W]` が `F[X]` 上整なので、核 ∩ `F[X]` = 0 から核 = 0 が出る
(`Ideal.eq_bot_of_comap_eq_bot`)。 -/
theorem translateHom_injective_of_transcendental (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀)
    (htr : ∀ p : Polynomial F, p ≠ 0 →
      Polynomial.eval₂ (algebraMap F W.FunctionField) (translateX W x₀ y₀) p ≠ 0) :
    Function.Injective (translateHom W hQ) := by
  rw [RingHom.injective_iff_ker_eq_bot]
  refine Ideal.eq_bot_of_comap_eq_bot (R := Polynomial F) ?_
  refine le_antisymm (fun p hp => ?_) bot_le
  rw [Ideal.mem_comap, RingHom.mem_ker, translateHom_polynomial] at hp
  by_contra hne
  exact htr p (by simpa using hne) hp

/-! ## ★出典の紐付け(`.src`) -/

def translateHom_injective_of_transcendental.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動の単射性を超越性に帰着する)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
