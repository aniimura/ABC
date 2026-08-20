import ABC3.Found.GaloisRep.WeilCharZero

/-!
# Galois (G5) 第 192 ブロック —— **★★★★★★係数を固定する自己同型の半線型輸送**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★Galois 同変性の土台

    σ(e_n(P, Q)) = e_n(σP, σQ)

★第 178 の `WeilSpec` は**データの存在**として書いてある。★★そこで `σ` で witness を
輸送すれば良い——`WeilSpec W n P Q c` の witness に `σ` を当てると
`WeilSpec W n (σP) (σQ) (σc)` の witness になる。

### ★★★★★まず座標環と関数体に `σ` を持ち上げる

`σ : L ≃+* L` が Weierstrass 係数 `a₁ … a₆` を固定するとする(`FixesCoeffs`)。
★このとき `σ` は `W.polynomial` を固定するので、`AdjoinRoot.lift` で

    Σ : W.CoordinateRing →+* W.CoordinateRing,   Σ(x) = x,  Σ(y) = y,  Σ(c) = σ(c)

が作れる。★★これは **`L`-線型ではない**(定数に `σ` が効く)——半線型である。
★★★`σ.symm` で作った側と合成すると恒等になる(第 119 の `coordinateRing_hom_ext` で
一意性を言う)ので、環同型 `semiCoordEquiv` になる。

### ★★★関数体へは局所化で延びる

`W.FunctionField = Frac(W.CoordinateRing)` なので
`IsLocalization.ringEquivOfRingEquiv` がそのまま使える。
★`Submonoid.map` の条件は `MulEquivClass.map_nonZeroDivisors` で埋まる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `FixesCoeffs` | `σ` が `a₁ … a₆` を固定すること |
| `map_polynomial_of_fixes` | `σ` は `W.polynomial` を固定する |
| `semiCoordEquiv` | ★★★★★★**座標環の半線型自己同型** |
| `semiFF` | ★★★★★★**関数体の半線型自己同型** |
| `semiFF_const` / `semiFF_coordX` / `semiFF_coordY` | 値の計算 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial nonZeroDivisors

variable {L : Type} [Field L] (W : WeierstrassCurve.Affine L)

/-! ## ★係数を固定するということ -/

/-- 体自己同型が Weierstrass 係数を固定すること。 -/
structure FixesCoeffs (σ : L ≃+* L) : Prop where
  a₁ : σ W.a₁ = W.a₁
  a₂ : σ W.a₂ = W.a₂
  a₃ : σ W.a₃ = W.a₃
  a₄ : σ W.a₄ = W.a₄
  a₆ : σ W.a₆ = W.a₆

theorem FixesCoeffs.symm {σ : L ≃+* L} (h : FixesCoeffs W σ) : FixesCoeffs W σ.symm :=
  ⟨σ.symm_apply_eq.2 h.a₁.symm, σ.symm_apply_eq.2 h.a₂.symm, σ.symm_apply_eq.2 h.a₃.symm,
   σ.symm_apply_eq.2 h.a₄.symm, σ.symm_apply_eq.2 h.a₆.symm⟩

/-- 係数を固定する体自己同型は Weierstrass 多項式を固定する。 -/
theorem map_polynomial_of_fixes {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    W.polynomial.map (Polynomial.mapRingHom (σ : L →+* L)) = W.polynomial := by
  simp only [WeierstrassCurve.Affine.polynomial, Polynomial.map_add, Polynomial.map_sub,
    Polynomial.map_mul, Polynomial.map_pow, Polynomial.map_X, Polynomial.map_C,
    Polynomial.coe_mapRingHom, RingHom.coe_coe]
  rw [h.a₁, h.a₂, h.a₃, h.a₄, h.a₆]

/-! ## ★★★★★★座標環への持ち上げ -/

/-- `σ` が座標環に誘導する半線型環準同型。 -/
noncomputable def semiCoord {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    W.CoordinateRing →+* W.CoordinateRing :=
  AdjoinRoot.lift ((AdjoinRoot.of W.polynomial).comp (Polynomial.mapRingHom (σ : L →+* L)))
    (AdjoinRoot.root W.polynomial)
    (by rw [← Polynomial.eval₂_map, map_polynomial_of_fixes W h]; exact AdjoinRoot.eval₂_root _)

theorem semiCoord_of {σ : L ≃+* L} (h : FixesCoeffs W σ) (p : Polynomial L) :
    semiCoord W h (AdjoinRoot.of W.polynomial p)
      = AdjoinRoot.of W.polynomial (p.map (σ : L →+* L)) := AdjoinRoot.lift_of _

theorem semiCoord_genY {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    semiCoord W h (genY W) = genY W := AdjoinRoot.lift_root _

theorem semiCoord_genX {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    semiCoord W h (genX W) = genX W := by
  show semiCoord W h (AdjoinRoot.of W.polynomial Polynomial.X)
    = AdjoinRoot.of W.polynomial Polynomial.X
  rw [semiCoord_of, Polynomial.map_X]

theorem semiCoord_algebraMap {σ : L ≃+* L} (h : FixesCoeffs W σ) (c : L) :
    semiCoord W h (algebraMap L W.CoordinateRing c) = algebraMap L W.CoordinateRing (σ c) := by
  show semiCoord W h (AdjoinRoot.of W.polynomial (Polynomial.C c))
    = AdjoinRoot.of W.polynomial (Polynomial.C (σ c))
  rw [semiCoord_of, Polynomial.map_C]
  rfl

/-- ★★★★★★**座標環の半線型自己同型**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`σ.symm` 側と合成すると恒等になる(第 119 の `coordinateRing_hom_ext` で一意性)。 -/
noncomputable def semiCoordEquiv {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    W.CoordinateRing ≃+* W.CoordinateRing := by
  refine RingEquiv.ofRingHom (semiCoord W h) (semiCoord W h.symm) ?_ ?_
  · refine coordinateRing_hom_ext _ _ (fun a => ?_) ?_ ?_
    · simp only [RingHom.comp_apply, semiCoord_algebraMap, RingHom.id_apply]
      rw [σ.apply_symm_apply]
    · simp only [RingHom.comp_apply, semiCoord_genX, RingHom.id_apply]
    · simp only [RingHom.comp_apply, semiCoord_genY, RingHom.id_apply]
  · refine coordinateRing_hom_ext _ _ (fun a => ?_) ?_ ?_
    · simp only [RingHom.comp_apply, semiCoord_algebraMap, RingHom.id_apply]
      rw [σ.symm_apply_apply]
    · simp only [RingHom.comp_apply, semiCoord_genX, RingHom.id_apply]
    · simp only [RingHom.comp_apply, semiCoord_genY, RingHom.id_apply]

theorem semiCoordEquiv_apply {σ : L ≃+* L} (h : FixesCoeffs W σ) (a : W.CoordinateRing) :
    semiCoordEquiv W h a = semiCoord W h a := rfl

theorem semiCoordEquiv_genX {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    semiCoordEquiv W h (genX W) = genX W := semiCoord_genX W h

theorem semiCoordEquiv_genY {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    semiCoordEquiv W h (genY W) = genY W := semiCoord_genY W h

theorem semiCoordEquiv_algebraMap {σ : L ≃+* L} (h : FixesCoeffs W σ) (c : L) :
    semiCoordEquiv W h (algebraMap L W.CoordinateRing c)
      = algebraMap L W.CoordinateRing (σ c) := semiCoord_algebraMap W h c

/-! ## ★★★★★★関数体への持ち上げ -/

/-- ★★★★★★**関数体の半線型自己同型**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`W.FunctionField = Frac(W.CoordinateRing)` なので局所化の同型でそのまま延びる。 -/
noncomputable def semiFF {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    W.FunctionField ≃+* W.FunctionField :=
  IsLocalization.ringEquivOfRingEquiv (M := W.CoordinateRing⁰) W.FunctionField W.FunctionField
    (semiCoordEquiv W h) (MulEquivClass.map_nonZeroDivisors (semiCoordEquiv W h))

theorem semiFF_algebraMap {σ : L ≃+* L} (h : FixesCoeffs W σ) (a : W.CoordinateRing) :
    semiFF W h (algebraMap W.CoordinateRing W.FunctionField a)
      = algebraMap W.CoordinateRing W.FunctionField (semiCoordEquiv W h a) :=
  IsLocalization.ringEquivOfRingEquiv_eq _ _

theorem semiFF_const {σ : L ≃+* L} (h : FixesCoeffs W σ) (c : L) :
    semiFF W h (algebraMap L W.FunctionField c) = algebraMap L W.FunctionField (σ c) := by
  rw [IsScalarTower.algebraMap_apply L W.CoordinateRing W.FunctionField c, semiFF_algebraMap,
    semiCoordEquiv_algebraMap, ← IsScalarTower.algebraMap_apply]

theorem semiFF_coordX {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    semiFF W h (coordX W) = coordX W := by
  show semiFF W h (algebraMap W.CoordinateRing W.FunctionField (genX W))
    = algebraMap W.CoordinateRing W.FunctionField (genX W)
  rw [semiFF_algebraMap, semiCoordEquiv_genX]

theorem semiFF_coordY {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    semiFF W h (coordY W) = coordY W := by
  show semiFF W h (algebraMap W.CoordinateRing W.FunctionField (genY W))
    = algebraMap W.CoordinateRing W.FunctionField (genY W)
  rw [semiFF_algebraMap, semiCoordEquiv_genY]

/-! ## ★出典の紐付け(`.src`) -/

def semiCoordEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の Galois 同変性——座標環への持ち上げ)",
    sectionId := "genell-thm-3-8" }

def semiFF.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の Galois 同変性——関数体への持ち上げ)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
