import ABC3.Found.GaloisRep.SemiLinear

/-!
# Galois (G5) 第 193 ブロック —— **★★★★★★★★★Weil 対の半線型同変性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★5 性質の 4 つ目

    σ(e_n(P, Q)) = e_n(σP, σQ)

★第 178 の `WeilSpec` は**データの存在**なので、witness を `σ` で輸送すれば良い。
★★`weilPairingVal` は `Classical.choose` で定めてあり、第 178 の `weilSpec_unique` が
一意性を保証するので、**述語の輸送だけで値の同変性が出る**。

### ★★★★★★★半線型な点の写像

`Point.map`(mathlib)は `f : F →ₐ[S] K` を要求する——**半線型写像は表せない**
(`σ` が底体 `L` に非自明に効くから)。★そこで加法則を直接輸送して `semiPoint` を作った:

| 補題 | 内容 |
|---|---|
| `equation_semi` / `nonsingular_semi` | 方程式・非特異性が保たれる |
| `negY_semi` / `addX_semi` / `addY_semi` / `slope_semi` | 加法公式が可換 |
| `semiPoint` | ★★★★★★★**加法準同型になる** |

★★`slope` の場合分け(`x₁ = x₂` と `y₁ = negY x₂ y₂`)は `σ` の単射性で対応が付く。

### ★★★★★★★★`Σ_F` は `n·生成点` を固定する

これが要である。★`Σ_F` は `coordX`・`coordY` と係数を固定するので、
`semiPoint` は生成点を固定し、加法準同型だから `n·生成点` も固定する。
★★したがって **`Σ_F xn = xn`・`Σ_F yn = yn`**、よって **`Σ_F ∘ μ = μ ∘ Σ_R`**
(第 119 の `coordinateRing_hom_ext` で一意性)。

### ★★★witness の各部品の行き先

| 部品 | 行き先 |
|---|---|
| `x, y, x₀, y₀` | `σx, σy, σx₀, σy₀` |
| `f_P` | `Σ_R f_P`(`XYIdeal` の輸送) |
| `μ`・`xn`・`yn`・`hns`・`hμP` | **そのまま**(`Σ_F` が固定するから) |
| `g` | `Σ_F g`(`(Σ_F g)^n = Σ_F(μ f_P) = μ(Σ_R f_P)`) |
| `τ` | `Σ_F ∘ τ ∘ Σ_F⁻¹`(`L`-代数同型になる) |
| `c` | `σ c` |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `semiPoint` | ★★★★★★★**半線型な点の写像** |
| `semiFF_fixes_mulByN` | ★★★★★★★★**`Σ_F` は `n·生成点` の座標を固定する** |
| `semiFF_comp_mu` | ★★★★★★★★**`Σ_F ∘ μ = μ ∘ Σ_R`** |
| `weilSpec_semi` | ★★★★★★★★★**`WeilSpec` の輸送** |
| `weilPairingVal_semi` | ★★★★★★★★★**Weil 対の半線型同変性** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial nonZeroDivisors

/-! ## ★★★★★★★加法公式の輸送 -/

section Formulas

variable {M : Type} [Field M] (V : WeierstrassCurve.Affine M) {Φ : M ≃+* M}

theorem equation_semi (h : FixesCoeffs V Φ) {x y : M} (he : V.Equation x y) :
    V.Equation (Φ x) (Φ y) := by
  rw [equation_iff] at he ⊢
  have hc := congrArg Φ he
  simp only [map_add, map_mul, map_pow] at hc
  rw [h.a₁, h.a₂, h.a₃, h.a₄, h.a₆] at hc
  exact hc

theorem nonsingular_semi (h : FixesCoeffs V Φ) {x y : M} (hns : V.Nonsingular x y) :
    V.Nonsingular (Φ x) (Φ y) := by
  rw [nonsingular_iff] at hns ⊢
  refine ⟨equation_semi V h hns.1, ?_⟩
  rcases hns.2 with hcc | hcc
  · left
    intro hcon
    apply hcc
    have hs := congrArg Φ.symm hcon
    simp only [map_add, map_mul, map_pow, RingEquiv.symm_apply_apply, map_ofNat] at hs
    rw [(h.symm).a₁, (h.symm).a₂, (h.symm).a₄] at hs
    exact hs
  · right
    intro hcon
    apply hcc
    have hs := congrArg Φ.symm hcon
    simp only [map_sub, map_neg, map_mul, RingEquiv.symm_apply_apply] at hs
    rw [(h.symm).a₁, (h.symm).a₃] at hs
    exact hs

theorem negY_semi (h : FixesCoeffs V Φ) (x y : M) :
    Φ (V.negY x y) = V.negY (Φ x) (Φ y) := by
  simp only [negY, map_sub, map_neg, map_mul, h.a₁, h.a₃]

theorem addX_semi (h : FixesCoeffs V Φ) (x₁ x₂ l : M) :
    Φ (V.addX x₁ x₂ l) = V.addX (Φ x₁) (Φ x₂) (Φ l) := by
  simp only [addX, map_sub, map_add, map_mul, map_pow, h.a₁, h.a₂]

theorem negAddY_semi (h : FixesCoeffs V Φ) (x₁ x₂ y₁ l : M) :
    Φ (V.negAddY x₁ x₂ y₁ l) = V.negAddY (Φ x₁) (Φ x₂) (Φ y₁) (Φ l) := by
  simp only [negAddY, map_sub, map_add, map_mul, addX_semi V h]

theorem addY_semi (h : FixesCoeffs V Φ) (x₁ x₂ y₁ l : M) :
    Φ (V.addY x₁ x₂ y₁ l) = V.addY (Φ x₁) (Φ x₂) (Φ y₁) (Φ l) := by
  rw [addY, addY, negY_semi V h, addX_semi V h, negAddY_semi V h]

theorem slope_semi [DecidableEq M] (h : FixesCoeffs V Φ) (x₁ x₂ y₁ y₂ : M) :
    Φ (V.slope x₁ x₂ y₁ y₂) = V.slope (Φ x₁) (Φ x₂) (Φ y₁) (Φ y₂) := by
  rw [slope, slope]
  by_cases hx : x₁ = x₂
  · by_cases hy : y₁ = V.negY x₂ y₂
    · rw [if_pos hx, if_pos hy, if_pos (congrArg Φ hx),
        if_pos (show Φ y₁ = V.negY (Φ x₂) (Φ y₂) by rw [hy, negY_semi V h]), map_zero]
    · rw [if_pos hx, if_neg hy, if_pos (congrArg Φ hx),
        if_neg (fun hcon => hy (Φ.injective (by rw [hcon, negY_semi V h])))]
      rw [map_div₀]
      simp only [map_sub, map_add, map_mul, map_pow, map_ofNat, h.a₁, h.a₂, h.a₄,
        negY_semi V h]
  · rw [if_neg hx, if_neg (fun hcon => hx (Φ.injective hcon)), map_div₀, map_sub, map_sub]

/-- ★★★★★★★**半線型な点の写像**——係数を固定する体自己同型が点に誘導する加法準同型。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `Point.map` は `F →ₐ[S] K` を要求するので半線型写像は表せない。 -/
noncomputable def semiPoint [DecidableEq M] (h : FixesCoeffs V Φ) : V.Point →+ V.Point where
  toFun P := match P with
    | 0 => 0
    | Point.some _ _ hns => Point.some _ _ (nonsingular_semi V h hns)
  map_zero' := rfl
  map_add' := by
    rintro (_ | ⟨x₁, y₁, h₁⟩) (_ | ⟨x₂, y₂, h₂⟩)
    any_goals rfl
    by_cases hxy : x₁ = x₂ ∧ y₁ = V.negY x₂ y₂
    · rw [Point.add_of_Y_eq hxy.1 hxy.2,
        Point.add_of_Y_eq (congrArg Φ hxy.1)
          (show Φ y₁ = V.negY (Φ x₂) (Φ y₂) by rw [hxy.2, negY_semi V h])]
    · rw [Point.add_some hxy,
        Point.add_some (fun hcon => hxy ⟨Φ.injective hcon.1,
          Φ.injective (by rw [negY_semi V h]; exact hcon.2)⟩)]
      show Point.some _ _ _ = Point.some _ _ _
      exact point_some_congr (by rw [addX_semi V h, slope_semi V h])
        (by rw [addY_semi V h, slope_semi V h]) _ _

theorem semiPoint_zero [DecidableEq M] (h : FixesCoeffs V Φ) : semiPoint V h 0 = 0 := rfl

theorem semiPoint_some [DecidableEq M] (h : FixesCoeffs V Φ) {x y : M} (hns : V.Nonsingular x y) :
    semiPoint V h (Point.some x y hns) = Point.some (Φ x) (Φ y) (nonsingular_semi V h hns) := rfl

theorem semiPoint_symm_semiPoint [DecidableEq M] (h : FixesCoeffs V Φ) (P : V.Point) :
    semiPoint V h.symm (semiPoint V h P) = P := by
  match P with
  | 0 => rw [map_zero, map_zero]
  | Point.some x y hns =>
      rw [semiPoint_some, semiPoint_some]
      exact point_some_congr (Φ.symm_apply_apply x) (Φ.symm_apply_apply y) _ _

end Formulas

/-! ## ★★★★★★★★`Σ_F` は `n·生成点` を固定する -/

section MulByN

variable {L : Type} [Field L] [DecidableEq L] (W : WeierstrassCurve.Affine L) [W.IsElliptic]

/-- `Σ_F` は基底変換した曲線の係数を固定する。 -/
theorem fixesCoeffs_FF {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    FixesCoeffs (W.map (algebraMap L W.FunctionField)) (semiFF W h) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · show semiFF W h (algebraMap L W.FunctionField W.a₁) = algebraMap L W.FunctionField W.a₁
    rw [semiFF_const, h.a₁]
  · show semiFF W h (algebraMap L W.FunctionField W.a₂) = algebraMap L W.FunctionField W.a₂
    rw [semiFF_const, h.a₂]
  · show semiFF W h (algebraMap L W.FunctionField W.a₃) = algebraMap L W.FunctionField W.a₃
    rw [semiFF_const, h.a₃]
  · show semiFF W h (algebraMap L W.FunctionField W.a₄) = algebraMap L W.FunctionField W.a₄
    rw [semiFF_const, h.a₄]
  · show semiFF W h (algebraMap L W.FunctionField W.a₆) = algebraMap L W.FunctionField W.a₆
    rw [semiFF_const, h.a₆]

/-- ★★★★★★★**`Σ_F` は生成点を固定する**。 -/
theorem semiPoint_generic {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    semiPoint (W.map (algebraMap L W.FunctionField)) (fixesCoeffs_FF W h)
      (ABC3.Found.GaloisRep.genericPoint W) = ABC3.Found.GaloisRep.genericPoint W := by
  rw [ABC3.Found.GaloisRep.genericPoint, semiPoint_some]
  exact point_some_congr (semiFF_coordX W h) (semiFF_coordY W h) _ _

/-- ★★★★★★★★**`Σ_F` は `n·生成点` の座標を固定する**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが半線型輸送の要である。 -/
theorem semiFF_fixes_mulByN {σ : L ≃+* L} (h : FixesCoeffs W σ) (n : ℕ)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap L W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns) :
    semiFF W h xn = xn ∧ semiFF W h yn = yn := by
  have hfix : semiPoint (W.map (algebraMap L W.FunctionField)) (fixesCoeffs_FF W h)
      (n • ABC3.Found.GaloisRep.genericPoint W) = n • ABC3.Found.GaloisRep.genericPoint W := by
    rw [map_nsmul, semiPoint_generic W h]
  rw [hμP, semiPoint_some] at hfix
  exact ⟨by injection hfix, by injection hfix⟩

/-- ★★★★★★★★**`Σ_F ∘ μ = μ ∘ Σ_R`**。 -/
theorem semiFF_comp_mu {σ : L ≃+* L} (h : FixesCoeffs W σ)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : L, μ (algebraMap L W.CoordinateRing d) = algebraMap L W.FunctionField d)
    (n : ℕ) {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap L W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn) (a : W.CoordinateRing) :
    semiFF W h (μ a) = μ (semiCoordEquiv W h a) := by
  obtain ⟨hxf, hyf⟩ := semiFF_fixes_mulByN W h n hns hμP
  have hring : ((semiFF W h : W.FunctionField →+* W.FunctionField).comp μ)
      = μ.comp (semiCoordEquiv W h : W.CoordinateRing →+* W.CoordinateRing) := by
    refine coordinateRing_hom_ext _ _ (fun d => ?_) ?_ ?_
    · show semiFF W h (μ (algebraMap L W.CoordinateRing d))
        = μ (semiCoordEquiv W h (algebraMap L W.CoordinateRing d))
      rw [hμF, semiFF_const, semiCoordEquiv_algebraMap, hμF]
    · show semiFF W h (μ (genX W)) = μ (semiCoordEquiv W h (genX W))
      rw [semiCoordEquiv_genX, hμx, hxf]
    · show semiFF W h (μ (genY W)) = μ (semiCoordEquiv W h (genY W))
      rw [semiCoordEquiv_genY, hμy, hyf]
  exact congrFun (congrArg (fun f : W.CoordinateRing →+* W.FunctionField => f.toFun) hring) a

theorem semiPoint_toFF {σ : L ≃+* L} (h : FixesCoeffs W σ) {x₀ y₀ : L}
    (hQ : W.Nonsingular x₀ y₀) :
    semiPoint (W.map (algebraMap L W.FunctionField)) (fixesCoeffs_FF W h)
        (toFF W (Point.some x₀ y₀ hQ))
      = toFF W (Point.some (σ x₀) (σ y₀) (nonsingular_semi W h hQ)) := by
  rw [toFF_some, semiPoint_some, toFF_some]
  exact point_some_congr (semiFF_const W h x₀) (semiFF_const W h y₀) _ _

/-- ★★★★★★★**`Σ_F` は平行移動座標を輸送する**。 -/
theorem semiFF_translate {σ : L ≃+* L} (h : FixesCoeffs W σ) {x₀ y₀ : L}
    (hQ : W.Nonsingular x₀ y₀) :
    semiFF W h (translateX W x₀ y₀) = translateX W (σ x₀) (σ y₀)
      ∧ semiFF W h (translateY W x₀ y₀) = translateY W (σ x₀) (σ y₀) := by
  have hbase := generic_add_toFF W hQ
  have himg := congrArg (semiPoint (W.map (algebraMap L W.FunctionField))
    (fixesCoeffs_FF W h)) hbase
  rw [map_add, semiPoint_generic W h, semiPoint_toFF W h hQ, semiPoint_some,
    generic_add_toFF W (nonsingular_semi W h hQ)] at himg
  exact ⟨by injection himg.symm, by injection himg.symm⟩

end MulByN

/-! ## ★イデアルと共役の輸送 -/

section Ideals

variable {L : Type} [Field L] (W : WeierstrassCurve.Affine L)

theorem semiCoord_XClass {σ : L ≃+* L} (h : FixesCoeffs W σ) (x : L) :
    semiCoordEquiv W h (CoordinateRing.XClass W x) = CoordinateRing.XClass W (σ x) := by
  show semiCoord W h (AdjoinRoot.of W.polynomial (Polynomial.X - Polynomial.C x))
    = AdjoinRoot.of W.polynomial (Polynomial.X - Polynomial.C (σ x))
  rw [semiCoord_of, Polynomial.map_sub, Polynomial.map_X, Polynomial.map_C]
  rfl

theorem semiCoord_YClass {σ : L ≃+* L} (h : FixesCoeffs W σ) (y : L) :
    semiCoordEquiv W h (CoordinateRing.YClass W (Polynomial.C y))
      = CoordinateRing.YClass W (Polynomial.C (σ y)) := by
  have hy : CoordinateRing.YClass W (Polynomial.C y)
      = genY W - algebraMap L W.CoordinateRing y := by
    show CoordinateRing.mk W (Polynomial.X - Polynomial.C (Polynomial.C y)) = _
    rw [map_sub]; rfl
  have hy' : CoordinateRing.YClass W (Polynomial.C (σ y))
      = genY W - algebraMap L W.CoordinateRing (σ y) := by
    show CoordinateRing.mk W (Polynomial.X - Polynomial.C (Polynomial.C (σ y))) = _
    rw [map_sub]; rfl
  rw [hy, hy', map_sub, semiCoordEquiv_genY, semiCoordEquiv_algebraMap]

theorem semiCoord_XYIdeal {σ : L ≃+* L} (h : FixesCoeffs W σ) (x y : L) :
    Ideal.map (semiCoordEquiv W h : W.CoordinateRing →+* W.CoordinateRing)
        (CoordinateRing.XYIdeal W x (Polynomial.C y))
      = CoordinateRing.XYIdeal W (σ x) (Polynomial.C (σ y)) := by
  rw [CoordinateRing.XYIdeal, CoordinateRing.XYIdeal, Ideal.map_span]
  congr 1
  rw [Set.image_insert_eq, Set.image_singleton]
  show {semiCoordEquiv W h (CoordinateRing.XClass W x),
    semiCoordEquiv W h (CoordinateRing.YClass W (Polynomial.C y))} = _
  rw [semiCoord_XClass W h, semiCoord_YClass W h]

theorem semiCoord_fP {σ : L ≃+* L} (h : FixesCoeffs W σ) {x y : L} (n : ℕ)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    (CoordinateRing.XYIdeal W (σ x) (Polynomial.C (σ y))) ^ n
      = Ideal.span {semiCoordEquiv W h fP} := by
  have hmap := congrArg (Ideal.map (semiCoordEquiv W h : W.CoordinateRing →+* W.CoordinateRing))
    hfP
  rw [Ideal.map_pow, semiCoord_XYIdeal W h, Ideal.map_span] at hmap
  rw [hmap, Set.image_singleton]
  rfl

theorem semiFF_symm_const {σ : L ≃+* L} (h : FixesCoeffs W σ) (r : L) :
    (semiFF W h).symm (algebraMap L W.FunctionField r)
      = algebraMap L W.FunctionField (σ.symm r) := by
  refine (semiFF W h).injective ?_
  rw [RingEquiv.apply_symm_apply, semiFF_const, σ.apply_symm_apply]

theorem semiFF_symm_coordX {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    (semiFF W h).symm (coordX W) = coordX W := by
  refine (semiFF W h).injective ?_
  rw [RingEquiv.apply_symm_apply, semiFF_coordX]

theorem semiFF_symm_coordY {σ : L ≃+* L} (h : FixesCoeffs W σ) :
    (semiFF W h).symm (coordY W) = coordY W := by
  refine (semiFF W h).injective ?_
  rw [RingEquiv.apply_symm_apply, semiFF_coordY]

/-- ★★★★★★`τ` の共役 `Σ_F ∘ τ ∘ Σ_F⁻¹`。 -/
noncomputable def semiConj {σ : L ≃+* L} (h : FixesCoeffs W σ)
    (τ : W.FunctionField ≃ₐ[L] W.FunctionField) : W.FunctionField ≃ₐ[L] W.FunctionField :=
  AlgEquiv.ofRingEquiv (f := (semiFF W h).symm.trans (τ.toRingEquiv.trans (semiFF W h)))
    (fun r => by
      show semiFF W h (τ ((semiFF W h).symm (algebraMap L W.FunctionField r)))
        = algebraMap L W.FunctionField r
      rw [semiFF_symm_const W h, τ.commutes, semiFF_const, σ.apply_symm_apply])

theorem semiConj_apply {σ : L ≃+* L} (h : FixesCoeffs W σ)
    (τ : W.FunctionField ≃ₐ[L] W.FunctionField) (z : W.FunctionField) :
    semiConj W h τ z = semiFF W h (τ ((semiFF W h).symm z)) := rfl

end Ideals

/-! ## ★★★★★★★★★`WeilSpec` の輸送と同変性 -/

section Transport

variable {L : Type} [Field L] [DecidableEq L] (W : WeierstrassCurve.Affine L) [W.IsElliptic]

/-- ★★★★★★★★★**`WeilSpec` は `σ` で輸送される**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`μ`・`xn`・`yn` はそのまま(`Σ_F` が固定する)。★★`f_P` は `Σ_R f_P`、`g` は `Σ_F g`、
`τ` は共役 `Σ_F ∘ τ ∘ Σ_F⁻¹`。 -/
theorem weilSpec_semi {σ : L ≃+* L} (h : FixesCoeffs W σ) (n : ℕ) {P Q : W.Point} {c : L}
    (hspec : WeilSpec W n P Q c) :
    WeilSpec W n (semiPoint W h P) (semiPoint W h Q) (σ c) := by
  obtain ⟨x, y, hP, x₀, y₀, hQ, fP, μ, g, τ, xn, yn, hns, hPeq, hQeq, hfP, hinj, hμF,
    hμx, hμy, hμP, hg, hτx, hτy, hd⟩ := hspec
  refine ⟨σ x, σ y, nonsingular_semi W h hP, σ x₀, σ y₀, nonsingular_semi W h hQ,
    semiCoordEquiv W h fP, μ, semiFF W h g, semiConj W h τ, xn, yn, hns,
    by rw [hPeq, semiPoint_some], by rw [hQeq, semiPoint_some],
    semiCoord_fP W h n fP hfP, hinj, hμF, hμx, hμy, hμP, ?_, ?_, ?_, ?_⟩
  · rw [← map_pow, hg, semiFF_comp_mu W h μ hμF n hns hμP hμx hμy]
  · rw [semiConj_apply, semiFF_symm_coordX, hτx, (semiFF_translate W h hQ).1]
  · rw [semiConj_apply, semiFF_symm_coordY, hτy, (semiFF_translate W h hQ).2]
  · rw [semiConj_apply, RingEquiv.symm_apply_apply, ← map_div₀, hd, semiFF_const]

/-- ★★★★★★★★★**Weil 対の半線型同変性**——`σ(e_n(P,Q)) = e_n(σP, σQ)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★witness が無い側は両方 1 になる(`σ.symm` で輸送して矛盾を出す)。 -/
theorem weilPairingVal_semi [IsAlgClosed L] (h2 : IsUnit (2 : L)) {σ : L ≃+* L}
    (h : FixesCoeffs W σ) {n : ℕ} (hn : 1 ≤ n) (P Q : W.Point) :
    σ (weilPairingVal W n P Q)
      = weilPairingVal W n (semiPoint W h P) (semiPoint W h Q) := by
  by_cases hex : ∃ c, WeilSpec W n P Q c
  · have hs := weilPairingVal_spec W n P Q hex
    exact (weilPairingVal_eq W h2 hn (weilSpec_semi W h n hs)).symm
  · rw [weilPairingVal_of_not W n P Q hex, map_one]
    refine (weilPairingVal_of_not W n _ _ ?_).symm
    rintro ⟨c', hc'⟩
    refine hex ⟨σ.symm c', ?_⟩
    have ht := weilSpec_semi W h.symm n hc'
    rwa [semiPoint_symm_semiPoint W h, semiPoint_symm_semiPoint W h] at ht

end Transport

/-! ## ★出典の紐付け(`.src`) -/

def semiPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の Galois 同変性——半線型な点の写像)",
    sectionId := "genell-thm-3-8" }

def weilPairingVal_semi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の Galois 同変性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
