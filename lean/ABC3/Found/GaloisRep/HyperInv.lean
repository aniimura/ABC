import ABC3.Found.GaloisRep.ValuationBasis

/-!
# Galois (G5) 第 147 ブロック —— **★★★★★★超楕円対合**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`f_P` と `f_{−P}` を結ぶ

D1' 場合 B の残りは `count(μ f_P) = count(μ f_{−P})` である。★これは

    ι(f_P) と f_{−P} は同伴

から出る。★★そのために**超楕円対合** `ι : F[W] → F[W]`(`y ↦ negY(x,y)`)を作る。

## ★★★★構成——`pointHom` の一般化

第 118 の `pointHom` は終域が `F(W)` に固定されていた。
★**任意の `F`-代数へ**一般化する(`curveHom`)。すると

    ι := curveHom W (x := x, y := negY(x,y) の生成点版)

が `F[W]` 自身への写像として得られる。仮説は mathlib の
`equation_neg`(`Equation x (negY x y) ↔ Equation x y`)と
第 114 の `equation_gen` で埋まる。

## ★★★★★イデアルの像

    ι(I_P) = I_{−P}

★機構: `ι(y) − y` と `y − negY(x,y)` は `a₁(x − x_P)` の差しかないので、
2 元生成のイデアルとしては同じものになる。

★★これを `n` 乗すると `ι((f_P)) = (f_{−P})`、すなわち **`ι(f_P)` と `f_{−P}` は同伴**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `curveHom` | ★★★★点から環準同型へ(任意の `F`-代数を終域に) |
| `hyperInv` | ★★★★★★**超楕円対合 `ι`** |
| `hyperInv_genX` / `hyperInv_genY_eq` | ★`ι(x) = x`、`ι(y) = −y − a₁x − a₃` |
| `hyperInv_algebraMap_poly` | ★★`ι` は `F[x]` を固定する |
| `hyperInv_xyIdeal` | ★★★★★**`ι(I_P) = I_{−P}`** |
| `hyperInv_fP_assoc` | ★★★★★**`ι(f_P)` と `f_{−P}` は同伴** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★★★点から環準同型へ(一般形) -/

/-- ★★★★**曲線の点は座標環からの環準同型を与える**(任意の `F`-代数を終域に)。

★第 118 の `pointHom` を終域について一般化したもの。 -/
noncomputable def curveHom (W : WeierstrassCurve.Affine F) {S : Type} [CommRing S] [Algebra F S]
    {x y : S} (h : (W.map (algebraMap F S)).Equation x y) : W.CoordinateRing →+* S :=
  AdjoinRoot.lift (Polynomial.eval₂RingHom (algebraMap F S) x) y (by
    rw [Polynomial.eval₂_eval₂RingHom_apply, ← WeierstrassCurve.Affine.map_polynomial]
    exact h)

theorem curveHom_genX (W : WeierstrassCurve.Affine F) {S : Type} [CommRing S] [Algebra F S]
    {x y : S} (h : (W.map (algebraMap F S)).Equation x y) : curveHom W h (genX W) = x := by
  rw [curveHom, genX, CoordinateRing.mk, AdjoinRoot.lift_mk]; simp

theorem curveHom_genY (W : WeierstrassCurve.Affine F) {S : Type} [CommRing S] [Algebra F S]
    {x y : S} (h : (W.map (algebraMap F S)).Equation x y) : curveHom W h (genY W) = y := by
  rw [curveHom, genY, CoordinateRing.mk, AdjoinRoot.lift_mk]; simp

theorem curveHom_algebraMap (W : WeierstrassCurve.Affine F) {S : Type} [CommRing S] [Algebra F S]
    {x y : S} (h : (W.map (algebraMap F S)).Equation x y) (c : F) :
    curveHom W h (algebraMap F W.CoordinateRing c) = algebraMap F S c := by
  have hc : algebraMap F W.CoordinateRing c
      = algebraMap (Polynomial F) W.CoordinateRing (Polynomial.C c) := rfl
  rw [hc, algebraMap_polynomial_coordinateRing, curveHom, CoordinateRing.mk, AdjoinRoot.lift_mk,
    Polynomial.eval₂_C, Polynomial.coe_eval₂RingHom, Polynomial.eval₂_C]

/-! ## ★★★★★★超楕円対合 -/

/-- ★★★★★★**超楕円対合**——`x ↦ x`、`y ↦ negY(x,y)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `equation_neg` と第 114 の `equation_gen` で仮説が埋まる。 -/
noncomputable def hyperInv (W : WeierstrassCurve.Affine F) :
    W.CoordinateRing →+* W.CoordinateRing :=
  curveHom W (S := W.CoordinateRing)
    ((equation_neg (W' := W.map (algebraMap F W.CoordinateRing)) (genX W) (genY W)).2
      (equation_gen W))

theorem hyperInv_genX (W : WeierstrassCurve.Affine F) : hyperInv W (genX W) = genX W :=
  curveHom_genX W _

theorem hyperInv_genY (W : WeierstrassCurve.Affine F) :
    hyperInv W (genY W)
      = (W.map (algebraMap F W.CoordinateRing)).negY (genX W) (genY W) :=
  curveHom_genY W _

theorem hyperInv_algebraMap (W : WeierstrassCurve.Affine F) (c : F) :
    hyperInv W (algebraMap F W.CoordinateRing c) = algebraMap F W.CoordinateRing c :=
  curveHom_algebraMap W _ c

/-- ★`ι(y) = −y − a₁x − a₃`。 -/
theorem hyperInv_genY_eq (W : WeierstrassCurve.Affine F) :
    hyperInv W (genY W) = -genY W - algebraMap F W.CoordinateRing W.a₁ * genX W
      - algebraMap F W.CoordinateRing W.a₃ := by
  rw [hyperInv_genY, WeierstrassCurve.Affine.negY, WeierstrassCurve.map_a₁,
    WeierstrassCurve.map_a₃]

/-- ★★`ι` は `F[x]` を固定する。 -/
theorem hyperInv_algebraMap_poly (W : WeierstrassCurve.Affine F) (p : Polynomial F) :
    hyperInv W (algebraMap (Polynomial F) W.CoordinateRing p)
      = algebraMap (Polynomial F) W.CoordinateRing p := by
  rw [← eval₂_genX_eq_algebraMap, Polynomial.hom_eval₂, hyperInv_genX]
  congr 1
  exact RingHom.ext (fun c => hyperInv_algebraMap W c)

/-! ## ★★★★★イデアルの像 -/

/-- ★★★★★**`ι(I_P) = I_{−P}`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ι(y) − y` と `y − negY(x,y)` は `a₁(x − x_P)` の差しかない。 -/
theorem hyperInv_xyIdeal (W : WeierstrassCurve.Affine F) (x y : F) :
    Ideal.map (hyperInv W) (CoordinateRing.XYIdeal W x (Polynomial.C y))
      = CoordinateRing.XYIdeal W x (Polynomial.C (W.negY x y)) := by
  have hneg : W.negY x y = -y - W.a₁ * x - W.a₃ := rfl
  have hkey : hyperInv W (genY W) - algebraMap F W.CoordinateRing y
      = -(algebraMap F W.CoordinateRing W.a₁)
          * (genX W - algebraMap F W.CoordinateRing x)
        - (genY W - algebraMap F W.CoordinateRing (W.negY x y)) := by
    rw [hyperInv_genY_eq, hneg, map_sub, map_sub, map_neg, map_mul]
    ring
  rw [xyIdeal_eq_span, xyIdeal_eq_span, Ideal.map_span]
  have himg : (hyperInv W) '' ({genX W - algebraMap F W.CoordinateRing x,
        genY W - algebraMap F W.CoordinateRing y} : Set W.CoordinateRing)
      = {genX W - algebraMap F W.CoordinateRing x,
        hyperInv W (genY W) - algebraMap F W.CoordinateRing y} := by
    rw [Set.image_insert_eq, Set.image_singleton, map_sub, hyperInv_genX, hyperInv_algebraMap,
      map_sub, hyperInv_algebraMap]
  rw [himg]
  refine le_antisymm (Ideal.span_le.2 ?_) (Ideal.span_le.2 ?_)
  · rintro z (rfl | rfl)
    · exact Ideal.subset_span (Set.mem_insert _ _)
    · rw [hkey]
      exact Ideal.sub_mem _
        (Ideal.mul_mem_left _ _ (Ideal.subset_span (Set.mem_insert _ _)))
        (Ideal.subset_span (Set.mem_insert_of_mem _ rfl))
  · rintro z (rfl | rfl)
    · exact Ideal.subset_span (Set.mem_insert _ _)
    · have hz : genY W - algebraMap F W.CoordinateRing (W.negY x y)
          = -(algebraMap F W.CoordinateRing W.a₁)
              * (genX W - algebraMap F W.CoordinateRing x)
            - (hyperInv W (genY W) - algebraMap F W.CoordinateRing y) := by
        rw [hkey]; ring
      rw [hz]
      exact Ideal.sub_mem _
        (Ideal.mul_mem_left _ _ (Ideal.subset_span (Set.mem_insert _ _)))
        (Ideal.subset_span (Set.mem_insert_of_mem _ rfl))

/-- ★★★★★**`ι(f_P)` と `f_{−P}` は同伴である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ι(I_P) = I_{−P}` を `n` 乗するだけ。 -/
theorem hyperInv_fP_assoc (W : WeierstrassCurve.Affine F) {x y : F} (n : ℕ)
    (fP fN : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (hfN : (CoordinateRing.XYIdeal W x (Polynomial.C (W.negY x y))) ^ n = Ideal.span {fN}) :
    Associated (hyperInv W fP) fN := by
  have h1 : Ideal.map (hyperInv W) (Ideal.span ({fP} : Set W.CoordinateRing))
      = Ideal.span {hyperInv W fP} := by
    rw [Ideal.map_span, Set.image_singleton]
  rw [← hfP, Ideal.map_pow, hyperInv_xyIdeal, hfN] at h1
  exact Ideal.span_singleton_eq_span_singleton.1 h1.symm

/-! ## ★出典の紐付け(`.src`) -/

def hyperInv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——超楕円対合と ι(I_P) = I_{−P})",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
