import ABC3.Found.GaloisRep.TranslateProd

/-!
# Galois (G5) 第 189 ブロック —— **★★★★★★★★`τ_T ∘ [n]^* = [n]^* ∘ τ_{nT}`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★一般の交換則

第 168 の `aut_comp_mulByN` は `n·T = O` の場合、すなわち **`τ_T ∘ [n]^* = [n]^*`** だった。
★交代性では `n·P' = P ≠ O` の点で平行移動するので、**一般形**が要る:

    τ_T ∘ [n]^*  =  [n]^* ∘ τ_{nT}

### ★★★★★点の言葉に翻訳すると 1 行

両辺は `F(W)` からの環準同型で、`coordX`・`coordY`・定数での値で決まる(`functionField_hom_ext`)。
★そこで生成点の像を 2 通りに計算する:

| 側 | 計算 |
|---|---|
| `τ_T ∘ [n]^*` | `autFF τ_T (n·generic) = n·(generic + T) = n·generic + nT` |
| `[n]^* ∘ τ_{nT}` | `muFF (generic + S) = n·generic + S`(`S = nT`) |

★★どちらも **`n·generic + toFF S`**。★★★`Point.some` の単射性で座標が一致する。

### ★★★`IsTranslate` は生成点の像で特徴づけられる

    IsTranslate W τ T  ⟺  autFF τ (generic) = generic + toFF T

★これで**平行移動の合成が平行移動**であることが `Point.map` の関手性から出る
(`τ_{P'}` を `n` 回合成して `τ_{nP'}` にする段で要る)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `functionField_hom_ext` | ★★★★★関数体からの環準同型は座標と定数で決まる |
| `muExt_const` / `muExtAlg` / `muFF` | `μ̃` を `F`-代数射・点の写像として見る |
| `muFF_toFF` / `muFF_generic` | `μ̃` は定数点を固定し、生成点を `n·generic` に送る |
| `autFF_generic_gen` | ★★★★★★`autFF τ (generic) = generic + T`(`T = O` 込み) |
| `isTranslate_iff` | ★★★★★★**`IsTranslate` の特徴づけ** |
| `isTranslate_trans` | ★★★★★★**平行移動の合成は平行移動** |
| `aut_comp_muExt_gen` | ★★★★★★★★**一般の交換則** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★関数体からの環準同型の一意性 -/

section Ext

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-- ★★★★★**関数体からの環準同型は座標と定数で決まる**。 -/
theorem functionField_hom_ext {S : Type} [CommRing S] (φ ψ : W.FunctionField →+* S)
    (hF : ∀ c : F, φ (algebraMap F W.FunctionField c) = ψ (algebraMap F W.FunctionField c))
    (hx : φ (coordX W) = ψ (coordX W)) (hy : φ (coordY W) = ψ (coordY W)) : φ = ψ := by
  refine IsLocalization.ringHom_ext (M := W.CoordinateRing⁰) ?_
  refine coordinateRing_hom_ext _ _ (fun a => ?_) hx hy
  simp only [RingHom.comp_apply, ← IsScalarTower.algebraMap_apply]
  exact hF a

end Ext

/-! ## ★`μ̃` を点の写像として見る -/

section MuFF

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]
variable {μ : W.CoordinateRing →+* W.FunctionField}

/-- `μ̃` は定数を固定する。 -/
theorem muExt_const (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (c : F) : muExt W hinj (algebraMap F W.FunctionField c) = algebraMap F W.FunctionField c := by
  have h1 : (algebraMap F W.FunctionField c)
      = algebraMap W.CoordinateRing W.FunctionField (algebraMap F W.CoordinateRing c) :=
    IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField c
  rw [h1, muExt_algebraMap, hμF c]
  exact h1

/-- `μ̃` を `F`-代数射として見る。 -/
noncomputable def muExtAlg (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c) :
    W.FunctionField →ₐ[F] W.FunctionField :=
  AlgHom.mk' (muExt W hinj) (fun c x => by
    rw [Algebra.smul_def, Algebra.smul_def, map_mul, muExt_const W hinj hμF c])

theorem muExtAlg_apply (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (z : W.FunctionField) : muExtAlg W hinj hμF z = muExt W hinj z := rfl

/-- `μ̃` が誘導する点の写像。 -/
noncomputable def muFF (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c) :
    (W.map (algebraMap F W.FunctionField)).Point →+
      (W.map (algebraMap F W.FunctionField)).Point :=
  WeierstrassCurve.Affine.Point.map (W' := W) (muExtAlg W hinj hμF)

theorem muFF_some (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {x y : W.FunctionField} (h : (W.map (algebraMap F W.FunctionField)).Nonsingular x y) :
    ∃ h' : (W.map (algebraMap F W.FunctionField)).Nonsingular
        (muExt W hinj x) (muExt W hinj y),
      muFF W hinj hμF (Point.some x y h) = Point.some (muExt W hinj x) (muExt W hinj y) h' :=
  ⟨_, WeierstrassCurve.Affine.Point.map_some _ _⟩

theorem muFF_toFF (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (S : W.Point) : muFF W hinj hμF (toFF W S) = toFF W S := by
  match S with
  | 0 => rw [map_zero, map_zero]
  | Point.some s₀ t₀ hS =>
      obtain ⟨h', heq⟩ := muFF_some W hinj hμF (mapNonsingular W hS)
      rw [toFF_some, heq]
      exact point_some_congr (muExt_const W hinj hμF s₀) (muExt_const W hinj hμF t₀) _ _

theorem muExt_coordX (hinj : Function.Injective μ) :
    muExt W hinj (coordX W) = μ (genX W) := muExt_algebraMap W hinj (genX W)

theorem muExt_coordY (hinj : Function.Injective μ) :
    muExt W hinj (coordY W) = μ (genY W) := muExt_algebraMap W hinj (genY W)

theorem muFF_generic (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn) :
    muFF W hinj hμF (ABC3.Found.GaloisRep.genericPoint W) = Point.some xn yn hns := by
  obtain ⟨h', heq⟩ := muFF_some W hinj hμF (nonsingular_coord W)
  have hg : muFF W hinj hμF (ABC3.Found.GaloisRep.genericPoint W)
      = Point.some (muExt W hinj (coordX W)) (muExt W hinj (coordY W)) h' := heq
  rw [hg]
  exact point_some_congr ((muExt_coordX W hinj).trans hμx)
    ((muExt_coordY W hinj).trans hμy) _ _

end MuFF

/-! ## ★★★★★★`IsTranslate` の特徴づけ -/

section Characterize

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-- ★★★★★★**平行移動の自己同型は生成点を `生成点 + T` に送る**(`T = O` 込み)。 -/
theorem autFF_generic_gen (τ : W.FunctionField ≃ₐ[F] W.FunctionField) (T : W.Point)
    (hT : IsTranslate W τ T) :
    autFF W τ (ABC3.Found.GaloisRep.genericPoint W)
      = ABC3.Found.GaloisRep.genericPoint W + toFF W T := by
  match T, hT with
  | 0, hT =>
      rw [show τ = AlgEquiv.refl from hT, map_zero, add_zero]
      obtain ⟨h', heq⟩ := autFF_some W (AlgEquiv.refl (A₁ := W.FunctionField) (R := F))
        (nonsingular_coord W)
      rw [ABC3.Found.GaloisRep.genericPoint, heq]
      exact point_some_congr rfl rfl _ _
  | Point.some x₀ y₀ hQ, hT => exact autFF_generic W τ hQ hT.1 hT.2

theorem autFF_trans (τ τ' : W.FunctionField ≃ₐ[F] W.FunctionField)
    (S : (W.map (algebraMap F W.FunctionField)).Point) :
    autFF W (τ'.trans τ) S = autFF W τ (autFF W τ' S) := by
  match S with
  | 0 => rw [map_zero, map_zero, map_zero]
  | Point.some a b h =>
      obtain ⟨h1, he1⟩ := autFF_some W (τ'.trans τ) h
      obtain ⟨h2, he2⟩ := autFF_some W τ' h
      obtain ⟨h3, he3⟩ := autFF_some W τ h2
      rw [he1, he2, he3]
      exact point_some_congr rfl rfl _ _

theorem autFF_toFF_gen (τ : W.FunctionField ≃ₐ[F] W.FunctionField) (S : W.Point) :
    autFF W τ (toFF W S) = toFF W S := by
  match S with
  | 0 => rw [map_zero, map_zero]
  | Point.some a b h =>
      obtain ⟨h1, he1⟩ := autFF_some W τ (mapNonsingular W h)
      rw [toFF_some, he1]
      exact point_some_congr (τ.commutes a) (τ.commutes b) _ _

/-- ★★★★★★**`IsTranslate` は生成点の像で特徴づけられる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem isTranslate_iff (τ : W.FunctionField ≃ₐ[F] W.FunctionField) (T : W.Point) :
    IsTranslate W τ T ↔ autFF W τ (ABC3.Found.GaloisRep.genericPoint W)
      = ABC3.Found.GaloisRep.genericPoint W + toFF W T := by
  constructor
  · exact autFF_generic_gen W τ T
  · intro hgen
    obtain ⟨h', heq⟩ := autFF_some W τ (nonsingular_coord W)
    have hg : autFF W τ (ABC3.Found.GaloisRep.genericPoint W)
        = Point.some (τ (coordX W)) (τ (coordY W)) h' := heq
    rw [hg] at hgen
    match T with
    | 0 =>
        rw [map_zero, add_zero] at hgen
        have hx : τ (coordX W) = coordX W := by injection hgen
        have hy : τ (coordY W) = coordY W := by injection hgen
        show τ = AlgEquiv.refl
        exact aut_ext W hx hy
    | Point.some x₀ y₀ hQ =>
        rw [generic_add_toFF W hQ] at hgen
        have hx : τ (coordX W) = translateX W x₀ y₀ := by injection hgen
        have hy : τ (coordY W) = translateY W x₀ y₀ := by injection hgen
        exact ⟨hx, hy⟩

/-- ★★★★★★**平行移動の合成は平行移動**。 -/
theorem isTranslate_trans {τ τ' : W.FunctionField ≃ₐ[F] W.FunctionField} {T T' : W.Point}
    (hT : IsTranslate W τ T) (hT' : IsTranslate W τ' T') :
    IsTranslate W (τ'.trans τ) (T + T') := by
  rw [isTranslate_iff, autFF_trans, (isTranslate_iff W τ' T').1 hT', map_add,
    (isTranslate_iff W τ T).1 hT, autFF_toFF_gen, map_add]
  abel

end Characterize

/-! ## ★★★★★★★★一般の交換則 -/

section Comm

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]
variable {μ : W.CoordinateRing →+* W.FunctionField}

/-- ★★★★★★★★**`τ_T ∘ [n]^* = [n]^* ∘ τ_{nT}`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 168 の `aut_comp_mulByN`(`n·T = O` の場合)の一般形。
★★交代性では `n·P' = P ≠ O` の点で平行移動するので、これが要る。 -/
theorem aut_comp_muExt_gen (hinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (n : ℕ) (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (τT τS : W.FunctionField ≃ₐ[F] W.FunctionField) (T S : W.Point)
    (hT : IsTranslate W τT T) (hS : IsTranslate W τS S) (hnT : S = n • T) (z : W.FunctionField) :
    τT (muExt W hinj z) = muExt W hinj (τS z) := by
  have hA : autFF W τT (Point.some xn yn hns)
      = n • ABC3.Found.GaloisRep.genericPoint W + toFF W S := by
    rw [← hμP, map_nsmul, autFF_generic_gen W τT T hT, nsmul_add, ← map_nsmul (toFF W), ← hnT]
  obtain ⟨hA', hAeq⟩ := autFF_some W τT hns
  rw [hAeq] at hA
  obtain ⟨hS', hSeq⟩ := autFF_some W τS (nonsingular_coord W)
  have hSgen : autFF W τS (ABC3.Found.GaloisRep.genericPoint W)
      = Point.some (τS (coordX W)) (τS (coordY W)) hS' := hSeq
  have hB : muFF W hinj hμF (Point.some (τS (coordX W)) (τS (coordY W)) hS')
      = n • ABC3.Found.GaloisRep.genericPoint W + toFF W S := by
    rw [← hSgen, autFF_generic_gen W τS S hS, map_add, muFF_toFF W hinj hμF,
      muFF_generic W hinj hμF hns hμx hμy, hμP]
  obtain ⟨hB', hBeq⟩ := muFF_some W hinj hμF hS'
  rw [hBeq] at hB
  have hcoord := hA.trans hB.symm
  have hx : τT xn = muExt W hinj (τS (coordX W)) := by injection hcoord
  have hy : τT yn = muExt W hinj (τS (coordY W)) := by injection hcoord
  have hring : ((τT : W.FunctionField →+* W.FunctionField).comp (muExt W hinj))
      = (muExt W hinj).comp (τS : W.FunctionField →+* W.FunctionField) := by
    refine functionField_hom_ext W _ _ (fun c => ?_) ?_ ?_
    · show τT (muExt W hinj (algebraMap F W.FunctionField c))
        = muExt W hinj (τS (algebraMap F W.FunctionField c))
      rw [muExt_const W hinj hμF c, τS.commutes, muExt_const W hinj hμF c]
      exact τT.commutes c
    · show τT (muExt W hinj (coordX W)) = muExt W hinj (τS (coordX W))
      rw [muExt_coordX W hinj, hμx]; exact hx
    · show τT (muExt W hinj (coordY W)) = muExt W hinj (τS (coordY W))
      rw [muExt_coordY W hinj, hμy]; exact hy
  exact congrFun (congrArg (fun f : W.FunctionField →+* W.FunctionField => f.toFun) hring) z

end Comm

/-! ## ★出典の紐付け(`.src`) -/

def isTranslate_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——IsTranslate の特徴づけ)",
    sectionId := "genell-thm-3-8" }

def aut_comp_muExt_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の交代性——τ_T ∘ [n]^* = [n]^* ∘ τ_{nT})",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
