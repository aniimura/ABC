/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeCartier
import ABC3.Found.FrdI.Sec6GaloisCat
import Mathlib.AlgebraicGeometry.Normalization
import Mathlib.AlgebraicGeometry.Stalk
import Mathlib.AlgebraicGeometry.Morphisms.QuasiSeparated
import Mathlib.AlgebraicGeometry.Morphisms.Integral

/-!
# `V[L]` —— `V` の `L` における正規化を `𝒟 = B(G)⁰` 上の関手として

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L]

## ★★`V[L]` の作り方

原文は `L ⊆ K̄` を `K = k(V)` の有限拡大として、`V[L]` を
「`V` の `L` における正規化」と定める。★mathlib には**相対正規化**があるので、

  `V[L] := (Spec L → V).normalization`

とすればよい。`Spec L → V` は `Spec.map (K → L)` と
`V.fromSpecStalk (genericPoint V)`(生成点からの射)の合成である。

## ★★★関手性は普遍性 1 本

`L → M` に対する `V[M] → V[L]` は mathlib の
**`Scheme.Hom.normalizationDesc`**(相対正規化の普遍性)そのもの。
★関手則(`map_id` / `map_comp`)は **`normalization.hom_ext`**(一意性)で出る。

★★**qcqs の側条件**は自動である —— 源が体のスペクトル(1 点、したがって Noether 空間、
かつアフィンなので quasi-separated)だから。

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `specToV` | `Spec L → V` |
| `normObj` | ★`V[L]` |
| `normMap` | ★`V[M] → V[L]`(普遍性) |
| `normFunctor` | ★★★**`L ↦ V[L]` は `(FinSub K K̄)ᵒᵖ ⥤ Scheme`** |
| `normObj_isIntegral` | `V[L]` は整スキーム |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-! ## ★1. `Spec L → V` -/

/-- ★**`Spec L → V`** —— 生成点からの射に体の拡大を合成したもの。 -/
noncomputable def specToV (L : FinSub V.functionField Kbar) :
    Spec (CommRingCat.of L.toIF) ⟶ V :=
  Spec.map (CommRingCat.ofHom (algebraMap V.functionField L.toIF))
    ≫ V.fromSpecStalk (genericPoint V)

instance quasiCompact_specToV (L : FinSub V.functionField Kbar) :
    QuasiCompact (specToV V L) :=
  quasiCompact_of_noetherianSpace_source _

instance quasiSeparated_specToV (L : FinSub V.functionField Kbar) :
    QuasiSeparated (specToV V L) :=
  QuasiSeparated.of_quasiSeparatedSpace _

/-- ★射 `f : L → M` が定める `Spec M → Spec L`。 -/
noncomputable def specMapOf {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    Spec (CommRingCat.of M.toIF) ⟶ Spec (CommRingCat.of L.toIF) :=
  Spec.map (CommRingCat.ofHom (FinSub.hom f).toRingHom)

/-- ★三角形の可換性 —— 中身は `AlgHom.commutes` 1 本。 -/
theorem specMapOf_specToV {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    specMapOf V f ≫ specToV V L = specToV V M := by
  rw [specToV, specToV, specMapOf, ← Category.assoc, ← Spec.map_comp]
  congr 2
  refine CommRingCat.hom_ext (RingHom.ext fun x => ?_)
  exact (FinSub.hom f).commutes x

@[simp] theorem specMapOf_id (L : FinSub V.functionField Kbar) :
    specMapOf V (𝟙 L) = 𝟙 _ := by
  rw [specMapOf]
  rw [show CommRingCat.ofHom (FinSub.hom (𝟙 L)).toRingHom = 𝟙 _ from rfl, Spec.map_id]

theorem specMapOf_comp {L M N : FinSub V.functionField Kbar} (f : L ⟶ M) (g : M ⟶ N) :
    specMapOf V (f ≫ g) = specMapOf V g ≫ specMapOf V f := by
  rw [specMapOf, specMapOf, specMapOf, ← Spec.map_comp]
  rfl

/-! ## ★2. `V[L]` -/

/-- ★★★**`V[L]`** —— `V` の `L` における正規化。 -/
noncomputable def normObj (L : FinSub V.functionField Kbar) : Scheme.{u} :=
  (specToV V L).normalization

/-- ★`V[L] → V`(整射)。 -/
noncomputable def normDown (L : FinSub V.functionField Kbar) : normObj V L ⟶ V :=
  Scheme.Hom.fromNormalization (specToV V L)

instance normDown_isIntegralHom (L : FinSub V.functionField Kbar) :
    IsIntegralHom (normDown V L) :=
  inferInstanceAs (IsIntegralHom (Scheme.Hom.fromNormalization (specToV V L)))

/-- ★`Spec L → V[L]`。 -/
noncomputable def normUp (L : FinSub V.functionField Kbar) :
    Spec (CommRingCat.of L.toIF) ⟶ normObj V L :=
  Scheme.Hom.toNormalization (specToV V L)

@[reassoc] theorem normUp_normDown (L : FinSub V.functionField Kbar) :
    normUp V L ≫ normDown V L = specToV V L :=
  Scheme.Hom.toNormalization_fromNormalization _

/-- ★`V[L]` は整スキーム。 -/
instance normObj_isIntegral (L : FinSub V.functionField Kbar) :
    IsIntegral (normObj V L) :=
  inferInstanceAs (IsIntegral (specToV V L).normalization)

/-! ## ★3. 関手性 -/

/-- ★★**`V[M] → V[L]`** —— 相対正規化の普遍性(`normalizationDesc`)。 -/
noncomputable def normMap {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    normObj V M ⟶ normObj V L :=
  Scheme.Hom.normalizationDesc (specToV V M) (specMapOf V f ≫ normUp V L) (normDown V L)
    (by rw [Category.assoc, normUp_normDown, specMapOf_specToV])

@[reassoc] theorem normMap_normDown {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    normMap V f ≫ normDown V L = normDown V M :=
  Scheme.Hom.normalizationDesc_comp _ _ _ _

@[reassoc] theorem normUp_normMap {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    normUp V M ≫ normMap V f = specMapOf V f ≫ normUp V L :=
  Scheme.Hom.toNormalization_normalizationDesc _ _ _ _

theorem normMap_id (L : FinSub V.functionField Kbar) :
    normMap V (𝟙 L) = 𝟙 (normObj V L) := by
  haveI : IsAffineHom (normDown V L) := inferInstance
  refine Scheme.Hom.normalization.hom_ext (f := specToV V L) (normMap V (𝟙 L))
    (𝟙 _) (normDown V L) ?_ ?_ ?_
  · show normUp V L ≫ normMap V (𝟙 L) = normUp V L ≫ 𝟙 (normObj V L)
    rw [normUp_normMap, specMapOf_id, Category.id_comp, Category.comp_id]
  · exact normMap_normDown V (𝟙 L)
  · show 𝟙 (normObj V L) ≫ normDown V L = normDown V L
    rw [Category.id_comp]

theorem normMap_comp {L M N : FinSub V.functionField Kbar} (f : L ⟶ M) (g : M ⟶ N) :
    normMap V (f ≫ g) = normMap V g ≫ normMap V f := by
  haveI : IsAffineHom (normDown V L) := inferInstance
  refine Scheme.Hom.normalization.hom_ext (f := specToV V N) (normMap V (f ≫ g))
    (normMap V g ≫ normMap V f) (normDown V L) ?_ ?_ ?_
  · show normUp V N ≫ normMap V (f ≫ g) = normUp V N ≫ (normMap V g ≫ normMap V f)
    rw [normUp_normMap, specMapOf_comp, Category.assoc, ← normUp_normMap, ← Category.assoc,
      ← normUp_normMap, Category.assoc]
  · exact normMap_normDown V (f ≫ g)
  · show (normMap V g ≫ normMap V f) ≫ normDown V L = normDown V N
    rw [Category.assoc, normMap_normDown, normMap_normDown]

/-- ★★★★★★**`L ↦ V[L]` は関手** —— `Example 6.1` の幾何の骨。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L] -/
noncomputable def normFunctor :
    (FinSub V.functionField Kbar)ᵒᵖ ⥤ Scheme.{u} where
  obj L := normObj V L.unop
  map {_ _} α := normMap V α.unop
  map_id L := normMap_id V L.unop
  map_comp {_ _ _} α β := normMap_comp V β.unop α.unop

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Example 6.1` の `V[L]`(相対正規化)とその関手性。 -/
def normFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — V[L](V の L における正規化)と L についての関手性",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
