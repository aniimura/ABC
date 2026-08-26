/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormNormal
import Mathlib.AlgebraicGeometry.Morphisms.UnderlyingMap

/-!
# `V[L]` の関数体の関手性(鎖 `normalize` の `B-functor` / `cartier-pullback` の共通の道具)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L]

## ★★何が要るか

`B(L)`(`NormB.lean`)も Cartier 因子の引き戻し(`cartier-pullback`)も、
`L → M` に対する**関数体の射 `K(V[L]) → K(V[M])`** を経由する。
★その射は `V[M] → V[L]` が**支配的**であることから出る:
支配的な射は生成点を生成点へ送るので、生成点の茎の間の `stalkMap` がそれである。

## ★★★支配性は 2 段

| 段 | 根拠 |
|---|---|
| `Spec M → Spec L` は全射 | 体のスペクトルは 1 点(`Unique (Spec (.of K))`) |
| `Spec L → V[L]` は支配的 | mathlib の `instance : IsDominant f.toNormalization` |
| ⟹ `V[M] → V[L]` は支配的 | `IsDominant.of_comp` を `normUp_normMap` に当てる |

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `genericPoint_eq_of_isDominant` | ★★**支配的な射は生成点を生成点へ送る** |
| `normMap_isDominant` | `V[M] → V[L]` は支配的 |
| `normMap_genericPoint` | 生成点は生成点へ |
| `normFF` | ★★★**`K(V[L]) ⟶ K(V[M])`** |
| `normFF_id` / `normFF_comp` | ★関手性 |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

/-! ## ★1. 支配的な射は生成点を生成点へ送る -/

/-- ★★**支配的な射は生成点を生成点へ送る**。

★`genericPoint_eq_of_isOpenImmersion`(mathlib)と同じ筋で、
開埋め込みの代わりに**像が稠密であること**だけを使う。 -/
theorem genericPoint_eq_of_isDominant {X Y : Scheme.{u}} (g : X ⟶ Y) [IsDominant g]
    [IrreducibleSpace X] [IrreducibleSpace Y] :
    g.base (genericPoint (X : Type u)) = genericPoint (Y : Type u) := by
  have hX : IsGenericPoint (genericPoint (X : Type u)) Set.univ := by
    simpa using genericPoint_spec (X : Type u)
  have hY : IsGenericPoint (genericPoint (Y : Type u)) Set.univ := by
    simpa using genericPoint_spec (Y : Type u)
  refine (hY.eq ?_).symm
  have h := hX.image g.continuous
  rwa [Set.image_univ, (Scheme.Hom.denseRange g).closure_range] at h

/-! ## ★2. `V[M] → V[L]` は支配的 -/

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

instance normUp_isDominant (L : FinSub V.functionField Kbar) :
    IsDominant (normUp V L) :=
  inferInstanceAs (IsDominant (Scheme.Hom.toNormalization (specToV V L)))

/-- ★体のスペクトルは 1 点なので、体の拡大が定める射は全射。 -/
instance specMapOf_surjective {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    Surjective (specMapOf V f) :=
  ⟨fun _ => ⟨default, Subsingleton.elim _ _⟩⟩

/-- ★★**`V[M] → V[L]` は支配的**。 -/
instance normMap_isDominant {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    IsDominant (normMap V f) := by
  have h : IsDominant (normUp V M ≫ normMap V f) := by
    rw [normUp_normMap]
    infer_instance
  exact IsDominant.of_comp (normUp V M) (normMap V f)

/-- ★★**`V[M] → V[L]` は整射** —— `normMap ≫ normDown = normDown` と
`IsIntegralHom.of_comp`(`normDown` はアフィン射なので分離的)。 -/
instance normMap_isIntegralHom {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    IsIntegralHom (normMap V f) := by
  haveI : IsIntegralHom (normMap V f ≫ normDown V L) := by
    rw [normMap_normDown]
    infer_instance
  exact IsIntegralHom.of_comp (normMap V f) (normDown V L)

/-- ★**生成点は生成点へ**。 -/
theorem normMap_genericPoint {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    (normMap V f).base (genericPoint (normObj V M : Type u))
      = genericPoint (normObj V L : Type u) :=
  genericPoint_eq_of_isDominant _

/-! ## ★3. 関数体の射 -/

/-- ★★★**`K(V[L]) ⟶ K(V[M])`** —— 生成点の茎の間の `stalkMap`。 -/
noncomputable def normFF {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    (normObj V L).functionField ⟶ (normObj V M).functionField :=
  ((normObj V L).presheaf.stalkCongr
      (Inseparable.of_eq (normMap_genericPoint V f).symm)).hom
    ≫ (normMap V f).stalkMap (genericPoint (normObj V M : Type u))

/-- ★関手性(恒等)。 -/
theorem normFF_id (L : FinSub V.functionField Kbar) :
    normFF V (𝟙 L) = 𝟙 ((normObj V L).functionField) := by
  rw [normFF, Scheme.Hom.stalkMap_congr_hom (normMap V (𝟙 L)) (𝟙 _) (normMap_id V L),
    Scheme.Hom.stalkMap_id]
  simp

/-- ★関手性(合成)。 -/
theorem normFF_comp {L M N : FinSub V.functionField Kbar} (f : L ⟶ M) (g : M ⟶ N) :
    normFF V (f ≫ g) = normFF V f ≫ normFF V g := by
  rw [normFF, normFF, normFF,
    Scheme.Hom.stalkMap_congr_hom (normMap V (f ≫ g)) (normMap V g ≫ normMap V f)
      (normMap_comp V f g),
    Scheme.Hom.stalkMap_comp]
  simp only [Category.assoc]
  rw [Scheme.Hom.stalkMap_congr_point_assoc (normMap V f)
    (genericPoint (normObj V M : Type u)) ((normMap V g).base (genericPoint (normObj V N : Type u)))
    (normMap_genericPoint V g).symm]
  simp

/-- ★★★★**`L ↦ K(V[L])` は関手**(共変)。 -/
noncomputable def normFFFunctor :
    FinSub V.functionField Kbar ⥤ CommRingCat.{u} where
  obj L := (normObj V L).functionField
  map f := normFF V f
  map_id L := normFF_id V L
  map_comp f g := normFF_comp V f g

/-- ★関数体の射は単射(体の射だから)。 -/
theorem normFF_injective {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    Function.Injective (normFF V f) :=
  (normFF V f).hom.injective

/-- ★★**単元群への誘導** —— `B(L)` の関手性で使う形。 -/
noncomputable def normFFUnits {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    ((normObj V L).functionField)ˣ →* ((normObj V M).functionField)ˣ :=
  Units.map (normFF V f).hom.toMonoidHom

@[simp] theorem normFFUnits_val {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (u : ((normObj V L).functionField)ˣ) :
    ((normFFUnits V f u : (normObj V M).functionField))
      = normFF V f (u : (normObj V L).functionField) := rfl

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `L ↦ K(V[L])` の関手性。 -/
def normFFFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — K(V[L]) の関手性(B(L) と Cartier 因子の引き戻しの共通の道具)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
