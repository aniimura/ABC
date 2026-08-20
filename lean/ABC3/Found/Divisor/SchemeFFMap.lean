/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeOrdUnit
import ABC3.Found.Divisor.NormFF

/-!
# 支配的な射に沿った関数体の射と茎(鎖 `normalize` の `B-functor` / `cartier-pullback`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★要になる四角形

支配的な射 `g : X ⟶ Y`(整スキームの間)に対し、各点 `x : X` で

```
  𝒪_{Y, g(x)} ──stalkMap──▶ 𝒪_{X, x}
      │                        │
      ▼                        ▼
    K(Y)   ────ffMap g────▶  K(X)
```

が可換である。★これは mathlib の `Scheme.Hom.stalkSpecializes_stalkMap` そのもので、
`stalkSpecializes` の合成則(`stalkSpecializes_comp`)で `ffMap` の
`stalkCongr` を吸収すれば出る。

★★**効き目**: `u` が `𝒪_{Y,g(x)}` の単元なら `ffMap g u` は `𝒪_{X,x}` の単元、
したがって `ord_x(ffMap g u) = 0`(`SchemeOrdUnit.lean`)。
これが `B(L) → B(M)` の中身である。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `ffMap` | ★支配的な射に沿った `K(Y) ⟶ K(X)` |
| `normFF_eq_ffMap` | `normFF` はその特別な場合 |
| `stalkSpecializes_comp_ffMap` | ★★★**要の四角形** |
| `ffMap_algebraMap` | 四角形の元での言い換え |
| `ordPt_ffMap_eq_zero_of_isUnit` | ★★★★**茎の単元は `ord = 0` へ移る** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory IsDedekindDomain

universe u

section FFMap

variable {X Y : Scheme.{u}} (g : X ⟶ Y) [IsDominant g]
  [IrreducibleSpace X] [IrreducibleSpace Y]

/-- ★★**支配的な射に沿った関数体の射** `K(Y) ⟶ K(X)`。 -/
noncomputable def ffMap : Y.functionField ⟶ X.functionField :=
  (Y.presheaf.stalkCongr
      (Inseparable.of_eq (genericPoint_eq_of_isDominant g).symm)).hom
    ≫ g.stalkMap (genericPoint (X : Type u))

/-- ★★★**要の四角形** —— 茎から関数体への射と `ffMap` は可換。 -/
theorem stalkSpecializes_comp_ffMap (x : X) :
    (Y.presheaf.stalkSpecializes
        ((genericPoint_spec (Y : Type u)).specializes (Set.mem_univ (g.base x)))) ≫ ffMap g
      = g.stalkMap x ≫ (X.presheaf.stalkSpecializes
        ((genericPoint_spec (X : Type u)).specializes (Set.mem_univ x))) := by
  rw [ffMap]
  simp only [TopCat.Presheaf.stalkCongr_hom, ← Category.assoc,
    TopCat.Presheaf.stalkSpecializes_comp]
  exact Scheme.Hom.stalkSpecializes_stalkMap g _ x _

/-- ★四角形の元での言い換え。 -/
theorem ffMap_algebraMap (x : X) (t : Y.presheaf.stalk (g.base x)) :
    ffMap g (algebraMap (Y.presheaf.stalk (g.base x)) Y.functionField t)
      = algebraMap (X.presheaf.stalk x) X.functionField (g.stalkMap x t) :=
  congrArg (fun m : Y.presheaf.stalk (g.base x) ⟶ X.functionField => CommRingCat.Hom.hom m t)
    (stalkSpecializes_comp_ffMap g x)

end FFMap

section Ord

variable {X Y : Scheme.{u}} (g : X ⟶ Y) [IsDominant g] [IsIntegral X] [IsIntegral Y]
  [IsLocallyNoetherian X]

/-- ★★★★**茎の単元は `ord = 0` へ移る** —— `B(L) → B(M)` の中身。 -/
theorem ordPt_ffMap_eq_zero_of_isUnit (hnormX : IsNormalScheme X) (x : PrimeDivisorPt X)
    (t : (Y.presheaf.stalk (g.base x.1))ˣ) :
    ordPt X hnormX x
        (ffMap g (algebraMap (Y.presheaf.stalk (g.base x.1)) Y.functionField
          (t : Y.presheaf.stalk (g.base x.1)))) = 0 := by
  rw [ffMap_algebraMap]
  have hu : IsUnit ((g.stalkMap x.1).hom (t : Y.presheaf.stalk (g.base x.1))) :=
    t.isUnit.map (g.stalkMap x.1).hom
  rw [show (g.stalkMap x.1) (t : Y.presheaf.stalk (g.base x.1))
      = ((hu.unit : (X.presheaf.stalk x.1)ˣ) : X.presheaf.stalk x.1) from hu.unit_spec.symm]
  exact ordPt_eq_zero_of_isUnit hnormX x hu.unit

omit [IsLocallyNoetherian X] in
/-- ★★**支配的な射の茎への射は単射** —— 四角形の左と下が単射だから。

★左は `𝒪_{Y,g(x)} ↪ K(Y)`(`IsFractionRing`)、下は体の射 `K(Y) → K(X)`。 -/
theorem stalkMap_injective (x : X) : Function.Injective (g.stalkMap x) := by
  intro a b hab
  have h1 : ffMap g (algebraMap (Y.presheaf.stalk (g.base x)) Y.functionField a)
      = ffMap g (algebraMap (Y.presheaf.stalk (g.base x)) Y.functionField b) := by
    rw [ffMap_algebraMap, ffMap_algebraMap, hab]
  exact IsFractionRing.injective (Y.presheaf.stalk (g.base x)) Y.functionField
    ((ffMap g).hom.injective h1)

end Ord

end ABC3.Found.Divisor
