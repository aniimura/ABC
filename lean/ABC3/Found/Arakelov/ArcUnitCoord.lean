import ABC3.Found.Arakelov.ArcScaleNorm

/-!
# Arakelov (C3) 第 287 ブロック —— **★★★★★単位加群のファイバーは 1 次元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`logMetric` を honest に定義するための土台

Green 関数 `logMetric` は**基準計量に対する比**でしか定まらない。
比が「どのベクトルで測っても同じ」であるためには、
★★**ファイバーが 1 次元**でなければならない。

## ★★スカラー作用の実測(2026-08-19)

mathlib の `modulesSpecToSheaf`(`Modules/Tilde.lean:44`)は

    … ⋙ sheafCompose _ (ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom)

★したがって `c • x` は **`(ΓSpecIso).inv.hom c • x`**(Γ 側の作用)である。
★★担い手の綴りが違うので `*` が直接書けない
——**逃げ道 7(型付き恒等関数の橋)** `specΓVal` を 1 本足して通した。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `specΓVal` | ★型付き恒等関数の橋 |
| `unitCoord` | ★★単位加群のファイバーの座標(ℂ の値) |
| `unitCoord_smul` | ★★★スカラー倍と両立 |
| `unitCoord_bijective` | ★★★★**全単射**——すなわち 1 次元 |
| `exists_smul_eq_unit` | ★★★★★**0 でない元は生成元である** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

/-- ★**型付き恒等関数の橋**——ファイバーの元を `Γ(Spec ℂ, ⊤)` の元と見る。 -/
def specΓVal (x : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
    (unitModules (Spec (CommRingCat.of ℂ))) : Type)) :
    ((Spec (CommRingCat.of ℂ)).presheaf.obj (op ⊤) : Type) := x

theorem specΓVal_bijective : Function.Bijective specΓVal :=
  ⟨fun _ _ h => h, fun y => ⟨y, rfl⟩⟩

/-- ★★スカラー作用は `ΓSpecIso.inv` を通した掛け算である(`rfl`)。 -/
theorem specΓVal_smul (c : (CommRingCat.of ℂ : Type))
    (x : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
      (unitModules (Spec (CommRingCat.of ℂ))) : Type)) :
    specΓVal (c • x) = (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom c * specΓVal x := rfl

/-- ★★**単位加群のファイバーの座標**。 -/
noncomputable def unitCoord (x : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
    (unitModules (Spec (CommRingCat.of ℂ))) : Type)) : (CommRingCat.of ℂ : Type) :=
  (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (specΓVal x)

/-- ★★★★**座標は全単射である**——ファイバーは 1 次元。 -/
theorem unitCoord_bijective : Function.Bijective unitCoord :=
  Function.Bijective.comp
    (ConcreteCategory.bijective_of_isIso (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom)
    specΓVal_bijective

/-- ★★★座標はスカラー倍と両立する。 -/
theorem unitCoord_smul (c : (CommRingCat.of ℂ : Type))
    (x : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
      (unitModules (Spec (CommRingCat.of ℂ))) : Type)) :
    unitCoord (c • x) = c * unitCoord x := by
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (specΓVal (c • x)) = _
  rw [specΓVal_smul, map_mul]
  congr 1
  exact congrArg (fun f : CommRingCat.of ℂ ⟶ CommRingCat.of ℂ => f.hom c)
    (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv_hom_id

theorem unitCoord_eq_zero_iff (x : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
    (unitModules (Spec (CommRingCat.of ℂ))) : Type)) : unitCoord x = 0 ↔ x = 0 := by
  constructor
  · intro h
    refine unitCoord_bijective.1 ?_
    rw [h]
    show (0 : (CommRingCat.of ℂ : Type)) = unitCoord 0
    show (0 : (CommRingCat.of ℂ : Type))
      = (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (specΓVal 0)
    rw [show specΓVal 0 = 0 from rfl, map_zero]
  · intro h
    rw [h]
    show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (specΓVal 0) = 0
    rw [show specΓVal 0 = 0 from rfl, map_zero]

/-- ★★★★★**0 でない元は生成元である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「計量の比がベクトルの選び方に依らない」ことの根拠である。 -/
theorem exists_smul_eq_unit (x y : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
    (unitModules (Spec (CommRingCat.of ℂ))) : Type)) (hx : x ≠ 0) :
    ∃ c : (CommRingCat.of ℂ : Type), c • x = y := by
  have hcx : unitCoord x ≠ 0 := fun h => hx ((unitCoord_eq_zero_iff x).1 h)
  refine ⟨unitCoord y * (unitCoord x)⁻¹, unitCoord_bijective.1 ?_⟩
  rw [unitCoord_smul, mul_assoc, inv_mul_cancel₀ hcx, mul_one]

/-! ## ★出典の紐付け(`.src`) -/

def exists_smul_eq_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——単位加群のファイバーが 1 次元であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
