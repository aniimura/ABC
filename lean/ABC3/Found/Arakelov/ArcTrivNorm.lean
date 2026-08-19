import ABC3.Found.Arakelov.ArcPointScheme

/-!
# Arakelov (C3) 第 248 ブロック —— ★★★★★**局所自明な層には計量が入る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★「等しいが綴りが違う」の新しい逃げ道

ファイバー `arcFiber p L` の `ℂ` 作用は `restrictScalars (ΓSpecIso).inv` を通しており、
係数環 `Γ(Spec ℂ, ⊤)` の作用と**定義的に等しいが綴りが違う**。
★`rw` も `show` も型上書き `(e : T)` も通らなかった
——インスタンス探索は**構文的な型**を見るからである。

★★★逃げ道は**ゴール側の綴りで `have` を立て、`exact` で定義的等しさに任せる**ことである:

    have h : topMap e (c • v) = ΓSpecIso.inv.hom c * topMap e v :=
      topMap_smul e (ΓSpecIso.inv.hom c) v

★`have` の**型は `ℂ` 作用**で書き、証明項は**係数作用**の補題である。
型検査は `isDefEq` なので通る——`rw` のような構文照合を経由しない。

★★★★これは「等しいが綴りが違う」の新しい逃げ道であり、
**「書き換えられないなら、書き換えずに主張する」**と要約できる。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `topMap` | ★自明化から `⊤` 切断の間の写像 |
| `topMap_bijective` / `topMap_smul` | ★★全単射・係数線形 |
| `fiberNorm` | ★★★ファイバーのノルム |
| `fiberNorm_smul` / `fiberNorm_eq_zero_iff` | ★★★★ノルムの 2 法則 |
| `arcMetricOf` | ★★★★★**局所自明 ⟹ 計量が入る**(連続性はまだ課さない) |
-/
universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {Z : Scheme.{u}} {M : Z.PresheafOfModules}

noncomputable def topMap
    (e : (restrictPresheafFunctor Z ⊤).obj M ≅ 𝟙_ (PresheafModulesOn Z ⊤))
    (v : (M.obj (op ⊤) : Type u)) : (Z.presheaf.obj (op ⊤) : Type u) :=
  (e.hom.app (op (Over.mk (𝟙 (⊤ : Z.Opens))))).hom (fVal M (op ⊤) v)

theorem topMap_bijective
    (e : (restrictPresheafFunctor Z ⊤).obj M ≅ 𝟙_ (PresheafModulesOn Z ⊤)) :
    Function.Bijective (topMap e) := by
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map e.hom) := inferInstance
  rw [NatTrans.isIso_iff_isIso_app] at this
  haveI := this (op (Over.mk (𝟙 (⊤ : Z.Opens))))
  have hb : Function.Bijective (((PresheafOfModules.toPresheaf _).map e.hom).app
      (op (Over.mk (𝟙 (⊤ : Z.Opens))))) := ConcreteCategory.bijective_of_isIso _
  exact hb

theorem topMap_smul
    (e : (restrictPresheafFunctor Z ⊤).obj M ≅ 𝟙_ (PresheafModulesOn Z ⊤))
    (c : (Z.presheaf.obj (op ⊤) : Type u)) (v : (M.obj (op ⊤) : Type u)) :
    topMap e (c • v) = c * topMap e v := by
  show (e.hom.app (op (Over.mk (𝟙 (⊤ : Z.Opens))))).hom (fVal M (op ⊤) (c • v)) = _
  rw [fVal_smul]
  exact map_smul (e.hom.app (op (Over.mk (𝟙 (⊤ : Z.Opens))))).hom (rVal (X := Z) (op ⊤) c)
    (fVal M (op ⊤) v)

theorem topMap_zero'
    (e : (restrictPresheafFunctor Z ⊤).obj M ≅ 𝟙_ (PresheafModulesOn Z ⊤)) :
    topMap e (0 : (M.obj (op ⊤) : Type u)) = 0 :=
  map_zero (e.hom.app (op (Over.mk (𝟙 (⊤ : Z.Opens))))).hom

section Fiber

variable {X : Scheme.{0}} (p : Spec (CommRingCat.of ℂ) ⟶ X) (L : X.Modules)
  (e : (restrictPresheafFunctor (Spec (CommRingCat.of ℂ)) ⊤).obj
      ((Scheme.Modules.pullback p).obj L).val
    ≅ 𝟙_ (PresheafModulesOn (Spec (CommRingCat.of ℂ)) ⊤))

noncomputable def fiberNorm (v : ↥(arcFiber p L)) : ℝ :=
  ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (topMap e v)‖

theorem fiberNorm_nonneg (v : ↥(arcFiber p L)) : 0 ≤ fiberNorm p L e v :=
  norm_nonneg _

theorem fiberNorm_smul (c : (CommRingCat.of ℂ : Type)) (v : ↥(arcFiber p L)) :
    fiberNorm p L e (c • v) = ‖c‖ * fiberNorm p L e v := by
  have h : topMap e (c • v)
      = (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom c * topMap e v :=
    topMap_smul e ((Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom c) v
  have hc : (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom c) = c :=
    congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) c)
      (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv_hom_id
  show ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (topMap e (c • v))‖
    = ‖c‖ * ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (topMap e v)‖
  rw [h, map_mul, hc]
  exact norm_mul _ _

theorem gammaIso_inv_hom (y : ((Spec (CommRingCat.of ℂ)).presheaf.obj (op ⊤) : Type)) :
    (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom
      ((Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom y) = y :=
  congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) y)
    (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom_inv_id

theorem fiberNorm_eq_zero_iff (v : ↥(arcFiber p L)) :
    fiberNorm p L e v = 0 ↔ v = 0 := by
  show ‖(Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom (topMap e v)‖ = 0 ↔ v = 0
  rw [norm_eq_zero]
  constructor
  · intro h
    have h0 : topMap e v = 0 := by
      have hh := congrArg (fun y => (Scheme.ΓSpecIso (CommRingCat.of ℂ)).inv.hom y) h
      simpa [gammaIso_inv_hom, map_zero] using hh
    refine (topMap_bijective e).1 ?_
    rw [h0]
    exact (topMap_zero' e).symm
  · intro h
    have h0 : topMap e v = 0 := by rw [h]; exact topMap_zero' e
    exact (congrArg (fun y => (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom y) h0).trans
      (map_zero _)

end Fiber

section Exists

variable {X : Scheme.{0}} (L : X.Modules)

/-- ★各複素点での自明化は存在する。 -/
theorem exists_triv (hL : IsLocallyTrivial X L.val) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    Nonempty ((restrictPresheafFunctor (Spec (CommRingCat.of ℂ)) ⊤).obj
        ((Scheme.Modules.pullback p).obj L).val
      ≅ 𝟙_ (PresheafModulesOn (Spec (CommRingCat.of ℂ)) ⊤)) :=
  trivial_of_locallyTrivial (isLocallyTrivial_pullbackModules p L hL)

/-- ★★★★★局所自明な層には計量が入る(連続性はまだ課さない)。 -/
noncomputable def arcMetricOf (hL : IsLocallyTrivial X L.val) : ArcMetric X L where
  nrm p := fiberNorm p L (Classical.choice (exists_triv L hL p))
  nonneg p := fiberNorm_nonneg p L _
  eq_zero_iff p := fiberNorm_eq_zero_iff p L _
  smul p := fiberNorm_smul p L _

end Exists

/-! ## ★出典の紐付け(`.src`) -/

def arcMetricOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——局所自明な層には計量が入る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
