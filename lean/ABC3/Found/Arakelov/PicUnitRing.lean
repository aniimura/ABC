import ABC3.Found.Arakelov.PicUnitSect

/-!
# Arakelov (B2) 第 166 ブロック —— ★★★★★★★**双対の加群構造**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★双対層への道が開いた

    unitEndEquiv : End (𝟙_ (PresheafModulesOn X U)) ≃+* Γ(X, U)
    dualModule   : Module Γ(X,U) ((restrictPresheafFunctor X U).obj F ⟶ 𝟙_)

★mathlib の `Preadditive.moduleEndRight : Module (End Y) (X ⟶ Y)` を
`Module.compHom` で `Γ(X,U)` へ移すだけである。

## ★★環同型の筋

| 段 | 内容 |
|---|---|
| `freeYonedaEquiv_apply'` | ★`freeYonedaEquiv m = m.app t (freeMk 𝟙)`(`symm_app` から) |
| `termIso_hom_gen` | ★`freeYonedaTermIso.hom` は生成元を `1` に送る |
| `unitEnd_unitMul` | ★第 165(全射) |
| `unitMul_unitEnd` | ★★**単射**——`hom ≫ −` で消去 |

★★`hom ≫ unitMul c = freeYonedaEquiv.symm c` が要点である
(`unitMul` の定義の `inv` が `hom` と消える)。

## ★★★逃げ道——**関手合成の括弧**

`freeYonedaEquiv_apply'` の `M` の型を
`(Over.forget U).op ⋙ X.presheaf ⋙ forget₂ …` と書くと `rw` が当たらない。
★`PresheafModulesOn X U` は **`((… ⋙ …) ⋙ …)`** なので括弧を合わせること。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (U : X.Opens)

/-- ★`freeYonedaEquiv` の値は生成元での値。 -/
theorem freeYonedaEquiv_apply'
    {M : PresheafOfModules.{u} (((Over.forget U).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})}
    (m : PresheafOfModules.freeObj (yoneda.obj (Over.mk (𝟙 U))) ⟶ M) :
    PresheafOfModules.freeYonedaEquiv m
      = (m.app (op (Over.mk (𝟙 U)))).hom (ModuleCat.freeMk (𝟙 (Over.mk (𝟙 U)))) := by
  have h := PresheafOfModules.freeYonedaEquiv_symm_app M (Over.mk (𝟙 U))
    (PresheafOfModules.freeYonedaEquiv m)
  rw [Equiv.symm_apply_apply] at h
  exact h.symm

/-- ★生成元での値。 -/
theorem termIso_hom_gen :
    ((freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom.app
      (op (Over.mk (𝟙 U)))).hom (ModuleCat.freeMk (𝟙 (Over.mk (𝟙 U)))) = unitOne U := by
  show (PresheafOfModules.freeObjDesc _).app _ _ = _
  erw [ModuleCat.freeDesc_apply]
  rfl

set_option maxHeartbeats 1000000 in
/-- ★★`unitMul ∘ unitEnd = id`。 -/
theorem unitMul_unitEnd (φ : End (𝟙_ (PresheafModulesOn X U))) :
    unitMul U (unitEnd U φ) = φ := by
  have key : ∀ c : (Γ(X, U) : Type u),
      (freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom ≫ unitMul U c
        = PresheafOfModules.freeYonedaEquiv.symm c := by
    intro c
    show _ ≫ ((freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).inv ≫ _) = _
    rw [← Category.assoc,
      (freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom_inv_id, Category.id_comp]
  have h3 : PresheafOfModules.freeYonedaEquiv
      ((freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom ≫ φ) = unitEnd U φ := by
    rw [PresheafOfModules.freeYonedaEquiv_comp, freeYonedaEquiv_apply' U, termIso_hom_gen]
    rfl
  have hcancel : (freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom
      ≫ unitMul U (unitEnd U φ)
    = (freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom ≫ φ := by
    rw [key, ← h3, Equiv.symm_apply_apply]
    rfl
  exact (cancel_epi _).1 hcancel

/-- ★★★`unitEnd` は環同型。 -/
noncomputable def unitEndEquiv :
    End (𝟙_ (PresheafModulesOn X U)) ≃+* (Γ(X, U) : Type u) :=
  RingEquiv.ofBijective (unitEnd U)
    ⟨fun a b h => by rw [← unitMul_unitEnd U a, ← unitMul_unitEnd U b, h],
     fun c => ⟨unitMul U c, unitEnd_unitMul U c⟩⟩

/-- ★★★★双対の `Γ(X,U)` 加群構造。 -/
noncomputable instance dualModule (F : X.PresheafOfModules) :
    Module (((X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op U)) : Type u)
      ((restrictPresheafFunctor X U).obj F ⟶ 𝟙_ (PresheafModulesOn X U)) :=
  Module.compHom _ (unitEndEquiv U).symm.toRingHom


/-- ★★★★双対の加群構造(**反対圏の索引版**)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★前層に組む段では係数環が `(…).obj W`(`W : (X.Opens)ᵒᵖ`)の形で現れる。
★★`(op W.unop) = W` は eta で等しいが **instance 探索は `op ?U =?= W` を解けない**ので、
索引を変えた同じ instance をもう 1 本置く。

★★★中身は同じ項に簡約されるので、経路が 2 本になる問題は起きない。 -/
noncomputable instance dualModuleOp (F : X.PresheafOfModules) (W : (X.Opens)ᵒᵖ) :
    Module (((X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj W) : Type u)
      ((restrictPresheafFunctor X W.unop).obj F ⟶ 𝟙_ (PresheafModulesOn X W.unop)) :=
  dualModule W.unop F

/-! ## ★出典の紐付け(`.src`) -/

def unitEndEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——双対の加群構造)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
