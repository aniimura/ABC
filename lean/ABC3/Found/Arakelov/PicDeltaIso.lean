import ABC3.Found.Arakelov.PicTensorMu
import ABC3.Found.Arakelov.PicDeltaLift

/-!
# Arakelov (B1) 第 40 ブロック —— ★★★★★★★★**`δ` は同型である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★引き戻しは strong monoidal である

    δ : f^*(P ⊗ Q) ⟶ f^*P ⊗ f^*Q

が**すべての `P, Q` で同型**であること——すなわち
**`f^*` が strong monoidal であること**が取れた。

★★★これは §9-41 で訂正した通り、平坦性を要しない事実(Stacks 01BJ)であり、
**mathlib には無い**(2026-08-18 実測)。

## ★★★★道のり(第 18–40 ブロック)

| 段 | ブロック | 内容 |
|---|---|---|
| モノイダル | 18–22 | 係数変換 lax → 押し出し lax → 引き戻し oplax → **比較射 δ** |
| 生成元 | 23–27 | 余極限保存・生成元の具体形・`⊓` の保存 |
| 持ち上げ | 28–30 | 余極限で同型を持ち上げ、**δ を生成元 1 点に帰着** |
| 計算 | 31–39 | 随伴関手定理で作られた**抽象な `unit` を生成元の上で計算し切る** |
| 完成 | 40 | ★本ブロック |

## ★★★★★★核心は 1 行

> **両辺とも `freeYonedaEquiv (unit)` の `f⁻¹(V ⊓ W)` への制限である。**

★左辺は `unit` の値(抽象)、右辺は第 24 ブロックの同型の `inv` の値(具体)。
★★第 36 ブロック(生成元の像は制限で決まる)+ 第 34 ブロック
(生成元での値が等しい)で一致する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `delta_free_eq` | ★★★★生成元の上で `δ` は同型の合成に等しい |
| `isIso_pullbackDelta_free` | ★★★★★生成元の上で `δ` は同型 |
| `isIso_pullbackDelta` | ★★★★★★★★**`δ` はすべての `P, Q` で同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★`freeY V ⊗ freeY W ≅ freeY (V ⊓ W)` の逆(第 26 ブロック)。 -/
noncomputable abbrev iotaY (V W : Y.Opens) :
    freeY (Y := Y) (V ⊓ W) ⟶ freeY (Y := Y) V ⊗ freeY (Y := Y) W :=
  (freeYonedaInfIso (R := Y.presheaf) V W).inv

/-- ★★`X` 側の交わりの同型。 -/
noncomputable abbrev targetIso (V W : Y.Opens) :=
  freeYonedaInfIso (R := X.presheaf)
    ((Opens.map f.base).obj V) ((Opens.map f.base).obj W)

/-! ## ★★★★生成元の上での等式 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★**生成元の上で `δ` は同型の合成に等しい**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★核心は「両辺とも `freeYonedaEquiv (unit)` の制限である」ことである
(第 34・36 ブロック)。 -/
theorem delta_free_eq (V W : Y.Opens) :
    (pullbackPre f).map (iotaY V W) ≫ pullbackDelta f (freeY V) (freeY W)
      = (pullbackFreeYonedaIso f (V ⊓ W)).hom
        ≫ (targetIso f V W).inv
        ≫ ((pullbackFreeYonedaIso f V).symm ⊗ᵢ (pullbackFreeYonedaIso f W).symm).hom := by
  refine ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv
    _ _).injective ?_
  erw [Adjunction.homEquiv_naturality_left]
  rw [homEquiv_pullbackDelta]
  erw [Adjunction.homEquiv_naturality_right]
  rw [homEquiv_pullbackFreeYonedaIso_hom]
  refine PresheafOfModules.freeYonedaEquiv.injective ?_
  simp only [PresheafOfModules.freeYonedaEquiv_comp]
  erw [freeYonedaEquiv_freeYonedaInfIso_inv]
  erw [pullbackIso_hom_unit_gen]
  erw [tensor_mu_app_tmul]
  erw [CategoryTheory.Functor.map_comp, CategoryTheory.comp_apply]
  have hpm : ∀ {A B : X.PresheafOfModules} (g : A ⟶ B) (Z : Y.Opens)
      (x : ((PresheafOfModules.pushforward (pullbackPhi f)).obj A).obj (op Z)),
      ((PresheafOfModules.pushforward (pullbackPhi f)).map g).app (op Z) x
        = g.app (op ((Opens.map f.base).obj Z)) x := fun _ _ _ => rfl
  rw [hpm, hpm]
  erw [← freeYonedaEquiv_apply, freeYonedaEquiv_freeYonedaInfIso_inv]
  erw [PresheafOfModules.Monoidal.tensorHom_app, ModuleCat.MonoidalCategory.tensorHom_tmul]
  erw [unit_app_eq_inv_app, unit_app_eq_inv_app]
  rfl

/-! ## ★★★★★生成元の上で同型 -/

/-- ★★★★★★**生成元の上で `δ` は同型である**。

★右辺は同型の合成であり、`f^*(iotaY)` も同型だから。 -/
theorem isIso_pullbackDelta_free (V W : Y.Opens) :
    IsIso (pullbackDelta f (freeY V) (freeY W)) := by
  haveI : IsIso ((pullbackPre f).map (iotaY V W)) :=
    (pullbackPre f).map_isIso _
  haveI : IsIso ((pullbackPre f).map (iotaY V W)
      ≫ pullbackDelta f (freeY V) (freeY W)) := by
    rw [delta_free_eq]
    infer_instance
  exact IsIso.of_isIso_comp_left ((pullbackPre f).map (iotaY V W)) _

/-! ## ★★★★★★★★すべての対象で同型 -/

/-- ★★★★★★★★**`δ` はすべての `P, Q` で同型である**——
すなわち **`f^*` は strong monoidal である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `f^*(L ⊗ M) ≅ f^*L ⊗ f^*M` が得られ、
`PicardData.pullback_mul` が書ける。

★機構は第 30 ブロックの**二重の持ち上げ**である。 -/
theorem isIso_pullbackDelta (P Q : Y.PresheafOfModules) :
    IsIso (pullbackDelta f P Q) :=
  isIso_pullbackDelta_of_free f (fun V W => isIso_pullbackDelta_free f V W) P Q

/-- ★★★★★★★★**引き戻しはテンソル積を保つ** `f^*(P ⊗ Q) ≅ f^*P ⊗ f^*Q`。 -/
noncomputable def pullbackTensorIso (P Q : Y.PresheafOfModules) :
    (pullbackPre f).obj (P ⊗ Q) ≅ (pullbackPre f).obj P ⊗ (pullbackPre f).obj Q :=
  haveI := isIso_pullbackDelta f P Q
  asIso (pullbackDelta f P Q)

/-! ## ★出典の紐付け(`.src`) -/

def delta_free_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元の上で δ が同型の合成に等しいこと)",
    sectionId := "genell-def-1-1-i" }

def isIso_pullbackDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻しがテンソル積を保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
