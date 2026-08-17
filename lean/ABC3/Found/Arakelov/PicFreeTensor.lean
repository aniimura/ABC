import ABC3.Found.Arakelov.PicPullFree
import Mathlib.Algebra.Category.ModuleCat.Adjunctions
import Mathlib.CategoryTheory.Monoidal.FunctorCategory

/-!
# Arakelov (B1) 第 25 ブロック —— **自由前層加群のテンソルは自由**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`Ring` インスタンスの二重路を抜けた

`↑((R ⋙ forget₂ CommRingCat RingCat).obj Z)` には `Ring` が **2 経路**で付く:

| 経路 | 出どころ | 使う側 |
|---|---|---|
| (a) | `RingCat` の構造 | `freeObj` / `restrictScalars` |
| (b) | `CommRing.toRing` | 前層加群のテンソル積 |

★★defeq だが構文が違うため、**`simp`/`rw` が発火しない**
——`CategoryTheory.comp_apply` は `@[simp]` なのに当たらない。

★★★★**抜け方は 2 つ組み合わせる**:

1. **加群の段で naturality を書かない**——`freeObjDesc` に渡せば
   mathlib が naturality を証明してくれる。★型前層の段なら等式は**要素の等式**である。
2. **`erw` を使う**——`rw` は `instances` 透明度で照合するので落ちるが、
   `erw` は `default` 透明度なので**二重路を越えられる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeTensorTypeHom` | ★型前層の段での比較射(naturality を要素で証明) |
| `freeTensorHom` | ★★`free (F ⊗ G) ⟶ free F ⊗ free G` |
| `freeTensorHom_app_freeMk` | ★★★生成元での値 |
| `freeTensorIso` | ★★★★★**`free (F ⊗ G) ≅ free F ⊗ free G`** |

★★★★★これで `δ` の同型性の第 3 段(生成元の上で同型)が揃う——
第 24 ブロック(`f^*(free (yoneda V)) ≅ free (yoneda (f⁻¹V))`)と
第 23 ブロック(逆像は `⊓` を保つ)と噛み合う。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite Limits

variable {C : Type u} [Category.{u} C] {R : Cᵒᵖ ⥤ CommRingCat.{u}}
  (F G : Cᵒᵖ ⥤ Type u)

/-! ## ★型前層の段での比較射 -/

/-- ★**型前層の段での比較射**。

★★★ここで作るのが要点である——加群の段で書くと naturality が
「`Ring` インスタンスの二重路」に阻まれるが、
型前層の段なら naturality は**要素の等式**であり、`erw` で抜ける。 -/
noncomputable def freeTensorTypeHom :
    (F ⊗ G) ⟶ (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) F
        ⊗ PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) G).presheaf
      ⋙ forget Ab where
  app Z := ↾fun (q : (F ⊗ G).obj Z) => show
      ((PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) F).obj Z
        ⊗ (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) G).obj Z
        : ModuleCat.{u} _) from
    (ModuleCat.freeMk (show F.obj Z × G.obj Z from q).1
      : (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) F).obj Z)
      ⊗ₜ (ModuleCat.freeMk (show F.obj Z × G.obj Z from q).2
      : (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) G).obj Z)
  naturality := by
    intro Z Z' φ
    ext q
    dsimp
    erw [PresheafOfModules.Monoidal.tensorObj_map_tmul]
    erw [PresheafOfModules.freeObj_map, PresheafOfModules.freeObj_map,
      ModuleCat.freeDesc_apply, ModuleCat.freeDesc_apply]
    rfl

/-! ## ★★加群の段へ上げる -/

/-- ★★**自由前層加群のテンソルへの比較射**。

★naturality は `freeObjDesc`(mathlib)が担う。 -/
noncomputable def freeTensorHom :
    PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) (F ⊗ G)
      ⟶ PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) F
        ⊗ PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) G :=
  PresheafOfModules.freeObjDesc (freeTensorTypeHom F G)

/-- ★★★**生成元での値**——これが同型性の根拠になる。 -/
theorem freeTensorHom_app_freeMk (Z : Cᵒᵖ) (q : (F ⊗ G).obj Z) :
    (freeTensorHom (R := R) F G).app Z (ModuleCat.freeMk q)
      = (ModuleCat.freeMk (show F.obj Z × G.obj Z from q).1
          : (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) F).obj Z)
        ⊗ₜ (ModuleCat.freeMk (show F.obj Z × G.obj Z from q).2
          : (PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) G).obj Z) :=
  ModuleCat.freeDesc_apply _ _

/-! ## ★★★★★同型であること -/

/-- ★★★★**比較射は各切断で同型である**——`ModuleCat.free` の `μIso` そのものだからである。 -/
theorem isIso_freeTensorHom_app (Z : Cᵒᵖ) :
    IsIso ((freeTensorHom (R := R) F G).app Z) := by
  have h : (freeTensorHom (R := R) F G).app Z
      = (ModuleCat.FreeMonoidal.μIso
          (((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj Z : RingCat.{u}) : Type u)
          (F.obj Z) (G.obj Z)).inv := by
    refine ModuleCat.free_hom_ext ?_
    intro q
    erw [freeTensorHom_app_freeMk, ModuleCat.FreeMonoidal.μIso_inv_freeMk]
    rfl
  rw [h]
  infer_instance

/-- ★★★★★**自由前層加群のテンソルは自由である** `free (F ⊗ G) ≅ free F ⊗ free G`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `δ` の同型性の第 3 段(生成元の上で同型)の残り半分である。 -/
noncomputable def freeTensorIso :
    PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) (F ⊗ G)
      ≅ PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) F
        ⊗ PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) G :=
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map (freeTensorHom (R := R) F G)) := by
    rw [NatTrans.isIso_iff_isIso_app]
    intro Z
    haveI := isIso_freeTensorHom_app F G (R := R) Z
    exact inferInstanceAs
      (IsIso ((forget₂ _ AddCommGrpCat).map ((freeTensorHom (R := R) F G).app Z)))
  haveI : IsIso (freeTensorHom (R := R) F G) :=
    isIso_of_reflects_iso _ (PresheafOfModules.toPresheaf _)
  asIso (freeTensorHom (R := R) F G)

/-! ## ★出典の紐付け(`.src`) -/

def freeTensorIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自由前層加群のテンソルが自由になること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
