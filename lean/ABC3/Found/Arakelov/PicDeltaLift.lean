import ABC3.Found.Arakelov.PicFreeYonedaLift
import ABC3.Found.Arakelov.PicPullTensorFree

/-!
# Arakelov (B1) 第 30 ブロック —— **`δ` の同型性を生成元の 1 点に帰着させる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★二重の持ち上げ

`δ` は 2 変数の比較射である。第 29 ブロックの器具は 1 変数用なので、**2 回**当てる:

| 回 | 固定するもの | 動かすもの | 自然性 |
|---|---|---|---|
| 1 | `P = free (yoneda V)` | `Q` | `δ_natural_right` |
| 2 | `Q` は任意 | `P` | `δ_natural_left` |

★★★これで

    「生成元 2 つの上で `δ` が同型」⟹「すべての `P, Q` で `δ` が同型」

が出る。★残るのは**生成元 2 つの 1 点**だけである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `deltaNatRight` | ★`Q` を動かしたときの自然変換 |
| `deltaNatLeft` | ★`P` を動かしたときの自然変換 |
| `isIso_pullbackDelta_of_free` | ★★★★★**生成元 2 つで同型なら全対象で同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-! ## ★自然変換としての梱包 -/

/-- ★**`Q` を動かしたときの `δ`** —— `f^*(P ⊗ -) ⟶ f^*P ⊗ f^*(-)`。 -/
noncomputable def deltaNatRight (P : Y.PresheafOfModules) :
    tensorLeft P ⋙ pullbackPre f ⟶ pullbackPre f ⋙ tensorLeft ((pullbackPre f).obj P) where
  app Q := pullbackDelta f P Q
  naturality _ _ g := (Functor.OplaxMonoidal.δ_natural_right (pullbackPre f) P g).symm

/-- ★**`P` を動かしたときの `δ`** —— `f^*(- ⊗ Q) ⟶ f^*(-) ⊗ f^*Q`。 -/
noncomputable def deltaNatLeft (Q : Y.PresheafOfModules) :
    tensorRight Q ⋙ pullbackPre f ⟶ pullbackPre f ⋙ tensorRight ((pullbackPre f).obj Q) where
  app P := pullbackDelta f P Q
  naturality _ _ g := (Functor.OplaxMonoidal.δ_natural_left (pullbackPre f) g Q).symm

/-! ## ★★★★★二重の持ち上げ -/

/-- ★★★★★★**生成元 2 つで同型なら、すべての `P, Q` で同型**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `δ` の同型性は**生成元 2 つの 1 点**だけに帰着した。 -/
theorem isIso_pullbackDelta_of_free
    (h : ∀ V W : Y.Opens, IsIso (pullbackDelta f
      ((PresheafOfModules.free (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
        (yoneda.obj V))
      ((PresheafOfModules.free (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
        (yoneda.obj W))))
    (P Q : Y.PresheafOfModules) : IsIso (pullbackDelta f P Q) := by
  haveI hpull : PreservesColimitsOfSize.{u, u} (pullbackPre f) := pullbackPre_preservesColimits f
  -- ★第 1 段: `P` を生成元に固定して `Q` を動かす
  have step1 : ∀ (V : Y.Opens) (Q : Y.PresheafOfModules),
      IsIso (pullbackDelta f
        ((PresheafOfModules.free (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
          (yoneda.obj V)) Q) := by
    intro V Q
    letI PV := (PresheafOfModules.free
      (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj (yoneda.obj V)
    exact isIso_app_of_freeYoneda (deltaNatRight f PV)
      (preservesColimitsOfSize_comp (tensorLeft PV) (pullbackPre f)
        (tensorLeft_preservesColimits (Y := Y) PV) hpull)
      (preservesColimitsOfSize_comp (pullbackPre f)
        (tensorLeft ((pullbackPre f).obj PV)) hpull
        (tensorLeft_preservesColimits (Y := X) ((pullbackPre f).obj PV)))
      (fun W => h V W) Q
  -- ★★第 2 段: `Q` を任意にして `P` を動かす
  exact isIso_app_of_freeYoneda (deltaNatLeft f Q)
    (preservesColimitsOfSize_comp (tensorRight Q) (pullbackPre f)
      (tensorRight_preservesColimits (Y := Y) Q) hpull)
    (preservesColimitsOfSize_comp (pullbackPre f)
      (tensorRight ((pullbackPre f).obj Q)) hpull
      (tensorRight_preservesColimits (Y := X) ((pullbackPre f).obj Q)))
    (fun V => step1 V Q) P

/-! ## ★出典の紐付け(`.src`) -/

def isIso_pullbackDelta_of_free.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——δ の同型性を生成元に帰着させること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
