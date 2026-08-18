import ABC3.Found.Arakelov.PicTildePre

/-!
# Arakelov (B1) 第 84 ブロック —— **層化の普遍性で降ろす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★層の射になった

第 83 ブロックの前層射

    tilde M ⊗ tilde N  ⟶  tilde (M ⊗_R N)      (前層)

の**終域は層**なので、層化の普遍性で

    tensorModules (tilde M) (tilde N)  ⟶  tilde (M ⊗_R N)      (層)

に降ろせる。★mathlib の `PresheafOfModules.sheafificationHomEquiv` を使う。

★★`restrictScalars (𝟙 _)` は**何もしない**(`rfl`、実測)ので余計な運搬は要らない。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `restrictScalarsId_obj` | ★`restrictScalars (𝟙 _)` は恒等(`rfl`) |
| `tildeTensorDesc` | ★★★★**層の射** |

## ★★★次

    茎で同型(第 77・78・79 の器具)  ⟹  `tildeTensorDesc` は同型

    茎(左) = M_p ⊗_{R_p} N_p     ← 層化は茎を変えない(第 77 の前提)
    茎(右) = (M ⊗_R N)_p          ← 第 78
    一致は第 79 の `localizedTensorEquiv`
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace TensorProduct

variable (R : CommRingCat.{u}) (M N : ModuleCat.{u} (R : Type u))

/-- ★**`restrictScalars (𝟙 _)` は何もしない**(`rfl`)。 -/
theorem restrictScalarsId_obj (F : (Spec R).Modules) :
    (PresheafOfModules.restrictScalars (𝟙 (Spec R).ringCatSheaf.obj)).obj
      ((SheafOfModules.forget _).obj F) = F.val := rfl

/-- ★★★★**層の射** `tensorModules (tilde M) (tilde N) ⟶ tilde (M ⊗_R N)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 83 ブロックの前層射を層化の普遍性で降ろしたものである。 -/
noncomputable def tildeTensorDesc :
    tensorModules (tilde M) (tilde N)
      ⟶ tilde (ModuleCat.of (R : Type u) (M ⊗[(R : Type u)] N)) :=
  (PresheafOfModules.sheafificationHomEquiv (𝟙 (Spec R).ringCatSheaf.obj)).symm
    (tildeTensorPre R M N)

/-! ## ★出典の紐付け(`.src`) -/

def tildeTensorDesc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層化の普遍性で降ろした層の射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
