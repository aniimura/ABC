import ABC3.Found.Arakelov.PicDualVal

/-!
# Arakelov (B2) 第 178 ブロック —— **貼り合わせた射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★双対が層であることの第 3 歩

第 177 で `dualGlueVal` が**加法的・線型・自然**であることが出た。
本ブロックはそれを束ねて**射**にする。

## ★★機構は `PresheafOfModules.homMk`

★§9-201 で測った: 構造体 `PresheafOfModules.Hom` を直接組むと
**`Module` の綴りが 2 通り**になって instance が見つからない
(`(Over.forget V).op ⋙ 𝒪 ⋙ forget₂` と `𝒪 ⋙ forget₂` の 2 経路)。
★★mathlib の `homMk`(**アーベル群の前層の射 + 線型性**から作る)を使えば、
`Module` は現れず**加法性だけ**で組める——[[ring-instance-two-paths]] の
**新しい逃げ道**である。

| 逃げ道 | 使う場面 |
|---|---|
| `homMk` | ★★前層加群の射を組むとき(本ブロック) |
| 型付き恒等関数 | 値の橋(第 173) |
| `letI` + `inferInstanceAs` | 局所化の係数(第 157) |

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `dual_eq_at` | ★共通部分の下では `φ i` と `φ j` は一致する |
| `dualHom_ext` | ★双対の切断は各対象での値で決まる |
| `dualGlueHom` | ★★★★**貼り合わせた射** |
| `dualGlueHom_app` | ★値の展開(`rfl`) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} {F : X.PresheafOfModules} {ι : Type u} {U : ι → X.Opens}
  (φ : ∀ i, ((dualPresheaf F).obj (op (U i)) : Type u))

/-- ★**共通部分の下では `φ i` と `φ j` は一致する**。 -/
theorem dual_eq_at (hφ : TopCat.Presheaf.IsCompatible (dualPresheaf F).presheaf U φ)
    {A : X.Opens} (i j : ι) (hAi : A ≤ U i) (hAj : A ≤ U j) (y : (F.obj (op A) : Type u)) :
    ((φ i).app (op (Over.mk (homOfLE hAi)))).hom y
      = ((φ j).app (op (Over.mk (homOfLE hAj)))).hom y :=
  congrArg (fun (ψ : ((dualPresheaf F).obj (op (U i ⊓ U j)) : Type u)) =>
    ((ψ.app (op (Over.mk (homOfLE (le_inf hAi hAj))))).hom y)) (hφ i j)

/-- ★**双対の切断は各対象での値で決まる**。 -/
theorem dualHom_ext {A : X.Opens} (a b : ((dualPresheaf F).obj (op A) : Type u))
    (h : ∀ (Z : Over A) (x : (F.obj (op Z.left) : Type u)),
      ((a.app (op Z))).hom x = ((b.app (op Z))).hom x) : a = b := by
  refine PresheafOfModules.Hom.ext ?_
  funext Z
  exact ModuleCat.hom_ext (LinearMap.ext (fun x => h Z.unop x))

variable (hφ : TopCat.Presheaf.IsCompatible (dualPresheaf F).presheaf U φ)

/-- ★★★★**貼り合わせた射**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`homMk` を使うので **`Module` の綴りの 2 経路**を踏まずに済む。 -/
noncomputable def dualGlueHom {V : X.Opens} (hV : V ≤ ⨆ i, U i) :
    ((dualPresheaf F).obj (op V) : Type u) :=
  PresheafOfModules.homMk
    { app := fun Z => AddCommGrpCat.ofHom
        (AddMonoidHom.mk' (fun x => dualGlueVal φ hφ ((leOfHom Z.unop.hom).trans hV) x)
          (fun x y => dualGlueVal_add φ hφ _ x y))
      naturality := by
        intro Z Z' f
        ext x
        exact (dualGlueVal_map φ hφ ((leOfHom Z.unop.hom).trans hV) (leOfHom f.unop.left) x).symm }
    (fun Z c x => dualGlueVal_smul φ hφ ((leOfHom Z.unop.hom).trans hV) c x)

/-- ★**貼り合わせた射の値**(`rfl`)。 -/
theorem dualGlueHom_app {V : X.Opens} (hV : V ≤ ⨆ i, U i) (Z : Over V)
    (x : (F.obj (op Z.left) : Type u)) :
    ((dualGlueHom φ hφ hV).app (op Z)).hom x
      = dualGlueVal φ hφ ((leOfHom Z.hom).trans hV) x := rfl

/-! ## ★出典の紐付け(`.src`) -/

def dualGlueHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——貼り合わせた射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
