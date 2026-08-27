/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.Modules.Sheaf
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Monoidal
import Mathlib.Algebra.Category.ModuleCat.Sheaf.LocallyFree
import ABC3.Meta.Claim

/-!
# スキーム上の**加群層のテンソル積**と**可逆層**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★`ample-and-projective-embedding` の段 A（回避経路）

台帳は段 A を「層のテンソル積（`SheafOfModules` の monoidal 構造）」としていた。
★2026-08-27 の実測で、mathlib には

* `PresheafOfModules.monoidalCategory` が**ある**（前層の側は完成）
* `Monoidal.transport` は**圏同値に沿ってのみ**
* `Localization.Monoidal` / `Monoidal.Reflective` / `Adjunction.monoidal` は**無い**

と判った。★★したがって「層化に沿って monoidal 構造が降りる」を作るのは
**mathlib への PR 規模**である。

★★★**回避経路**: 全 monoidal 構造（結合子・単位子・五角形・三角形）を作らず、
**`tensorObj` だけを定義して `L^{⊗n}` を再帰で作る**。
`ample` の定義（[Stacks] 28.27.1）にはそれで足りる。本ファイルはそれを取る。

## ★★★★★★機構は 3 本の合成

    X.Modules --forget--> X.PresheafOfModules --⊗--> X.PresheafOfModules
      --sheafification--> X.Modules

★`X.ringCatSheaf.obj = X.presheaf ⋙ forget₂ CommRingCat RingCat` は **`rfl`** なので
（2026-08-27 実測）、`PresheafOfModules.Monoidal.tensorObj` がそのまま当たる。
★★層化は `PresheafOfModules.sheafification (𝟙 _)` である
——構造前層は既に層なので `α = 𝟙` が局所全単射になる。

## ★★★可逆層は「逆元を持つ」で定義する

[Stacks] 01CR の定義（`M ⊗ N ≅ 𝒪_X` なる `N` が存在する）を採る。
★階数 1 の局所自由層という定義（mathlib の `IsLocallyFree` を絞る道）もあるが、
**逆元の存在のほうが `ample` の消費に近い**。

## ★残っている段（明示）

★**結合子・単位子は作っていない**ので、`L^{⊗n}` は右結合の固定した再帰である。
★★`ample` の定義には `X_s`（切断の非消失軌跡）が要り、それは本ファイルには無い。
★★★段 B（`Pic(X)`）には monoidal 構造の全体が要る。
-/

namespace ABC3.Found.GenEll

open CategoryTheory AlgebraicGeometry

/-! ## ★★★★★★★★加群層のテンソル積 -/

/-- ★★★★★★★★**スキーム上の加群層のテンソル積**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★前層のテンソル積（mathlib の `PresheafOfModules.Monoidal.tensorObj`）を層化する。
★★`X.ringCatSheaf.obj = X.presheaf ⋙ forget₂ CommRingCat RingCat` が `rfl` なので
型はそのまま合う。 -/
noncomputable def schemeTensorObj (X : Scheme.{0}) (M N : X.Modules) : X.Modules :=
  (PresheafOfModules.sheafification (R := X.ringCatSheaf) (𝟙 X.ringCatSheaf.obj)).obj
    (PresheafOfModules.Monoidal.tensorObj
      ((Scheme.Modules.toPresheafOfModules X).obj M)
      ((Scheme.Modules.toPresheafOfModules X).obj N))

/-- ★★★★★★**テンソル冪 `L^{⊗n}`**（右結合の再帰）。

★★`ample` の定義（[Stacks] 28.27.1）が要求するのはこれだけである。
★結合子を作っていないので **`L^{⊗(m+n)} ≅ L^{⊗m} ⊗ L^{⊗n}` は本ファイルには無い**。 -/
noncomputable def schemeTensorPow (X : Scheme.{0}) (M : X.Modules) : ℕ → X.Modules
  | 0 => SheafOfModules.unit X.ringCatSheaf
  | n + 1 => schemeTensorObj X M (schemeTensorPow X M n)

@[simp] theorem schemeTensorPow_zero (X : Scheme.{0}) (M : X.Modules) :
    schemeTensorPow X M 0 = SheafOfModules.unit X.ringCatSheaf := rfl

@[simp] theorem schemeTensorPow_succ (X : Scheme.{0}) (M : X.Modules) (n : ℕ) :
    schemeTensorPow X M (n + 1) = schemeTensorObj X M (schemeTensorPow X M n) := rfl

/-! ## ★★★可逆層 -/

/-- ★★★★★**可逆層** —— 逆元を持つ加群層（[Stacks] 01CR の定義）。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★階数 1 の局所自由層という定義（mathlib の `SheafOfModules.IsLocallyFree` を絞る道）
もあるが、**逆元の存在のほうが `ample` の消費に近い**。 -/
def IsInvertibleModule (X : Scheme.{0}) (M : X.Modules) : Prop :=
  ∃ N : X.Modules, Nonempty (schemeTensorObj X M N ≅ SheafOfModules.unit X.ringCatSheaf)

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (iv)` の証明が使う語彙（ample な直線束）を作るための**配管**である。
`ample` の定義そのものではない。 -/

def schemeTensorObj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(語彙——スキーム上の加群層のテンソル積)",
    sectionId := "genell-prop-1-4" }

def schemeTensorPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(語彙——テンソル冪 L^{⊗n})",
    sectionId := "genell-prop-1-4" }

def IsInvertibleModule.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(語彙——可逆層。ample の定義は含まない)",
    sectionId := "genell-prop-1-4" }

def schemeTensorObj.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "PresheafOfModules.Monoidal.tensorObj(前層のテンソル積)"
      (.inMathlib "PresheafOfModules.Monoidal.tensorObj") 7,
    .citation "[mathlib]" "PresheafOfModules.sheafification(層化)"
      (.inMathlib "PresheafOfModules.sheafification") 7,
    .citation "[mathlib]" "SheafOfModules.unit(構造層を加群層と見る)"
      (.inMathlib "SheafOfModules.unit") 7,
    .implicitStep
      ("★★★★台帳の段 A は「層化に沿って monoidal 構造が降りる」であり " ++
       "mathlib への PR 規模である(Localization.Monoidal は無い、2026-08-27 実測)。" ++
       "★本ファイルは**回避経路**——tensorObj だけを定義して L^{⊗n} を再帰で作る") 7,
    .implicitStep
      ("★★残っているもの: 結合子・単位子を作っていないので " ++
       "L^{⊗(m+n)} ≅ L^{⊗m} ⊗ L^{⊗n} は無い。" ++
       "★ample の定義には X_s(切断の非消失軌跡)も要る。" ++
       "★★段 B(Pic(X))には monoidal 構造の全体が要る") 7 ]

end ABC3.Found.GenEll
