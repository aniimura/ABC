/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.BaseChangeSetup
import ABC3.Found.Arakelov.DegArithPre
import ABC3.Found.Arakelov.PicBaseChange
import ABC3.Meta.Claim

/-!
# `Γ(L)` は**係数環 `R` の上でも可逆**である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★これは何か——`R` と `Γ(Spec R, ⊤)` の摩擦を渡る

`§9-779` の `invertible_gammaPre` は `Γ(L)` が **`Γ(Spec R, ⊤)`-加群として**可逆であると言う。
★底変換では `Γ(L)` を **`R`-加群**（`gammaModPre`、`§9-786`）として扱うので、
そちらでも可逆であることが要る。

★★機構は在庫の `invertible_of_surjective_algebraMap`（`PicBaseChange.lean`、第 197 ブロック）である
——`algebraMap` が**全射**なら係数環を取り替えても可逆。
`Γ-Spec` 同型は全射なので、そのまま渡せる。

## ★配管の記録（`lean-idioms` 向け）

`gammaModPre` は `def` なので**インスタンス探索が中身を見ない**
——`Module ↑Γ(Spec R,⊤) ↑(gammaModPre R L)` が見つからない。
★そこで `ModuleCat.restrictScalars … |>.obj …` と**書き下した形**で証明し
（そこでは `ModuleCat.instModuleCarrierObjRestrictScalars` が効く）、
`gammaModPre` 版は `rfl` 相当で受け直す。
★★`attribute [reducible]` を後付けすることはできない
（`failed to set [local reducible]`、2026-08-28 実測）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`Γ(L)` は `R` の上でも可逆である**（書き下した形）。

★在庫の `invertible_of_surjective_algebraMap` に `Γ-Spec` 同型の全射性を渡す。 -/
theorem invertible_gammaRestrict (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Invertible (R : Type)
      (((ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj
        (L.carrier.sheaf.obj (op ⊤))) : Type) := by
  letI algI : Algebra (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type) (R : Type) :=
    (Scheme.ΓSpecIso R).hom.hom.toAlgebra
  haveI hinv : Module.Invertible (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type)
      (((ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj
        (L.carrier.sheaf.obj (op ⊤))) : Type) := invertible_gammaPre L
  haveI tower : IsScalarTower (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type) (R : Type)
      (((ModuleCat.restrictScalars (Scheme.ΓSpecIso R).inv.hom).obj
        (L.carrier.sheaf.obj (op ⊤))) : Type) := by
    constructor
    intro c t x
    show ((Scheme.ΓSpecIso R).inv.hom ((Scheme.ΓSpecIso R).hom.hom c * t)) • x
      = c • (((Scheme.ΓSpecIso R).inv.hom t) • x)
    rw [map_mul, mul_smul]
    congr 1
    exact congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m c) (Scheme.ΓSpecIso R).hom_inv_id
  exact invertible_of_surjective_algebraMap (gammaSpecRingEquiv R).surjective

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**`gammaModPre R L` は可逆 `R`-加群である**。

★★これで底変換の議論（`§9-785` の `surjective_of_map_surjective`）に渡せる
——そこで要る「可逆な因子」がこれである。 -/
theorem invertible_gammaModPre (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Invertible (R : Type) (gammaModPre R L.carrier.sheaf : Type) :=
  invertible_gammaRestrict R L

/-! ### ★出典の紐付け(`.src`) -/

def invertible_gammaModPre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(Γ(L) は係数環 R の上でも可逆である)",
    sectionId := "genell-def-1-1-ii" }

def invertible_gammaModPre.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "invertible_gammaPre(前層の大域切断は可逆、§9-779)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.invertible_gammaPre") 4,
    .citation "[ABC3]" "invertible_of_surjective_algebraMap(係数環を取り替えても可逆、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.invertible_of_surjective_algebraMap") 3,
    .implicitStep
      ("★配管の記録: gammaModPre は def なのでインスタンス探索が中身を見ない。" ++
       "ModuleCat.restrictScalars … |>.obj … と書き下した形で証明し、" ++
       "gammaModPre 版は rfl 相当で受け直した。" ++
       "★★attribute [reducible] の後付けはできない(2026-08-28 実測)") 4 ]

end ABC3.Found.Arakelov
