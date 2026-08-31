/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackFunctorial
import ABC3.Found.Arakelov.PicSheafifyTrivial
import ABC3.Meta.Claim

/-!
# **層化に沿った自明化の輸送**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★台帳 `arakelov-degF-finite-places` の**段 A** への一歩

無条件の `deg_F`（[Szp] を仮定しない形）に要るのは

    `ADiv(F)/APrc(F) ≅ APic(Spec 𝒪_F)`

だけになった（`§9-746`・`§9-747`）。その入口が台帳の**段 A（アフィンの橋）**である:

> `Spec R` 上の局所自明な**前層**加群 `L` から**可逆 `R`-加群** `Γ(L,⊤)` を得る

★★`APicM` は**前層の水準**なので `Γ(L,⊤)` は何も保証されない
（局所自明でも層とは限らない）。★★★したがって段 A の中身は
**層化に沿って計量を移す**ことである。

| 段 | 内容 | 状態 |
|---|---|---|
| A1 | 層化は局所自明性を保つ | ✅ 在庫（`isLocallyTrivial_sheafify`） |
| A2 | **計量も移る**（自明化・遷移単元の輸送） | ★★**本ファイルはその道具** |
| A3 | `Γ(L^sh,⊤)` が可逆 `R`-加群（`equivPicRingSheaf` へ） | ⬜ |

## ★★★★★★★★輸送そのものは在庫の証明の中にあった

`isLocallyTrivial_sheafify` の証明の最後の行が

    `(asIso ((restrictPresheafFunctor X V).map (sheafifyUnit X P))).symm ≪≫ e`

である。★本ファイルはそれに **`sheafifyTriv` と名前を付ける**
——`PullbackTriv.lean` が `isLocallyTrivial_pullbackPre` の中身に名前を付けたのと同じ手である。

## ★★★★★★遷移単元は層化で**変わらない**

    `transUnit (P^sh) V (sheafifyTriv P e) (sheafifyTriv P e') = transUnit P V e e'`

★機構は本ファイルの一般補題 `transUnit_isoComp`:

    `transUnit A V (ψ ≪≫ e) (ψ ≪≫ e') = transUnit B V e e'`

——**共通の `ψ` は遷移単元から消える**。★★これは `transUnit_pullTriv`
（`AMetricIso.lean`）の一般化であり、`pullTriv` が `mapIso φ ≪≫ e` の形だったのに対し、
層化では `ψ` が `V` ごとの `asIso` なので一般形が要る。

## ★`IsIso` インスタンスの扱い

`sheafifyTriv P e` の `asIso` は `e` から導出したインスタンスに依存するので、
`e` と `e'` で**項が異なる**。★`IsIso` は `Subsingleton` なので
`sheafifyTrivOf`（インスタンスを引数に取る形）へ落として揃える
（`sheafifyTriv_eq_of`）。★★これが `tools\lean-idioms.md` 向けの新しい失敗形である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

/-! ## ★★★★★★共通の同型は遷移単元から消える -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★自明化に共通の同型を前置しても、座標の読みは前置分だけずれる。★**`rfl`**。 -/
theorem trivEquiv_isoComp {X : Scheme.{0}} {A B : X.PresheafOfModules} {V : X.Opens}
    (ψ : (restrictPresheafFunctor X V).obj A ≅ (restrictPresheafFunctor X V).obj B)
    (e : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (x : (A.obj (op V) : Type)) :
    trivEquiv A V (ψ ≪≫ e) x = trivEquiv B V e (ψ.hom.app (op (Over.mk (𝟙 V))) x) := rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★★同上——生成切断は `ψ.inv` で戻る。 -/
theorem trivGen_isoComp {X : Scheme.{0}} {A B : X.PresheafOfModules} {V : X.Opens}
    (ψ : (restrictPresheafFunctor X V).obj A ≅ (restrictPresheafFunctor X V).obj B)
    (e : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivGen A V (ψ ≪≫ e) = ψ.inv.app (op (Over.mk (𝟙 V))) (trivGen B V e) := by
  apply (trivEquiv A V (ψ ≪≫ e)).injective
  show (trivEquiv A V (ψ ≪≫ e)) ((trivEquiv A V (ψ ≪≫ e)).symm 1) = _
  rw [LinearEquiv.apply_symm_apply,
    trivEquiv_isoComp ψ e (ψ.inv.app (op (Over.mk (𝟙 V))) (trivGen B V e))]
  have hiv : ψ.hom.app (op (Over.mk (𝟙 V))) (ψ.inv.app (op (Over.mk (𝟙 V))) (trivGen B V e))
      = trivGen B V e :=
    congrArg (fun (χ : (restrictPresheafFunctor X V).obj B
      ⟶ (restrictPresheafFunctor X V).obj B) => χ.app (op (Over.mk (𝟙 V))) (trivGen B V e))
      ψ.inv_hom_id
  rw [hiv]
  show (1 : (Γ(X, V) : Type)) = (trivEquiv B V e) ((trivEquiv B V e).symm 1)
  rw [LinearEquiv.apply_symm_apply]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**共通の同型は遷移単元から消える**。

★これが `transUnit_pullTriv`（`AMetricIso.lean`）の一般化である
——`pullTriv` は `mapIso φ ≪≫ e` の形だったが、層化では `ψ` が
`V` ごとの `asIso` なので一般形が要る。 -/
theorem transUnit_isoComp {X : Scheme.{0}} {A B : X.PresheafOfModules} {V : X.Opens}
    (ψ : (restrictPresheafFunctor X V).obj A ≅ (restrictPresheafFunctor X V).obj B)
    (e e' : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit A V (ψ ≪≫ e) (ψ ≪≫ e') = transUnit B V e e' := by
  rw [transUnit_eq_trivEquiv_trivGen, transUnit_eq_trivEquiv_trivGen, trivGen_isoComp,
    trivEquiv_isoComp ψ e' (ψ.inv.app (op (Over.mk (𝟙 V))) (trivGen B V e))]
  congr 1
  exact congrArg (fun (χ : (restrictPresheafFunctor X V).obj B
    ⟶ (restrictPresheafFunctor X V).obj B) => χ.app (op (Over.mk (𝟙 V))) (trivGen B V e))
    ψ.inv_hom_id

/-! ## ★★★★★層化に沿った自明化の輸送 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★**自明化があれば制限は層である**（`isLocallyTrivial_sheafify` の中身に名前を付けた）。 -/
theorem isSheaf_restrict_of_triv (X : Scheme.{0}) (V : X.Opens) (P : X.PresheafOfModules)
    (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙ PresheafOfModules.presheaf P) := by
  have hunit : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙
        (PresheafOfModules.unit X.ringCatSheaf.obj).presheaf) :=
    isSheaf_restrict X V _
      ((sheafCompose (Opens.grothendieckTopology X)
        (forget₂ RingCat.{0} AddCommGrpCat.{0})).obj X.ringCatSheaf).property
  exact (Presheaf.isSheaf_of_iso_iff ((PresheafOfModules.toPresheaf _).mapIso e)).2 hunit

set_option backward.isDefEq.respectTransparency false in
/-- ★**インスタンスを引数に取る形**の輸送（`IsIso` を揃えるため）。 -/
noncomputable def sheafifyTrivOf {X : Scheme.{0}} (P : X.PresheafOfModules) {V : X.Opens}
    (inst : IsIso ((restrictPresheafFunctor X V).map (sheafifyUnit X P)))
    (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    (restrictPresheafFunctor X V).obj ((sheafifyFunctor X).obj P).val
      ≅ 𝟙_ (PresheafModulesOn X V) :=
  (@asIso _ _ _ _ _ inst).symm ≪≫ e

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**層化に沿った自明化の輸送** `(P^sh)|_V ≅ 𝟙_`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★在庫の `isLocallyTrivial_sheafify` の証明の最後の行に名前を付けたものである
（`PullbackTriv.lean` が `isLocallyTrivial_pullbackPre` にしたのと同じ手）。 -/
noncomputable def sheafifyTriv {X : Scheme.{0}} (P : X.PresheafOfModules) {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    (restrictPresheafFunctor X V).obj ((sheafifyFunctor X).obj P).val
      ≅ 𝟙_ (PresheafModulesOn X V) :=
  sheafifyTrivOf P (isIso_restrictMap_sheafifyUnit X V P (isSheaf_restrict_of_triv X V P e)) e

set_option backward.isDefEq.respectTransparency false in
/-- ★`IsIso` は `Subsingleton` なので、どのインスタンスで作っても同じ。 -/
theorem sheafifyTriv_eq_of {X : Scheme.{0}} (P : X.PresheafOfModules) {V : X.Opens}
    (inst : IsIso ((restrictPresheafFunctor X V).map (sheafifyUnit X P)))
    (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    sheafifyTriv P e = sheafifyTrivOf P inst e := by
  show sheafifyTrivOf P (isIso_restrictMap_sheafifyUnit X V P
    (isSheaf_restrict_of_triv X V P e)) e = _
  congr 1

/-! ## ★★★★★★遷移単元は層化で変わらない -/

set_option backward.isDefEq.respectTransparency false in
theorem transUnit_sheafifyTrivOf {X : Scheme.{0}} (P : X.PresheafOfModules) {V : X.Opens}
    (inst : IsIso ((restrictPresheafFunctor X V).map (sheafifyUnit X P)))
    (e e' : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit ((sheafifyFunctor X).obj P).val V (sheafifyTrivOf P inst e)
        (sheafifyTrivOf P inst e')
      = transUnit P V e e' :=
  transUnit_isoComp (@asIso _ _ _ _ _ inst).symm e e'

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★**遷移単元は層化で変わらない**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これが「計量を層化に沿って移す」段（台帳の段 A2）の核である
——計量が使うのは遷移単元だけだからである。 -/
theorem transUnit_sheafifyTriv {X : Scheme.{0}} (P : X.PresheafOfModules) {V : X.Opens}
    (e e' : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit ((sheafifyFunctor X).obj P).val V (sheafifyTriv P e) (sheafifyTriv P e')
      = transUnit P V e e' := by
  haveI inst := isIso_restrictMap_sheafifyUnit X V P (isSheaf_restrict_of_triv X V P e)
  rw [sheafifyTriv_eq_of P inst e, sheafifyTriv_eq_of P inst e']
  exact transUnit_sheafifyTrivOf P inst e e'

/-! ## ★★★★輸送は制限と可換である -/

set_option backward.isDefEq.respectTransparency false in
/-- ★制限は層化の unit と可換である。★**`rfl`**。 -/
theorem restrictOn_map_sheafifyUnit {X : Scheme.{0}} (P : X.PresheafOfModules) {V W : X.Opens}
    (hWV : W ≤ V) :
    (restrictOnFunctor hWV).map ((restrictPresheafFunctor X V).map (sheafifyUnit X P))
      = (restrictPresheafFunctor X W).map (sheafifyUnit X P) := rfl

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★**層化に沿った輸送は制限と可換である**。

★`rfl` では**ない**（`asIso` のインスタンスが `V` と `W` で別に導出されるから）。
`Functor.map_inv` と `restrictOn_map_sheafifyUnit` で潰す。 -/
theorem sheafifyTriv_restrict {X : Scheme.{0}} (P : X.PresheafOfModules) {V W : X.Opens}
    (hWV : W ≤ V) (e : (restrictPresheafFunctor X V).obj P ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivialOfLe ((sheafifyFunctor X).obj P).val hWV (sheafifyTriv P e)
      = sheafifyTriv P (trivialOfLe P hWV e) := by
  haveI instV := isIso_restrictMap_sheafifyUnit X V P (isSheaf_restrict_of_triv X V P e)
  haveI instW := isIso_restrictMap_sheafifyUnit X W P
    (isSheaf_restrict_of_triv X W P (trivialOfLe P hWV e))
  rw [sheafifyTriv_eq_of P instV e, sheafifyTriv_eq_of P instW (trivialOfLe P hWV e)]
  refine Iso.ext ?_
  show (restrictOnFunctor hWV).map
      (CategoryTheory.inv ((restrictPresheafFunctor X V).map (sheafifyUnit X P)) ≫ e.hom)
    = CategoryTheory.inv ((restrictPresheafFunctor X W).map (sheafifyUnit X P))
      ≫ (restrictOnFunctor hWV).map e.hom
  rw [Functor.map_comp, Functor.map_inv]
  congr 1

/-! ### ★出典の紐付け(`.src`) -/

def transUnit_isoComp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(共通の同型は遷移単元から消えること)",
    sectionId := "genell-def-1-1-i" }

def sheafifyTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層化に沿った自明化の輸送)",
    sectionId := "genell-def-1-1-i" }

def transUnit_sheafifyTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(遷移単元は層化で変わらないこと)",
    sectionId := "genell-def-1-1-i" }

def sheafifyTriv_restrict.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層化に沿った輸送が制限と可換であること)",
    sectionId := "genell-def-1-1-i" }

def transUnit_sheafifyTriv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isLocallyTrivial_sheafify(層化は局所自明性を保つ＝台帳の段 A1)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isLocallyTrivial_sheafify") 3,
    .citation "[ABC3]" "isIso_restrictMap_sheafifyUnit(自明化があれば層化の unit は同型)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isIso_restrictMap_sheafifyUnit") 3,
    .implicitStep
      ("★台帳 arakelov-degF-finite-places の**段 A2**(計量を層化に沿って移す)の道具である。" ++
       "段 A1(局所自明性の保存)は在庫、段 A3(Γ(L^sh,⊤) が可逆 R-加群)が残る") 3,
    .implicitStep
      ("★★`IsIso` インスタンスの扱い: sheafifyTriv の asIso は e から導出した" ++
       "インスタンスに依存するので e と e' で**項が異なる**。IsIso は Subsingleton なので " ++
       "sheafifyTrivOf(インスタンスを引数に取る形)へ落として揃える(sheafifyTriv_eq_of)") 3 ]

end ABC3.Found.Arakelov
