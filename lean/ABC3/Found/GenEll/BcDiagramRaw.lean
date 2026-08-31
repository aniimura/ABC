/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DescentUnique

/-!
# [GenEll] Remark 1.5.1 —— **生の `pullback` で書いた底変換図式**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★なぜ書き直すのか —— 配管の設計判断

`BaseChangeRatTower.lean` の `baseChangeRatTowerDiagram` は

    overRatTowerDiagram ⋙ Over.pullback f ⋙ Over.forget X

という**関手の合成**である。★随伴で極限を押すにはこれが正しい形だが、
**下流の計算では `rw`/`simp` が一切効かない**（2026-08-27 に 3 ブロック費やして判明）:

`(D f).obj i` は `pullback (over.obj i).hom f` と **defeq だが構文的に違う**ので、
`pullback.hom_ext` が生む目標は「外側は関手合成の型・内側は生の型」の混在になり、
★`Category.assoc` も `Category.comp_id` も `pullback.lift_fst` も発火しない。
（エラー末尾の `Note: The target expression is not type-correct under the
instances transparency level` が目印。`tools/lean-idioms.md` にも記録した。）

★★**したがって下流用には生の形の図式を別に持つ。**

## ★★★★★★次の設計 —— `Over` を経由しない直接構成

★★★`bcDiagram` の `IsLimit` は `Over` を経由せず**直に**取れる:

* 錐の頂点は `pullback (Spec ℚ ⟶ Spec ℤ) f`（＝ `X_ℚ`）
* 任意の錐 `(T, τ_i)` に対し
  - `T ⟶ X` は `τ_i ≫ pullback.snd`（★段に依らない）
  - `T ⟶ Spec ℚ` は `τ_i ≫ pullback.fst` が `specRatTowerDiagram` 上の錐をなすので
    `specRatTowerIsLimit.lift` で得る
  - 2 つを `pullback.lift` で束ねる
* 一意性は `pullback.hom_ext` ＋ `specRatTowerIsLimit.uniq`

★これなら**全部が生の `pullback` の言葉**なので `rw`/`simp` が効く。
★★`BaseChangeRatTower.lean` の随伴による `IsLimit` は**残す**——
2 つが同じ極限であることの確認に使える。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-- ★★**生の `pullback` で書いた底変換図式** `n ↦ X ×_ℤ Spec ℤ[1/n!]`。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★`Functor.ofOpSequence` を使うので `map_id` / `map_comp` は mathlib が持つ
（`tools/lean-idioms.md` の「`Subring.closure` を台にした型を圏論の構造に載せると
核が発散する」と同じ理由で、構造体を素朴に書かない）。

★★★**`obj` が構文的に `pullback` である**ことが要点である——
`BaseChangeRatTower.lean` の合成版と defeq だが、下流の `rw`/`simp` はこちらでしか効かない。 -/
noncomputable def bcDiagram {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) : ℕᵒᵖ ⥤ Scheme.{0} :=
  Functor.ofOpSequence
    (X := fun n => pullback (overRatTowerDiagram.obj (Opposite.op n)).hom f)
    (fun n => pullback.map _ _ _ _
      ((overRatTowerDiagram.map (homOfLE (Nat.le_succ n)).op).left) (𝟙 X) (𝟙 _)
      (by simpa using (Over.w (overRatTowerDiagram.map (homOfLE (Nat.le_succ n)).op)).symm)
      (by simp))

/-- ★合成版と生版の対象は同じもの（`rfl`）。 -/
theorem bcDiagram_obj {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (i : ℕᵒᵖ) :
    (bcDiagram f).obj i = pullback (overRatTowerDiagram.obj i).hom f := rfl

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** 本ファイルは下流の計算のための**器**であり、
`IsLimit` の直接構成（上の設計）と、その先の同型の降下・因子の降下・
conductor の一致が残っている。 -/

def bcDiagram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(生の pullback で書いた底変換図式——IsLimit は含まない)",
    sectionId := "genell-rem-1-5-1" }

def bcDiagram.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "baseChangeRatTowerIsLimit(合成版の極限——随伴で押したもの)"
      (.inProject "ABC3" "ABC3.Found.GenEll.baseChangeRatTowerIsLimit") 9,
    .citation "[mathlib]" "Functor.ofOpSequence(ℕᵒᵖ 上の図式を後続の射だけで作る)"
      (.inMathlib "CategoryTheory.Functor.ofOpSequence") 9,
    .implicitStep
      ("★★★★★配管の設計判断: baseChangeRatTowerDiagram は関手の合成なので " ++
       "(D f).obj i が pullback (over.obj i).hom f と defeq だが構文的に違い、" ++
       "下流で rw/simp が一切効かない(2026-08-27 に 3 ブロック費やして判明、" ++
       "tools/lean-idioms.md に記録)。そこで下流用に生の形を別に持つ") 9,
    .implicitStep
      ("★★次の設計: bcDiagram の IsLimit は Over を経由せず直に取れる。" ++
       "頂点は pullback (Spec ℚ ⟶ Spec ℤ) f。任意の錐 (T, τ_i) に対し " ++
       "T ⟶ X は τ_i ≫ pullback.snd(段に依らない)、T ⟶ Spec ℚ は " ++
       "τ_i ≫ pullback.fst が錐をなすので specRatTowerIsLimit.lift で得て、" ++
       "2 つを pullback.lift で束ねる。一意性は pullback.hom_ext ＋ uniq") 9 ]

end ABC3.Found.GenEll
