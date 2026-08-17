/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Otri
import ABC3.Found.FrdI.Prop44Gp

/-!
# [FrdI] `Proposition 4.4, (ii)` —— Ore の四角形と零因子準同型

★★**2026-08-18。「壁」と呼んでいた `Proposition 4.4, (ii)` を分解して、
最初の 2 つの葉をここで閉じる。**

## ★測り違いの訂正

2026-08-17 の `Gap/FrdI/Prop44.lean` には

> `preStepSpan` は `X ⟶ A`・`X ⟶ E` の span を与えるだけで、
> 必要な `A ⟶ E` 向きの co-angular pre-step は与えない

と書いた。★**見る条を間違えていた。** 効くのは (i)(b) の `preStepSpan` ではなく
**(iii)(d) の `coaPreOverEquiv`**、すなわち

  `(𝒞^coa-pre)_A ≃ Order(Φ(A))^opp`

である。★`Order(Φ(A))` は**前順序圏**なので射は高々 1 本、そして `Φ(A)` は
`δ ≼ δ + ε` により**有向**である。ゆえに `(𝒞^coa-pre)_A` は**余フィルター**で、
`δ_c + δ_s` に対応する対象から `(E,c)`・`(E,s)` の双方へ射が降りる。

★★★**したがって Ore の四角形は一般の Frobenioid で作れる。**
`Thm52Frob.lean` の `exists_ore_square` は model について因子を明示して作ったものだったが、
**model であることは要らなかった**。

## ★このファイルが閉じる葉(`ResearchPaper/frdi-decomposition.json` の `otricomm`)

| 葉 | 宣言 |
|---|---|
| `ore-square` | `exists_ore_square_coaPre` / `exists_ore_common` |
| `divgp-hom` | `otriDivGp_mul` / `otriDivGpHom` |

★残るのは `ker-eq-image`(層1)・`central-ext`(層2)・**`pairing-vanishes`(層3、本体)**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★1. Ore の四角形 —— (iii)(d) から

原文 (FrdI p.25):
> categories.
-/

include P G in
/-- ★★★**Ore の四角形**(一般の Frobenioid)。

`c`・`s : E ⟶ A` が co-angular pre-step なら、co-angular pre-step
`m`・`p : X ⟶ E` があって `p ≫ s = m ≫ c` となる。

★**道具は (iii)(d) `coaPreOverEquiv` ひとつ**である:
`(𝒞^coa-pre)_A ≃ Order(Φ(A))^opp` の**本質的全射性**で `δ_c + δ_s` に対応する
対象 `W` を取り、**充満性**で `W ⟶ (E,c)`・`W ⟶ (E,s)` を得る。
`Φ(A)` が `δ ≼ δ + ε` で有向であることがすべてを動かしている。

★これは `Thm52Frob.lean` の `exists_ore_square`(model 版、因子を明示して構成)の
**一般化**である —— model であることは要らなかった。 -/
theorem exists_ore_square_coaPre {E A : C} (c s : E ⟶ A)
    (hc : IsCoAngular P c) (hcp : IsPreStep P c)
    (hs : IsCoAngular P s) (hsp : IsPreStep P s) :
    ∃ (X : C) (m p : X ⟶ E), IsCoAngular P m ∧ IsPreStep P m ∧
      IsCoAngular P p ∧ IsPreStep P p ∧ p ≫ s = m ≫ c := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  haveI hcb : IsIso (P.Base c) := hcp.2
  haveI hsb : IsIso (P.Base s) := hsp.2
  -- ★`(𝒞^coa-pre)_A` の対象として見る
  let Zc : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨E⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨c, ⟨hc, hcp⟩⟩)
  let Zs : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨E⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨s, ⟨hs, hsp⟩⟩)
  -- ★対応する `Φ(A)` の元
  let dc : Φ.val (P.toElem.obj A).base := Φ.map (inv (P.Base c)) (P.Div c)
  let ds : Φ.val (P.toElem.obj A).base := Φ.map (inv (P.Base s)) (P.Div s)
  -- ★本質的全射性で `dc + ds` の対象を取る
  let W := (coaPreOverFunctor P A).objPreimage (Opposite.op (toOrderCat (dc + ds)))
  have hiso : (coaPreOverFunctor P A).obj W ≅ Opposite.op (toOrderCat (dc + ds)) :=
    (coaPreOverFunctor P A).objObjPreimageIso _
  -- ★`Φ(A)` の有向性
  have hle1 : (toOrderCat dc : OrderCat (Φ.val (P.toElem.obj A).base)) ≤ toOrderCat (dc + ds) :=
    ⟨ds, rfl⟩
  have hle2 : (toOrderCat ds : OrderCat (Φ.val (P.toElem.obj A).base)) ≤ toOrderCat (dc + ds) :=
    ⟨dc, add_comm ds dc⟩
  -- ★充満性で射を降ろす
  obtain ⟨f1, -⟩ := (coaPreOverFunctor P A).map_surjective
    (show (coaPreOverFunctor P A).obj W ⟶ (coaPreOverFunctor P A).obj Zc from
      hiso.hom ≫ (homOfLE hle1).op)
  obtain ⟨f2, -⟩ := (coaPreOverFunctor P A).map_surjective
    (show (coaPreOverFunctor P A).obj W ⟶ (coaPreOverFunctor P A).obj Zs from
      hiso.hom ≫ (homOfLE hle2).op)
  have hw1 : f1.left.hom ≫ c = W.hom.hom :=
    congrArg InducedWideCategory.Hom.hom (Over.w f1)
  have hw2 : f2.left.hom ≫ s = W.hom.hom :=
    congrArg InducedWideCategory.Hom.hom (Over.w f2)
  -- ★`rw` を使わない: `WideSubcategory` の型が構文的に一致しないため(既知の罠)
  exact ⟨W.left.obj, f1.left.hom, f2.left.hom, f1.left.property.1, f1.left.property.2,
    f2.left.property.1, f2.left.property.2, hw2.trans hw1.symm⟩

include P G in
/-- ★★**Ore の四角形の有限族版** —— `ι` 個の co-angular pre-step `f i : E ⟶ A` を
**同時に**共通の `w : X ⟶ A` へ持ち上げる。

★交換子の計算では `c`・`s`・`t` の 3 本を同じ `X` の上に載せる必要があるので、
2 本版では足りない。`δ := ∑ i, δ_i` に取れば `δ_i ≼ δ` がすべて同時に成り立つ。 -/
theorem exists_ore_common {ι : Type*} [Fintype ι] [DecidableEq ι] {E A : C}
    (f : ι → (E ⟶ A)) (hc : ∀ i, IsCoAngular P (f i)) (hp : ∀ i, IsPreStep P (f i)) :
    ∃ (X : C) (m : ι → (X ⟶ E)),
      (∀ i, IsCoAngular P (m i) ∧ IsPreStep P (m i)) ∧
      ∃ w : X ⟶ A, ∀ i, m i ≫ f i = w := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  haveI hb : ∀ i, IsIso (P.Base (f i)) := fun i => (hp i).2
  let Z : ι → Over (⟨A⟩ : WideSubcategory (coaPreProp P)) := fun i =>
    Over.mk (show (⟨E⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨f i, ⟨hc i, hp i⟩⟩)
  let d : ι → Φ.val (P.toElem.obj A).base := fun i => Φ.map (inv (P.Base (f i))) (P.Div (f i))
  let W := (coaPreOverFunctor P A).objPreimage (Opposite.op (toOrderCat (∑ i, d i)))
  have hiso : (coaPreOverFunctor P A).obj W ≅ Opposite.op (toOrderCat (∑ i, d i)) :=
    (coaPreOverFunctor P A).objObjPreimageIso _
  have hle : ∀ i, (toOrderCat (d i) : OrderCat (Φ.val (P.toElem.obj A).base))
      ≤ toOrderCat (∑ j, d j) := fun i =>
    ⟨∑ j ∈ Finset.univ.erase i, d j, Finset.add_sum_erase _ _ (Finset.mem_univ i)⟩
  choose g hg using fun i => (coaPreOverFunctor P A).map_surjective
    (show (coaPreOverFunctor P A).obj W ⟶ (coaPreOverFunctor P A).obj (Z i) from
      hiso.hom ≫ (homOfLE (hle i)).op)
  refine ⟨W.left.obj, fun i => (g i).left.hom,
    fun i => ⟨(g i).left.property.1, (g i).left.property.2⟩, W.hom.hom, fun i => ?_⟩
  exact congrArg InducedWideCategory.Hom.hom (Over.w (g i))

/-! ## ★2. `𝒪^▷(A^birat)` の零因子は準同型

★`biratDivGp` は `𝒞^birat` の射に `Φ^gp` の元を割り当てる(`Prop44Gp.lean`)。
`𝒪^▷(A^birat)` の元は底が恒等・次数が 1 なので、合成の公式

  `Div(f ≫ g) = Φ(Base f)(Div g) + deg(g) · Div(f)`

の 2 つの補正項がどちらも消え、**加法的**になる。 -/

variable {P G}

/-- ★**`𝒪^▷(A^birat)` 上で零因子は加法的**。

★`End` の積は `x * y = y ≫ x` である(mathlib の規約)ことに注意。 -/
theorem otriDivGp_mul (X : BiratCat P G) (x y : OTri (biratPre P G) X) :
    biratDivGp ((x * y : End X) : X ⟶ X)
      = biratDivGp ((x : End X) : X ⟶ X) + biratDivGp ((y : End X) : X ⟶ X) := by
  show biratDivGp (((y : End X) : X ⟶ X) ≫ ((x : End X) : X ⟶ X)) = _
  rw [biratDivGp_comp', gpMap_biratBase_of_baseIdentity y.2.1,
    show biratDeg ((x : End X) : X ⟶ X) = 1 from x.2.2]
  simp

/-- ★★**`𝒪^▷(A^birat) → Φ^gp(A)` はモノイド準同型**。

★★これが `Proposition 4.4, (ii)` の残る 1 条を測るための座標である ——
可換性の問題は、この準同型の**核**(`𝒪^×(A)` の像、これは中心に入る)と
**像**(`Φ^gp` の部分群、これはアーベル)の間の**中心拡大**の問題になる。 -/
noncomputable def otriDivGpHom (X : BiratCat P G) :
    OTri (biratPre P G) X →*
      Multiplicative (Gp (Φ.val (P.toElem.obj (biratDown P G X)).base)) where
  toFun x := Multiplicative.ofAdd (biratDivGp ((x : End X) : X ⟶ X))
  map_one' := by
    show Multiplicative.ofAdd (biratDivGp ((1 : End X) : X ⟶ X)) = 1
    rw [show ((1 : End X) : X ⟶ X) = 𝟙 X from rfl, biratDivGp_id]
    rfl
  map_mul' x y := by
    show Multiplicative.ofAdd (biratDivGp ((x * y : End X) : X ⟶ X))
      = Multiplicative.ofAdd _ * Multiplicative.ofAdd _
    rw [otriDivGp_mul]
    rfl

@[simp] theorem otriDivGpHom_apply (X : BiratCat P G) (x : OTri (biratPre P G) X) :
    otriDivGpHom X x = Multiplicative.ofAdd (biratDivGp ((x : End X) : X ⟶ X)) := rfl

/-! ## ★出典の紐付け

★★**条つきの `.src` である。** `Proposition 4.4, (ii)` はまだ閉じていない
(`pairing-vanishes` が残っている)ので、条なしの `.src` は書かない。
`tools/frdi-progress.mjs` はこれを進捗として数えない。 -/

def exists_ore_square_coaPre.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 85, item := "Proposition 4.4, (ii)",
    sectionId := "frdi-prop-4-4" }

def otriDivGpHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 85, item := "Proposition 4.4, (ii)",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
