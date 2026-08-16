import ABC3.Found.FrdI.Frobenioid

/-!
# [FrdI] pull-back 射を右へ寄せる —— 分解を合成するための道具

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.24。

原文 (FrdI p.24):
> equivalence of categories [cf. §0].

## ★この段の目的

`Definition 1.3, (iv), (a)` の分解 `φ = δ ≫ γ ≫ β ≫ α`(`α` が pull-back)は
**pull-back を最後に置く**。★したがって 2 つの射の分解を**合成**しようとすると

  `φ ≫ φ' = δ ≫ γ ≫ β ≫ (α ≫ δ' ≫ γ' ≫ β') ≫ α'`

の真ん中で **pull-back `α` が後ろの因子とぶつかる**。これを解くには
**`α` を右へ通す**必要がある。

## ★★道具は 2 つだけ

1. **pull-back の普遍性**(`IsPullBack` の定義そのもの) ——
   `Hom(X, Y) ≅ {(f : X ⟶ Z, b : X_𝒟 ⟶ Y_𝒟) | Base f = b ≫ Base α}`
2. **基底変換**(`Definition 1.3, (i), (c)` の `plBkEquiv`) ——
   `Over(A)^{plbk} ≃ Over(A_𝒟)`、すなわち**どの底の射に沿っても pull-back が取れる**

★★**これだけで押し出しが作れる**: `α ≫ ρ` の底 `Base (α ≫ ρ)` に沿って
`Z` を引き戻して `α̃ : Ỹ ⟶ Z` を作り、`α̃` の普遍性に `(f, b) := (α ≫ ρ, k)` を
入れれば `ρ̃ : Y ⟶ Ỹ` が**一意に**得られる。

★★★**しかも `Base ρ̃` は同型になる** —— 引き戻す底を `Base (α ≫ ρ)` に取ったから。
`degFr` と `Div` も落ちてくる:

| 量 | 値 |
|---|---|
| `Base ρ̃` | 同型(`plBkEquiv` の与える `k`) |
| `degFr ρ̃` | `degFr ρ` |
| `Div ρ̃` | `Φ.map (Base α) (Div ρ)` |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★手 1 —— pull-back の普遍性を使いやすい形に -/

include P in
/-- ★★**pull-back の普遍性**(持ち上げ) —— `IsPullBack` の定義の言い換え。

`f : X ⟶ Z` と底の射 `b` が `Base f = b ≫ Base α` で両立するなら、
**`g ≫ α = f` かつ `Base g = b` なる `g` が一意に存在する**。 -/
theorem IsPullBack.lift {Y Z : C} {α : Y ⟶ Z} (hα : IsPullBack P α) (X : C)
    (f : X ⟶ Z) (b : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base)
    (hb : P.Base f = b ≫ P.Base α) :
    ∃! g : X ⟶ Y, g ≫ α = f ∧ P.Base g = b := by
  obtain ⟨g, hg⟩ := (hα X).2 ⟨(f, b), hb⟩
  have hg' : (g ≫ α, P.Base g) = (f, b) := congrArg Subtype.val hg
  refine ⟨g, ⟨congrArg Prod.fst hg', congrArg Prod.snd hg'⟩, ?_⟩
  intro g' hg'2
  refine (hα X).1 ?_
  apply Subtype.ext
  show (g' ≫ α, P.Base g') = (g ≫ α, P.Base g)
  rw [hg'2.1, hg'2.2]
  exact hg'.symm

include P in
/-- ★**pull-back は底の上で忠実** —— 普遍性の**単射性**の側。 -/
theorem IsPullBack.hom_ext {Y Z : C} {α : Y ⟶ Z} (hα : IsPullBack P α) {X : C}
    (g g' : X ⟶ Y) (h1 : g ≫ α = g' ≫ α) (h2 : P.Base g = P.Base g') : g = g' := by
  refine (hα X).1 ?_
  apply Subtype.ext
  show (g ≫ α, P.Base g) = (g' ≫ α, P.Base g')
  rw [h1, h2]

include P in
/-- ★★**底が同型な pull-back は同型** —— 普遍性から直ちに逆射が作れる。

★★分解 `δ ≫ γ ≫ β ≫ α` において**全体の底が同型**なら
(`δ`・`γ`・`β` の底はいつも同型なので)`α` の底も同型になり、
したがって **`α` 自身が同型**になる。★合成した分解を正規形に戻すときに効く。 -/
theorem isIso_of_isPullBack_of_baseIso {Y Z : C} {α : Y ⟶ Z} (hα : IsPullBack P α)
    (hb : IsIso (P.Base α)) : IsIso α := by
  haveI := hb
  obtain ⟨g, ⟨hg1, hg2⟩, -⟩ :=
    IsPullBack.lift P hα Z (𝟙 Z) (inv (P.Base α)) (by rw [P.Base_id, IsIso.inv_hom_id])
  refine ⟨g, ?_, hg1⟩
  refine IsPullBack.hom_ext P hα _ _ ?_ ?_
  · rw [Category.assoc, hg1, Category.comp_id, Category.id_comp]
  · rw [P.Base_comp, hg2, IsIso.hom_inv_id, P.Base_id]

/-! ## ★手 2 —— どの底の射に沿っても pull-back が取れる -/

include P in
/-- ★★**基底変換** —— `Definition 1.3, (i), (c)` の圏同値 `plBkEquiv` の
**本質的全射性**だけを取り出したもの。

底の射 `w : W ⟶ Z_𝒟` に対し、pull-back `α̃ : Ỹ ⟶ Z` と**同型** `k : Ỹ_𝒟 ≅ W`
があって `Base α̃ = k.hom ≫ w`。 -/
theorem plBk_baseChange (F : FrobenioidCore P) (Z : C) {W : D}
    (w : W ⟶ (P.toElem.obj Z).base) :
    ∃ (Yt : C) (αt : Yt ⟶ Z) (k : (P.toElem.obj Yt).base ≅ W),
      IsPullBack P αt ∧ P.Base αt = k.hom ≫ w := by
  haveI := F.plBkEquiv Z
  set O : Over ((P.toElem.obj Z).base) := Over.mk w with hO
  set V : Over (⟨Z⟩ : PlBk P) := (plBkOverFunctor P Z).objPreimage O with hV
  obtain ⟨e⟩ : Nonempty ((plBkOverFunctor P Z).obj V ≅ O) :=
    ⟨(plBkOverFunctor P Z).objObjPreimageIso O⟩
  refine ⟨V.left.obj, V.hom.hom, ⟨e.hom.left, e.inv.left, ?_, ?_⟩, V.hom.property, ?_⟩
  · exact congrArg CommaMorphism.left e.hom_inv_id
  · exact congrArg CommaMorphism.left e.inv_hom_id
  · exact (Over.w e.hom).symm

/-! ## ★★★本体 —— pull-back を右へ寄せる -/

include P in
/-- ★★★**pull-back の押し出し** —— pull-back `α` の後ろに任意の `ρ` が来ても、
**`ρ̃ ≫ α̃` の形に組み替えられる**(`α̃` は pull-back、`Base ρ̃` は同型)。

★これが分解どうしの合成を可能にする。 -/
theorem plBk_shuffle (F : FrobenioidCore P) {Y B Z : C} {α : Y ⟶ B}
    (hα : IsPullBack P α) (ρ : B ⟶ Z) :
    ∃ (Yt : C) (ρt : Y ⟶ Yt) (αt : Yt ⟶ Z),
      α ≫ ρ = ρt ≫ αt ∧ IsPullBack P αt ∧ IsIso (P.Base ρt) := by
  obtain ⟨Yt, αt, k, hαt, hbase⟩ := plBk_baseChange P F Z (P.Base (α ≫ ρ))
  obtain ⟨ρt, ⟨hcomp, hbρ⟩, -⟩ :=
    IsPullBack.lift P hαt Y (α ≫ ρ) k.inv (by rw [hbase, ← Category.assoc, k.inv_hom_id,
      Category.id_comp])
  refine ⟨Yt, ρt, αt, hcomp.symm, hαt, ?_⟩
  rw [hbρ]
  infer_instance

include P in
/-- ★**押し出しの `degFr`** —— `α` は linear なので次数は変わらない。 -/
theorem plBk_shuffle_degFr (F : FrobenioidCore P) {Y Yt B Z : C} {α : Y ⟶ B}
    (hα : IsPullBack P α) {ρ : B ⟶ Z} {ρt : Y ⟶ Yt} {αt : Yt ⟶ Z}
    (hcomp : α ≫ ρ = ρt ≫ αt) (hαt : IsPullBack P αt) :
    P.degFr ρt = P.degFr ρ := by
  have h := congrArg P.degFr hcomp
  rw [P.degFr_comp, P.degFr_comp, (F.pullBackLB α hα).2, (F.pullBackLB αt hαt).2,
    mul_one, one_mul] at h
  exact h.symm

include P in
/-- ★★**押し出しの `Div`** —— `Div ρ̃ = Φ.map (Base α) (Div ρ)`。

★`α`・`α̃` はどちらも pull-back なので `Div = 0`・`degFr = 1`、
`Base ρ̃` は同型だが `Div α̃ = 0` なので効かない。 -/
theorem plBk_shuffle_div (F : FrobenioidCore P) {Y Yt B Z : C} {α : Y ⟶ B}
    (hα : IsPullBack P α) {ρ : B ⟶ Z} {ρt : Y ⟶ Yt} {αt : Yt ⟶ Z}
    (hcomp : α ≫ ρ = ρt ≫ αt) (hαt : IsPullBack P αt) :
    P.Div ρt = Φ.map (P.Base α) (P.Div ρ) := by
  have h := congrArg P.Div hcomp
  rw [P.Div_comp, P.Div_comp, (F.pullBackLB α hα).1.2, (F.pullBackLB αt hαt).1.2,
    (F.pullBackLB αt hαt).2] at h
  simpa using h.symm

end ABC3.Found.FrdI
