/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.CategoryVocabulary
import Mathlib.FieldTheory.IntermediateField.Adjoin.Basic
import Mathlib.FieldTheory.Galois.Basic

/-!
# [FrdI] Example 6.1 の底の圏 `𝒟 = B(G)⁰`

原文 (FrdI p.109):
> connected objects of the Galois category B(G) [cf. §0] may be thought of as schemes

★★`Example 6.1` の底の圏 `𝒟` は **`B(G)⁰`**、すなわち
`G = Gal(K̄/K)` の連続作用を持つ**連結**な有限集合の圏である。
原文自身が「`Spec(L)`(`L ⊆ K̄` は `K` の有限拡大)と思ってよい」と書くので、
★**本ファイルはその形で作る** —— `K̄/K` の有限部分拡大を対象、`K`-代数射を射とする
圏 `FinSub` の**反対圏**である(`Spec` を取るので反対になる)。

## ★これで何が済むか

`Theorem 5.2, (ii)`(実装済み)が要求する底の圏の条件は
**連結**と **totally epimorphic** の 2 つである。本ファイルはその 2 つを与える:

* `finSubOp_isConnected` —— `⊥`(= `K` 自身)を経由するジグザグ
* `finSubOp_totallyEpimorphic` —— `FinSub` の射は体の射なので**すべて単射**、
  したがって mono。反対圏では epi になる。

## ★残り

`𝒟` が **of FSM-type** であることは別途要る(`Φ^birat` を monoid on `𝒟` に
するために足した仮定の、幾何での根拠)。★中身は
「`FinSub` の epi は同型」——`M/f(L)` が非自明な有限分離拡大なら
相異なる 2 つの埋め込みが取れる、という主張であり、分離性の在庫を要する。
★節点 `galois-cat-B0` に残す。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe u

/-- ★★**`B(G)⁰` の対象** —— `K̄/K` の有限部分拡大。

★原文 (FrdI p.109) の「`Spec(L)`, where `L ⊆ K̄` is a finite extension of `K`」。
★`K̄/K` が分離的なら `L/K` は自動的に分離的なので、構造には持たせない。 -/
structure FinSub (K : Type u) (Kbar : Type u) [Field K] [Field Kbar] [Algebra K Kbar] where
  /-- 中間体 -/
  toIF : IntermediateField K Kbar
  /-- 有限次であること -/
  fin : FiniteDimensional K toIF

variable {K Kbar : Type u} [Field K] [Field Kbar] [Algebra K Kbar]

/-- ★`FinSub` の圏構造(射は `K`-代数射)。 -/
instance finSubCategory : Category (FinSub K Kbar) where
  Hom L M := L.toIF →ₐ[K] M.toIF
  id L := AlgHom.id K L.toIF
  comp f g := g.comp f

/-- ★圏の射を `K`-代数射として読む。 -/
def FinSub.hom {L M : FinSub K Kbar} (f : L ⟶ M) : L.toIF →ₐ[K] M.toIF := f

@[simp] theorem FinSub.hom_id (L : FinSub K Kbar) :
    FinSub.hom (𝟙 L) = AlgHom.id K L.toIF := rfl

@[simp] theorem FinSub.hom_comp {L M N : FinSub K Kbar} (f : L ⟶ M) (g : M ⟶ N) :
    FinSub.hom (f ≫ g) = (FinSub.hom g).comp (FinSub.hom f) := rfl

theorem FinSub.hom_ext {L M : FinSub K Kbar} {f g : L ⟶ M}
    (h : FinSub.hom f = FinSub.hom g) : f = g := h

/-- ★体の射は単射。 -/
theorem FinSub.hom_injective {L M : FinSub K Kbar} (f : L ⟶ M) :
    Function.Injective (FinSub.hom f) :=
  (FinSub.hom f).toRingHom.injective

/-- ★★**`FinSub` の射はすべて mono**(体の射は単射だから)。 -/
theorem finSub_mono {L M : FinSub K Kbar} (f : L ⟶ M) : Mono f := by
  refine ⟨fun {Z} g h e => ?_⟩
  have he : FinSub.hom (g ≫ f) = FinSub.hom (h ≫ f) := congrArg FinSub.hom e
  rw [FinSub.hom_comp, FinSub.hom_comp] at he
  refine FinSub.hom_ext (AlgHom.ext fun x => FinSub.hom_injective f ?_)
  exact congrArg (fun t : Z.toIF →ₐ[K] M.toIF => t x) he

/-- ★`K` 自身(`⊥`)は `FinSub` の対象。 -/
def botSub (K Kbar : Type u) [Field K] [Field Kbar] [Algebra K Kbar] : FinSub K Kbar :=
  ⟨⊥, inferInstance⟩

/-- ★`⊥` からはどこへでも射がある(`⊥ ≃ₐ K` と `Algebra.ofId`)。 -/
noncomputable def botHom (L : FinSub K Kbar) : botSub K Kbar ⟶ L :=
  show (⊥ : IntermediateField K Kbar) →ₐ[K] L.toIF from
    (Algebra.ofId K L.toIF).comp (IntermediateField.botEquiv K Kbar).toAlgHom

/-- ★★`FinSub` は連結(`⊥` を経由するジグザグ)。 -/
instance finSub_isConnected : IsConnected (FinSub K Kbar) := by
  haveI : Nonempty (FinSub K Kbar) := ⟨botSub K Kbar⟩
  refine zigzag_isConnected ?_
  intro L M
  exact (Zigzag.of_inv (botHom L)).trans (Zigzag.of_hom (botHom M))

/-- ★★★**`𝒟 = B(G)⁰` は連結**。 -/
instance finSubOp_isConnected : IsConnected (FinSub K Kbar)ᵒᵖ := by
  haveI : Nonempty (FinSub K Kbar)ᵒᵖ := ⟨Opposite.op (botSub K Kbar)⟩
  refine zigzag_isConnected ?_
  intro L M
  exact (Zigzag.of_hom ((botHom L.unop).op)).trans (Zigzag.of_inv ((botHom M.unop).op))

/-- ★★★**`𝒟 = B(G)⁰` は totally epimorphic** —— `FinSub` の射がすべて mono だから。 -/
theorem finSubOp_totallyEpimorphic :
    IsTotallyEpimorphic (FinSub K Kbar)ᵒᵖ := by
  intro A B f
  haveI := finSub_mono f.unop
  refine ⟨fun {Z} g h e => ?_⟩
  have h1 : g.unop ≫ f.unop = h.unop ≫ f.unop := congrArg Quiver.Hom.unop e
  exact Quiver.Hom.unop_inj ((cancel_mono f.unop).mp h1)

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の底の圏 `𝒟 = B(G)⁰`(連結・totally epimorphic)。 -/
def finSubOp_totallyEpimorphic.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 底の圏 𝒟 = B(G)⁰ は連結かつ totally epimorphic",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.FrdI
