/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.CategoryVocabulary
import Mathlib.FieldTheory.IntermediateField.Adjoin.Basic
import Mathlib.FieldTheory.Galois.Basic
import Mathlib.FieldTheory.Galois.Infinite

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

## ★of FSM-type も揃った

`𝒟` が **of FSM-type** であることも本ファイルで示す(`finSubOp_isOfFSMType`)。★これは `Φ^birat` を monoid on `𝒟` にするために
足した仮定((B) として `index.html` に開示)が、**幾何の応用では本当に成り立つ**ことの根拠である。
★中身は「`FinSub` の epi は同型」で、無限 Galois 理論の `fixedField (fixingSubgroup L) = L` 1 本で出る。
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

/-! ## ★2. `of FSM-type`

★★`Φ^birat` を `monoid on 𝒟` にするために足した仮定(`IsOfFSMType 𝒟`、
`index.html` に (B) として開示)が、**幾何の応用では本当に成り立つ**ことの根拠である。

★`B(G)⁰` の言葉で言えば「推移的 `G`-集合の間の単射な `G`-写像は全単射」——
反対圏の mono は `FinSub` の epi であり、それが同型であることを示せばよい。 -/

section FSM

variable [IsGalois K Kbar]

/-- ★★★★**`FinSub` の epi は同型**。

★★`M/f(L)` が非自明なら、`K̄/K` が Galois なので `f(L)` を固定して `x` を動かす
`σ ∈ Gal(K̄/K)` が取れる(無限 Galois 理論の
`InfiniteGalois.fixedField_fixingSubgroup`)。
★`M` と `σ(M)` の合成体 `N` は有限次なので `FinSub` の対象であり、
`M ⟶ N` の相異なる 2 射が `f` の後で一致してしまう —— epi に反する。 -/
theorem finSub_isIso_of_epi {L M : FinSub K Kbar} (f : L ⟶ M) (hepi : Epi f) : IsIso f := by
  haveI := hepi
  set L₁ : IntermediateField K Kbar := ((M.toIF.val).comp (FinSub.hom f)).fieldRange with hL₁
  have hge : M.toIF ≤ L₁ := by
    by_contra hcon
    obtain ⟨x, hxM, hxL⟩ := SetLike.not_le_iff_exists.mp hcon
    have hfix : IntermediateField.fixedField L₁.fixingSubgroup = L₁ :=
      InfiniteGalois.fixedField_fixingSubgroup L₁
    have hx : x ∉ IntermediateField.fixedField L₁.fixingSubgroup := by rw [hfix]; exact hxL
    rw [IntermediateField.mem_fixedField_iff] at hx
    push Not at hx
    obtain ⟨σ, hσ, hne⟩ := hx
    have hfixy : ∀ y ∈ L₁, σ y = y := fun y hy => by simpa using hσ ⟨y, hy⟩
    haveI := M.fin
    haveI : FiniteDimensional K (M.toIF.map σ.toAlgHom) :=
      (IntermediateField.intermediateFieldMap σ M.toIF).toLinearEquiv.finiteDimensional
    set N : FinSub K Kbar :=
      ⟨M.toIF ⊔ M.toIF.map σ.toAlgHom, IntermediateField.finiteDimensional_sup _ _⟩ with hN
    have hgmem : ∀ y : M.toIF, σ (y : Kbar) ∈ N.toIF := by
      intro y
      refine le_sup_right (α := IntermediateField K Kbar) ?_
      exact ⟨(y : Kbar), y.2, rfl⟩
    set g' : M ⟶ N :=
      show M.toIF →ₐ[K] N.toIF from
        (σ.toAlgHom.comp M.toIF.val).codRestrict N.toIF.toSubalgebra (fun y => hgmem y) with hg'
    set h' : M ⟶ N :=
      show M.toIF →ₐ[K] N.toIF from IntermediateField.inclusion le_sup_left with hh'
    have hcomp : f ≫ g' = f ≫ h' := by
      refine FinSub.hom_ext (AlgHom.ext fun y => ?_)
      refine Subtype.ext ?_
      show σ (((FinSub.hom f) y : M.toIF) : Kbar) = (((FinSub.hom f) y : M.toIF) : Kbar)
      exact hfixy _ ⟨y, rfl⟩
    have hgh : g' = h' := (cancel_epi f).mp hcomp
    have hfixx : σ x = x := by
      have h1 := congrArg FinSub.hom hgh
      exact congrArg (fun t : M.toIF →ₐ[K] N.toIF => ((t ⟨x, hxM⟩ : N.toIF) : Kbar)) h1
    exact hne hfixx
  have hsurj : Function.Surjective (FinSub.hom f) := by
    intro z
    obtain ⟨y, hy⟩ := hge z.2
    exact ⟨y, Subtype.ext hy⟩
  let e : L.toIF ≃ₐ[K] M.toIF :=
    AlgEquiv.ofBijective (FinSub.hom f) ⟨FinSub.hom_injective f, hsurj⟩
  refine ⟨show M ⟶ L from (e.symm : M.toIF →ₐ[K] L.toIF), ?_, ?_⟩
  · exact FinSub.hom_ext (AlgHom.ext fun y => e.symm_apply_apply y)
  · exact FinSub.hom_ext (AlgHom.ext fun y => e.apply_symm_apply y)

/-- ★★★★★**`𝒟 = B(G)⁰` は of FSM-type**。

★FSM = fiberwise-surjective ∧ mono だが、**mono だけで足りる**。 -/
theorem finSubOp_isOfFSMType : IsOfFSMType (FinSub K Kbar)ᵒᵖ := by
  intro B A β hβ
  have hmono : Mono β := hβ.2
  have hepi : Epi β.unop := by
    refine ⟨fun {Z} g h e => ?_⟩
    have h1 : g.op ≫ β = h.op ≫ β := Quiver.Hom.unop_inj e
    exact Quiver.Hom.op_inj ((cancel_mono β).mp h1)
  haveI := finSub_isIso_of_epi β.unop hepi
  exact ⟨(inv β.unop).op,
    Quiver.Hom.unop_inj (IsIso.inv_hom_id β.unop),
    Quiver.Hom.unop_inj (IsIso.hom_inv_id β.unop)⟩

end FSM

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の底の圏 `𝒟 = B(G)⁰`(連結・totally epimorphic)。 -/
def finSubOp_totallyEpimorphic.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 底の圏 𝒟 = B(G)⁰ は連結かつ totally epimorphic",
    sectionId := "frdi-example-6-1" }

/-- ★locator —— `Example 6.1` の底の圧が of FSM-type であること。 -/
def finSubOp_isOfFSMType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 底の圧 B(G)⁰ は of FSM-type",
    sectionId := "frdi-example-6-1" }

/-- ★★**`FinSub` の自己射は同型**(有限次元だから全単射)。

★これが `Theorem 6.2, (iii)` / `Theorem 6.4, (i)` の non-dilating に効く。 -/
theorem finSub_isIso_of_endo {L : FinSub K Kbar} (σ : L ⟶ L) : IsIso σ := by
  haveI := L.fin
  have hb : Function.Bijective (FinSub.hom σ) := AlgHom.bijective (FinSub.hom σ)
  refine ⟨(AlgEquiv.ofBijective (FinSub.hom σ) hb).symm.toAlgHom, ?_, ?_⟩
  · refine FinSub.hom_ext ?_
    refine AlgHom.ext fun x => ?_
    exact (AlgEquiv.ofBijective (FinSub.hom σ) hb).symm_apply_apply x
  · refine FinSub.hom_ext ?_
    refine AlgHom.ext fun x => ?_
    exact (AlgEquiv.ofBijective (FinSub.hom σ) hb).apply_symm_apply x

/-- ★反対圏でも同じ。 -/
theorem finSubOp_isIso_of_endo {A : (FinSub K Kbar)ᵒᵖ} (e : A ⟶ A) : IsIso e := by
  haveI := finSub_isIso_of_endo e.unop
  refine ⟨(inv e.unop).op, ?_, ?_⟩
  · exact Quiver.Hom.unop_inj (IsIso.inv_hom_id e.unop)
  · exact Quiver.Hom.unop_inj (IsIso.hom_inv_id e.unop)

end ABC3.Found.FrdI
