import ABC3.Found.Arakelov.UltraValuation
import Mathlib.AlgebraicGeometry.ValuativeCriterion

/-!
# Arakelov (C2) 第 81–82 ブロック —— **★★★★超フィルターを付値判定法に載せる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★超フィルターから超積体の点を作る

`Arc X` 上の超フィルター `𝒰` を取る。★`X` は擬コンパクト(固有だから)なので
**有限個のアフィン開で覆える**——超フィルターの性質から、
**あるアフィン chart `U` に 𝒰-ほとんどすべての点が入る**。

★★その族は環準同型の族 `Γ(X,U) → ℂ` を与え、超積を取ると

    Γ(X,U) → *ℂ = ℂ^{Arc X}/𝒰

★★★これが `Spec *ℂ ⟶ X` を与える。

## ★★★★★付値判定法が持ち上げる

    Spec *ℂ ⟶ X
        ↓         ↘
    Spec O   ⟶   Spec ℤ

★`O` は有限元の付値環(第 80)、`*ℂ` はその分数体。
★★mathlib の `ValuativeCriterion`(固有射の判定法)がこの四角形を**持ち上げる**。
★★★持ち上げ `Spec O ⟶ X` の**閉点の像**が、求める極限点である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `arcPt` | ★ℂ-点の像 |
| `exists_lift_of_mem` | ★アフィン開に入る点は chart から来る |
| `exists_affineOpen_mem_ultrafilter` | ★★**超フィルターは chart に集中する** |
| `starHom` / `starPoint` | ★★★**超積体の点** |
| `exists_ultraLift` | ★★★★**付値判定法による持ち上げ** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Filter

variable {X : Scheme.{0}}

/-! ## ★chart への集中 -/

/-- ★ℂ-点の像。 -/
noncomputable def arcPt (p : Spec (CommRingCat.of ℂ) ⟶ X) : X := p.base default

/-- ★アフィン開の中に像がある点は、その chart から来る。 -/
theorem exists_lift_of_mem (U : X.Opens)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : arcPt p ∈ U) :
    ∃ q : Spec (CommRingCat.of ℂ) ⟶ U.toScheme, q ≫ U.ι = p := by
  refine ⟨IsOpenImmersion.lift U.ι p ?_, IsOpenImmersion.lift_fac _ _ _⟩
  rintro _ ⟨z, rfl⟩
  have hz : z = default := Subsingleton.elim _ _
  subst hz
  rw [Scheme.Opens.range_ι]
  exact hp

/-- ★★**超フィルターはあるアフィン chart に集中する**。

★`X` が擬コンパクトなので有限被覆が取れ、超フィルターは有限和の一つを含む。 -/
theorem exists_affineOpen_mem_ultrafilter [CompactSpace X]
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X)) :
    ∃ U : X.affineOpens, {p : Spec (CommRingCat.of ℂ) ⟶ X | arcPt p ∈ U.1} ∈ 𝒰 := by
  have hcov : (Set.univ : Set X) ⊆ ⋃ U : X.affineOpens, (U.1 : Set X) := by
    intro x _
    obtain ⟨_, ⟨V, hV, rfl⟩, hxV, -⟩ :=
      (X.isBasis_affineOpens).exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
    exact Set.mem_iUnion.2 ⟨⟨V, hV⟩, hxV⟩
  obtain ⟨t, ht⟩ := isCompact_univ.elim_finite_subcover (fun U : X.affineOpens => (U.1 : Set X))
    (fun U => U.1.2) hcov
  have hunion : (Set.univ : Set (Spec (CommRingCat.of ℂ) ⟶ X))
      = ⋃ U ∈ (t : Set X.affineOpens), {p : Spec (CommRingCat.of ℂ) ⟶ X | arcPt p ∈ U.1} := by
    ext p
    simp only [Set.mem_univ, true_iff, Set.mem_iUnion, Set.mem_setOf_eq]
    have hp := ht (Set.mem_univ (arcPt p))
    simp only [Set.mem_iUnion] at hp
    obtain ⟨U, hU, hmem⟩ := hp
    exact ⟨U, hU, hmem⟩
  have hmem : (⋃ U ∈ (t : Set X.affineOpens),
      {p : Spec (CommRingCat.of ℂ) ⟶ X | arcPt p ∈ U.1}) ∈ 𝒰 := by
    rw [← hunion]
    exact Filter.univ_mem
  obtain ⟨U, -, hU⟩ := (Ultrafilter.finite_biUnion_mem_iff t.finite_toSet).1 hmem
  exact ⟨U, hU⟩

/-! ## ★★★超積体の点 -/

open scoped Classical in
/-- ★chart 内の点から得られる環準同型(chart 外では既定値 `d` を使う)。 -/
noncomputable def arcHomOf (U : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : Γ(X, U.1) →+* ℂ :=
  (Spec.preimage
    ((if h : arcPt p ∈ U.1 then Classical.choose (exists_lift_of_mem U.1 p h) else d)
      ≫ (U.2.isoSpec).hom)).hom

open scoped Classical in
theorem arcHomOf_of_mem (U : X.affineOpens) (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : arcPt p ∈ U.1) :
    ∃ q : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme, q ≫ U.1.ι = p ∧
      arcHomOf U d p = (Spec.preimage (q ≫ (U.2.isoSpec).hom)).hom := by
  refine ⟨Classical.choose (exists_lift_of_mem U.1 p hp),
    Classical.choose_spec (exists_lift_of_mem U.1 p hp), ?_⟩
  rw [arcHomOf, dif_pos hp]

/-- ★★**超積の点**——`Γ(X,U) → *ℂ`。 -/
noncomputable def starHom (U : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X)) :
    Γ(X, U.1) →+* Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ :=
  (Germ.coeRingHom _).comp (RingHom.pi (fun p => arcHomOf U d p))

theorem starHom_apply (U : X.affineOpens) (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X)) (x : Γ(X, U.1)) :
    starHom U d 𝒰 x
      = (((fun p => arcHomOf U d p x) : (Spec (CommRingCat.of ℂ) ⟶ X) → ℂ) :
          Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ) := rfl

/-- ★★★**超積体の点** `Spec *ℂ ⟶ X`。 -/
noncomputable def starPoint (U : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X)) :
    Spec (CommRingCat.of (Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ)) ⟶ X :=
  Spec.map (CommRingCat.ofHom (starHom U d 𝒰)) ≫ (U.2.isoSpec).inv ≫ U.1.ι

/-! ## ★★★★付値判定法 -/

/-- ★★★**付値可換四角形**——超積体の点と有限元の環。 -/
noncomputable def ultraSq (U : X.affineOpens)
    (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X)) :
    ValuativeCommSq (specZIsTerminal.from X) where
  R := ↥(finGermSub 𝒰)
  K := Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ
  i₁ := starPoint U d 𝒰
  i₂ := specZIsTerminal.from _
  commSq := ⟨specZIsTerminal.hom_ext _ _⟩

/-- ★★★★**付値判定法で持ち上がる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★これが Chow の補題の代わりに固有性を使う場所である。 -/
theorem exists_ultraLift (hval : ValuativeCriterion (specZIsTerminal.from X))
    (U : X.affineOpens) (d : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme)
    (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X)) :
    ∃ l : Spec (CommRingCat.of ↥(finGermSub 𝒰)) ⟶ X,
      Spec.map (CommRingCat.ofHom
          (algebraMap ↥(finGermSub 𝒰)
            (Germ (𝒰 : Filter (Spec (CommRingCat.of ℂ) ⟶ X)) ℂ))) ≫ l
        = starPoint U d 𝒰 := by
  obtain ⟨inst⟩ := hval (ultraSq U d 𝒰)
  exact ⟨inst.default.l, inst.default.fac_left⟩

/-- ★★**極限点の候補**——持ち上げを閉点で見る。 -/
noncomputable def ultraLimit (𝒰 : Ultrafilter (Spec (CommRingCat.of ℂ) ⟶ X))
    (l : Spec (CommRingCat.of ↥(finGermSub 𝒰)) ⟶ X) : Spec (CommRingCat.of ℂ) ⟶ X :=
  Spec.map (CommRingCat.ofHom (stHom 𝒰)) ≫ l

/-! ## ★出典の紐付け(`.src`) -/

def exists_ultraLift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——固有性から超フィルターの極限を作る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
