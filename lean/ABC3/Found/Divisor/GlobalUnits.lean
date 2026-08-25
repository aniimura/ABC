/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormFinite
import ABC3.Meta.Claim

/-!
# proper normal なら大域函数は `k_L`(鎖 `normalize` の `global-units` の 2 段目)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> that [since V [L] is a proper normal variety] for A ∈Ob(CV,

原文 (FrdI p.110):
> where kL denotes the algebraic closure of k in L [so kL is a finite separable extension

## ★★★★測定の訂正(2026-08-25)—— `proper-global-sections` は mathlib に**あった**

2026-08-21 の測定は「`IsProper` と大域切断を結ぶ補題は grep 0 件」だったが、
`Mathlib/AlgebraicGeometry/Morphisms/Proper.lean` の `## Main results` に

| 宣言 | 中身 |
|---|---|
| `isField_of_universallyClosed` | `X` が整で `Spec k` 上 universally closed なら `Γ(X,⊤)` は**体** |
| `finite_appTop_of_universallyClosed` | さらに locally of finite type なら `Γ(X,⊤)` は `k` 上**有限** |

が**そのものある**。`IsProper = IsSeparated ∧ UniversallyClosed ∧ LocallyOfFiniteType`
なので proper ならそのまま当たる。★探索先が浅かった。

## ★★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `exists_globalSection` | ★★★**貼り合わせ** —— 各開集合で `u` を実現する切断があれば大域切断がある |
| `mem_range_germ_top_of_forall_affine` | ★アフィン開の族に対する系 |
| `globalSections_isField` / `globalSections_finite` | mathlib の言い換え |
| `exists_unit_globalSection` | ★★★★**`u ≠ 0` が大域切断なら `Γ(X,⊤)` の単元** |

## ★★★貼り合わせが効くのは `X` が既約だから

`X` が整なら空でない開集合は必ず生成点を含むので、**2 つの空でない開集合は必ず交わる**。
したがって「各 `U i` の切断が同じ `u ∈ K(X)` を与える」なら、
`Γ(X, U i ⊓ U j) ↪ K(X)` の単射性(`germToFunctionField_injective`)から
重なりでの一致が**自動**で出る。★これが層の貼り合わせ条件そのものになる。

## ★残るもの(`global-units` の 1 段目)

「`ord_v(u) ≥ 0` がすべての余次元 1 の点で成り立つ ⟹ 各アフィン開の切断に入る」
——**代数的 Hartogs のスキーム版**。環の側(`Found/Divisor/Hartogs.lean`)は
閉じているので、残るのは `Γ(X,U)` の高さ 1 素イデアルと `X` の余次元 1 の点の
対応を通す配線だけである。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory TopologicalSpace ABC3.Meta

universe u

/-! ## ★1. 生成点と空でない開集合 -/

/-- ★**整スキームでは空でない開集合は生成点を含む**。 -/
theorem genericPoint_mem_of_nonempty_opens {X : Scheme.{u}} [IsIntegral X] {U : X.Opens}
    (hU : Nonempty U) : genericPoint (X : Type u) ∈ U := by
  have hg : IsGenericPoint (genericPoint (X : Type u)) Set.univ := by
    simpa using genericPoint_spec (X : Type u)
  refine (hg.mem_open_set_iff U.2).mpr ?_
  have hne : ((U : Set X)).Nonempty := ⟨hU.some.1, hU.some.2⟩
  simpa using hne

instance nonemptyTopOpens {X : Scheme.{u}} [IsIntegral X] : Nonempty ((⊤ : X.Opens)) :=
  ⟨⟨(IsIntegral.nonempty (X := X)).some, trivial⟩⟩

/-! ## ★2. 貼り合わせ -/

/-- ★★★★**貼り合わせ** —— 空でない開被覆の各成分で同じ `u ∈ K(X)` を実現する
切断があるなら、`u` を実現する大域切断がある。

★重なりでの一致は `X` の既約性から**自動**である
(2 つの空でない開集合は生成点で交わり、`Γ(X,−) ↪ K(X)` は単射)。 -/
theorem exists_globalSection {X : Scheme.{u}} [IsIntegral X]
    {ι : Type*} (U : ι → X.Opens) (hcover : iSup U = ⊤) [hne : ∀ i, Nonempty (U i)]
    (u : X.functionField) (s : ∀ i, Γ(X, U i))
    (hs : ∀ i, (ConcreteCategory.hom (X.germToFunctionField (U i))) (s i) = u) :
    ∃ t : Γ(X, ⊤), (ConcreteCategory.hom (X.germToFunctionField ⊤)) t = u := by
  classical
  haveI : Nonempty ι := by
    by_contra h
    simp only [not_nonempty_iff] at h
    have h1 : (⊤ : X.Opens) = ⊥ := by rw [← hcover, iSup_of_empty]
    have h2 : genericPoint (X : Type u) ∈ (⊥ : X.Opens) := by
      rw [← h1]; exact genericPoint_mem_of_nonempty_opens (nonemptyTopOpens (X := X))
    simp at h2
  have hcompat : TopCat.Presheaf.IsCompatible X.presheaf U s := by
    intro i j
    haveI : Nonempty ((U i ⊓ U j : X.Opens)) :=
      ⟨⟨genericPoint (X : Type u),
        ⟨genericPoint_mem_of_nonempty_opens (hne i), genericPoint_mem_of_nonempty_opens (hne j)⟩⟩⟩
    refine X.germToFunctionField_injective (U i ⊓ U j) ?_
    show (ConcreteCategory.hom (X.presheaf.germ _ (genericPoint (X : Type u)) _)) _ = _
    rw [TopCat.Presheaf.germ_res_apply, TopCat.Presheaf.germ_res_apply]
    exact (hs i).trans (hs j).symm
  obtain ⟨t, ht, -⟩ := X.sheaf.existsUnique_gluing' U ⊤ (fun i => homOfLE le_top)
    (by rw [hcover]) s hcompat
  refine ⟨t, ?_⟩
  have hi := Classical.arbitrary ι
  have ht' : (ConcreteCategory.hom
      (X.presheaf.map (homOfLE (le_top) : U hi ⟶ (⊤ : X.Opens)).op)) t = s hi := ht hi
  show (ConcreteCategory.hom (X.presheaf.germ ⊤ (genericPoint (X : Type u)) _)) t = u
  rw [← TopCat.Presheaf.germ_res_apply X.presheaf (homOfLE le_top : U hi ⟶ ⊤)
    (genericPoint (X : Type u)) (genericPoint_mem_of_nonempty_opens (hne hi)), ht']
  exact hs hi

/-- ★**空でないアフィン開だけで被覆できる**。 -/
theorem iSup_nonempty_affineOpens_eq_top (X : Scheme.{u}) :
    ⨆ (U : {U : X.affineOpens // Nonempty ((U : X.Opens))}), ((U.1 : X.Opens)) = ⊤ := by
  refine le_antisymm le_top ?_
  rintro x -
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  exact TopologicalSpace.Opens.mem_iSup.mpr ⟨⟨⟨U, hU⟩, ⟨⟨x, hxU⟩⟩⟩, hxU⟩

/-- ★★★**すべての空でないアフィン開の切断に入るなら大域切断である**。 -/
theorem mem_range_germ_top_of_forall_affine {X : Scheme.{u}} [IsIntegral X]
    (u : X.functionField)
    (h : ∀ (U : X.Opens) (_ : IsAffineOpen U) (hne : Nonempty (U)),
      u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField U))) :
    u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField ⊤)) := by
  classical
  set ι := {U : X.affineOpens // Nonempty ((U : X.Opens))} with hι
  haveI hne : ∀ i : ι, Nonempty (((i.1 : X.Opens))) := fun i => i.2
  have hsec : ∀ i : ι, ∃ t : Γ(X, (i.1 : X.Opens)),
      (ConcreteCategory.hom (X.germToFunctionField (i.1 : X.Opens))) t = u := fun i =>
    h (i.1 : X.Opens) i.1.2 i.2
  choose s hs using hsec
  obtain ⟨t, ht⟩ := exists_globalSection (fun i : ι => (i.1 : X.Opens))
    (iSup_nonempty_affineOpens_eq_top X) u s hs
  exact ⟨t, ht⟩

/-! ## ★3. proper なら大域切断は `k` 上有限次の体 -/

/-- ★★**`Γ(X,⊤)` は体**(整スキームが `Spec k` 上 proper なら)。 -/
theorem globalSections_isField {X : Scheme.{u}} [IsIntegral X] {k : Type u} [Field k]
    (g : X ⟶ Spec (CommRingCat.of k)) (hg : IsProper g) : IsField Γ(X, ⊤) := by
  haveI := hg
  exact isField_of_universallyClosed k g

/-- ★★**`Γ(X,⊤)` は `k` 上有限**。★これが原文の `k_L`(「`k` の `L` における代数閉包」)。 -/
theorem globalSections_finite {X : Scheme.{u}} [IsIntegral X] {k : Type u} [Field k]
    (g : X ⟶ Spec (CommRingCat.of k)) (hg : IsProper g) : g.appTop.hom.Finite := by
  haveI := hg
  exact finite_appTop_of_universallyClosed k g

/-- ★★★★★**`u ≠ 0` が大域切断で実現されるなら、それは `Γ(X,⊤)` の単元**。

★`Γ(X,⊤)` が体だから。これが原文の `𝒪^×(A) = k_L^×` の右半分である。 -/
theorem exists_unit_globalSection {X : Scheme.{u}} [IsIntegral X] {k : Type u} [Field k]
    (g : X ⟶ Spec (CommRingCat.of k)) (hg : IsProper g)
    {u : X.functionField} (hu : u ≠ 0)
    (h : u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField ⊤))) :
    ∃ t : Γ(X, ⊤), IsUnit t ∧ (ConcreteCategory.hom (X.germToFunctionField ⊤)) t = u := by
  obtain ⟨t, ht⟩ := h
  letI := (globalSections_isField g hg).toField
  refine ⟨t, ?_, ht⟩
  rw [isUnit_iff_ne_zero]
  intro h0
  apply hu
  rw [← ht, h0, map_zero]

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_globalSection.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — 大域函数の貼り合わせ(整スキームでは重なりの一致が自動)",
    sectionId := "frdi-example-6-1" }

def exists_globalSection.needs : List ProofObligation :=
  [ .citation "[mathlib]" "TopCat.Sheaf.existsUnique_gluing'(層の貼り合わせ)"
      (.inMathlib "TopCat.Sheaf.existsUnique_gluing'") 110,
    .citation "[mathlib]" "Scheme.germToFunctionField_injective(整スキームの切断は函数体に単射)"
      (.inMathlib "AlgebraicGeometry.Scheme.germToFunctionField_injective") 110,
    .derivation "既約なので 2 つの空でない開集合は生成点で交わる。重なりでの一致は単射性から自動" 110 ]

def globalSections_isField.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — proper normal なら Γ(X,⊤) は体",
    sectionId := "frdi-example-6-1" }

def globalSections_isField.needs : List ProofObligation :=
  [ .citation "[mathlib]" "isField_of_universallyClosed"
      (.inMathlib "AlgebraicGeometry.isField_of_universallyClosed") 110,
    .implicitStep "★原文は角括弧で「[since V [L] is a proper normal variety]」と畳む" 110 ]

def globalSections_finite.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — k_L は k の有限拡大",
    sectionId := "frdi-example-6-1" }

def globalSections_finite.needs : List ProofObligation :=
  [ .citation "[mathlib]" "finite_appTop_of_universallyClosed"
      (.inMathlib "AlgebraicGeometry.finite_appTop_of_universallyClosed") 110,
    .implicitStep "★原文は「[so kL is a finite separable extension of k]」と角括弧で置く" 110 ]

def exists_unit_globalSection.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — O^×(A) = k_L^× の右半分",
    sectionId := "frdi-example-6-1" }

def exists_unit_globalSection.needs : List ProofObligation :=
  [ .citation "[ABC3]" "globalSections_isField"
      (.inProject "ABC3" "ABC3.Found.Divisor.globalSections_isField") 110,
    .derivation "体の 0 でない元は単元" 110 ]

def mem_range_germ_top_of_forall_affine.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — アフィン開ごとの切断から大域切断へ",
    sectionId := "frdi-example-6-1" }

def mem_range_germ_top_of_forall_affine.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_globalSection"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_globalSection") 110,
    .derivation "空でないアフィン開だけで被覆できる" 110 ]

end ABC3.Found.Divisor
