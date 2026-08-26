/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.GlobalUnits

/-!
# 正規スキームでは「整な有理函数は切断」(鎖 `normalize` の `normalization-universal-normal` の核)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism

## ★★これが「正規化の普遍性」の中身である

原文 `Theorem 6.2, (i)` が要求するのは
「**正規整スキームからの支配射は相対正規化を経由する**」——
`W` が正規で `g : W ⟶ Y` が支配的、`L ⊆ K(W)` なら `W ⟶ Y[L]` が作れる、である。

★★その**数学の中身**は次の 1 本に尽きる:

```
u ∈ K(W) が Ω 内のどのアフィン開 V でも Γ(W,V) 上整  ⟹  u ∈ Γ(W, Ω)
```

`Y[L]|_U = Spec(integralClosure Γ(Y,U) L)` なので、その切断を `Γ(W, g⁻¹U)` へ送る
環準同型がこれで作れる。★残りは**射の貼り合わせ**(`Scheme.OpenCover.glueMorphisms`)
という配管である。

## ★★★中身は 2 段

| 段 | 根拠 |
|---|---|
| アフィン開ごと: 整 ⟹ 切断 | `IsIntegrallyClosed.isIntegral_iff` ＋ `functionField_isFractionRing_of_isAffineOpen` |
| 開集合 `Ω` へ貼り合わせ | ★`exists_sectionOn`(本ファイル、`GlobalUnits.lean` の `⊤` 版の一般化) |

★★貼り合わせの重なりの一致は **`X` が既約であること**から自動である
(空でない開集合は必ず生成点を含むので必ず交わり、`Γ(X,−) ↪ K(X)` は単射)。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory TopologicalSpace ABC3.Meta

universe u

/-! ## ★1. 開集合の上での貼り合わせ -/

/-- ★**空でないアフィン開だけで開集合 `Ω` を被覆できる**。 -/
theorem iSup_affineOpens_le_eq (X : Scheme.{u}) (Ω : X.Opens) :
    ⨆ (V : {V : X.affineOpens // Nonempty ((V : X.Opens)) ∧ (V : X.Opens) ≤ Ω}),
      ((V.1 : X.Opens)) = Ω := by
  refine le_antisymm (iSup_le fun V => V.2.2) ?_
  intro x hx
  obtain ⟨_, ⟨V, hV, rfl⟩, hxV, hsub⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open hx Ω.2
  exact TopologicalSpace.Opens.mem_iSup.mpr ⟨⟨⟨V, hV⟩, ⟨⟨⟨x, hxV⟩⟩, hsub⟩⟩, hxV⟩

/-- ★★★★**開集合 `Ω` の上での貼り合わせ** ——
`Ω` を覆う空でない開集合の各成分で同じ `u ∈ K(X)` を実現する切断があるなら、
`Ω` の上に `u` を実現する切断がある。

★`GlobalUnits.lean` の `exists_globalSection`(`Ω = ⊤` の場合)の一般化。 -/
theorem exists_sectionOn {X : Scheme.{u}} [IsIntegral X]
    {ι : Type*} (U : ι → X.Opens) (Ω : X.Opens) (hcover : Ω ≤ iSup U)
    [hne : ∀ i, Nonempty (U i)] [Nonempty Ω] (hle : ∀ i, U i ≤ Ω)
    (u : X.functionField) (s : ∀ i, Γ(X, U i))
    (hs : ∀ i, (ConcreteCategory.hom (X.germToFunctionField (U i))) (s i) = u) :
    ∃ t : Γ(X, Ω), (ConcreteCategory.hom (X.germToFunctionField Ω)) t = u := by
  classical
  haveI : Nonempty ι := by
    by_contra h
    simp only [not_nonempty_iff] at h
    have h1 : Ω ≤ ⊥ := by rw [← iSup_of_empty U]; exact hcover
    have h2 : genericPoint (X : Type u) ∈ (⊥ : X.Opens) :=
      h1 (genericPoint_mem_of_nonempty_opens (X := X) (U := Ω) inferInstance)
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
  obtain ⟨t, ht, -⟩ := X.sheaf.existsUnique_gluing' U Ω (fun i => homOfLE (hle i)) hcover s hcompat
  refine ⟨t, ?_⟩
  have hi := Classical.arbitrary ι
  have ht' : (ConcreteCategory.hom
      (X.presheaf.map (homOfLE (hle hi) : U hi ⟶ Ω).op)) t = s hi := ht hi
  show (ConcreteCategory.hom (X.presheaf.germ Ω (genericPoint (X : Type u)) _)) t = u
  rw [← TopCat.Presheaf.germ_res_apply X.presheaf (homOfLE (hle hi) : U hi ⟶ Ω)
    (genericPoint (X : Type u)) (genericPoint_mem_of_nonempty_opens (hne hi)), ht']
  exact hs hi

/-- ★★★**アフィン開ごとに切断なら `Ω` の上で切断**。 -/
theorem mem_range_germ_of_forall_affine_le {X : Scheme.{u}} [IsIntegral X]
    (Ω : X.Opens) [Nonempty Ω] (u : X.functionField)
    (h : ∀ (V : X.Opens) (_ : IsAffineOpen V) (hne : Nonempty V) (_ : V ≤ Ω),
      u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField V))) :
    u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField Ω)) := by
  classical
  set ι := {V : X.affineOpens // Nonempty ((V : X.Opens)) ∧ (V : X.Opens) ≤ Ω} with hι
  haveI hne : ∀ i : ι, Nonempty (((i.1 : X.Opens))) := fun i => i.2.1
  have hsec : ∀ i : ι, ∃ t : Γ(X, (i.1 : X.Opens)),
      (ConcreteCategory.hom (X.germToFunctionField (i.1 : X.Opens))) t = u := fun i =>
    h (i.1 : X.Opens) i.1.2 i.2.1 i.2.2
  choose s hs using hsec
  exact exists_sectionOn (fun i : ι => (i.1 : X.Opens)) Ω
    (le_of_eq (iSup_affineOpens_le_eq X Ω).symm) (fun i => i.2.2) u s hs

/-! ## ★2. 整な有理函数は切断 -/

/-- ★★**アフィン開の上: 整閉なら整な元は切断**。 -/
theorem mem_range_germ_of_isIntegral_affine {X : Scheme.{u}} [IsIntegral X]
    {V : X.Opens} (hV : IsAffineOpen V) (hne : Nonempty V)
    (hic : IsIntegrallyClosed Γ(X, V)) (u : X.functionField)
    (hint : IsIntegral Γ(X, V) u) :
    u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField V)) := by
  haveI := hne
  haveI := hic
  haveI := functionField_isFractionRing_of_isAffineOpen (X := X) V hV
  exact (IsIntegrallyClosed.isIntegral_iff).mp hint

/-- ★★★★★★**正規スキームでは「整な有理函数は切断」** ——
`u ∈ K(X)` が `Ω` 内のどのアフィン開 `V` でも `Γ(X,V)` 上整なら `u ∈ Γ(X,Ω)`。

★★これが `Theorem 6.2, (i)` の「正規整スキームからの支配射は相対正規化を経由する」の
**数学の中身**である。`Y[L]|_U = Spec(integralClosure Γ(Y,U) L)` の切断を
`Γ(W, g⁻¹U)` へ送る環準同型がこれで作れる。 -/
theorem mem_range_germ_of_forall_isIntegral {X : Scheme.{u}} [IsIntegral X]
    (Ω : X.Opens) [Nonempty Ω]
    (hic : ∀ (V : X.Opens) (_ : IsAffineOpen V) (hne : Nonempty V) (_ : V ≤ Ω),
      IsIntegrallyClosed Γ(X, V))
    (u : X.functionField)
    (hint : ∀ (V : X.Opens) (_ : IsAffineOpen V) (hne : Nonempty V) (_ : V ≤ Ω),
      IsIntegral Γ(X, V) u) :
    u ∈ Set.range (ConcreteCategory.hom (X.germToFunctionField Ω)) :=
  mem_range_germ_of_forall_affine_le Ω u
    (fun V hV hne hVΩ => mem_range_germ_of_isIntegral_affine hV hne (hic V hV hne hVΩ) u
      (hint V hV hne hVΩ))

/-- ★★★**底の切断上で整なら、上の切断上でも整** —— `g` に沿った移送。

★`Y[L]` の切断は `Γ(Y,U)` 上整な `L` の元なので、これで
`Γ(W,V)` 上整に読み替えられる。 -/
theorem isIntegral_sections_of_appLE {X Y : Scheme.{u}} [IsIntegral X]
    (g : X ⟶ Y) {U : Y.Opens} {V : X.Opens} (hVU : V ≤ g ⁻¹ᵁ U) [Nonempty V]
    (u : X.functionField)
    (h : ((g.appLE U V hVU) ≫ X.germToFunctionField V).hom.IsIntegralElem u) :
    IsIntegral Γ(X, V) u := by
  letI : Algebra Γ(Y, U) Γ(X, V) := (g.appLE U V hVU).hom.toAlgebra
  letI : Algebra Γ(Y, U) X.functionField :=
    ((g.appLE U V hVU) ≫ X.germToFunctionField V).hom.toAlgebra
  haveI : IsScalarTower Γ(Y, U) Γ(X, V) X.functionField :=
    IsScalarTower.of_algebraMap_eq' rfl
  have h1 : IsIntegral Γ(Y, U) u := h
  exact h1.tower_top

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def mem_range_germ_of_forall_isIntegral.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 正規スキームでは整な有理函数は切断",
    sectionId := "frdi-thm-6-2" }

def mem_range_germ_of_forall_isIntegral.needs : List ProofObligation :=
  [ .citation "[mathlib]" "IsIntegrallyClosed.isIntegral_iff"
      (.inMathlib "IsIntegrallyClosed.isIntegral_iff") 110,
    .citation "[mathlib]" "functionField_isFractionRing_of_isAffineOpen"
      (.inMathlib "AlgebraicGeometry.functionField_isFractionRing_of_isAffineOpen") 110,
    .citation "[ABC3]" "開集合の上での貼り合わせ"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_sectionOn") 110,
    .derivation "既約なので重なりの一致は Γ(X,−) ↪ K(X) の単射性から自動" 110 ]

def exists_sectionOn.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 開集合の上での切断の貼り合わせ",
    sectionId := "frdi-thm-6-2" }

def exists_sectionOn.needs : List ProofObligation :=
  [ .citation "[mathlib]" "TopCat.Sheaf.existsUnique_gluing'"
      (.inMathlib "TopCat.Sheaf.existsUnique_gluing'") 110,
    .citation "[ABC3]" "germToFunctionField_injective 経由の重なりの一致"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_globalSection") 110 ]

def isIntegral_sections_of_appLE.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 底の切断上で整なら上の切断上でも整",
    sectionId := "frdi-thm-6-2" }

def isIntegral_sections_of_appLE.needs : List ProofObligation :=
  [ .citation "[mathlib]" "IsIntegral.tower_top"
      (.inMathlib "IsIntegral.tower_top") 110 ]

end ABC3.Found.Divisor
