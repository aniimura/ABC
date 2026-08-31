/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.EnlargedFamily
import ABC3.Found.GenEll.GlobalChartSurjective
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★段 E3 が閉じた —— チャート写像は全射である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★★これは何か —— 段 E3 の到達点

`§9-918` の族の拡大と `§9-847` の全射判定を繋ぐ:

    分母の族 `{s_i}`（アフィンなチャート）⟹ 拡大した族 `s'` について
    **どの分母 `ρ i` でもチャート写像 `A⁰_{x_{ρ i}} → Γ(X, X_{s'_{ρ i}})` は全射**

★★これが `§9-913`／`§9-916` の `hsurj` そのものである。

## ★★★機構 —— `subst` で開集合を潰す

★全射判定（`exists_finset_surjective_criterion`）は**開集合 `V` だけに依る**ので、
拡大の**前**に `V_i ≔ X_{s_i}` について取っておける。

★★拡大後のチャート環は `Γ(X, X_{s'_{ρ i}})` だが、
`nonVanishing_unit_secPow` により `X_{s'_{ρ i}} = X_{s_i}` なので**同じ環**である。
★★★型の上でそれを合わせるのが `surjective_globalAwayHom_of_ratios` で、
`V` を変数にしておいて **`subst`** すれば制限写像が恒等になる（`res_self`）。

## ★配管の記録

★★★★**証明無関係（proof irrelevance）が効く**——
`homOfLE p` と `homOfLE q`（`p q : U ≤ V` の 2 つの証明）は**定義上等しい**ので、
`≤` の証明の食い違いは気にしなくてよい。
★食い違うのは**開集合そのもの**（`X_{s'_{ρ i}}` と `X_{(s_i)^{⊗(n+1)}}`）の方であり、
そちらは `∀ u, s' (ρ i) = u → …` と汎化してから `rintro rfl` で潰す。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}} {A : CommRingCat.{0}}

/-! ## ★★★★★★★★★比が揃えば全射である -/

/-- ★★★★★★★★★**試験元が全部比なら、チャート写像は全射である**。

★全射判定 `hT` は**開集合 `V` だけに依る**ので、拡大の前に取っておける。
★★`V` を変数にしておいて `subst` すれば制限写像が恒等になり、型が合う。 -/
theorem surjective_globalAwayHom_of_ratios (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    (M' : X.PresheafOfModules) (hM' : IsLocallyTrivial X M')
    {N' : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s' : Fin (N' + 1) → (M'.obj (op ⊤) : Type)) (i₀ : Fin (N' + 1))
    {V : X.Opens} (hVeq : nonVanishing M' (s' i₀) = V)
    (T : Finset ((Γ(X, V) : Type)))
    (hT : ∀ {R' : Type} [CommRing R'] (ψ : R' →+* (Γ(X, V) : Type)),
      Set.range (f.appLE ⊤ V (by rw [← hVeq]; simp)).hom ⊆ Set.range ψ →
      (T : Set ((Γ(X, V) : Type))) ⊆ Set.range ψ → Function.Surjective ψ)
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ ⊤ (by simp) r ∈ Set.range φ)
    (hratio : ∀ g ∈ T, ∃ j : Fin (N' + 1),
      X.presheaf.map (homOfLE (le_of_eq hVeq)).op g = globalRatio M' hM' (s' j) (s' i₀)) :
    Function.Surjective (globalAwayHom M' hM' φ s' i₀) := by
  subst hVeq
  refine hT (globalAwayHom M' hM' φ s' i₀)
    (range_appLE_subset_global f M' hM' φ s' i₀ hφ) ?_
  intro g hg
  obtain ⟨j, hj⟩ := hratio g hg
  rw [res_self] at hj
  exact mem_range_globalAwayHom M' hM' φ s' i₀ j g hj

/-! ## ★★★★★★★★★★★★★★★段 E3 の到達点 -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★**拡大した族では分母のチャート写像が全射である** —— 段 E3。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    アフィンなチャート `{X_{s_i}}` ⟹ ある `n` と拡大した族 `s'` があって
    * `X_{s'_{ρ i}} = X_{s_i}`（開集合は変わらない）
    * ★**チャート写像 `A⁰_{x_{ρ i}} → Γ(X, X_{s'_{ρ i}})` は全射**

★★これが `§9-916` の `hsurj`（部分族についてのもの）そのものである。 -/
theorem exists_enlarged_family_surjective (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ ⊤ (by simp) r ∈ Set.range φ) :
    ∃ (n N' : ℕ)
      (s' : Fin (N' + 1) → (((sheafifyFunctor X).obj
        (presheafTensorPow M (n + 1))).val.obj (op (⊤ : X.Opens)) : Type))
      (ρ : Fin (N + 1) → Fin (N' + 1)),
      (∀ i, nonVanishing ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
          (s' (ρ i)) = nonVanishing M (s i)) ∧
      (∀ i, Function.Surjective (globalAwayHom
        ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
        (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
        φ s' (ρ i))) := by
  choose Tfam hTfam using fun i : Fin (N + 1) => exists_finset_surjective_criterion f (haff i)
  obtain ⟨n, N', s', ρ, hden, hratio⟩ :=
    exists_enlarged_family M hM U hcov hU hUij e s Tfam
  have hVeq : ∀ i, nonVanishing ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
      (s' (ρ i)) = nonVanishing M (s i) := by
    intro i
    rw [hden i]
    exact nonVanishing_unit_secPow M hM (s i) n
  refine ⟨n, N', s', ρ, hVeq, ?_⟩
  intro i
  refine surjective_globalAwayHom_of_ratios f _ _ φ s' (ρ i) (hVeq i)
    (Tfam i) (hTfam i) hφ ?_
  intro g hg
  obtain ⟨j, hj⟩ := hratio i g hg
  refine ⟨j, ?_⟩
  have key : ∀ (u : ((((sheafifyFunctor X).obj
        (presheafTensorPow M (n + 1))).val.obj (op (⊤ : X.Opens))) : Type)),
      s' (ρ i) = u →
      (∀ (hle : nonVanishing ((sheafifyFunctor X).obj
          (presheafTensorPow M (n + 1))).val u ≤ nonVanishing M (s i)),
        X.presheaf.map (homOfLE hle).op g
          = globalRatio ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
              (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
              (s' j) u) →
      X.presheaf.map (homOfLE (le_of_eq (hVeq i))).op g
        = globalRatio ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
            (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
            (s' j) (s' (ρ i)) := by
    rintro u rfl h
    exact h _
  exact key _ (hden i) (fun _ => hj)

/-! ## ★被覆の移送 -/

/-- ★**部分族が覆えば全体も覆う**（拡大した族の被覆条件）。 -/
theorem iSup_nonVanishing_of_subfamily {M' : X.PresheafOfModules} {N' N : ℕ}
    (s' : Fin (N' + 1) → (M'.obj (op ⊤) : Type)) (ρ : Fin (N + 1) → Fin (N' + 1))
    (W : Fin (N + 1) → X.Opens)
    (hVeq : ∀ i, nonVanishing M' (s' (ρ i)) = W i)
    (hcovs : (⨆ i, W i) = ⊤) :
    (⨆ k, nonVanishing M' (s' k)) = ⊤ ∧
    (⨆ j : {j // ∃ i, ρ i = j}, nonVanishing M' (s' j.1)) = ⊤ := by
  have hsub : (⨆ j : {j // ∃ i, ρ i = j}, nonVanishing M' (s' j.1)) = ⊤ := by
    refine eq_top_iff.2 ?_
    rw [← hcovs]
    refine iSup_le (fun i => ?_)
    refine le_trans (le_of_eq (hVeq i).symm) ?_
    exact le_iSup (fun j : {j // ∃ i, ρ i = j} => nonVanishing M' (s' j.1)) ⟨ρ i, ⟨i, rfl⟩⟩
  refine ⟨?_, hsub⟩
  refine eq_top_iff.2 ?_
  rw [← hsub]
  exact iSup_le (fun j => le_iSup (fun k => nonVanishing M' (s' k)) j.1)

/-! ## ★出典の紐付け(`.src`) -/

def surjective_globalAwayHom_of_ratios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(試験元が全部比ならチャート写像は全射である)",
    sectionId := "genell-prop-1-4" }

def exists_enlarged_family_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3——拡大した族では分母のチャート写像が全射である)",
    sectionId := "genell-prop-1-4" }

def iSup_nonVanishing_of_subfamily.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(部分族が覆えば全体も覆う)",
    sectionId := "genell-prop-1-4" }

def exists_enlarged_family_surjective.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_enlarged_family(族の拡大、§9-918)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_enlarged_family") 3,
    .citation "[ABC3]" "exists_finset_surjective_criterion(有限個の試験元で全射判定、§9-833)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finset_surjective_criterion") 2,
    .citation "[ABC3]" "nonVanishing_unit_secPow(冪を上げても非消失軌跡は変わらない)"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_unit_secPow") 1,
    .implicitStep
      ("★★全射判定は**開集合 V だけに依る**ので、拡大の前に V_i ≔ X_{s_i} について" ++
       "取っておける。拡大後のチャート環は Γ(X, X_{s'_{ρ i}}) だが " ++
       "nonVanishing_unit_secPow により X_{s'_{ρ i}} = X_{s_i} なので同じ環である") 3,
    .implicitStep
      ("★★★★配管: 証明無関係(proof irrelevance)が効くので " ++
       "homOfLE p と homOfLE q(≤ の 2 つの証明)は**定義上等しい**。" ++
       "食い違うのは**開集合そのもの**の方であり、" ++
       "∀ u, s' (ρ i) = u → … と汎化してから rintro rfl で潰す") 2,
    .implicitStep
      ("★★★これで §9-916 の hsurj が揃った。" ++
       "残るのは有限アフィン自明化被覆(U・hU・hUij・e)を ample から出すことである") 4 ]

end ABC3.Found.GenEll
