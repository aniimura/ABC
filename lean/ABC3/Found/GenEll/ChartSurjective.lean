/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Globalize
import ABC3.Found.GenEll.ChartFiniteType
import ABC3.Meta.Claim

/-!
# ★★★★★★★段 E3c の還元 —— 有限個の生成元が像に入れば全射（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★これは何か —— 段 E3c を「有限個の確認」に落とす

段 E3c は「環準同型 `A⁰_{x_i} → Γ(X, X_{s_i})` が**全射**である」ことを言う。
★本ファイルはそれを

> **有限個**の生成元 `T`（段 E3b、`§9-832`）が像に入れば全射である

に落とす。★★これで「無限個の元を全部当てる」問題が「有限個の確認」になった。

## ★★★もう一方の入口 —— `X_s` の上の大域な関数

段 E3a-3（`§9-831`）はチャートごとの関数の族から出発する形だったが、
段 E3c が実際に扱うのは **`Γ(X, X_s)` の元 1 個**である。
★`exists_global_section_of_globalFun`（本ファイル）がその橋である
——`nonVanishing_inf`（段 D2）で `X_s ⊓ U_i = D(f_i)` なので、`g` を各チャートへ制限すれば
族になり、**重なりでの一致は自動**である（どれも同じ `g` の制限だから）。

## ★★残っている段（明示）

★★★**残るのは「`a/s^n` が `A⁰_{x_i}` の像に入る」ことだけ**である。
それには射影埋め込みの**座標の選び方**が要る——`M^{⊗n}` の大域切断のうち、
生成元それぞれを延ばして得た `a_g` を**座標に含める**こと。
★そこは `§9-816`（`chartToProj`）と `§9-815`（`awayToSectionHom`）の側の仕事であり、
本ファイルには無い。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

/-! ## ★★★★★純代数 —— 生成元が像に入れば全射 -/

/-- ★★★★★**生成元が像に入れば全射である**（純代数）。

`S` が `R` 上 `T` で生成され、`ψ : R' → S` の像が `R` の像と `T` を含むなら、
★`ψ` は全射である。

★★機構は「`ψ.range` は `R`-部分代数である」＋ `Algebra.adjoin_le` だけである。 -/
theorem surjective_of_adjoin_le_range {R S R' : Type} [CommRing R] [CommRing S] [CommRing R']
    (φ : R →+* S) (ψ : R' →+* S) (T : Set S)
    (hT : @Algebra.adjoin R S _ _ φ.toAlgebra T = ⊤)
    (hφ : Set.range φ ⊆ Set.range ψ) (hTψ : T ⊆ Set.range ψ) :
    Function.Surjective ψ := by
  letI := φ.toAlgebra
  let B : Subalgebra R S :=
    { ψ.range with algebraMap_mem' := fun r => hφ ⟨r, rfl⟩ }
  have hsub : Algebra.adjoin R T ≤ B := Algebra.adjoin_le hTψ
  rw [hT] at hsub
  intro y
  exact hsub (Algebra.mem_top (R := R) (x := y))

variable {X : Scheme.{0}} {A : CommRingCat.{0}}

/-! ## ★★★★★★★段 E3c の還元 -/

/-- ★★★★★★★**段 E3c の還元** —— 有限個の生成元が像に入れば全射である。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`f : X ⟶ Spec A` が有限型で `V` がアフィン開なら、★**ある有限集合 `T ⊆ Γ(X, V)`** があって、
任意の環準同型 `ψ : R' → Γ(X, V)` について

* `ψ` の像が `Γ(Spec A, ⊤) → Γ(X, V)` の像を含み、
* `ψ` の像が `T` を含む

なら `ψ` は**全射**である。

★★これで段 E3c が「無限個の元を全部当てる」問題から
「**有限個の確認**」に落ちた。 -/
theorem exists_finset_surjective_criterion (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    {V : X.Opens} (hV : IsAffineOpen V) :
    ∃ T : Finset (Γ(X, V) : Type),
      ∀ {R' : Type} [CommRing R'] (ψ : R' →+* (Γ(X, V) : Type)),
        Set.range (f.appLE ⊤ V (by simp)).hom ⊆ Set.range ψ →
        (T : Set (Γ(X, V) : Type)) ⊆ Set.range ψ →
          Function.Surjective ψ := by
  obtain ⟨T, hT⟩ := exists_finset_generating f hV
  exact ⟨T, fun ψ hφ hTψ => surjective_of_adjoin_le_range _ ψ _ hT hφ hTψ⟩

/-! ## ★★★★★★★`X_s` の上の大域な関数から出発する形 -/

/-- ★★★★★★★**`X_s` の上の大域な関数は、冪を上げれば大域切断に延びる**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

段 E3a-3（`§9-831`）はチャートごとの族から出発する形だったが、
段 E3c が扱うのは **`Γ(X, X_s)` の元 1 個**である。★本補題がその橋である。

★★**重なりでの一致は自動**である——どのチャートの関数も同じ `g` の制限だからである。
★機構は `basicOpen_trivValue_le`（段 D2——`D(f_i) ≤ X_s`）だけである。 -/
theorem exists_global_section_of_globalFun {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (g : (Γ(X, nonVanishing M s) : Type)) :
    ∃ (n : ℕ) (a : ∀ i, (Γ(X, U i) : Type))
      (t : (((sheafifyFunctor X).obj (presheafTensorPow M n)).val.obj
        (op (⊤ : X.Opens)) : Type)),
      (∀ i, X.presheaf.map
            (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op g
          * (algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type)
            (trivValue M (U i) (e i) s)) ^ n
          = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type) (a i)) ∧
      (∀ i, ((sheafifyFunctor X).obj (presheafTensorPow M n)).val.map
          (homOfLE (le_top : U i ≤ ⊤)).op t
        = ((sheafifyUnit X (presheafTensorPow M n)).app (op (U i))).hom
            (secOfFun M (U i) (e i) n (a i))) := by
  refine exists_global_section_of_localData M U hcov hU hUij e s
    (fun i => X.presheaf.map (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op g) ?_
  intro i j
  rw [res_trans, res_trans]

/-! ## ★出典の紐付け(`.src`) -/

def surjective_of_adjoin_le_range.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(生成元が像に入れば全射——純代数)",
    sectionId := "genell-prop-1-4" }

def exists_finset_surjective_criterion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(段 E3c の還元——有限個の確認に落とす)",
    sectionId := "genell-prop-1-4" }

def exists_global_section_of_globalFun.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(X_s の上の大域な関数は冪を上げれば大域切断に延びる)",
    sectionId := "genell-prop-1-4" }

def exists_finset_surjective_criterion.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_finset_generating(段 E3b、§9-832)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finset_generating") 2,
    .citation "[ABC3]" "exists_global_section_of_localData(段 E3a-3、§9-831)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_global_section_of_localData") 3,
    .citation "[mathlib]" "Algebra.adjoin_le(生成部分代数の最小性)"
      (.inMathlib "Algebra.adjoin_le") 2,
    .implicitStep
      ("★★★**残るのは「a/s^n が A⁰_{x_i} の像に入る」ことだけ**である。" ++
       "それには射影埋め込みの**座標の選び方**が要る——M^{⊗n} の大域切断のうち、" ++
       "生成元それぞれを延ばして得た a_g を**座標に含める**こと。" ++
       "★そこは §9-816(chartToProj)と §9-815(awayToSectionHom)の側の仕事である") 8,
    .implicitStep
      ("★重なりでの一致が自動になるのは「どのチャートの関数も同じ g の制限だから」である" ++
       "——nonVanishing_inf(段 D2)で X_s ⊓ U_i = D(f_i) となることが効いている") 3 ]

end ABC3.Found.GenEll
