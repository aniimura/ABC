/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LocalPointChart
import ABC3.Found.GenEll.ImmersionGlobalToProj
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★紐の付け替え —— `ℙᴺ` の側から `X` の側へ（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★これは何か —— チャートの添字が指定できる

`§9-939` は「局所化した点は**どれかの** `X_{s_i}` を通る」を与えた。
★しかし `§9-940` の整合（`𝔞_Q = (x_i)`）は**特定の添字** `i = j`
（`§9-943` の最小割り切り成分）で必要である。

★★★★本ファイルはその**添字を指定した形**を与える:

    `Spec (𝓞_F)_Q ⟶ ℙᴺ` が `D₊(x_j)` を通る  ⟹  `Spec (𝓞_F)_Q ⟶ X` が `X_{s_j}` を通る

★機構は `§9-913`（`ψ⁻¹(D₊(x_j)) = X_{s_j}`）だけである。

## ★★★これで鎖が繋がった

| 段 | 内容 |
|---|---|
| `§9-943` | 有限素点では座標の 1 つ `x_j` が他を全部割る（`𝔞_Q = (x_j)`） |
| `§9-944` | 比の組 `x_k/x_j` から `Spec (𝓞_F)_Q ⟶ D₊(x_j)` を構成 |
| `§9-945` | 分離的なら生成点で決まる（付値判定法） |
| `§9-946` | 比の組が点を決める（体の上で） |
| `§9-947` | ★局所化した点は `D₊(x_j)` を通る |
| **本ファイル** | ★★**したがって `X_{s_j}` を通る** |

## ★残っている段（明示）

★★残るのは `§9-940` の `hw`——`(s_0/s_j)(y_Q)·x_j = x_0`——である。
`§9-943` は `g·x_j = x_0` なる `g` を与えているので、
★**その `g` がチャート射の値 `(s_0/s_j)(y_Q)` と一致する**ことを言えばよい。
機構は `§9-886`（チャート射の値は比の切断の値）である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial NumberField ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★★★★添字を指定した局所チャート -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★**`ℙᴺ` の側で `D₊(x_j)` を通るなら、
`X` の側で `X_{s_j}` を通る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は `§9-913`（`ψ⁻¹(D₊(x_j)) = X_{s_j}`）と `IsOpenImmersion.lift` だけである。
★★`§9-939` は「どれかのチャート」だったが、本定理は**添字を指定できる**
——`§9-940` の整合は特定の添字（`§9-943` の最小割り切り成分）で必要だからである。 -/
theorem exists_localChart_at_index (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (F : Type) [Field F] [NumberField F]
    (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal)
    (xF : specRingOfIntegers F ⟶ X) (j : Fin (N + 1))
    (hrange : haveI := hQ.isPrime; Set.range
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF
        ≫ globalToProj M hM φ s hcov).base
      ⊆ Set.range (chartA N ℤ j).base) :
    haveI := hQ.isPrime
    ∃ y : Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶ (nonVanishing M (s j)).toScheme,
      Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF
        = y ≫ (nonVanishing M (s j)).ι := by
  haveI := hQ.isPrime
  have hsub : Set.range
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF).base
      ⊆ Set.range ((nonVanishing M (s j)).ι).base := by
    rw [Scheme.Opens.range_ι]
    rintro _ ⟨z, rfl⟩
    have hz2 := hrange ⟨z, rfl⟩
    have hz3 : (Spec.map (CommRingCat.ofHom
          (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF ≫ globalToProj M hM φ s hcov).base z
        ∈ Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)
            (MvPolynomial.X j) := by
      rw [← Proj.opensRange_awayι (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)
        (MvPolynomial.X j) (MvPolynomial.isHomogeneous_X ℤ j) one_pos]
      exact hz2
    have hmem : (Spec.map (CommRingCat.ofHom
          (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF).base z
        ∈ globalToProj M hM φ s hcov ⁻¹ᵁ
          (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)
            (MvPolynomial.X j)) := hz3
    rwa [globalToProj_preimage_basicOpen M hM φ s hcov j] at hmem
  exact ⟨IsOpenImmersion.lift (nonVanishing M (s j)).ι _ hsub,
    (IsOpenImmersion.lift_fac _ _ _).symm⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_localChart_at_index.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ℙᴺ の側で D₊(x_j) を通るなら X の側で X_{s_j} を通る)",
    sectionId := "genell-prop-1-4" }

def exists_localChart_at_index.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "globalToProj_preimage_basicOpen(ψ⁻¹(D₊(x_j)) = X_{s_j}、§9-913)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalToProj_preimage_basicOpen") 3,
    .citation "[ABC3]" "range_localized_subset_chart(局所化した点は D₊(x_j) を通る、§9-947)"
      (.inProject "ABC3" "ABC3.Found.GenEll.range_localized_subset_chart") 3,
    .citation "[mathlib]" "IsOpenImmersion.lift(像が入っていれば射は持ち上がる)"
      (.inMathlib "AlgebraicGeometry.IsOpenImmersion.lift") 1,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-939 は「どれかのチャート」だったが、" ++
       "§9-940 の整合(𝔞_Q = (x_i))は**特定の添字**" ++
       "(§9-943 の最小割り切り成分 j)で必要である。" ++
       "★本定理は添字を指定した形を与え、§9-947 と繋がる") 5,
    .implicitStep
      ("★残るのは §9-940 の hw——(s_0/s_j)(y_Q)·x_j = x_0——である。" ++
       "§9-943 は g·x_j = x_0 なる g を与えているので、" ++
       "★その g がチャート射の値 (s_0/s_j)(y_Q) と一致することを言えばよい" ++
       "(機構は §9-886——チャート射の値は比の切断の値)") 4 ]

end ABC3.Found.GenEll
