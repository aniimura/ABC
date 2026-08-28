/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.EnlargedSurjective
import ABC3.Found.GenEll.ImmersionSubfamily
import ABC3.Found.GenEll.FinCover
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★ample なら射影空間へ埋め込める（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★★★★★★★これは何か —— 原文が 1 行で引く定理

原文は `Proposition 1.4, (iv)`（Northcott）の証明で

> that [some positive tensor power of] the ample line bundle `L_ℚ` **yields an embedding**

と 1 行で書く。★これは Serre の定理であり、本ファイルがその**射の側**を型にする:

    `IsAmple M`  ⟹  ある `L`・`n`・`N'` と切断族 `s'` があって
                    `ψ ≔ globalToProj … : X ⟶ ℙᴺ'` は **`IsImmersion`**

## ★★★★★★これまでの鎖（12 ブロック）

| 段 | 内容 | ブロック |
|---|---|---|
| E1d 前段 | 重なりの上で比は単元 | `§9-907` |
| E1d 前段 | 比のコサイクル則 | `§9-908` |
| E1d 前段 | 重なりの環準同型 | `§9-909` |
| E1d 前段 | それが `i`・`j` について対称 | `§9-910` |
| **E1d** | **`ψ : X ⟶ ℙᴺ_R` が貼れた** | `§9-911` |
| — | 比が消えない所は `X_s ∩ X_t` | `§9-912` |
| — | **`ψ⁻¹(D₊(x_i)) = X_{s_i}`・`ψ` は埋め込み** | `§9-913` |
| E2 | 被覆を `Fin (N+1)` に並べ直す | `§9-914` |
| — | アフィン性を `IsAmple` から運ぶ | `§9-915` |
| — | 埋め込み性は**部分族**で足りる | `§9-916` |
| E3 | 分母が複数でも単一の指数 | `§9-917` |
| E3 | 族の拡大 | `§9-918` |
| **E3** | **チャート写像は全射** | `§9-919` |
| **到達** | **本ファイル** | `§9-920` |

## ★★残っている仮定（明示）

★★★**幾何の側の仮定はすべて「有限アフィン自明化被覆」に集約された**:

* `U : ι → X.Opens`（`Fintype ι`）、`⨆ U i = ⊤`
* `hU`・`hUij`——`U i` と `U i ⊓ U j` がアフィン（＝分離的＋準コンパクト）
* `e`——`M` が `U i` の上で自明（＝直線束であること）

★これらは「`X` が `Spec A` 上分離的・準コンパクトで `M` が可逆」から出るはずのもので、
**Arakelov 理論ではなくスキーム論の一般論**である。
★★`f : X ⟶ Spec A` の `LocallyOfFiniteType` と `hφ`（底環の像が `φ` の像に入る）は
`§9-833` の全射判定がそのまま要求するものである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}} {A : CommRingCat.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★チャートのデータから -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★**チャートのデータから埋め込みが出る**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

`X` を覆うアフィンな非消失軌跡 `{X_{s_i}}` があれば、
★族を拡大して（`§9-918`・`§9-919`）**埋め込み** `X ⟶ ℙᴺ'` が作れる。 -/
theorem exists_immersion_of_chart_data (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcovs : (⨆ i, nonVanishing M (s i)) = ⊤)
    (haff : ∀ i, IsAffineOpen (nonVanishing M (s i)))
    {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ ⊤ (by simp) r ∈ Set.range φ) :
    ∃ (n N' : ℕ)
      (s' : Fin (N' + 1) → (((sheafifyFunctor X).obj
        (presheafTensorPow M (n + 1))).val.obj (op (⊤ : X.Opens)) : Type))
      (hcov' : (⨆ k, nonVanishing ((sheafifyFunctor X).obj
        (presheafTensorPow M (n + 1))).val (s' k)) = ⊤),
      IsImmersion (globalToProj ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
        (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
        φ s' hcov') := by
  obtain ⟨n, N', s', ρ, hVeq, hsurj⟩ :=
    exists_enlarged_family_surjective f M hM U hcov hU hUij e s haff φ hφ
  obtain ⟨hcov', hcov₀⟩ :=
    iSup_nonVanishing_of_subfamily s' ρ (fun i => nonVanishing M (s i)) hVeq hcovs
  refine ⟨n, N', s', hcov', ?_⟩
  refine isImmersion_globalToProj_of_subfamily _ _ φ s' hcov'
    (fun j => ∃ i, ρ i = j) hcov₀ ?_ ?_
  · rintro ⟨j, i, rfl⟩
    rw [hVeq i]
    exact haff i
  · rintro ⟨j, i, rfl⟩
    exact hsurj i

/-! ## ★★★★★★★★★★★★★★★★★★★★`IsAmple` から -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**ample なら射影空間へ埋め込める**
—— 原文の "yields an embedding"。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `IsAmple M` ⟹ ある `L`・`n`・`N'` と切断族 `s'` があって
                  `ψ ≔ globalToProj … : X ⟶ ℙᴺ'` は **`IsImmersion`**

★★`L` は「共通次数」（`§9-CommonDegree`）、`n` は「分母を払う指数」（`§9-917`）であり、
どちらも原文の「some positive tensor power」の中身である。 -/
theorem exists_immersion_of_isAmple (f : X ⟶ Spec A) [LocallyOfFiniteType f]
    {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (hample : IsAmple M) (x₀ : (X : Type))
    {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (hφ : ∀ r : (Γ(Spec A, (⊤ : (Spec A).Opens)) : Type),
      f.appLE ⊤ ⊤ (by simp) r ∈ Set.range φ) :
    ∃ (L n N' : ℕ)
      (s' : Fin (N' + 1) → (((sheafifyFunctor X).obj
        (presheafTensorPow (presheafTensorPow M L) (n + 1))).val.obj
          (op (⊤ : X.Opens)) : Type))
      (hcov' : (⨆ k, nonVanishing ((sheafifyFunctor X).obj
        (presheafTensorPow (presheafTensorPow M L) (n + 1))).val (s' k)) = ⊤),
      IsImmersion (globalToProj
        ((sheafifyFunctor X).obj (presheafTensorPow (presheafTensorPow M L) (n + 1))).val
        (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow
          (isLocallyTrivial_presheafTensorPow hM L) (n + 1)))
        φ s' hcov') := by
  obtain ⟨L, hL, N, s, hcovs, haffs⟩ := exists_fin_cover_of_isAmple M hM hample x₀
  obtain ⟨n, N', s', hcov', himm⟩ := exists_immersion_of_chart_data f
    (presheafTensorPow M L) (isLocallyTrivial_presheafTensorPow hM L)
    U hcov hU hUij (fun i => tensorPowTriv (e i) L) s hcovs haffs φ hφ
  exact ⟨L, n, N', s', hcov', himm⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_immersion_of_chart_data.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートのデータから埋め込みが出る)",
    sectionId := "genell-prop-1-4" }

def exists_immersion_of_isAmple.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(ample なら射影空間へ埋め込める——原文の yields an embedding)",
    sectionId := "genell-prop-1-4" }

def exists_immersion_of_isAmple.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_fin_cover_of_isAmple(IsAmple から被覆を作る、§9-914・915)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_fin_cover_of_isAmple") 2,
    .citation "[ABC3]" "exists_enlarged_family_surjective(段 E3、§9-919)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_enlarged_family_surjective") 4,
    .citation "[ABC3]" "isImmersion_globalToProj_of_subfamily(埋め込み性は部分族で足りる、§9-916)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isImmersion_globalToProj_of_subfamily") 3,
    .citation "[ABC3]" "globalToProj(貼った射、§9-911)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalToProj") 3,
    .implicitStep
      ("★★★幾何の側の仮定はすべて「有限アフィン自明化被覆」に集約された: " ++
       "U : ι → X.Opens(Fintype ι)・⨆ U i = ⊤・U i と U i ⊓ U j がアフィン・" ++
       "M が U i の上で自明。これらは「X が Spec A 上分離的・準コンパクトで M が可逆」" ++
       "から出るはずのもので、**Arakelov 理論ではなくスキーム論の一般論**である") 4,
    .implicitStep
      ("★★L は「共通次数」(§9-CommonDegree)、n は「分母を払う指数」(§9-917)であり、" ++
       "どちらも原文の「some positive tensor power」の中身である") 2,
    .implicitStep
      ("★残るのは**高さの側**——D^{⊗n} = ψ^*(超平面)を言えば " ++
       "northcott_of_veryAmple(§9-882)に繋がり、Proposition 1.4, (iv) が閉じる") 6 ]

end ABC3.Found.GenEll
