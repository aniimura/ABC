/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottFinal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★局所チャートは**常に取れる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★これは何か —— 局所チャートは仮定でなく定理である

`§9-935`〜`§9-938` の鎖は点ごとに

* 有限素点 `Q`——`Spec (𝓞_F)_Q ⟶ X` が `X_{s_i}` を通ること（`chartOf`・`y`・`hy`）
* 無限素点 `v`——`Spec ℂ ⟶ ℙᴺ` が `D₊(x_i)` を通ること（`archChart`・`ρ`・`hfac`）

を仮定として受けていた。★★★★**本ファイルはその 2 つを定理にする**。

## ★★★機構

### 有限素点 —— 局所環のスペクトルは閉点の近傍に収まる

`Spec R`（`R` 局所環）の**すべての点は閉点に特殊化する**
（`IsLocalRing.specializes_closedPoint`）。
★開集合は**生成化で閉じている**（`Specializes.mem_open`）から、
閉点の像を含む開集合は**像全体を含む**。
★★あとは `IsOpenImmersion.lift` で射が持ち上がる。

### 無限素点 —— `ℙᴺ` の複素点は必ずどれかのチャートを通る

`§9-870` の `chartIndexOf`・`specMap_projChartHom` がそのまま与える。

## ★これで `Proposition 1.4, (iv)` に残った点ごとの仮定

★★★残るのは**座標とチャートの整合**だけである:

* `hspan`——`𝔞·(𝓞_F)_Q = (x_{i(Q)})`
* `hw`——`(s_0/s_{i(Q)})(y_Q)·x_{i(Q)} = x_0`
* `hcv`——`σ_v(x_k) = ρ_v(x_k/x_{i(v)})·σ_v(x_{i(v)})`

★これらは「座標 `x` を点から**構成する**」段と一体である
——`x` を点から作れば、整合は構成の定義から出る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite ABC3.Found.Arakelov NumberField MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★有限素点 —— 局所環の点はチャートを通る -/

/-- ★★★★★★★★★★★★**局所環のスペクトルからの射は、被覆のどれかを通る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は 2 つだけ:
* `Spec R`（`R` 局所環）の**すべての点は閉点に特殊化する**
* 開集合は**生成化で閉じている**——だから閉点の像を含む開集合は像全体を含む

★★これが「素点ごとのチャートは常に取れる」ことの中身である。 -/
theorem exists_chart_of_isLocalRing (M : X.PresheafOfModules) {N : ℕ}
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (R : Type) [CommRing R] [IsLocalRing R]
    (g : Spec (CommRingCat.of R) ⟶ X) :
    ∃ (i : Fin (N + 1)) (y : Spec (CommRingCat.of R) ⟶ (nonVanishing M (s i)).toScheme),
      g = y ≫ (nonVanishing M (s i)).ι := by
  have hcl : g.base (IsLocalRing.closedPoint R) ∈ (⊤ : X.Opens) := trivial
  rw [← hcov] at hcl
  obtain ⟨i, hi⟩ := TopologicalSpace.Opens.mem_iSup.mp hcl
  refine ⟨i, IsOpenImmersion.lift (nonVanishing M (s i)).ι g ?_,
    (IsOpenImmersion.lift_fac _ _ _).symm⟩
  rw [Scheme.Opens.range_ι]
  rintro _ ⟨z, rfl⟩
  have hsp : z ⤳ IsLocalRing.closedPoint R := IsLocalRing.specializes_closedPoint z
  exact (hsp.map g.base.hom.continuous).mem_open (nonVanishing M (s i)).2 hi

/-- ★★★★★★★★★★★★★**素点ごとのチャートの族は常に取れる**。

★`exists_chart_of_isLocalRing` を各極大イデアルで選ぶだけである。 -/
theorem exists_localChart_family (M : X.PresheafOfModules) {N : ℕ}
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X) :
    ∃ (chartOf : ∀ (Q : Ideal (𝓞 F)), Q.IsMaximal → Fin (N + 1))
      (y : ∀ (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal),
        Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶
          (nonVanishing M (s (chartOf Q hQ))).toScheme),
      ∀ (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal),
        Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF
          = y Q hQ ≫ (nonVanishing M (s (chartOf Q hQ))).ι := by
  have h : ∀ (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal),
      ∃ (i : Fin (N + 1)) (y : Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶
        (nonVanishing M (s i)).toScheme),
        Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF
          = y ≫ (nonVanishing M (s i)).ι := by
    intro Q hQ
    exact exists_chart_of_isLocalRing M s hcov (Localization.AtPrime Q) _
  choose chartOf y hy using h
  exact ⟨chartOf, y, hy⟩

/-! ## ★★★★★★★無限素点 —— `ℙᴺ` の複素点は必ずチャートを通る -/

/-- ★★★★★★★**`ℙᴺ` の複素点は必ずどれかのチャートを通る**。

★`§9-870` の `chartIndexOf`・`specMap_projChartHom` がそのまま与える。 -/
theorem exists_chartA_factorization (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :
    ∃ (i : Fin (N + 1)) (ρ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ), p = Spec.map ρ ≫ chartA N ℤ i :=
  ⟨chartIndexOf N p, projChartHom N ℤ ℂ p (chartIndexOf N p) (chartIndexOf_spec N p),
    (specMap_projChartHom N ℤ ℂ p (chartIndexOf N p) (chartIndexOf_spec N p)).symm⟩

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★**無限素点ごとのチャートの族は常に取れる**。 -/
theorem exists_archChart_family (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X) :
    ∃ (archChart : InfinitePlace F → Fin (N + 1))
      (ρ : ∀ v : InfinitePlace F, CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X (archChart v)))
        ⟶ CommRingCat.of ℂ),
      ∀ v : InfinitePlace F, archPoint xF v ≫ globalToProj M hM φ s hcov
        = Spec.map (ρ v) ≫ chartA N ℤ (archChart v) := by
  have h : ∀ v : InfinitePlace F, ∃ (i : Fin (N + 1))
      (ρ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ),
      archPoint xF v ≫ globalToProj M hM φ s hcov = Spec.map ρ ≫ chartA N ℤ i :=
    fun v => exists_chartA_factorization N (archPoint xF v ≫ globalToProj M hM φ s hcov)
  choose archChart ρ hfac using h
  exact ⟨archChart, ρ, hfac⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_chart_of_isLocalRing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所環のスペクトルからの射は被覆のどれかを通る)",
    sectionId := "genell-prop-1-4" }

def exists_localChart_family.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(素点ごとのチャートの族は常に取れる)",
    sectionId := "genell-prop-1-4" }

def exists_chartA_factorization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ℙᴺ の複素点は必ずどれかのチャートを通る)",
    sectionId := "genell-prop-1-4" }

def exists_archChart_family.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(無限素点ごとのチャートの族は常に取れる)",
    sectionId := "genell-prop-1-4" }

def exists_localChart_family.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsLocalRing.specializes_closedPoint(すべての点は閉点に特殊化する)"
      (.inMathlib "IsLocalRing.specializes_closedPoint") 2,
    .citation "[mathlib]" "Specializes.mem_open(開集合は生成化で閉じている)"
      (.inMathlib "Specializes.mem_open") 2,
    .citation "[mathlib]" "IsOpenImmersion.lift(像が入っていれば射は持ち上がる)"
      (.inMathlib "AlgebraicGeometry.IsOpenImmersion.lift") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-935〜938 が仮定として受けていた" ++
       "「素点ごとのチャート」は**定理である**。" ++
       "Spec R(R 局所環)のすべての点は閉点に特殊化し、開集合は生成化で閉じているので、" ++
       "閉点の像を含む開集合は像全体を含むからである") 5,
    .implicitStep
      ("★★残る点ごとの仮定は**座標とチャートの整合**だけである" ++
       "(hspan・hw・hcv)——これらは「座標 x を点から構成する」段と一体であり、" ++
       "x を点から作れば整合は構成の定義から出る") 4 ]

end ABC3.Found.GenEll
