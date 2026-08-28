/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjCoordProportional
import ABC3.Found.GenEll.ArcModel
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★`ℙᴺ_ℤ` の射影モデル（`ArcModel`）（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

## ★★★★★★★★★★★★★★これは何か —— `ArcModel` を実際に 1 つ作る

`Found/GenEll/ArcModel.lean` は「射影モデルは**与えられたものとして受ける**」と書き、
**`X` ごとに構成することは主張しない**としていた。★本ファイルは
`X = ℙᴺ_ℤ` について**実際に構成する**。

    `projArcModel N : ArcModel (ℙᴺ_ℤ) (Fin (N+1) → ℂ)`

★★これで `Proposition 1.6` のアルキメデス側と `Proposition 1.4, (iii)` が、
`ℙᴺ` については**仮定なしで**使えるようになる。

## ★★★機構 —— 4 つの欄はすべて既存の部品

| 欄 | 部品 |
|---|---|
| `emb` | `projPointCoord`（`§9-C2b`）——射影的に well-defined（`§9-874`） |
| `cone` | `univ`（閉錐であることは自明） |
| `emb_injective` | ★`ext_of_projCoord`（`§9-850`）＋ 座標の比較 |
| `emb_range` | ★★`projPointOfCoords`（`§9-873`）で全射 |

★★★本ファイルの新しい段は 2 つ:
`range_of_projCoord_ne_zero`（`§9-871` の**逆向き**——座標が `0` でなければチャートに入る）と、
`eq_of_embVec_smul`（座標ベクトルが定数倍なら点は等しい）。

## ★★★★測定の記録

★`ℤ` からの環準同型は一意（`Subsingleton (ℤ →+* R)`）なので、
`ext_of_projCoord` の「定数での一致」は**自動**である（`hom_awayConst_eq`）。
★★これは基礎環を `ℤ` に取っていることの利得で、一般の `R` では別に要る（2026-08-28 実測）。
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory

open scoped LinearAlgebra.Projectivization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★座標が `0` でなければチャートに入る（`§9-871` の逆向き） -/

/-- ★★★**座標が `0` でなければチャートに入る** —— `§9-871` の `projCoord_ne_zero_of_range`
の逆向き。★機構は `specMap_mem_basicOpen_iff` が**同値**であることだけ。 -/
theorem range_of_projCoord_ne_zero (N : ℕ) (i i' : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ)
    (hne : φ.hom (projCoord N ℤ i i') ≠ 0) :
    Set.range (Spec.map φ ≫ chartA N ℤ i).base ⊆ Set.range (chartA N ℤ i').base := by
  haveI : Subsingleton (Spec (CommRingCat.of ℂ)) :=
    inferInstanceAs (Subsingleton (PrimeSpectrum ℂ))
  have hpre : (Spec.map φ).base (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ))
      ∈ PrimeSpectrum.basicOpen (projCoord N ℤ i i') :=
    (specMap_mem_basicOpen_iff φ _).2 hne
  rw [← isLocalizationElem_eq_projCoord,
    ← Proj.awayι_preimage_basicOpen _ (MvPolynomial.isHomogeneous_X ℤ i) one_pos
      (MvPolynomial.isHomogeneous_X ℤ i') one_pos] at hpre
  have hmem : (Spec.map φ ≫ chartA N ℤ i).base
      (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ))
      ∈ Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i') := hpre
  rintro _ ⟨y, rfl⟩
  have hy : y = (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ)) :=
    Subsingleton.elim _ _
  subst hy
  rw [← Proj.opensRange_awayι (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
    (MvPolynomial.X i') (MvPolynomial.isHomogeneous_X ℤ i') one_pos] at hmem
  exact hmem

/-! ## ★★チャートの射は一意である -/

/-- ★★**`projChartHom` は分解を与える唯一の射である**。 -/
theorem projChartHom_unique (N : ℕ) (R : Type) [CommRing R] (F : Type) [Field F]
    (x : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))
    (i : Fin (N + 1))
    (hx : Set.range x.base ⊆ Set.range (chartA N R i).base)
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i))
      ⟶ CommRingCat.of F)
    (hφ : Spec.map φ ≫ chartA N R i = x) :
    projChartHom N R F x i hx = φ := by
  apply Spec.map_injective
  rw [← cancel_mono (chartA N R i), specMap_projChartHom, hφ]

theorem projPointCoord_of_hom (N : ℕ) (i : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ)
    (hx : Set.range (Spec.map φ ≫ chartA N ℤ i).base ⊆ Set.range (chartA N ℤ i).base)
    (k : Fin (N+1)) :
    projPointCoord N ℤ ℂ (Spec.map φ ≫ chartA N ℤ i) i hx k = φ.hom (projCoord N ℤ i k) := by
  rw [projPointCoord, projChartHom_unique N ℤ ℂ _ i hx φ rfl]

/-! ## ★★★定数での一致は自動である -/

/-- ★**`awayConst` は整数のキャストである** —— `ℤ` からの環準同型は一意だから。 -/
theorem awayConst_eq_intCast (N : ℕ) (j : Fin (N+1)) (c : ℤ) :
    awayConst N ℤ j c = ((c : ℤ) : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X j)) := by
  have h1 : (awayConstHom N ℤ j : ℤ →+* _) = Int.castRingHom _ := Subsingleton.elim _ _
  exact congrArg (fun f : ℤ →+* _ => f c) h1

/-- ★★**したがって定数での一致は自動である** —— `ext_of_projCoord` の第 1 仮定。 -/
theorem hom_awayConst_eq (N : ℕ) (j : Fin (N+1))
    (ψ χ : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X j) →+* ℂ)
    (c : ℤ) : ψ (awayConst N ℤ j c) = χ (awayConst N ℤ j c) := by
  rw [awayConst_eq_intCast, map_intCast, map_intCast]

/-! ## ★★★★★★★★★★座標ベクトルと単射性 -/

/-- ★**複素点の座標ベクトル**（選んだチャートで正規化したもの）。 -/
noncomputable def embVec (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :
    Fin (N+1) → ℂ :=
  projPointCoord N ℤ ℂ p (chartIndexOf N p) (chartIndexOf_spec N p)

theorem embVec_self (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :
    embVec N p (chartIndexOf N p) = 1 :=
  projPointCoord_self N ℤ ℂ p (chartIndexOf N p) (chartIndexOf_spec N p)

theorem embVec_ne_zero (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :
    embVec N p ≠ 0 := by
  intro h
  have h1 := embVec_self N p
  rw [h] at h1
  exact one_ne_zero h1.symm

/-- ★★★★★★★★★★**座標ベクトルが定数倍なら点は等しい**。

★機構は `§9-850` の `ext_of_projCoord`（座標と定数で環準同型が決まる）で、
定数の側は `ℤ` からの環準同型の一意性で自動である。 -/
theorem eq_of_embVec_smul (N : ℕ)
    (p q : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (c : ℂ) (hc : c ≠ 0) (h : ∀ k, embVec N p k = c * embVec N q k) : p = q := by
  have hpj : embVec N p (chartIndexOf N q) = c := by
    rw [h (chartIndexOf N q), embVec_self, mul_one]
  have hpjne : embVec N p (chartIndexOf N q) ≠ 0 := by rw [hpj]; exact hc
  have hxp : Set.range p.base ⊆ Set.range (chartA N ℤ (chartIndexOf N q)).base := by
    have hfac := specMap_projChartHom N ℤ ℂ p (chartIndexOf N p) (chartIndexOf_spec N p)
    have hsub := range_of_projCoord_ne_zero N (chartIndexOf N p) (chartIndexOf N q)
      (projChartHom N ℤ ℂ p (chartIndexOf N p) (chartIndexOf_spec N p)) hpjne
    rwa [hfac] at hsub
  have hcoord : ∀ k, projPointCoord N ℤ ℂ p (chartIndexOf N q) hxp k = embVec N q k := by
    intro k
    have hcg : embVec N p k
        = projPointCoord N ℤ ℂ p (chartIndexOf N q) hxp k * embVec N p (chartIndexOf N q) :=
      projPointCoord_congr N p (chartIndexOf N p) (chartIndexOf N q)
        (chartIndexOf_spec N p) hxp k
    rw [hpj, h k, mul_comm c (embVec N q k)] at hcg
    exact (mul_right_cancel₀ hc hcg).symm
  have hhom : (projChartHom N ℤ ℂ p (chartIndexOf N q) hxp).hom
      = (projChartHom N ℤ ℂ q (chartIndexOf N q) (chartIndexOf_spec N q)).hom :=
    ext_of_projCoord (chartIndexOf N q) _ _
      (fun c => hom_awayConst_eq N (chartIndexOf N q) _ _ c) (fun k => hcoord k)
  calc p = Spec.map (projChartHom N ℤ ℂ p (chartIndexOf N q) hxp)
        ≫ chartA N ℤ (chartIndexOf N q) :=
      (specMap_projChartHom N ℤ ℂ p (chartIndexOf N q) hxp).symm
    _ = Spec.map (projChartHom N ℤ ℂ q (chartIndexOf N q) (chartIndexOf_spec N q))
        ≫ chartA N ℤ (chartIndexOf N q) := by rw [CommRingCat.hom_ext hhom]
    _ = q := specMap_projChartHom N ℤ ℂ q (chartIndexOf N q) (chartIndexOf_spec N q)

/-! ## ★★★★★★★★★★★★★★射影モデル -/

/-- ★★**複素点を `ℙ(ℂ^{N+1})` へ送る**。 -/
noncomputable def projEmb (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :
    ℙ ℂ (Fin (N+1) → ℂ) :=
  Projectivization.mk ℂ (embVec N p) (embVec_ne_zero N p)

theorem projEmb_injective (N : ℕ) : Function.Injective (projEmb N) := by
  intro p q hpq
  rw [projEmb, projEmb, Projectivization.mk_eq_mk_iff] at hpq
  obtain ⟨a, ha⟩ := hpq
  refine eq_of_embVec_smul N p q (a : ℂ) (Units.ne_zero a) (fun k => ?_)
  rw [← ha]
  rfl

theorem projEmb_surjective (N : ℕ) : Function.Surjective (projEmb N) := by
  intro P
  obtain ⟨i, hi⟩ := Function.ne_iff.1 (Projectivization.rep_nonzero P)
  refine ⟨projPointOfCoords N P.rep i hi, ?_⟩
  have hxi : Set.range (projPointOfCoords N P.rep i hi).base
      ⊆ Set.range (chartA N ℤ i).base := by
    rintro _ ⟨y, rfl⟩
    exact ⟨(Spec.map (CommRingCat.ofHom (awayHomOfCoords N P.rep i hi))).base y, rfl⟩
  have hcoord : ∀ k, projPointCoord N ℤ ℂ (projPointOfCoords N P.rep i hi) i hxi k
      = P.rep k / P.rep i := by
    intro k
    have hk := projPointCoord_of_hom N i
      (CommRingCat.ofHom (awayHomOfCoords N P.rep i hi)) hxi k
    rw [CommRingCat.hom_ofHom, awayHomOfCoords_projCoord] at hk
    exact hk
  have hne : embVec N (projPointOfCoords N P.rep i hi) i ≠ 0 :=
    projPointCoord_ne_zero N (projPointOfCoords N P.rep i hi)
      (chartIndexOf N (projPointOfCoords N P.rep i hi)) i
      (chartIndexOf_spec N (projPointOfCoords N P.rep i hi)) hxi
  have hrel : ∀ k, embVec N (projPointOfCoords N P.rep i hi) k
      = (embVec N (projPointOfCoords N P.rep i hi) i / P.rep i) * P.rep k := by
    intro k
    have hcg : embVec N (projPointOfCoords N P.rep i hi) k
        = projPointCoord N ℤ ℂ (projPointOfCoords N P.rep i hi) i hxi k
          * embVec N (projPointOfCoords N P.rep i hi) i :=
      projPointCoord_congr N (projPointOfCoords N P.rep i hi)
        (chartIndexOf N (projPointOfCoords N P.rep i hi)) i
        (chartIndexOf_spec N (projPointOfCoords N P.rep i hi)) hxi k
    rw [hcg, hcoord k]
    field_simp
  refine Eq.trans ?_ (Projectivization.mk_rep P)
  rw [projEmb, Projectivization.mk_eq_mk_iff]
  refine ⟨Units.mk0 (embVec N (projPointOfCoords N P.rep i hi) i / P.rep i)
    (div_ne_zero hne hi), ?_⟩
  funext k
  exact (hrel k).symm

/-- ★★★★★★★★★★★★★★**`ℙᴺ_ℤ` の射影モデル**。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★`ArcModel.lean` が「与えられたものとして受ける」としていたデータを、
`X = ℙᴺ_ℤ` について**実際に構成した**ものである。
★★これで `Proposition 1.6` のアルキメデス側と `Proposition 1.4, (iii)` が
`ℙᴺ` については仮定なしで使える。 -/
noncomputable def projArcModel (N : ℕ) :
    ArcModel (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) (Fin (N+1) → ℂ) where
  emb := projEmb N
  cone := Set.univ
  cone_closed := isClosed_univ
  cone_isCone := fun _ _ _ => Set.mem_univ _
  emb_injective := projEmb_injective N
  emb_range := by
    rw [Set.range_eq_univ.2 (projEmb_surjective N)]
    ext P
    simp

/-! ## ★出典の紐付け(`.src`) -/

def range_of_projCoord_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(座標が 0 でなければチャートに入る)",
    sectionId := "genell-def-1-1" }

def eq_of_embVec_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(座標ベクトルが定数倍なら点は等しい)",
    sectionId := "genell-def-1-1" }

def projEmb.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(複素点を ℙ(ℂ^{N+1}) へ送る)",
    sectionId := "genell-def-1-1" }

def projArcModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(ℙᴺ_ℤ の射影モデル)",
    sectionId := "genell-def-1-1" }

def projArcModel.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ext_of_projCoord(座標と定数で環準同型が決まる、§9-850)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ext_of_projCoord") 2,
    .citation "[ABC3]" "projPointCoord_congr(座標はチャートを変えると定数倍、§9-874)"
      (.inProject "ABC3" "ABC3.Found.GenEll.projPointCoord_congr") 2,
    .citation "[ABC3]" "projPointOfCoords(座標の組から複素点、§9-873)"
      (.inProject "ABC3" "ABC3.Found.GenEll.projPointOfCoords") 2,
    .implicitStep
      ("★ArcModel.lean は「射影モデルは**与えられたものとして受ける**」と書き、" ++
       "X ごとに構成することは主張しない としていた。本ファイルは X = ℙᴺ_ℤ について" ++
       "**実際に構成する**") 3,
    .implicitStep
      ("★★測定: ℤ からの環準同型は一意(Subsingleton (ℤ →+* R))なので、" ++
       "ext_of_projCoord の「定数での一致」は**自動**である(hom_awayConst_eq)。" ++
       "これは基礎環を ℤ に取っていることの利得で、一般の R では別に要る(2026-08-28 実測)") 3,
    .implicitStep
      ("★★★これで段 C2c-3 に残るのは、greenFS と与えられた計量の**差**が " ++
       "projArcModel の位相で連続であることの確認だけである") 4 ]

end ABC3.Found.GenEll
