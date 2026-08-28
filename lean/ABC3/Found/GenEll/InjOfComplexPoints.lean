/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottLocal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★`hinj` は複素点の単射性から出る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★これは何か —— 点の側に残った最後の仮定

`§9-935` まで来て、`Proposition 1.4, (iv)` の点の側に残った仮定は

    `hinj : Function.Injective (fun p k => ((x p k / x p idx : F_p) : ℂ))`

——**正規化座標が点を分けること**——だけである。

★★★★本ファイルはそれを**幾何の言葉**に翻訳する:

    `p ↦ (点 p の定める ℙᴺ の複素点)` が単射  ⟹  `hinj`

★★★★★これは `ψ` が埋め込み（`§9-920`）であることの点版であり、
「点の族 `P` が本当に相異なる点の族である」という**当たり前の条件**になった。

## ★★★機構 —— 座標ベクトルは射影的に決まっている

`§9-874`（`projPointCoord_congr`）と `§9-875`（`eq_of_embVec_smul`）で
**座標ベクトルが定数倍なら複素点は等しい**ことが取れている。

★局所チャートの分解 `archPoint(x_p) ≫ ψ = Spec ρ_p ≫ chartA i_p` と
`σ_v(x_p k) = ρ_p(x_k/x_{i_p})·σ_v(x_p i_p)` から

    `embVec(点 p) = a_p · σ_v(x_p ·)`   （`a_p ≠ 0`）

が出るので、正規化座標が一致すれば `embVec` が定数倍で一致し、複素点が一致する。

## ★残っている仮定（明示）

★★`hemb`——選んだ無限素点 `v_p` が `F_p ↪ ℂ` の包含が定めるものであること
（`archRingHom (F_p) (v_p) = (包含)`）。
★これは `InfinitePlace.mk` の埋め込みが元の埋め込みか**その複素共役**かの差であり、
共役の場合も（共役が単射だから）結論は変わらない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★無限素点での値は 0 でない -/

/-- ★**`x ≠ 0` なら `σ_v(x) ≠ 0`**。 -/
theorem archRingHom_ne_zero (F : Type) [Field F] [NumberField F] (v : InfinitePlace F)
    {a : 𝓞 F} (ha : a ≠ 0) : archRingHom F v a ≠ 0 := by
  intro hc
  have hpos : (0:ℝ) < v ((a : F)) := (AbsoluteValue.pos_iff _).2 (by simpa using ha)
  have h1 : ‖archRingHom F v a‖ = v ((a : F)) := norm_archRingHom F v a
  rw [hc, norm_zero] at h1
  exact hpos.ne' h1.symm

/-! ## ★★★★★★★★★★★座標ベクトルが比例していれば単射性が移る -/

/-- ★★★★★★★★★★★**座標ベクトルが定数倍で与えられていれば、
正規化座標の単射性は複素点の単射性から出る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は `§9-875`（`eq_of_embVec_smul`——座標ベクトルが定数倍なら点は等しい）だけである。 -/
theorem injective_of_embVec_prop (N : ℕ) {P : Type}
    (A : P → complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (c : P → Fin (N+1) → ℂ) (idx : Fin (N+1)) (hidx : ∀ p, c p idx ≠ 0)
    (a : P → ℂ) (ha : ∀ p, a p ≠ 0)
    (hprop : ∀ p k, embVec N (A p) k = a p * c p k)
    (hA : Function.Injective A) :
    Function.Injective (fun (p : P) (k : Fin (N+1)) => c p k / c p idx) := by
  intro p q h
  refine hA ?_
  have hk : ∀ k, c p k / c p idx = c q k / c q idx := fun k => congrFun h k
  have hne : a q * c q idx ≠ 0 := mul_ne_zero (ha q) (hidx q)
  refine eq_of_embVec_smul N (A p) (A q) ((a p * c p idx) / (a q * c q idx))
    (div_ne_zero (mul_ne_zero (ha p) (hidx p)) hne) (fun k => ?_)
  rw [hprop, hprop]
  have h1 : c p k = c p idx * (c q k / c q idx) := by
    rw [← hk k, mul_div_cancel₀ _ (hidx p)]
  rw [h1]
  field_simp [ha q, ha p, hidx p, hidx q]

/-! ## ★★★★★★★★★★★★チャートから座標ベクトルの比例を作る -/

/-- ★★★★★★★★★★★★**チャートを通る複素点の座標ベクトルは `σ_v(x_·)` に比例する**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-874`（座標はチャートを変えると定数倍）で `embVec` をチャート `i` の座標に直し、
仮定 `hcv`（`σ_v(x_k) = ρ(x_k/x_i)·σ_v(x_i)`）を当てる。 -/
theorem embVec_chart_prop (N : ℕ) (F : Type) [Field F] [NumberField F]
    (x : Fin (N + 1) → 𝓞 F) (v : InfinitePlace F) (i : Fin (N + 1))
    (ρ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ)
    (hcv : ∀ k, archRingHom F v (x k) = ρ.hom (projCoord N ℤ i k) * archRingHom F v (x i))
    (hi : x i ≠ 0) :
    ∃ a : ℂ, a ≠ 0 ∧ ∀ k, embVec N (Spec.map ρ ≫ chartA N ℤ i) k
      = a * archRingHom F v (x k) := by
  set Apt := Spec.map ρ ≫ chartA N ℤ i with hApt
  have hx : Set.range Apt.base ⊆ Set.range (chartA N ℤ i).base := by
    rintro _ ⟨y, rfl⟩
    exact ⟨(Spec.map ρ).base y, rfl⟩
  have hco : ∀ k, projPointCoord N ℤ ℂ Apt i hx k = ρ.hom (projCoord N ℤ i k) :=
    projPointCoord_of_hom N i ρ hx
  have hcg : ∀ k, embVec N Apt k = projPointCoord N ℤ ℂ Apt i hx k * embVec N Apt i :=
    fun k => projPointCoord_congr N Apt (chartIndexOf N Apt) i (chartIndexOf_spec N Apt) hx k
  have hine : embVec N Apt i ≠ 0 :=
    projPointCoord_ne_zero N Apt (chartIndexOf N Apt) i (chartIndexOf_spec N Apt) hx
  have hσi : archRingHom F v (x i) ≠ 0 := archRingHom_ne_zero F v hi
  refine ⟨embVec N Apt i / archRingHom F v (x i), div_ne_zero hine hσi, fun k => ?_⟩
  rw [hcv k, hcg k, hco k]
  field_simp

/-! ## ★★★★★★★★★★★★★★★★到達 —— `hinj` は複素点の単射性である -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★**`hinj` は「点の族が相異なる複素点を与えること」である**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★`Proposition 1.4, (iv)` の点の側に残っていた `hinj`（正規化座標が点を分けること）は、
**`p ↦ (点 p の定める ℙᴺ の複素点)` の単射性**に他ならない
——`ψ` が埋め込み（`§9-920`）であることの点版である。 -/
theorem hinj_of_injective_archPoint (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
    (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
    (v : ∀ p, haveI := hnf p; InfinitePlace (fld p))
    (archChart : P → Fin (N + 1))
    (ρ : ∀ p, CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ) (MvPolynomial.X (archChart p)))
      ⟶ CommRingCat.of ℂ)
    (hfac : ∀ p, haveI := hnf p; archPoint (xF p) (v p) ≫ globalToProj M hM φ s hcov
      = Spec.map (ρ p) ≫ chartA N ℤ (archChart p))
    (hcv : ∀ p, haveI := hnf p; ∀ k, archRingHom (fld p) (v p) (x p k)
      = (ρ p).hom (projCoord N ℤ (archChart p) k)
        * archRingHom (fld p) (v p) (x p (archChart p)))
    (hiv : ∀ p, x p (archChart p) ≠ 0)
    (idx : Fin (N + 1)) (hidx : ∀ p, x p idx ≠ 0)
    (hemb : ∀ p, haveI := hnf p; ∀ k,
      archRingHom (fld p) (v p) (x p k) = ((x p k : fld p) : ℂ))
    (hpt : Function.Injective (fun p => haveI := hnf p;
      archPoint (xF p) (v p) ≫ globalToProj M hM φ s hcov)) :
    Function.Injective (fun (p : P) (k : Fin (N + 1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)) := by
  choose a ha hprop using fun p => haveI := hnf p;
    embVec_chart_prop N (fld p) (x p) (v p) (archChart p) (ρ p) (hcv p) (hiv p)
  have hcidx : ∀ p, haveI := hnf p; (((x p idx : fld p)) : ℂ) ≠ 0 := by
    intro p
    haveI := hnf p
    rw [← hemb p idx]
    exact archRingHom_ne_zero (fld p) (v p) (hidx p)
  have key := injective_of_embVec_prop N
    (fun p => haveI := hnf p; archPoint (xF p) (v p) ≫ globalToProj M hM φ s hcov)
    (fun p k => haveI := hnf p; (((x p k : fld p)) : ℂ)) idx hcidx a ha
    (fun p k => by
      haveI := hnf p
      rw [hfac p, hprop p k, hemb p k]) hpt
  intro p q hpq
  refine key ?_
  funext k
  have h1 := congrFun hpq k
  haveI := hnf p
  haveI := hnf q
  show (((x p k : fld p)) : ℂ) / (((x p idx : fld p)) : ℂ)
    = (((x q k : fld q)) : ℂ) / (((x q idx : fld q)) : ℂ)
  push_cast at h1
  exact h1

/-! ## ★出典の紐付け(`.src`) -/

def archRingHom_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(x ≠ 0 なら無限素点での値は 0 でない)",
    sectionId := "genell-prop-1-4" }

def injective_of_embVec_prop.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(座標ベクトルが定数倍なら正規化座標の単射性が移る)",
    sectionId := "genell-prop-1-4" }

def embVec_chart_prop.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートを通る複素点の座標ベクトルは σ_v(x) に比例する)",
    sectionId := "genell-prop-1-4" }

def hinj_of_injective_archPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(hinj は点の族が相異なる複素点を与えることである)",
    sectionId := "genell-prop-1-4" }

def hinj_of_injective_archPoint.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "eq_of_embVec_smul(座標ベクトルが定数倍なら点は等しい、§9-875)"
      (.inProject "ABC3" "ABC3.Found.GenEll.eq_of_embVec_smul") 3,
    .citation "[ABC3]" "projPointCoord_congr(座標はチャートを変えると定数倍、§9-874)"
      (.inProject "ABC3" "ABC3.Found.GenEll.projPointCoord_congr") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): Proposition 1.4, (iv) の点の側に残っていた hinj" ++
       "(正規化座標が点を分けること)は**幾何の言葉に翻訳できる**" ++
       "——p ↦ (点 p の定める ℙᴺ の複素点)が単射であること。" ++
       "★これは ψ が埋め込み(§9-920)であることの点版であり、" ++
       "「点の族 P が本当に相異なる点の族である」という当たり前の条件である") 5,
    .implicitStep
      ("★hemb(選んだ無限素点 v_p が F_p ↪ ℂ の包含が定めるものであること)は" ++
       "InfinitePlace.mk の埋め込みが元の埋め込みかその複素共役かの差であり、" ++
       "共役の場合も(共役が単射だから)結論は変わらない") 2 ]

end ABC3.Found.GenEll
