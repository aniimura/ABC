/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottDivisorOfSection
import ABC3.Found.GenEll.FinCover
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★`IsAmple` から Northcott まで一本（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★これは何か —— 鎖を 1 本にする

`§9-914`（`IsAmple` から被覆）から `§9-927`（切断の零点因子についての Northcott）までを
**1 つの定理**にまとめる:

    `IsAmple M`  ⟹  ある `L > 0`・`N`・切断族 `s` があって
                    「`D` が `s₀` の零点因子＋引き戻した計量なら `ht_D` について Northcott」

★これが原文 `Proposition 1.4, (iv)` の**構成の側での到達点**である。

## ★★★残っている仮定（明示）—— 2 つだけ

### (1) 幾何の側 —— 有限アフィン自明化被覆

`U`・`hcovU`・`hU`・`eU`。★`X` が `Spec A` 上分離的・準コンパクトで `M` が可逆なら
出るはずのものであり、**Arakelov 理論ではなくスキーム論の一般論**である。

### (2) 点の側 —— `chart` と `hinj`

* **`chart`・`xF`**——点 `Spec 𝓞_F ⟶ X` が**どれかのチャートを通る**こと。
  ★★★★**測定（2026-08-29）**: これは単なる技術的仮定ではない。
  `Spec 𝓞_F ⟶ ℙᴺ_ℤ` がチャート `D₊(x_i)` を通るのは
  「座標の生成するイデアル `(x_0,…,x_N)` が `x_i` で生成される」ときに限る
  ——**イデアル類が自明でないと通らない**。
  ★原文の高さの計算（素点ごとの寄与の和）はこの制限を受けないので、
  外すには **`pullbackIdeal` を素点ごとに読む**段が要る。
* **`hinj`**——座標が点を分けること。★`ψ` が埋め込み（`§9-920`）であることの点版であり、
  点の族 `P` を座標で添字づけるための道具である。

★★★★★どちらも `Proposition 1.4, (iv)` の**最後の 2 つ**であり、
幾何の側（Serre 道 56 段）は `A2`（そもそも塡がりでない）を除いて閉じている。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★**`IsAmple` から Northcott まで一本**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `IsAmple M` ⟹ ある `L > 0`・`N`・切断族 `s` があって
                  「`D` が `s₀` の零点因子＋引き戻した計量なら `ht_D` について Northcott」

★★`§9-914`（`IsAmple` から被覆）→ `§9-920`（埋め込み）→ `§9-926`（段 E0）→
`§9-927`（Northcott）の鎖を 1 本にしたものである。 -/
theorem northcott_of_isAmple (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {ι : Type} (U : ι → X.Opens) (hcovU : (⨆ j, U j) = ⊤)
    (hU : ∀ j, IsAffineOpen (U j))
    (eU : ∀ j, (restrictPresheafFunctor X (U j)).obj M ≅ 𝟙_ (PresheafModulesOn X (U j)))
    (hample : IsAmple M) (x₀ : (X : Type))
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type)) (d : ℕ) :
    ∃ (L : ℕ) (_ : 0 < L) (N : ℕ)
      (s : Fin (N + 1) → ((presheafTensorPow M L).obj (op ⊤) : Type))
      (hcov : (⨆ k, nonVanishing (presheafTensorPow M L) (s k)) = ⊤),
      (∀ i, IsAffineOpen (nonVanishing (presheafTensorPow M L) (s i))) ∧
      ∀ (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) ℤ)))
        (_ : E.divisor = hyperplaneIdeal N ℤ)
        (_ : @Continuous _ _ (projArcModel N).topology _
          (fun p => E.green p - greenFS N p))
        (D : ArithCartier X)
        (_ : D.divisor = divisorOfSection (presheafTensorPow M L) (s 0))
        (_ : D.green = fun p => E.green (p ≫ globalToProj (presheafTensorPow M L)
          (isLocallyTrivial_presheafTensorPow hM L) φ s hcov))
        {P : Type}
        (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
        (_ : ∀ p, Module.finrank ℚ (fld p) ≤ d)
        (chart : P → Fin (N + 1))
        (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶
          (nonVanishing (presheafTensorPow M L) (s (chart p))).toScheme)
        (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
        (_ : ∀ p, x p ≠ 0)
        (_ : ∀ p, haveI := hnf p; ∀ k, x p k =
          (Spec.preimage (xF p ≫ globalChartMorphism (presheafTensorPow M L)
            (isLocallyTrivial_presheafTensorPow hM L) φ s (chart p))).hom
            (projCoord N ℤ (chart p) k) * x p (chart p))
        (_ : ∀ p, x p 0 ≠ 0)
        (idx : Fin (N + 1))
        (_ : Function.Injective (fun (p : P) (k : Fin (N + 1)) =>
          ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
        (C : ℝ),
        {p : P | haveI := hnf p;
          htArith (fld p) D
            (xF p ≫ (nonVanishing (presheafTensorPow M L) (s (chart p))).ι) ≤ C}.Finite := by
  obtain ⟨L, hL, N, s, hcov, haffs⟩ := exists_fin_cover_of_isAmple M hM hample x₀
  refine ⟨L, hL, N, s, hcov, haffs, ?_⟩
  intro E hdiv hcont D hDdiv hDgreen P fld hnf hdeg chart xF x hx hw h0 idx hinj C
  exact northcott_of_divisorOfSection (presheafTensorPow M L)
    (isLocallyTrivial_presheafTensorPow hM L) d φ s hcov haffs
    U hcovU hU (fun j => tensorPowTriv (eU j) L) E hdiv hcont D hDdiv hDgreen
    fld hnf hdeg chart xF x hx hw h0 idx hinj C

/-! ## ★出典の紐付け(`.src`) -/

def northcott_of_isAmple.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(IsAmple から Northcott まで一本)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_fin_cover_of_isAmple(IsAmple から被覆、§9-914・915)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_fin_cover_of_isAmple") 2,
    .citation "[ABC3]" "northcott_of_divisorOfSection(切断の零点因子についての Northcott、§9-927)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_divisorOfSection") 6,
    .implicitStep
      ("★★★幾何の側の仮定は「有限アフィン自明化被覆」(U・hcovU・hU・eU)に集約されている" ++
       "——X が Spec A 上分離的・準コンパクトで M が可逆なら出るはずのもので、" ++
       "**Arakelov 理論ではなくスキーム論の一般論**である") 4,
    .implicitStep
      ("★★★★測定(2026-08-29): 点の側の仮定 chart(点がどれかのチャートを通ること)は" ++
       "**単なる技術的仮定ではない**。Spec 𝓞_F ⟶ ℙᴺ_ℤ がチャート D₊(x_i) を通るのは" ++
       "「座標の生成するイデアル (x_0,…,x_N) が x_i で生成される」ときに限る" ++
       "——イデアル類が自明でないと通らない。" ++
       "★原文の高さの計算(素点ごとの寄与の和)はこの制限を受けないので、" ++
       "外すには pullbackIdeal を**素点ごとに読む**段が要る") 6,
    .implicitStep
      ("★★hinj(座標が点を分けること)は ψ が埋め込み(§9-920)であることの点版であり、" ++
       "点の族 P を座標で添字づけるための道具である") 3 ]

end ABC3.Found.GenEll
