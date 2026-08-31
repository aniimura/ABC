/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottChartUnion
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★[GenEll] Proposition 1.4, (iv) —— チャートを跨いだ組み立て（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★これは何か —— Serre の道の組み立て

★`§9-888`（段 E3h）: チャート 1 枚（`X_{s_i}`）の点について Northcott
★★`§9-889`（段 E3i）: チャートごとに有限なら全体でも有限

★★★本ファイルはその 2 つを**実際に組む**。点の集合を

    `Σ i, P i`   （`P i` ＝ チャート `i` を通る点の集合）

の形で受けると、依存型の輸送なしにきれいに組める。

    `{ q : Σ i, P i | ht(q) ≤ C }` は有限

## ★★★測定の記録 —— なぜ `Σ` 型か

★★★★点ごとにチャートが違うと、`xF p : Spec 𝓞 ⟶ X_{s_{chart p}}` の**型が点に依存する**。
`chart p = i` に沿った輸送（`▸`）を仮定の中で書くことになり、読めなくなる。
★`Σ i, P i` で受けると、各 `i` の成分は**型が固定**され、
`§9-888` をそのまま当てられる（2026-08-28 実測）。

## ★残っている段（明示）

★★仮定は `§9-888` と同じ 3 つ（`hdeg`・整数表示・点の分離）で、
いずれも出どころが判っている（`§9-854`・`§9-849`・`§9-851`）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★`Σ` 型の各成分ごとに有限なら全体でも有限 -/

/-- ★★★**`Σ` 型の各成分ごとに有限なら全体でも有限**。 -/
theorem finite_sigma_of_finite_fibers {ι : Type} [Finite ι] {P : ι → Type}
    (ht : (Σ i, P i) → ℝ) (C : ℝ)
    (hfin : ∀ i : ι, {p : P i | ht ⟨i, p⟩ ≤ C}.Finite) :
    {q : (Σ i, P i) | ht q ≤ C}.Finite := by
  refine finite_of_finite_charts ht C Sigma.fst (fun i => ?_)
  refine Set.Finite.subset ((hfin i).image (Sigma.mk i)) ?_
  rintro ⟨j, p⟩ ⟨hj, hp⟩
  cases hj
  exact ⟨p, hp, rfl⟩

/-! ## ★★★★★★★★★★★★★★★★★チャートを跨いだ Northcott -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★**[GenEll] Proposition 1.4, (iv)**——チャートを跨いで。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★点の集合を `Σ i, P i`（`P i` ＝ チャート `i` を通る点）の形で受けると、
依存型の輸送なしに `§9-888` を各成分へ当てられる。 -/
theorem northcott_globalChart_sigma (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (d : ℕ)
    {P : Fin (N+1) → Type}
    (fld : ∀ i, P i → IntermediateField ℚ ℂ) (hnf : ∀ i p, NumberField (fld i p))
    (hdeg : ∀ i p, Module.finrank ℚ (fld i p) ≤ d)
    (xF : ∀ i (p : P i), haveI := hnf i p;
      specRingOfIntegers (fld i p) ⟶ (nonVanishing M (s i)).toScheme)
    (x : ∀ i (p : P i), Fin (N+1) → NumberField.RingOfIntegers (fld i p))
    (hx : ∀ i p, x i p ≠ 0)
    (hw : ∀ i p, haveI := hnf i p; ∀ k, x i p k
      = (Spec.preimage (xF i p ≫ globalChartMorphism M hM φ s i)).hom
          (projCoord N ℤ i k) * x i p i)
    (h0 : ∀ i p, x i p 0 ≠ 0)
    (idx : Fin (N+1))
    (hinj : ∀ i, Function.Injective (fun (p : P i) (k : Fin (N+1)) =>
      ((((x i p k : fld i p)) / ((x i p idx : fld i p)) : fld i p) : ℂ)))
    (C : ℝ) :
    {q : (Σ i, P i) | haveI := hnf q.1 q.2;
      htArith (fld q.1 q.2) ((hyperplaneArith N).comap (globalChartToProj M hM φ s q.1))
        (xF q.1 q.2) ≤ C}.Finite :=
  finite_sigma_of_finite_fibers _ C (fun i =>
    northcott_globalChart M hM φ s i d (fld i) (hnf i) (hdeg i) (xF i) (x i)
      (hx i) (hw i) (h0 i) idx (hinj i) C)

/-! ## ★出典の紐付け(`.src`) -/

def finite_sigma_of_finite_fibers.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(Σ 型の各成分ごとに有限なら全体でも有限)",
    sectionId := "genell-prop-1-4" }

def northcott_globalChart_sigma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートを跨いだ Northcott)",
    sectionId := "genell-prop-1-4" }

def northcott_globalChart_sigma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_globalChart(チャート 1 枚分、段 E3h、§9-888)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_globalChart") 2,
    .citation "[ABC3]" "finite_of_finite_charts(有限個の和、段 E3i、§9-889)"
      (.inProject "ABC3" "ABC3.Found.GenEll.finite_of_finite_charts") 1,
    .implicitStep
      ("★★★★測定: 点ごとにチャートが違うと xF p : Spec 𝓞 ⟶ X_{s_{chart p}} の" ++
       "**型が点に依存する**ので、chart p = i に沿った輸送(▸)を仮定の中で書くことになり" ++
       "読めなくなる。★Σ i, P i で受けると各 i の成分は**型が固定**され、" ++
       "§9-888 をそのまま当てられる(2026-08-28 実測)") 3,
    .implicitStep
      ("★★仮定は §9-888 と同じ 3 つ(hdeg・整数表示・点の分離)で、" ++
       "いずれも出どころが判っている(§9-854・§9-849・§9-851)") 3 ]

end ABC3.Found.GenEll
