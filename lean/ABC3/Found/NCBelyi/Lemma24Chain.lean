/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.DescendBeta
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★[NCBelyi] `Lemma 2.4` が閉じた（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.4–p.5。

原文 (NCBelyi p.4):
> Lemma 2.4.

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`Lemma 2.4`（Reduction to the Rational Case）は、`S ⊆ ℙ¹(ℚ̄)` と `β ∈ ℚ̄∖S` に対し
非定数有理関数 `f ∈ ℚ(x)` と `S_φ ⊆ ℙ¹(ℚ)` を出して

* **(a)** `φ(S) ⊆ S_φ`
* **(b)** `φ(β) ∈ ℙ¹(ℚ)∖S_φ`
* **(c)** `φ` は `ℙ¹_ℚ∖S_φ` 上不分岐

を言う。★★★**本ファイルはこれを閉じる**（`lemma_2_4_chain`）。

## ★有理関数の表し方 —— `Chain`

原文の `f(x) ∈ ℚ(x)` は、帰納の各段で作られる
**Möbius `x ↦ c/(x−λ)` と多項式 `f₀` の交互合成**である。
★本ファイルはそれを `Chain`（そのリスト）として表す:

    `Chain.nil` … 恒等
    `Chain.cons λ c f rest` … `x ↦ rest(f(c/(x−λ)))`

* `Chain.eval` … 評価
* `Chain.IsCrit` … 臨界点（★**極 `x = λ` は除く**——極は `∞` へ写るので
  `S_φ ⊆ ℙ¹(ℚ)` の `∞` が受ける。原文が `S_φ ⊆ ℙ¹(ℚ)` と書き `𝔸¹` と書かない理由がこれである）

★★これは `Found/NCBelyi/BelyiPoly.lean` や `Thm25Step3.lean` が
「多項式 ＋ `∞ ↦ ∞`」で書いているのと**同じ流儀**である。

## ★★★組み立て

| 段 | 出どころ |
|---|---|
| `β` を運ぶ入れ子帰納法 | `§9-986` `nested_induction_descend_beta` |
| 降下段（正規化 → `f₀` → 分離 → 置き換え） | `§9-987` ＋ 本ファイルの `exists_descend_chain_step` |
| 基底段（`m(S) = 0`、すなわち `S ⊆ ℙ¹(ℚ)`） | `Chain.nil` でよい（臨界点を持たない） |

★基底段が易しいのは **`Chain.nil` が臨界点を持たない**からである
——`Thm25Step3`／`Lemma22` の「1 次式は臨界点を持たない」と同じ理由。

## ★★★★★逸脱（明示）

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| `f(x) ∈ ℚ(x)` | 有理関数そのもの | **`Chain`（Möbius と多項式の交互合成）** | `ℙ¹` 上の有理写像の一般論を作らずに済む。★中身は同じものである |
| `S_φ ⊆ ℙ¹(ℚ)` | `∞` を含みうる | **`Finset ℚ` ＋「極は除く」** | 極は `∞` へ写り、`∞ ∈ S_φ` が受ける |
| `β ∈ ℚ̄∖S` | `ℚ̄` の点 | **`β ∈ ℚ`** | 原文が最初に `x ↦ x − β` で `β ∈ ℙ¹(ℚ)` に帰着させている |
| `S` の `Gal`-安定性 | 「without loss of generality」 | **仮定として受ける** | `conjClosure`（第 410）で作れる。測度も変わらない |
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★★★★★★★有理関数の表し方 -/

/-- ★★★★★★★★**Möbius `x ↦ c/(x−λ)` と多項式の交互合成**。

原文 (NCBelyi p.4):
> Lemma 2.4.

★原文の `f(x) ∈ ℚ(x)` は、帰納の各段で作られるこの形の合成である。 -/
inductive Chain : Type
  | nil : Chain
  | cons : ℚ → ℚ → ℚ[X] → Chain → Chain

/-- 1 段ぶんの Möbius `x ↦ c/(x−λ)`。 -/
noncomputable def mob (lam c : ℚ) (x : ℂ) : ℂ := (c : ℂ) * (x - (lam : ℂ))⁻¹

/-- 鎖の評価。 -/
noncomputable def Chain.eval : Chain → ℂ → ℂ
  | .nil, x => x
  | .cons lam c f rest, x => Chain.eval rest (aeval (mob lam c x) f)

/-- 鎖の臨界点。★**極 `x = λ` は除く**——極は `∞` へ写るので `S_φ` の `∞` が受ける。 -/
def Chain.IsCrit : Chain → ℂ → Prop
  | .nil, _ => False
  | .cons lam c f rest, x =>
      x ≠ (lam : ℂ) ∧
        (aeval (mob lam c x) (derivative f) = 0 ∨ Chain.IsCrit rest (aeval (mob lam c x) f))

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★降下段（鎖のデータつき） -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**降下段（鎖を作るのに要る全部）**。

原文 (NCBelyi p.5):
> hypothesis, and composing the resulting morphisms yields

`§9-987` の `exists_descend_data_beta` に、**鎖を組み立てるための出力**
（`λ`, `c`, `f₀` と、像・`β′`・臨界点の対応）を足したものである。 -/
theorem exists_descend_chain_step (S : Finset ℂ) (hint : ∀ x ∈ S, IsIntegral ℚ x)
    (hstab : IsConjStable S) (hS : 0 < maxRedDeg S) (β : ℚ) (hβ : ((β : ℂ)) ∉ S) :
    ∃ (T : Finset ℂ) (β' lam c : ℚ) (f₀ : ℚ[X]),
      (∀ x ∈ T, IsIntegral ℚ x) ∧ IsConjStable T ∧ ((β' : ℂ)) ∉ T
      ∧ (maxRedDeg T < maxRedDeg S ∨ (maxRedDeg T = maxRedDeg S ∧ dSum T < dSum S))
      ∧ (∀ x ∈ S, x ≠ (lam : ℂ)) ∧ ((β : ℂ)) ≠ (lam : ℂ)
      ∧ (∀ x ∈ S, aeval (mob lam c x) f₀ ∈ T)
      ∧ ((β' : ℂ)) = aeval (mob lam c ((β : ℂ))) f₀
      ∧ (∀ w : ℂ, aeval (mob lam c w) (derivative f₀) = 0 → aeval (mob lam c w) f₀ ∈ T) := by
  classical
  obtain ⟨lam, c, hc0, hlamS, hlamβ, hnorm1, hnorm4⟩ := exists_normalization_complex S β hβ
  have hcne : c ≠ 0 := ne_of_gt hc0
  set σ : ℂ → ℂ := fun x => (c : ℂ) * (x - (lam : ℂ))⁻¹ + (((0:ℚ)) : ℂ) with hσ
  have hσmob : ∀ x : ℂ, σ x = mob lam c x := by
    intro x; rw [hσ, mob]; simp
  set S₁ : Finset ℂ := S.image σ with hS1
  have hlamne : ∀ x ∈ S, x ≠ (lam : ℂ) := fun x hx hc => hlamS (hc ▸ hx)
  have hS1int : ∀ y ∈ S₁, IsIntegral ℚ y := by
    intro y hy
    obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 hy
    exact isIntegral_mobius lam c 0 (hint x hxS)
  have hS1stab : IsConjStable S₁ := isConjStable_mobius lam c 0 hcne S hint hlamne hstab
  have hS1max : maxRedDeg S₁ = maxRedDeg S := maxRedDeg_mobius lam c 0 hcne S hlamS
  have hS1d : dSum S₁ = dSum S := dSum_mobius lam c 0 hcne S hlamS
  have hS1norm : ∀ y ∈ S₁, ‖y‖ ≤ 1 := by
    intro y hy
    obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 hy
    simpa [hσ] using hnorm1 x hxS
  have hS1pos : 0 < maxRedDeg S₁ := by rw [hS1max]; exact hS
  obtain ⟨E, x₀, hx₀S, hx₀top, hno, hE, hx₀drop, hEdef⟩ := exists_descend_data S₁ hS1int hS1pos
  have hx₀int : IsIntegral ℚ x₀ := hS1int x₀ hx₀S
  have hdeg : 0 < (minpoly ℚ x₀).natDegree := by
    rw [natDegree_minpoly_eq_redDeg_succ hx₀int]; omega
  have hd0 : derivative (minpoly ℚ x₀) ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := minpoly ℚ x₀)).1 hc
    omega
  have hmapne : ((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)) ≠ 0 :=
    (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
  set β₁ : ℚ := c / (β - lam) with hβ1
  have hβlam : (β : ℚ) - lam ≠ 0 := sub_ne_zero.mpr (fun h => hlamβ h.symm)
  have hβ1C : ((β₁ : ℂ)) = σ ((β : ℂ)) := by
    rw [hβ1, hσ]
    push_cast
    rw [div_eq_mul_inv, add_zero]
  have hβ1norm : 3 < ‖((β₁ : ℂ))‖ := by
    rw [hβ1C, hσ]
    simp only [Rat.cast_zero, add_zero]
    linarith [hnorm4]
  have hRint : ∀ w ∈ ((((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)).roots).toFinset),
      IsIntegral ℚ w := by
    intro w hw
    rw [Multiset.mem_toFinset] at hw
    have hwr : aeval w (derivative (minpoly ℚ x₀)) = 0 := by
      have h2 := (Polynomial.mem_roots hmapne).1 hw
      simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using h2
    exact isIntegral_of_aeval_eq_zero hd0 hwr
  have hcast : ((((minpoly ℚ x₀).eval β₁ : ℚ)) : ℂ) = aeval ((β₁ : ℂ)) (minpoly ℚ x₀) := by
    have h : ((β₁ : ℂ)) = algebraMap ℚ ℂ β₁ := by simp
    rw [h, Polynomial.aeval_algebraMap_apply]
    simp
  refine ⟨S₁.image (fun y => aeval y (minpoly ℚ x₀)) ∪ E, (minpoly ℚ x₀).eval β₁, lam, c,
    minpoly ℚ x₀, ?_, ?_, ?_, ?_, hlamne, ?_, ?_, ?_, ?_⟩
  · intro x hx
    rcases Finset.mem_union.1 hx with h | h
    · obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 h
      exact isIntegral_aeval (hS1int y hy) _
    · rw [hEdef] at h
      obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 h
      exact isIntegral_aeval (hRint w hw) _
  · rw [hEdef]
    exact isConjStable_union (isConjStable_image_aeval _ S₁ hS1int hS1stab)
      (isConjStable_image_aeval _ _ hRint
        (isConjStable_rootSet (derivative (minpoly ℚ x₀)) hd0))
  · rw [hcast, hEdef]
    exact separation_step_not_mem S₁ hS1stab hS1norm hx₀S hx₀int hdeg hβ1norm
  · have hmeas := measure_lt S₁ E (fun y => aeval y (minpoly ℚ x₀)) hno hE hx₀S hx₀top
      hx₀drop hS1pos
    rw [hS1max, hS1d] at hmeas
    exact hmeas
  · intro h
    exact hlamβ (by exact_mod_cast h.symm)
  · intro x hxS
    refine Finset.mem_union_left _ (Finset.mem_image.2 ⟨mob lam c x, ?_, rfl⟩)
    rw [← hσmob]
    exact Finset.mem_image.2 ⟨x, hxS, rfl⟩
  · rw [hcast, hβ1C, hσmob]
  · intro w hw
    refine Finset.mem_union_right _ ?_
    rw [hEdef]
    refine Finset.mem_image.2 ⟨mob lam c w, ?_, rfl⟩
    rw [Multiset.mem_toFinset]
    refine (Polynomial.mem_roots hmapne).2 ?_
    simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using hw

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Lemma 2.4` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[NCBelyi] Lemma 2.4**
(Reduction to the Rational Case)。

原文 (NCBelyi p.4):
> Lemma 2.4.

`S ⊆ ℂ`（有限・代数的・`Gal`-安定）と `β ∈ ℚ∖S` に対し、鎖 `ch`（＝ `ℚ(x)` の有理関数）と
`S_φ ⊆ ℚ` が取れて:

* **(a)** `∀ x ∈ S, ch(x) ∈ S_φ`
* **(b)** `ch(β) = b ∈ ℚ` かつ `b ∉ S_φ`
* **(c)** `∀ w`（極でない臨界点）, `ch(w) ∈ S_φ`
  ——すなわち `ch` は `ℙ¹∖S_φ` 上不分岐

★★★組み立ては `§9-986`（`β` を運ぶ帰納）と `exists_descend_chain_step`（降下段）だけである。
★基底段（`m(S) = 0`）は `Chain.nil`——**臨界点を持たない**から易しい。 -/
theorem lemma_2_4_chain (S : Finset ℂ) (β : ℚ)
    (hint : ∀ x ∈ S, IsIntegral ℚ x) (hstab : IsConjStable S) (hβ : ((β : ℂ)) ∉ S) :
    ∃ (ch : Chain) (Sphi : Finset ℚ) (b : ℚ),
      (∀ x ∈ S, ∃ q ∈ Sphi, ch.eval x = (q : ℂ))
      ∧ ch.eval ((β : ℂ)) = (b : ℂ) ∧ b ∉ Sphi
      ∧ (∀ w : ℂ, ch.IsCrit w → ∃ q ∈ Sphi, ch.eval w = (q : ℂ)) := by
  classical
  refine nested_induction_descend_beta
    (P := fun S β => ∃ (ch : Chain) (Sphi : Finset ℚ) (b : ℚ),
      (∀ x ∈ S, ∃ q ∈ Sphi, ch.eval x = (q : ℂ))
      ∧ ch.eval ((β : ℂ)) = (b : ℂ) ∧ b ∉ Sphi
      ∧ (∀ w : ℂ, ch.IsCrit w → ∃ q ∈ Sphi, ch.eval w = (q : ℂ)))
    (Q := fun S β => (∀ x ∈ S, IsIntegral ℚ x) ∧ IsConjStable S ∧ ((β : ℂ)) ∉ S)
    ?_ ?_ S β ⟨hint, hstab, hβ⟩
  · -- ★基底段: `m(S) = 0` すなわち `S ⊆ ℙ¹(ℚ)`
    rintro T b hT ⟨hTint, -, hbT⟩
    have hrat : ∀ x ∈ T, ∃ q : ℚ, x = (q : ℂ) := by
      intro x hx
      exact (redDeg_eq_zero_iff (hTint x hx)).1 (maxRedDeg_eq_zero_iff.1 hT x hx)
    choose qf hqf using hrat
    refine ⟨Chain.nil, T.attach.image (fun p => qf p.1 p.2), b, ?_, rfl, ?_, ?_⟩
    · intro x hx
      exact ⟨qf x hx, Finset.mem_image.2 ⟨⟨x, hx⟩, Finset.mem_attach _ _, rfl⟩,
        by rw [Chain.eval]; exact hqf x hx⟩
    · intro hmem
      obtain ⟨p, -, hp⟩ := Finset.mem_image.1 hmem
      exact hbT (by rw [← hp] at *; rw [← hqf p.1 p.2]; exact p.2)
    · intro w hw
      exact absurd hw (by rw [Chain.IsCrit]; exact id)
  · -- ★★降下段
    rintro U b hU ⟨hUint, hUstab, hbU⟩
    obtain ⟨T, b', lam, c, f₀, hTint, hTstab, hb'T, hmeas, hlamne, hblam, hmemT, hb'eq,
      hcritT⟩ := exists_descend_chain_step U hUint hUstab hU b hbU
    refine ⟨T, b', hmeas, ⟨hTint, hTstab, hb'T⟩, ?_⟩
    rintro ⟨ch', Sphi, bb, hval', hbeta', hbb, hcrit'⟩
    refine ⟨Chain.cons lam c f₀ ch', Sphi, bb, ?_, ?_, hbb, ?_⟩
    · intro x hx
      rw [Chain.eval]
      exact hval' _ (hmemT x hx)
    · rw [Chain.eval, ← hb'eq]
      exact hbeta'
    · intro w hw
      rw [Chain.IsCrit] at hw
      obtain ⟨-, hcase⟩ := hw
      rw [Chain.eval]
      rcases hcase with h | h
      · exact hval' _ (hcritT w h)
      · exact hcrit' _ h

/-! ## ★出典の紐付け(`.src`) -/

def Chain.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.4(有理関数の表し方——Möbius と多項式の交互合成)",
    sectionId := "ncbelyi-lemma-2-4" }

def exists_descend_chain_step.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(降下段——鎖のデータつき)",
    sectionId := "ncbelyi-lemma-2-4" }

def lemma_2_4_chain.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4, item := "Lemma 2.4",
    sectionId := "ncbelyi-lemma-2-4" }

def lemma_2_4_chain.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "nested_induction_descend_beta(β を運ぶ帰納、§9-986)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.nested_induction_descend_beta") 3,
    .citation "[ABC3]" "exists_descend_data_beta(降下段、§9-987)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.exists_descend_data_beta") 5,
    .citation "[ABC3]" "separation_step_not_mem(f₀(β) ∉ S′、§9-984)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.separation_step_not_mem") 4,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 526 ブロック): " ++
       "**[NCBelyi] Lemma 2.4 が閉じた**。" ++
       "★組み立ては §9-986(β を運ぶ帰納)と exists_descend_chain_step(降下段)だけである。" ++
       "★★基底段(m(S) = 0)は Chain.nil——**臨界点を持たない**から易しい" ++
       "(Lemma 2.2 の『1 次式は臨界点を持たない』と同じ理由)") 8,
    .implicitStep
      ("★★★★★逸脱(明示): 原文の f(x) ∈ ℚ(x) を **Chain**" ++
       "(Möbius x ↦ c/(x−λ) と多項式の交互合成)として表している" ++
       "——ℙ¹ 上の有理写像の一般論を作らずに済ませるためであり、中身は同じものである。" ++
       "★臨界点からは**極 x = λ を除いてある**——極は ∞ へ写り、" ++
       "原文の S_φ ⊆ ℙ¹(ℚ) の ∞ が受ける" ++
       "(原文が 𝔸¹ ではなく ℙ¹ と書く理由がこれである)。" ++
       "★★また β ∈ ℚ とし、S の Gal-安定性は仮定として受けている" ++
       "(前者は原文が x ↦ x − β で帰着させる段、後者は conjClosure(第 410)で作れる)") 7,
    .implicitStep
      ("★★★これで [NCBelyi] Theorem 2.5 に残るのは**段 1 だけ**である" ++
       "——Riemann–Roch / Serre 双対性で ψ : X → ℙ¹ を作り T = {β} に帰着させる段" ++
       "([Stacks] 53.4 / 53.5 / 48.27、第 419 で手元にあることを実測)。" ++
       "★段 2(Lemma 2.4)は本ファイル、段 3 は §9-981、段 4 は第 398-404 である") 8 ]

end ABC3.Found.NCBelyi
