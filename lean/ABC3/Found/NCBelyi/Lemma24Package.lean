/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.DescendData
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★`Lemma 2.4` の (a)(c) と、残る (b) の形（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.4–p.5。

原文 (NCBelyi p.4):
> Lemma 2.4.

## ★★★★★★★★★★★★★★★★★これは何か

原文の `Lemma 2.4` は、`S ⊆ ℙ¹(ℚ̄)` と `β ∈ ℚ̄∖S` に対し、非定数有理関数
`f ∈ ℚ(x)` と `S_φ ⊆ ℙ¹(ℚ)` を出して

* **(a)** `φ(S) ⊆ S_φ`
* **(b)** `φ(β) ∈ ℙ¹(ℚ)∖S_φ`
* **(c)** `φ` は `ℙ¹_ℚ∖S_φ` 上不分岐

を言う。★★★**本ファイルは (a)(c) を取り、(b) が要求するものを型で 1 つに絞る。**

## ★機構 —— (a)(c) は `exists_poly_image_rat_crit` の**言い換え**である

`§第 417-418` の `exists_poly_image_rat_crit` は

    `∃ h ∈ ℚ[X], 0 < deg h ∧ h(S) ⊆ ℚ ∧ (h の臨界値) ⊆ ℚ`

を与える。★`S_φ ≔ h(S) ∪ (臨界値)`（どちらも `ℚ` の有限集合）と置けば
**(a) と (c) はそのまま出る**——(c) は「臨界値が `S_φ` に入っている」ことそのものである。

★★これが `lemma_2_4_ac` である。

## ★★★★★★残るのは (b) だけ ——**分離**

(b) `h(β) ∉ S_φ` は、`h` が `β` を `S` からも臨界点からも**分離する**ことを要求する:

    `∀ x ∈ S,  h(x) ≠ h(β)`   かつ   `∀ w, h′(w) = 0 → h(w) ≠ h(β)`

★★`lemma_2_4_of_separation` はこれを仮定として受け、**(a)(b)(c) をすべて出す**。
★★★したがって `Lemma 2.4` に残っているのは**この分離だけ**である。

## ★原文はこれをどう作るか（次に着手する段の設計図）

原文 p.5:

> multiplying by some positive rational number, we may assume that |α| ≤1, for all

★すなわち `Lemma 2.3` の自己同型と正の有理数倍で
**`|α| ≤ 1`（`∀ α ∈ S`）かつ `|β| ≥ C`（`C` は `d₀` に応じて十分大きい）**に正規化する。
★★すると最小多項式 `f₀`（次数 `d₀`）の係数は `≤ d₀^{d₀}` で押さえられ、
`|f₀(α)|` と `|f₀(S₀)|` は `d₀` だけで決まる量で抑えられる一方、
`|f₀(β)|` は `C^{d₀}` 程度になる——★★★**`C` を十分大きく取れば `f₀(β) ∉ S′`**。

★これに要る部品は本プロジェクトに**ある**（第 405-410、第 416-418）:

| 部品 | 場所 |
|---|---|
| 正規化（`|α| ≤ 1`、`|β| ≥ C`） | `Found/NCBelyi/RatSeparation.lean` の `exists_rat_normalization` |
| 係数の限界 `≤ d₀^{d₀}` | `Found/NCBelyi/CoeffBound.lean` の `norm_coeff_le_choose` |
| 入れ子帰納法 | `Found/NCBelyi/NestedInduction.lean`／`DescendData.lean` |

★★**残る仕事は、`nested_induction_descend'` を `β` を運ぶ形に書き直すこと**である
——`exists_descend_data` が返す `g` に「`g β ∉ S′`」を足す。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★★★★★★★★★★★(a) と (c) -/

/-- ★★★★★★★★★★★★**`Lemma 2.4` の (a) と (c)**。

原文 (NCBelyi p.4):
> Lemma 2.4.

`S_φ ⊆ ℙ¹(ℚ)` が取れて `φ(S) ⊆ S_φ`（(a)）かつ `φ` は `ℙ¹∖S_φ` 上不分岐（(c)）。

★`exists_poly_image_rat_crit`（第 417-418）を `S_φ ≔ h(S) ∪ (臨界値)` の形に言い換えただけ。
★★(c) は「**臨界値が `S_φ` に入っている**」ことそのものである。 -/
theorem lemma_2_4_ac (S : Finset ℂ) (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    ∃ (h : ℚ[X]) (Sphi : Finset ℚ), 0 < h.natDegree
      ∧ (∀ x ∈ S, ∃ q ∈ Sphi, aeval x h = (q : ℂ))
      ∧ (∀ w : ℂ, aeval w (derivative h) = 0 → ∃ q ∈ Sphi, aeval w h = (q : ℂ)) := by
  classical
  obtain ⟨h, hdeg, hval, hcrit⟩ := exists_poly_image_rat_crit S hint
  have hd0 : derivative h ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := h)).1 hc
    omega
  have hmapne : ((derivative h).map (algebraMap ℚ ℂ)) ≠ 0 :=
    (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
  set R : Finset ℂ := (((derivative h).map (algebraMap ℚ ℂ)).roots).toFinset with hR
  have hRiff : ∀ w : ℂ, w ∈ R ↔ aeval w (derivative h) = 0 := by
    intro w
    rw [hR, Multiset.mem_toFinset, Polynomial.mem_roots hmapne]
    simp [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]
  have hRmem : ∀ w ∈ R, aeval w (derivative h) = 0 := fun w hw => (hRiff w).1 hw
  choose qS hqS using hval
  choose qC hqC using hcrit
  refine ⟨h, (S.attach.image (fun p => qS p.1 p.2))
      ∪ (R.attach.image (fun p => qC p.1 (hRmem p.1 p.2))), hdeg, ?_, ?_⟩
  · intro x hx
    refine ⟨qS x hx, Finset.mem_union_left _ ?_, hqS x hx⟩
    exact Finset.mem_image.2 ⟨⟨x, hx⟩, Finset.mem_attach _ _, rfl⟩
  · intro w hw
    have hwR : w ∈ R := (hRiff w).2 hw
    refine ⟨qC w hw, Finset.mem_union_right _ ?_, hqC w hw⟩
    exact Finset.mem_image.2 ⟨⟨w, hwR⟩, Finset.mem_attach _ _, rfl⟩

/-! ## ★★★★★★★★★★★★★★★★★(b) を仮定すれば (a)(b)(c) がすべて出る -/

/-- ★★★★★★★★★★★★★★★★★**分離を仮定すれば `Lemma 2.4` の (a)(b)(c) が出る**。

原文 (NCBelyi p.4):
> Lemma 2.4.

★**仮定は 2 行だけ**である:

    `∀ x ∈ S,  h(x) ≠ h(β)`   かつ   `∀ w, h′(w) = 0 → h(w) ≠ h(β)`

★★これが原文の「`for a suitable choice of C, it follows that f₀(β) ∉ S′`」にあたる。
★★★**`Lemma 2.4` に残っているのはこの分離だけ**である
——(a)(c) は `lemma_2_4_ac` で、`h` の存在は `exists_poly_image_rat_crit` で取れている。 -/
theorem lemma_2_4_of_separation (S : Finset ℂ) (β : ℚ)
    (hsep : ∃ h : ℚ[X], 0 < h.natDegree
       ∧ (∀ x ∈ S, ∃ q : ℚ, aeval x h = (q : ℂ))
       ∧ (∀ w : ℂ, aeval w (derivative h) = 0 → ∃ q : ℚ, aeval w h = (q : ℂ))
       ∧ (∀ x ∈ S, aeval x h ≠ ((h.eval β : ℚ) : ℂ))
       ∧ (∀ w : ℂ, aeval w (derivative h) = 0 → aeval w h ≠ ((h.eval β : ℚ) : ℂ))) :
    ∃ (h : ℚ[X]) (Sphi : Finset ℚ), 0 < h.natDegree
      ∧ (∀ x ∈ S, ∃ q ∈ Sphi, aeval x h = (q : ℂ))
      ∧ h.eval β ∉ Sphi
      ∧ (∀ w : ℂ, aeval w (derivative h) = 0 → ∃ q ∈ Sphi, aeval w h = (q : ℂ)) := by
  classical
  obtain ⟨h, hdeg, hval, hcrit, hsepS, hsepC⟩ := hsep
  have hd0 : derivative h ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := h)).1 hc
    omega
  have hmapne : ((derivative h).map (algebraMap ℚ ℂ)) ≠ 0 :=
    (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
  set R : Finset ℂ := (((derivative h).map (algebraMap ℚ ℂ)).roots).toFinset with hR
  have hRiff : ∀ w : ℂ, w ∈ R ↔ aeval w (derivative h) = 0 := by
    intro w
    rw [hR, Multiset.mem_toFinset, Polynomial.mem_roots hmapne]
    simp [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]
  have hRmem : ∀ w ∈ R, aeval w (derivative h) = 0 := fun w hw => (hRiff w).1 hw
  choose qS hqS using hval
  choose qC hqC using hcrit
  refine ⟨h, (S.attach.image (fun p => qS p.1 p.2))
      ∪ (R.attach.image (fun p => qC p.1 (hRmem p.1 p.2))), hdeg, ?_, ?_, ?_⟩
  · intro x hx
    refine ⟨qS x hx, Finset.mem_union_left _ ?_, hqS x hx⟩
    exact Finset.mem_image.2 ⟨⟨x, hx⟩, Finset.mem_attach _ _, rfl⟩
  · intro hmem
    rcases Finset.mem_union.1 hmem with hm | hm
    · obtain ⟨p, -, hp⟩ := Finset.mem_image.1 hm
      exact hsepS p.1 p.2 (by rw [hqS p.1 p.2, hp])
    · obtain ⟨p, -, hp⟩ := Finset.mem_image.1 hm
      exact hsepC p.1 (hRmem p.1 p.2) (by rw [hqC p.1 (hRmem p.1 p.2), hp])
  · intro w hw
    have hwR : w ∈ R := (hRiff w).2 hw
    refine ⟨qC w hw, Finset.mem_union_right _ ?_, hqC w hw⟩
    exact Finset.mem_image.2 ⟨⟨w, hwR⟩, Finset.mem_attach _ _, rfl⟩

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def lemma_2_4_ac.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.4((a) φ(S) ⊆ S_φ と (c) ℙ¹∖S_φ 上不分岐)",
    sectionId := "ncbelyi-lemma-2-4" }

def lemma_2_4_of_separation.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.4(分離を仮定すれば (a)(b)(c) が出る)",
    sectionId := "ncbelyi-lemma-2-4" }

def lemma_2_4_of_separation.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_poly_image_rat_crit(多項式の段、第 417-418)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.exists_poly_image_rat_crit") 4,
    .implicitStep
      ("★★★★★★測定(2026-08-29): Lemma 2.4 の (a)(c) は " ++
       "exists_poly_image_rat_crit の**言い換え**である" ++
       "(S_φ ≔ h(S) ∪ 臨界値 と置くだけ。(c) は『臨界値が S_φ に入っている』ことそのもの)。" ++
       "★★残るのは (b) すなわち**分離**だけである: " ++
       "∀ x ∈ S, h(x) ≠ h(β) かつ ∀ w, h′(w) = 0 → h(w) ≠ h(β)") 6,
    .implicitStep
      ("★★★次に着手する段の設計図(原文 p.5): Lemma 2.3 の自己同型と正の有理数倍で " ++
       "|α| ≤ 1(∀ α ∈ S)かつ |β| ≥ C に正規化すると、" ++
       "最小多項式 f₀(次数 d₀)の係数は ≤ d₀^{d₀} で押さえられ、" ++
       "|f₀(α)| と |f₀(S₀)| は d₀ だけで決まる量で抑えられる一方 |f₀(β)| は C^{d₀} 程度になる。" ++
       "★C を十分大きく取れば f₀(β) ∉ S′。" ++
       "★★部品は本プロジェクトにある(exists_rat_normalization、norm_coeff_le_choose、" ++
       "nested_induction_descend′)——残る仕事は " ++
       "**nested_induction_descend′ を β を運ぶ形に書き直すこと**である") 7 ]

end ABC3.Found.NCBelyi
