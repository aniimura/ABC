/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.Lemma22
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★`Theorem 2.5` の段 3（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.7。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★★★★★★★★★★★★★★★★★★★★★これは何か

`Theorem 2.5` の証明（原文 p.6–p.7）は 4 段である:

| 段 | 中身 | 状態 |
|---|---|---|
| 1 | Riemann–Roch / Serre 双対性で `ψ : X → ℙ¹` を作り `T = {β}` に帰着 | ☆残る |
| 2 | `Lemma 2.4` で `S ⊆ ℙ¹(ℚ)` に帰着 | △ 多項式の段は済（第 417–418） |
| ★3 | `Lemma 2.3`（`C = 4`）と `x ↦ ν·x + μ` で `Lemma 2.2` の仮定を作る | ★★**本ファイル** |
| 4 | `Lemma 2.2` で閉じる | ✅ 済（第 398–404） |

原文 p.7:

> Finally, by applying an automorphism as in Lemma 2.3 [for, say, C = 4], followed by
> a suitable automorphism of the form x ↦ν · x + μ, where ν ∈{±1} and μ ∈Q,
> gives rise to a situation in which the hypotheses of Lemma 2.2 are valid.

★★★**本ファイルはこの 1 文を型にして証明する。**

## ★機構 —— なぜ `C = 4` でちょうど足りるか

`λ` を `Lemma 2.3`（`C = 4`）で取り、`b ≔ 1/(β−λ)`、`a_α ≔ 1/(α−λ)` と置くと
`4·|a_α| ≤ |b|`。`ν ≔ sign(b)`、`B ≔ |b| = ν·b > 0` とすれば

* `T ≔ {0} ∪ {ν·a_α}` は `|x| ≤ B/4` を満たす（`0` は `∞` の像）
* `m ≔ min T ∈ [−B/4, 0]`、`μ ≔ −m ∈ [0, B/4]`
* `S₀ ≔ T + μ` は `min = 0`、`0 ≤ x ≤ B/4 + B/4 = B/2`
* `β₀ ≔ B + μ ≥ B`

★したがって `2x ≤ B ≤ β₀`——**`Lemma 2.2` の `2α ≤ β` がちょうど出る**。
★★`μ ≥ 0` が効くので、原文の `C = 4` で足りる。

## ★★`x ↦ 1/(x−λ)` は `ℙ¹` の自己同型だから分岐を増やさない

★合成 `f ∘ g` の臨界点は `g⁻¹`（`f` の臨界点）だけである。
★★本ファイルは `g` の単射性（`mobius_inj`）を取り、その事実を型で確かめる。

## ★★★★★逸脱（明示）

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| `ℙ¹(ℚ)` の `∞` | 点として扱う | **`∞` の像 `μ` を別に出す** | `Lemma 2.2` が多項式（＝`∞ ↦ ∞`）で書かれているため |
| 段 1・段 2 | 証明する | **含めない** | 本ファイルは段 3 だけである |
| `Lemma 2.3` | `absInvShift` の形（`Found/NCBelyi/RatSeparation.lean`） | **`ℚ` の中で書き直した**（`exists_sep_rat`） | `S ⊆ ℙ¹(ℚ)` に落ちた後は `ℂ` を経由する必要が無い |
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★★★★★★★★`Lemma 2.3` の `ℚ` 版 -/

/-- ★★★★★★★★★**`λ = β + ε` を小さく取れば分離できる**（`Lemma 2.3` の `ℚ` 版）。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

    `C·|1/(α−λ)| ≤ |1/(β−λ)|`   （`∀ α ∈ S`）

★`δ ≔ min_{α∈S}|α−β| > 0` として `ε ≔ min(δ/2, δ/(2C))` と取ればよい。
★★`|β−λ| = ε` は小さく、`|α−λ| ≥ δ/2` は下から押さえられる。 -/
theorem exists_sep_rat (S : Finset ℚ) (β : ℚ) (hβ : β ∉ S) (C : ℚ) (hC : 0 < C) :
    ∃ lam : ℚ, lam ∉ S ∧ lam ≠ β ∧
      ∀ α ∈ S, C * |1 / (α - lam)| ≤ |1 / (β - lam)| := by
  classical
  by_cases hS : S.Nonempty
  · set D : Finset ℚ := S.image (fun α => |α - β|) with hD
    have hDne : D.Nonempty := hS.image _
    set delta : ℚ := D.min' hDne with hdelta
    have hdpos : 0 < delta := by
      rw [hdelta]
      refine (Finset.lt_min'_iff D hDne).mpr ?_
      intro y hy
      rw [hD, Finset.mem_image] at hy
      obtain ⟨α, hα, rfl⟩ := hy
      have h : α - β ≠ 0 := sub_ne_zero.mpr (fun h => hβ (h ▸ hα))
      exact abs_pos.mpr h
    have hdle : ∀ α ∈ S, delta ≤ |α - β| := by
      intro α hα
      rw [hdelta]
      exact D.min'_le _ (by rw [hD]; exact Finset.mem_image_of_mem _ hα)
    set eps : ℚ := min (delta / 2) (delta / (2 * C)) with heps
    have hepos : 0 < eps := lt_min (by positivity) (by positivity)
    have heps1 : eps ≤ delta / 2 := min_le_left _ _
    have heps2 : eps ≤ delta / (2 * C) := min_le_right _ _
    have heps2' : 2 * C * eps ≤ delta := by
      rw [le_div_iff₀ (by positivity)] at heps2
      linarith
    refine ⟨β + eps, ?_, ?_, ?_⟩
    · intro hmem
      have h := hdle _ hmem
      rw [show β + eps - β = eps by ring, abs_of_pos hepos] at h
      linarith
    · intro h; simp at h; linarith
    · intro α hα
      have hb : β - (β + eps) = -eps := by ring
      have habs : |1 / (β - (β + eps))| = 1 / eps := by
        rw [hb, abs_div, abs_neg, abs_one, abs_of_pos hepos]
      have hlow : delta / 2 ≤ |α - (β + eps)| := by
        have h1 : |α - β| - eps ≤ |α - (β + eps)| := by
          have h2 := abs_sub_abs_le_abs_sub (α - β) eps
          rw [show α - β - eps = α - (β + eps) by ring, abs_of_pos hepos] at h2
          linarith
        have h3 := hdle α hα
        linarith
      have hpos2 : (0:ℚ) < |α - (β + eps)| := lt_of_lt_of_le (by positivity) hlow
      rw [habs, abs_div, abs_one]
      have h1 : C * (1 / |α - (β + eps)|) ≤ C * (1 / (delta / 2)) := by
        refine mul_le_mul_of_nonneg_left ?_ hC.le
        exact one_div_le_one_div_of_le (by positivity) hlow
      have h2 : C * (1 / (delta / 2)) ≤ 1 / eps := by
        have h : C * (1 / (delta / 2)) = 2 * C / delta := by field_simp
        rw [h, div_le_div_iff₀ hdpos hepos]
        nlinarith
      linarith
  · rw [Finset.not_nonempty_iff_eq_empty] at hS
    exact ⟨β + 1, by simp [hS], by intro h; simp at h, by simp [hS]⟩

/-! ## ★★★`x ↦ ν/(x−λ) + μ` は単射 -/

/-- ★★★**Möbius 変換 `x ↦ ν/(x−λ) + μ` は（`λ` を除いて）単射である**。

★`ℙ¹` の自己同型だから**分岐を増やさない**——合成 `f ∘ g` の臨界点は
`g⁻¹`（`f` の臨界点）だけである。 -/
theorem mobius_inj (lam nu mu : ℚ) (hnu : nu ≠ 0) :
    ∀ x y : ℚ, x ≠ lam → y ≠ lam →
      nu * (1 / (x - lam)) + mu = nu * (1 / (y - lam)) + mu → x = y := by
  intro x y hx hy h
  have hx' : x - lam ≠ 0 := sub_ne_zero.mpr hx
  have hy' : y - lam ≠ 0 := sub_ne_zero.mpr hy
  have h2 : nu * (1 / (x - lam)) = nu * (1 / (y - lam)) := by linarith
  have h3 : (1 / (x - lam)) = (1 / (y - lam)) := mul_left_cancel₀ hnu h2
  field_simp at h3
  linarith

/-! ## ★★★★★★★★★★★★★★★★★★段 3 —— `Lemma 2.2` の仮定を作る -/

/-- ★★★★★★★★★★★★★★★★★★**原文 p.7 の段 3**
——`Lemma 2.3`（`C = 4`）と `x ↦ ν·x + μ` で `Lemma 2.2` の仮定を作る。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

★`ν ∈ {±1}`、`μ ∈ ℚ` は原文の指定どおりである。
★★`S₀` は `S ∪ {∞}` の像であり、`μ` が `∞` の像である。
★★★`C = 4` でちょうど足りる理由は本ファイルの冒頭を参照。 -/
theorem exists_lemma22_setup (S : Finset ℚ) (β : ℚ) (hβ : β ∉ S) :
    ∃ (lam nu mu : ℚ) (S₀ : Finset ℚ),
      (nu = 1 ∨ nu = -1) ∧ lam ∉ S ∧ lam ≠ β
      ∧ (∀ α ∈ S, nu * (1 / (α - lam)) + mu ∈ S₀)
      ∧ mu ∈ S₀
      ∧ (0 : ℚ) ∈ S₀
      ∧ (∀ x ∈ S₀, x ≠ 0 → 0 < x)
      ∧ nu * (1 / (β - lam)) + mu ∉ S₀
      ∧ nu * (1 / (β - lam)) + mu ≠ 0
      ∧ (∀ x ∈ S₀, x ≠ 0 → 2 * x ≤ nu * (1 / (β - lam)) + mu) := by
  classical
  obtain ⟨lam, hlamS, hlamβ, hsep⟩ := exists_sep_rat S β hβ 4 (by norm_num)
  set b : ℚ := 1 / (β - lam) with hb
  have hβlam : β - lam ≠ 0 := sub_ne_zero.mpr (fun h => hlamβ h.symm)
  have hbne : b ≠ 0 := by rw [hb]; exact one_div_ne_zero hβlam
  set nu : ℚ := if 0 < b then 1 else -1 with hnudef
  have hnu : nu = 1 ∨ nu = -1 := by rw [hnudef]; split_ifs <;> simp
  have hnuabs : |nu| = 1 := by rcases hnu with h | h <;> rw [h] <;> norm_num
  have hnub : nu * b = |b| := by
    rw [hnudef]
    split_ifs with h
    · rw [one_mul, abs_of_pos h]
    · rw [abs_of_neg (lt_of_le_of_ne (not_lt.mp h) hbne)]; ring
  set B : ℚ := |b| with hB
  have hBpos : 0 < B := abs_pos.mpr hbne
  set T : Finset ℚ := insert (0:ℚ) (S.image (fun α => nu * (1 / (α - lam)))) with hT
  have hTne : T.Nonempty := ⟨0, Finset.mem_insert_self _ _⟩
  have hT0 : (0:ℚ) ∈ T := Finset.mem_insert_self _ _
  have hTbd : ∀ x ∈ T, |x| ≤ B / 4 := by
    intro x hx
    rw [hT, Finset.mem_insert] at hx
    rcases hx with rfl | hx
    · rw [abs_zero]; positivity
    · rw [Finset.mem_image] at hx
      obtain ⟨α, hα, rfl⟩ := hx
      rw [abs_mul, hnuabs, one_mul]
      have := hsep α hα
      linarith
  set m : ℚ := T.min' hTne with hm
  have hmT : m ∈ T := T.min'_mem hTne
  have hm0 : m ≤ 0 := T.min'_le 0 hT0
  have hmlb : -(B/4) ≤ m := by
    have h := hTbd m hmT
    rw [abs_le] at h
    exact h.1
  set mu : ℚ := -m with hmu
  have hmu0 : 0 ≤ mu := by rw [hmu]; linarith
  have hmuub : mu ≤ B / 4 := by rw [hmu]; linarith
  set S₀ : Finset ℚ := T.image (fun x => x + mu) with hS0
  have hmemS0 : ∀ y ∈ S₀, ∃ t ∈ T, y = t + mu := by
    intro y hy
    rw [hS0, Finset.mem_image] at hy
    obtain ⟨t, ht, rfl⟩ := hy
    exact ⟨t, ht, rfl⟩
  have hnn : ∀ y ∈ S₀, 0 ≤ y ∧ y ≤ B / 2 := by
    intro y hy
    obtain ⟨t, ht, rfl⟩ := hmemS0 y hy
    have h1 : m ≤ t := T.min'_le t ht
    have h2 := hTbd t ht
    rw [abs_le] at h2
    refine ⟨by rw [hmu]; linarith, by linarith⟩
  have hbeta : nu * (1 / (β - lam)) + mu = B + mu := by rw [← hb, hnub]
  refine ⟨lam, nu, mu, S₀, hnu, hlamS, hlamβ, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro α hα
    rw [hS0]
    exact Finset.mem_image_of_mem _ (Finset.mem_insert_of_mem
      (Finset.mem_image_of_mem _ hα))
  · rw [hS0]
    have h := Finset.mem_image_of_mem (fun x => x + mu) hT0
    simpa using h
  · rw [hS0]
    have h := Finset.mem_image_of_mem (fun x => x + mu) hmT
    simpa [hmu] using h
  · intro x hx hx0
    exact lt_of_le_of_ne (hnn x hx).1 (Ne.symm hx0)
  · rw [hbeta]
    intro hmem
    have := (hnn _ hmem).2
    linarith
  · rw [hbeta]; intro h; linarith
  · intro x hx _
    rw [hbeta]
    have := (hnn x hx).2
    linarith

/-! ## ★★★★★★★★★★★★★★★★★★★★★★段 3 ＋ 段 4 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**`Theorem 2.5` の `X = ℙ¹_ℚ`、`S ⊆ ℙ¹(ℚ)`、
`T = {β}` の場合**（原文 p.7 の最後の 2 段）。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

`g(x) ≔ ν/(x−λ) + μ`（`ν ∈ {±1}`、`λ, μ ∈ ℚ`）と多項式 `f ∈ ℚ[x]` が取れて:

* (b) `f(g(α)) ∈ {0,1}`（`∀ α ∈ S`）、`f(g(∞)) = f(μ) ∈ {0,1}`
* (c) `f(g(β)) ∉ {0,1}`
* (a) `f` は `ℙ¹_ℚ∖{0,1,∞}` 上不分岐

★`g` は `ℙ¹` の自己同型（`mobius_inj`）だから、合成 `f ∘ g` の分岐は `f` の分岐だけである。
★★**残るのは段 1（Riemann–Roch）と段 2（`Lemma 2.4`）だけ**である。 -/
theorem theorem_2_5_p1_rat (S : Finset ℚ) (β : ℚ) (hβ : β ∉ S) :
    ∃ (lam nu mu : ℚ) (f : Polynomial ℚ),
      (nu = 1 ∨ nu = -1) ∧ lam ∉ S ∧ lam ≠ β ∧ 0 < f.natDegree
      ∧ (∀ α ∈ S, f.eval (nu * (1 / (α - lam)) + mu) = 0
          ∨ f.eval (nu * (1 / (α - lam)) + mu) = 1)
      ∧ (f.eval mu = 0 ∨ f.eval mu = 1)
      ∧ f.eval (nu * (1 / (β - lam)) + mu) ≠ 0
      ∧ f.eval (nu * (1 / (β - lam)) + mu) ≠ 1
      ∧ (∀ x : ℂ, (derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
          (f.map (algebraMap ℚ ℂ)).eval x = 0 ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1) := by
  obtain ⟨lam, nu, mu, S₀, hnu, hlamS, hlamβ, hmemS, hmuS, h0S, hposS, hβS, hβ0, hratio⟩ :=
    exists_lemma22_setup S β hβ
  obtain ⟨f, hdeg, hval, hβne0, hβne1, hcrit⟩ :=
    lemma_2_2 S₀ (nu * (1 / (β - lam)) + mu) h0S hposS hβS hβ0 hratio
  exact ⟨lam, nu, mu, f, hnu, hlamS, hlamβ, hdeg,
    fun α hα => hval _ (hmemS α hα), hval _ hmuS, hβne0, hβne1, hcrit⟩

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def exists_sep_rat.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Theorem 2.5(段 3——Lemma 2.3 の ℚ 版)",
    sectionId := "ncbelyi-thm-2-5" }

def mobius_inj.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Theorem 2.5(段 3——x ↦ ν/(x−λ)+μ は単射)",
    sectionId := "ncbelyi-thm-2-5" }

def exists_lemma22_setup.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Theorem 2.5(段 3——Lemma 2.2 の仮定を作る)",
    sectionId := "ncbelyi-thm-2-5" }

def theorem_2_5_p1_rat.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Theorem 2.5(X = ℙ¹_ℚ、S ⊆ ℙ¹(ℚ)、T = {β} の場合)",
    sectionId := "ncbelyi-thm-2-5" }

def theorem_2_5_p1_rat.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_2_2(段 4、第 398-404)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.lemma_2_2") 3,
    .otherPaper "[NCBelyi]"
      ("Lemma 2.4(段 2——S ⊆ ℙ¹(ℚ) への帰着)。" ++
       "★多項式の段は第 417-418 で閉じた(exists_poly_image_rat)。" ++
       "残るのは ℙ¹ 上の有理関数としての組み立て(不分岐性の帳簿)だけである") 5,
    .otherPaper "[Stacks]"
      "53.4 / 53.5 / 48.27(段 1——Riemann–Roch / Serre 双対性。第 419 で手元にあることを実測)" 4441,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 原文 p.7 の『by applying an automorphism as in " ++
       "Lemma 2.3 [for, say, C = 4], followed by a suitable automorphism of the form " ++
       "x ↦ ν·x + μ』の 1 文を型にして証明した。" ++
       "★**C = 4 でちょうど足りる**——μ ≥ 0(= −min T)が効くので " ++
       "2x ≤ 2·(B/4 + B/4) = B ≤ B + μ = β₀ になる") 6,
    .implicitStep
      ("★★★逸脱: ℙ¹(ℚ) の ∞ は点として扱わず、その像 μ を別に出している" ++
       "(Lemma 2.2 が多項式＝∞ ↦ ∞ で書かれているため)。" ++
       "★また Lemma 2.3 は absInvShift の形(Found/NCBelyi/RatSeparation.lean)ではなく" ++
       "**ℚ の中で書き直した**——S ⊆ ℙ¹(ℚ) に落ちた後は ℂ を経由する必要が無い") 5,
    .implicitStep
      ("★★★★★これで [NCBelyi] Theorem 2.5 に残るのは**段 1(Riemann–Roch)と" ++
       "段 2 の最後(Lemma 2.4 の不分岐性の帳簿)だけ**である。" ++
       "★Theorem 2.5 は [GenEll] Theorem 2.1 の (ii) ⟹ (i) が要求する 2 点のうちの 1 つ" ++
       "(ResearchPaper/mathlib-gap.json の belyi-noncritical)") 7 ]

end ABC3.Found.NCBelyi
