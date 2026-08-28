/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.Lemma24Chain
import ABC3.Found.NCBelyi.Thm25Step3
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★[NCBelyi] `Theorem 2.5` の `X = ℙ¹` の場合（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5–p.7。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`Theorem 2.5` の証明 4 段のうち、**段 2・3・4 が閉じた**:

| 段 | 中身 | 出どころ |
|---|---|---|
| 1 | Riemann–Roch / Serre 双対性で `ψ : X → ℙ¹` を作り `T = {β}` に帰着 | ☆残る |
| 2 | `Lemma 2.4` で `S ⊆ ℙ¹(ℚ)` に帰着 | ★`§9-988` |
| 3 | `Lemma 2.3`(C=4) ＋ `x ↦ νx+μ` で `Lemma 2.2` の仮定を作る | ★`§9-981` |
| 4 | `Lemma 2.2` で閉じる | ★第 398-404 |

★★★**本ファイルは段 2・3・4 を合成して、`X = ℙ¹` の場合の `Theorem 2.5` を取る**:

    `S ⊆ ℂ`（有限・代数的・`Gal`-安定）と `β ∈ ℚ∖S` に対し、鎖 `ch`（＝ `ℚ(x)` の有理関数）で
      (b) `ch(S) ⊆ {0,1}`
      (c) `ch(β) ∉ {0,1}`
      (a) `ch` は `ℙ¹∖{0,1,∞}` 上不分岐（極でない臨界点の値が `{0,1}` に入る）

★★★★これが原文の `Theorem 2.5` の結論そのものである（`X = ℙ¹`、`T = {β}`）。

## ★合成の仕方 —— `Chain.snoc`

`§9-988` の `Chain` は「Möbius → 多項式」を**先頭に**積む形だったので、
末尾に 1 段足す `Chain.snoc` を入れた。★`§9-981` の `x ↦ ν/(x−λ) + μ` の `+μ` は
多項式側へ吸収する（`f.comp (X + C μ)`）。

## ★★残るのは段 1 だけ

☆Riemann–Roch / Serre 双対性（`deg(ω_X ⊗ L^{-1}(x)) ≤ −2 ⟹ H¹(X, L(−x)) = 0
⟹ Γ(X,L) ↠ L ⊗ k(x)`）。★mathlib には曲線の因子・種数・Riemann–Roch が無い
（2026-08-29 に `#check` で確認）。

★`.src` は条つき——`Theorem 2.5` を項目として閉じるには段 1 が要る。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★★★★鎖の末尾に 1 段足す -/

/-- ★★★★★**鎖の末尾に 1 段足す**（`§9-988` の `Chain` は先頭に積む形だった）。 -/
def Chain.snoc : Chain → ℚ → ℚ → ℚ[X] → Chain
  | .nil, lam, c, f => .cons lam c f .nil
  | .cons lam₀ c₀ f₀ rest, lam, c, f => .cons lam₀ c₀ f₀ (rest.snoc lam c f)

/-- ★`snoc` の評価は「もとの鎖の値に 1 段当てる」。 -/
theorem Chain.eval_snoc (ch : Chain) (lam c : ℚ) (f : ℚ[X]) (x : ℂ) :
    (ch.snoc lam c f).eval x = aeval (mob lam c (ch.eval x)) f := by
  induction ch generalizing x with
  | nil => rfl
  | cons lam₀ c₀ f₀ rest ih => simp only [Chain.snoc, Chain.eval, ih]

/-- ★★`snoc` の臨界点は「もとの鎖の臨界点」か「最後の段の臨界点」。 -/
theorem Chain.isCrit_snoc {ch : Chain} {lam c : ℚ} {f : ℚ[X]} {x : ℂ}
    (h : (ch.snoc lam c f).IsCrit x) :
    ch.IsCrit x ∨ (ch.eval x ≠ (lam : ℂ) ∧
      aeval (mob lam c (ch.eval x)) (derivative f) = 0) := by
  induction ch generalizing x with
  | nil =>
      obtain ⟨hne, hcase⟩ := h
      rcases hcase with hd | hc
      · exact Or.inr ⟨hne, hd⟩
      · exact absurd hc (by rw [Chain.IsCrit]; exact id)
  | cons lam₀ c₀ f₀ rest ih =>
      obtain ⟨hne, hcase⟩ := h
      rcases hcase with hd | hc
      · exact Or.inl ⟨hne, Or.inl hd⟩
      · rcases ih hc with h1 | h2
        · exact Or.inl ⟨hne, Or.inr h1⟩
        · exact Or.inr h2

/-- ★★★**有理点での最後の段の値**——`§9-981` の `x ↦ ν/(x−λ)+μ` を
多項式側へ吸収した形（`f.comp (X + C μ)`）。 -/
theorem eval_snoc_rat (lam nu mu : ℚ) (f : ℚ[X]) (q : ℚ) (hq : q ≠ lam) :
    aeval (mob lam nu ((q : ℂ))) (f.comp (X + C mu))
      = ((f.eval (nu * (1 / (q - lam)) + mu) : ℚ) : ℂ) := by
  have hqC : ((q : ℂ)) - (lam : ℂ) ≠ 0 := by
    intro hc
    exact hq (by exact_mod_cast sub_eq_zero.mp hc)
  have hcast : mob lam nu ((q : ℂ)) + (mu : ℂ)
      = ((nu * (1 / (q - lam)) + mu : ℚ) : ℂ) := by
    rw [mob]; push_cast; rw [one_div]
  rw [Polynomial.aeval_comp]
  simp only [map_add, aeval_X, aeval_C]
  rw [show (algebraMap ℚ ℂ) mu = ((mu : ℚ) : ℂ) by simp, hcast]
  have h : (((nu * (1 / (q - lam)) + mu : ℚ)) : ℂ)
      = algebraMap ℚ ℂ (nu * (1 / (q - lam)) + mu) := by simp
  rw [h, Polynomial.aeval_algebraMap_apply]
  simp

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Theorem 2.5`（`X = ℙ¹`） -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[NCBelyi] Theorem 2.5 の
`X = ℙ¹_ℚ̄`、`T = {β}` の場合**。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

`S ⊆ ℂ`（有限・代数的・`Gal`-安定）と `β ∈ ℚ∖S` に対し、鎖 `ch` が取れて:

* **(b)** `ch(S) ⊆ {0,1}`
* **(c)** `ch(β) ∉ {0,1}`
* **(a)** `ch` は `ℙ¹∖{0,1,∞}` 上不分岐（極でない臨界点の値が `{0,1}`）

★★★段 2（`§9-988` の `Lemma 2.4`）・段 3（`§9-981`）・段 4（第 398-404 の `Lemma 2.2`）の
**合成**である。☆残るのは段 1（Riemann–Roch）だけ。 -/
theorem theorem_2_5_p1 (S : Finset ℂ) (β : ℚ)
    (hint : ∀ x ∈ S, IsIntegral ℚ x) (hstab : IsConjStable S) (hβ : ((β : ℂ)) ∉ S) :
    ∃ ch : Chain,
      (∀ x ∈ S, ch.eval x = 0 ∨ ch.eval x = 1)
      ∧ (ch.eval ((β : ℂ)) ≠ 0 ∧ ch.eval ((β : ℂ)) ≠ 1)
      ∧ (∀ w : ℂ, ch.IsCrit w → ch.eval w = 0 ∨ ch.eval w = 1) := by
  classical
  obtain ⟨ch, Sphi, b, hval, hbeq, hbnot, hcrit⟩ := lemma_2_4_chain S β hint hstab hβ
  obtain ⟨lam, nu, mu, f, hnu, hlamS, hlamb, hdeg, hf1, hf2, hf3, hf4, hfcrit⟩ :=
    theorem_2_5_p1_rat Sphi b hbnot
  have hqne : ∀ q ∈ Sphi, q ≠ lam := fun q hq h => hlamS (h ▸ hq)
  refine ⟨ch.snoc lam nu (f.comp (X + C mu)), ?_, ?_, ?_⟩
  · intro x hx
    obtain ⟨q, hqS, hqe⟩ := hval x hx
    rw [Chain.eval_snoc, hqe, eval_snoc_rat lam nu mu f q (hqne q hqS)]
    rcases hf1 q hqS with h | h <;> rw [h] <;> simp
  · constructor
    · rw [Chain.eval_snoc, hbeq, eval_snoc_rat lam nu mu f b (fun h => hlamb h.symm)]
      intro hc
      exact hf3 (by exact_mod_cast hc)
    · rw [Chain.eval_snoc, hbeq, eval_snoc_rat lam nu mu f b (fun h => hlamb h.symm)]
      intro hc
      exact hf4 (by exact_mod_cast hc)
  · intro w hw
    rcases Chain.isCrit_snoc hw with h1 | ⟨hne, hd⟩
    · obtain ⟨q, hqS, hqe⟩ := hcrit w h1
      rw [Chain.eval_snoc, hqe, eval_snoc_rat lam nu mu f q (hqne q hqS)]
      rcases hf1 q hqS with h | h <;> rw [h] <;> simp
    · rw [Chain.eval_snoc]
      set y : ℂ := mob lam nu (ch.eval w) with hy
      have hdC : eval (y + (mu : ℂ)) (derivative (f.map (algebraMap ℚ ℂ))) = 0 := by
        rw [Polynomial.derivative_map, Polynomial.eval_map, ← Polynomial.aeval_def]
        rw [Polynomial.derivative_comp] at hd
        simpa [Polynomial.aeval_comp] using hd
      have hres := hfcrit (y + (mu : ℂ)) hdC
      rw [Polynomial.eval_map, ← Polynomial.aeval_def] at hres
      rw [Polynomial.aeval_comp]
      simp only [map_add, aeval_X, aeval_C]
      rw [show (algebraMap ℚ ℂ) mu = ((mu : ℚ) : ℂ) by simp]
      exact hres

/-! ## ★出典の紐付け(`.src`)——★**条つきである。段 1 が残る** -/

def Chain.snoc.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Theorem 2.5(鎖の末尾に 1 段足す)",
    sectionId := "ncbelyi-thm-2-5" }

def theorem_2_5_p1.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Theorem 2.5(X = ℙ¹、T = {β} の場合——段 2・3・4 の合成)",
    sectionId := "ncbelyi-thm-2-5" }

def theorem_2_5_p1.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "lemma_2_4_chain(段 2、§9-988)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.lemma_2_4_chain") 5,
    .citation "[ABC3]" "theorem_2_5_p1_rat(段 3 ＋ 段 4、§9-981)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.theorem_2_5_p1_rat") 4,
    .otherPaper "[Stacks]"
      ("53.4 / 53.5 / 48.27(段 1——Riemann–Roch / Serre 双対性)。" ++
       "★mathlib には曲線の因子・種数・Riemann–Roch が**無い**" ++
       "(2026-08-29 に #check で確認)。**残る**") 4441,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 550 ブロック): " ++
       "[NCBelyi] Theorem 2.5 の証明 4 段のうち**段 2・3・4 が閉じ、本ファイルで合成した**" ++
       "——`X = ℙ¹_ℚ̄`、`T = {β}` の場合の Theorem 2.5 である。" ++
       "★残るのは段 1(Riemann–Roch / Serre 双対性で ψ : X → ℙ¹ を作り T = {β} に帰着)だけ") 8,
    .implicitStep
      ("★★逸脱: 有理関数を Chain(Möbius と多項式の交互合成)で表している(§9-988 と同じ)。" ++
       "★極は臨界点から除いてある——極は ∞ へ写り、{0,1,∞} の ∞ が受ける。" ++
       "★★β ∈ ℚ とし S の Gal-安定性は仮定として受けている" ++
       "(前者は原文が x ↦ x − β で帰着させる段、後者は conjClosure で作れる)") 6 ]

end ABC3.Found.NCBelyi
