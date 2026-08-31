/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Elementary

/-!
# [DelSB616] Exemple 1.4 —— **段 5・6（初等評価）**（`Found`）

原典: P. Deligne, *Preuve des conjectures de Tate et de Shafarevitch (d'après G. Faltings)*,
Séminaire Bourbaki 616 (1983/84)、物理 p.5。

原文 (DelSB616 p.5):
> Exemple 1.4 Le cas k = Q , g = 1 .- Soient E une courbe elliptique sur Q et w

## ★★★★★★★★★`Skeleton/GenEll/DeligneHeight.lean` が測った 6 段のうち 2 段

`Skeleton/GenEll/DeligneHeight.lean` は `[DelSB616] Exemple 1.4` を 6 段に割り、
**段 5 と段 6 を「初等」**と見積もっていた。★本ファイルはその 2 段を取る。

| 段 | 内容 | 本ファイル |
|---|---|---|
| 5 | `Δ` は整数なので `\|j\| ≤ H*^{12}`、ゆえに `H* ≤ C₁·H·√(max 1 (log H*))` | `sqrt_bound_of_j_bound` |
| 6 | ★**`H* ≤ C₂·H·√(max 1 (log H))`** | `exists_bound_of_self_referential_sqrt_log` |

★★段 6 が要である——**自己参照的な評価を右辺から消す**。

## ★★★★★★機構（段 6）

`A ≤ C₁·B·√(max 1 (log A))` の両辺の対数を取ると

    log A ≤ log C₁ + log B + (1/2)·log(log A)

★`log x ≤ x − 1`（`Real.log_le_sub_one_of_pos`）から `log(log A) ≤ log A` なので、

    (1/2)·log A ≤ log C₁ + log B   すなわち   log A ≤ 2·log C₁ + 2·log B

★★これで右辺の `log A` が `log B` に置き換わる。定数は `C₁·√(2 log C₁ + 2)`。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★**楕円曲線の言葉を一切使わない**——`H*`・`H`・`|j|`・`Im τ` を
**実数の変数**として述べる。
★★段 1〜4（複素解析の側：`(E,ω) ≅ (ℂ/Λ, dz)`、`Im τ ∼ sup(1, log|j|)/2π`）が
揃ったときに、そのまま代入できる形である。
★★★これは `Found/GenEll/Elementary.lean` の `lemma_3_6` と同じ立場である
（原文の初等評価を実数の言葉で取る）。
-/

namespace ABC3.Found.GenEll

open Real

/-! ## ★★★★★段 6 —— 自己参照的な評価を消す -/

/-- ★★★★★★★★**自己参照的な評価から右辺の `log A` を消す**。

原文 (DelSB616 p.5):
> Exemple 1.4 Le cas k = Q , g = 1 .- Soient E une courbe elliptique sur Q et w

> `A ≤ C₁·B·√(max 1 (log A))`  ⟹  `A ≤ C₂·B·√(max 1 (log B))`

★★機構は対数を取って `log(log A) ≤ log A` を使うだけである
（`Real.log_le_sub_one_of_pos`）。★定数は `C₁·√(2 log C₁ + 2)`。

★★★これが `[DelSB616] Exemple 1.4` の**段 6**（`H*` の評価を `H` の評価に直す段）である。 -/
theorem exists_bound_of_self_referential_sqrt_log {C₁ : ℝ} (hC₁ : 1 ≤ C₁) :
    ∃ C₂ : ℝ, 0 < C₂ ∧ ∀ A B : ℝ, 1 ≤ A → 1 ≤ B →
      A ≤ C₁ * B * Real.sqrt (max 1 (Real.log A)) →
      A ≤ C₂ * B * Real.sqrt (max 1 (Real.log B)) := by
  have hC₁pos : (0:ℝ) < C₁ := lt_of_lt_of_le one_pos hC₁
  have hlogC₁ : (0:ℝ) ≤ Real.log C₁ := Real.log_nonneg hC₁
  set K : ℝ := 2 * Real.log C₁ + 2 with hK
  have hKpos : (0:ℝ) < K := by rw [hK]; linarith
  refine ⟨C₁ * Real.sqrt K, by positivity, ?_⟩
  intro A B hA hB hle
  have hApos : (0:ℝ) < A := lt_of_lt_of_le one_pos hA
  have hBpos : (0:ℝ) < B := lt_of_lt_of_le one_pos hB
  have hlogA : (0:ℝ) ≤ Real.log A := Real.log_nonneg hA
  have hlogB : (0:ℝ) ≤ Real.log B := Real.log_nonneg hB
  have hmB1 : (1:ℝ) ≤ max 1 (Real.log B) := le_max_left _ _
  have hkey : max 1 (Real.log A) ≤ K * max 1 (Real.log B) := by
    rcases le_or_gt (Real.log A) 1 with h | h
    · have hm : max 1 (Real.log A) = 1 := max_eq_left h
      rw [hm]
      have hK1 : (1:ℝ) ≤ K := by rw [hK]; linarith
      calc (1:ℝ) = 1 * 1 := by ring
        _ ≤ K * max 1 (Real.log B) := mul_le_mul hK1 hmB1 zero_le_one (le_of_lt hKpos)
    · have hmA : max 1 (Real.log A) = Real.log A := max_eq_right (le_of_lt h)
      have hLpos : (0:ℝ) < Real.log A := lt_trans one_pos h
      have hsqL : Real.sqrt (max 1 (Real.log A)) = Real.sqrt (Real.log A) := by rw [hmA]
      have hlog := Real.log_le_log hApos hle
      rw [Real.log_mul (by positivity) (by positivity),
        Real.log_mul (by positivity) (by positivity),
        hsqL, Real.log_sqrt (le_of_lt hLpos)] at hlog
      have hloglog : Real.log (Real.log A) ≤ Real.log A := by
        have := Real.log_le_sub_one_of_pos hLpos
        linarith
      have hhalf : Real.log A ≤ 2 * Real.log C₁ + 2 * Real.log B := by linarith
      rw [hmA]
      have h1 : Real.log B ≤ max 1 (Real.log B) := le_max_right _ _
      rw [hK]
      nlinarith [hlogC₁, h1, hmB1]
  have hsq := Real.sqrt_le_sqrt hkey
  calc A ≤ C₁ * B * Real.sqrt (max 1 (Real.log A)) := hle
    _ ≤ C₁ * B * Real.sqrt (K * max 1 (Real.log B)) :=
        mul_le_mul_of_nonneg_left hsq (by positivity)
    _ = C₁ * Real.sqrt K * B * Real.sqrt (max 1 (Real.log B)) := by
        rw [Real.sqrt_mul (le_of_lt hKpos)]; ring

/-! ## ★★★★段 5 —— `j` の評価を `H*` の評価に直す -/

/-- ★★★★★★**`|j| ≤ H*^{12}` と `Im τ ≤ c·max 1 (log|j|)` から
`H* ≤ C₁·H·√(max 1 (log H*))`**。

原文 (DelSB616 p.5):
> Exemple 1.4 Le cas k = Q , g = 1 .- Soient E une courbe elliptique sur Q et w

★`Δ` が整数であることから `|j| = |c₄|³/|Δ| ≤ |c₄|³ ≤ H*^{12}` が出る
（本ファイルはその不等式を**仮説として受ける**——段 1〜4 の側である）。
★★あとは `log` を取って `12` を外に出すだけ。 -/
theorem sqrt_bound_of_j_bound {c C₀ : ℝ} (hc : 0 < c) (hC₀ : 0 < C₀) :
    ∃ C₁ : ℝ, 0 < C₁ ∧ ∀ (Hstar H jabs imTau : ℝ),
      1 ≤ Hstar → 0 < H → 0 ≤ imTau → 1 ≤ jabs →
      jabs ≤ Hstar ^ (12:ℕ) →
      imTau ≤ c * max 1 (Real.log jabs) →
      Hstar ≤ C₀ * H * Real.sqrt imTau →
      Hstar ≤ C₁ * H * Real.sqrt (max 1 (Real.log Hstar)) := by
  refine ⟨C₀ * Real.sqrt (12 * c), by positivity, ?_⟩
  intro Hstar H jabs imTau hHs hH himτ hj hjle hτ hle
  have hHspos : (0:ℝ) < Hstar := lt_of_lt_of_le one_pos hHs
  have hlogHs : (0:ℝ) ≤ Real.log Hstar := Real.log_nonneg hHs
  have hlogj : Real.log jabs ≤ 12 * Real.log Hstar := by
    have h1 : Real.log jabs ≤ Real.log (Hstar ^ (12:ℕ)) :=
      Real.log_le_log (lt_of_lt_of_le one_pos hj) hjle
    rwa [Real.log_pow, Nat.cast_ofNat] at h1
  have hmax : max 1 (Real.log jabs) ≤ 12 * max 1 (Real.log Hstar) := by
    have h1 : (1:ℝ) ≤ max 1 (Real.log Hstar) := le_max_left _ _
    have h2 : Real.log Hstar ≤ max 1 (Real.log Hstar) := le_max_right _ _
    exact max_le (by linarith) (by linarith)
  have hτ2 : imTau ≤ 12 * c * max 1 (Real.log Hstar) := by nlinarith [le_of_lt hc]
  have hsq : Real.sqrt imTau ≤ Real.sqrt (12 * c * max 1 (Real.log Hstar)) :=
    Real.sqrt_le_sqrt hτ2
  calc Hstar ≤ C₀ * H * Real.sqrt imTau := hle
    _ ≤ C₀ * H * Real.sqrt (12 * c * max 1 (Real.log Hstar)) :=
        mul_le_mul_of_nonneg_left hsq (by positivity)
    _ = C₀ * Real.sqrt (12 * c) * H * Real.sqrt (max 1 (Real.log Hstar)) := by
        rw [Real.sqrt_mul (by positivity : (0:ℝ) ≤ 12 * c)]
        ring

/-! ## ★★★★★★★★★段 5 → 段 6 の合成 -/

/-- ★★★★★★★★★**[DelSB616] Exemple 1.4 の段 5・6 の到達点**。

原文 (DelSB616 p.5):
> Exemple 1.4 Le cas k = Q , g = 1 .- Soient E une courbe elliptique sur Q et w

> `|j| ≤ H*^{12}`・`Im τ ≤ c·max 1 (log|j|)`・`H* ≤ C₀·H·√(Im τ)`
>   ⟹ **`H* ≤ C₂·H·√(max 1 (log H))`**

★★これが原文の「素朴な高さは Faltings 高さで（対数の補正つきで）抑えられる」の
**初等的な側**である。

★★★残るのは段 1〜4（複素解析の側）:
`(E,ω) ≅ (ℂ/Λ, dz)`・`H(E) = (1/|ω₁|)·|Im τ/π|^{-1/2}`・
`H* ≤ C₀·H·√(Im τ)`・`Im τ ∼ sup(1, log|j|)/2π`。
★並行セッションが `LatticeCurve.lean` ほかで構築中である。 -/
theorem naive_bound_of_chain {c C₀ : ℝ} (hc : 0 < c) (hC₀ : 0 < C₀) :
    ∃ C₂ : ℝ, 0 < C₂ ∧ ∀ (Hstar H jabs imTau : ℝ),
      1 ≤ Hstar → 1 ≤ H → 0 ≤ imTau → 1 ≤ jabs →
      jabs ≤ Hstar ^ (12:ℕ) →
      imTau ≤ c * max 1 (Real.log jabs) →
      Hstar ≤ C₀ * H * Real.sqrt imTau →
      Hstar ≤ C₂ * H * Real.sqrt (max 1 (Real.log H)) := by
  obtain ⟨C₁, hC₁pos, hC₁⟩ := sqrt_bound_of_j_bound hc hC₀
  obtain ⟨C₂, hC₂pos, hC₂⟩ :=
    exists_bound_of_self_referential_sqrt_log (C₁ := max C₁ 1) (le_max_right _ _)
  refine ⟨C₂, hC₂pos, ?_⟩
  intro Hstar H jabs imTau hHs hH himτ hj hjle hτ hle
  have hstep1 := hC₁ Hstar H jabs imTau hHs (lt_of_lt_of_le one_pos hH) himτ hj hjle hτ hle
  have hstep1' : Hstar ≤ max C₁ 1 * H * Real.sqrt (max 1 (Real.log Hstar)) := by
    refine le_trans hstep1 ?_
    have hle1 : C₁ ≤ max C₁ 1 := le_max_left _ _
    have hpos : (0:ℝ) ≤ H * Real.sqrt (max 1 (Real.log Hstar)) := by positivity
    nlinarith [hpos]
  exact hC₂ Hstar H hHs hH hstep1'

/-! ### ★出典の紐付け(`.src`)

★★`[DelSB616] Exemple 1.4` の**段 5・6 だけ**である。段 1〜4（複素解析）は含まない。 -/

def exists_bound_of_self_referential_sqrt_log.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 5,
    item := "Exemple 1.4(段 6——自己参照的な評価から右辺の log A を消す)",
    sectionId := "delsb616-ex-1-4" }

def sqrt_bound_of_j_bound.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 5,
    item := "Exemple 1.4(段 5——|j| ≤ H*^12 から H* の評価へ)",
    sectionId := "delsb616-ex-1-4" }

def naive_bound_of_chain.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 5,
    item := "Exemple 1.4(段 5・6 の合成。段 1〜4 の複素解析は含まない)",
    sectionId := "delsb616-ex-1-4" }

def naive_bound_of_chain.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Real.log_le_sub_one_of_pos(log x ≤ x − 1)"
      (.inMathlib "Real.log_le_sub_one_of_pos") 5,
    .implicitStep
      ("★★★★段 6 の要は「自己参照的な評価を右辺から消す」ことである" ++
       "——log を取ると log(log A) ≤ log A で (1/2)log A が吸収される") 5,
    .implicitStep
      ("★★残る段 1〜4 は複素解析の側である: (E,ω) ≅ (ℂ/Λ, dz)・" ++
       "H(E) = (1/|ω₁|)·|Im τ/π|^{-1/2}・H* ≤ C₀·H·√(Im τ)・" ++
       "Im τ ∼ sup(1, log|j|)/2π。★並行セッションが LatticeCurve.lean ほかで構築中") 5,
    .implicitStep
      ("★逸脱: 楕円曲線の言葉を使わず、H*・H・|j|・Im τ を実数の変数として述べた。" ++
       "★★段 1〜4 が揃ったときにそのまま代入できる形である") 5 ]

end ABC3.Found.GenEll
