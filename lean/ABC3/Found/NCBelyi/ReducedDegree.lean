import ABC3.Meta.Claim
import Mathlib.FieldTheory.IntermediateField.Adjoin.Basic
import Mathlib.FieldTheory.Minpoly.Field
import Mathlib.LinearAlgebra.Complex.Module
import Mathlib.Algebra.Algebra.Rat

/-!
# [NCBelyi] Lemma 2.4 —— reduced degree と次数の降下(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> the reduced degree of F. Write

## ★★本ファイルが取るもの

`Lemma 2.4` の入れ子帰納法が回るのは、`S′ ≝ f₀(S) ∪ f₀(S₀)` の各点の
**reduced degree(`[F:ℚ] − 1`)が上がらない**からである。原文はそれを 1 文で言う:

原文 (NCBelyi p.5):
> every α ∈f0(S0) [since f0(x),f

その中身は 2 つの不等式である。

| 主張 | ここでの形 |
|---|---|
| `α ∈ S` の像は次数を上げない | `finrank_adjoin_aeval_le` |
| `f₀′` の根 `w` は `[ℚ(w):ℚ] < d₀` | `finrank_adjoin_root_derivative_lt` |

★★どちらも `ℚ(g(x)) ⊆ ℚ(x)` と「最小多項式の次数 ≤ 消す多項式の次数」だけで出る。

## ★★★`ℚ̄` ではなく `ℂ` の中で書く

原文は `ℚ̄ ⊆ ℂ` の中の点を扱う。★ここでは**代数的とは仮定せずに `ℂ` の点として書き**、
代数性が要る補題だけ `IsIntegral ℚ x` を仮定に置く。
★★`Lemma 2.4` の後段(ノルムの評価)が `ℂ` の中で回るので、
`ℚ̄` という別の型を持ち込まないほうが繋ぎやすい。

## ★★★★`redDeg` は `ℕ` の引き算である

原文の reduced degree は `[F:ℚ] − 1` である。`[F:ℚ] ≥ 1` なので `ℕ` の切り捨てで正しいが、
★不等式を扱うときは `finrank` の側で書き、`redDeg` へ落とすのは最後にする。
-/

namespace ABC3.Found.NCBelyi

open Polynomial IntermediateField

/-! ## ★整性と次数 -/

/-- ★`g ≠ 0` を消すなら整である(`g` をモニックに直すだけ)。 -/
theorem isIntegral_of_aeval_eq_zero {x : ℂ} {g : ℚ[X]} (hg : g ≠ 0)
    (hgx : Polynomial.aeval x g = 0) : IsIntegral ℚ x := by
  refine ⟨g * Polynomial.C g.leadingCoeff⁻¹, monic_mul_leadingCoeff_inv hg, ?_⟩
  have h0 : Polynomial.eval₂ (algebraMap ℚ ℂ) x g = 0 := by
    simpa [Polynomial.aeval_def] using hgx
  simp only [Polynomial.eval₂_mul, Polynomial.eval₂_C, h0, zero_mul]

/-- ★★**`[ℚ(x):ℚ] ≤ deg g`**(`g(x) = 0`、`g ≠ 0`)。

原文 (NCBelyi p.5):
> every α ∈f0(S0) [since f0(x),f

★原文の角括弧の中身がこれである——『`f₀′` の次数は `≤ d₀ − 1`』。 -/
theorem finrank_adjoin_le_natDegree {x : ℂ} {g : ℚ[X]} (hg : g ≠ 0)
    (hgx : Polynomial.aeval x g = 0) :
    Module.finrank ℚ ℚ⟮x⟯ ≤ g.natDegree := by
  rw [IntermediateField.adjoin.finrank (isIntegral_of_aeval_eq_zero hg hgx)]
  exact natDegree_le_natDegree (minpoly.degree_le_of_ne_zero ℚ x hg hgx)

/-- ★★★**多項式の像は次数を上げない** —— `ℚ(g(x)) ⊆ ℚ(x)` だから。

★原文はこれを明示的には書かない(『`m(S′) ≤ m(S)`』を当然としている)。
★★しかし帰納法の測度が下がることの**半分はこれ**である。 -/
theorem finrank_adjoin_aeval_le {x : ℂ} (hx : IsIntegral ℚ x) (g : ℚ[X]) :
    Module.finrank ℚ ℚ⟮Polynomial.aeval x g⟯ ≤ Module.finrank ℚ ℚ⟮x⟯ := by
  haveI : FiniteDimensional ℚ ℚ⟮x⟯ := IntermediateField.adjoin.finiteDimensional hx
  have hmem : Polynomial.aeval x g ∈ ℚ⟮x⟯ :=
    (IntermediateField.mem_adjoin_simple_iff ℚ _).2 ⟨g, 1, by simp⟩
  exact IntermediateField.finrank_le_of_le_right (IntermediateField.adjoin_simple_le_iff.2 hmem)

/-! ## ★★★★`f₀′` の根の次数 -/

/-- ★★★★**`f₀′` の根は `d₀` 未満の次数を持つ**。

原文 (NCBelyi p.5):
> every α ∈f0(S0) [since f0(x),f

★`f₀′ ∈ ℚ[x]` で `deg f₀′ ≤ d₀ − 1` だから、そのまま。 -/
theorem finrank_adjoin_root_derivative_lt {f : ℚ[X]} (hdeg : 0 < f.natDegree) {w : ℂ}
    (hw : Polynomial.aeval w (derivative f) = 0) :
    Module.finrank ℚ ℚ⟮w⟯ < f.natDegree := by
  have hd0 : derivative f ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := f)).1 hc
    omega
  have h1 := finrank_adjoin_le_natDegree hd0 hw
  have h2 : (derivative f).natDegree < f.natDegree :=
    Polynomial.natDegree_derivative_lt (by omega)
  omega

/-- ★★★★★**`f₀(S₀)` の点は `d₀` 未満の次数を持つ** —— 上の 2 つを繋いだもの。

原文 (NCBelyi p.5):
> every α ∈f0(S0) [since f0(x),f

★★これが `m(S′) < m(S)` または `d(S′) < d(S)` の**片側**である
——`S₀` から来る点は最大次数の層に入らない。 -/
theorem finrank_adjoin_eval_root_derivative_lt {f : ℚ[X]} (hdeg : 0 < f.natDegree) {w : ℂ}
    (hw : Polynomial.aeval w (derivative f) = 0) :
    Module.finrank ℚ ℚ⟮Polynomial.aeval w f⟯ < f.natDegree := by
  have hd0 : derivative f ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := f)).1 hc
    omega
  exact lt_of_le_of_lt
    (finrank_adjoin_aeval_le (isIntegral_of_aeval_eq_zero hd0 hw) f)
    (finrank_adjoin_root_derivative_lt hdeg hw)

/-! ## ★★★★★★reduced degree -/

/-- **reduced degree** —— 原文の `[F : ℚ] − 1`。

原文 (NCBelyi p.5):
> the reduced degree of F. Write -/
noncomputable def redDeg (x : ℂ) : ℕ := Module.finrank ℚ ℚ⟮x⟯ - 1

/-- ★`[ℚ(x):ℚ] ≥ 1` —— `redDeg` の `ℕ` の引き算が切り捨てにならないこと。 -/
theorem one_le_finrank_adjoin (x : ℂ) (hx : IsIntegral ℚ x) :
    1 ≤ Module.finrank ℚ ℚ⟮x⟯ := by
  haveI : FiniteDimensional ℚ ℚ⟮x⟯ := IntermediateField.adjoin.finiteDimensional hx
  exact Module.finrank_pos

/-- ★**整な点の多項式の像も整**。 -/
theorem isIntegral_aeval {w : ℂ} (hint : IsIntegral ℚ w) (f : ℚ[X]) :
    IsIntegral ℚ (Polynomial.aeval w f) := by
  haveI : FiniteDimensional ℚ ℚ⟮w⟯ := IntermediateField.adjoin.finiteDimensional hint
  have hmem : Polynomial.aeval w f ∈ ℚ⟮w⟯ :=
    (IntermediateField.mem_adjoin_simple_iff ℚ _).2 ⟨f, 1, by simp⟩
  have h1 : IsIntegral ℚ (⟨Polynomial.aeval w f, hmem⟩ : ℚ⟮w⟯) := IsIntegral.of_finite ℚ _
  simpa using h1.map (IsScalarTower.toAlgHom ℚ ℚ⟮w⟯ ℂ)

/-- ★★**有理点の reduced degree は `0`**。 -/
theorem redDeg_ratCast (q : ℚ) : redDeg (q : ℂ) = 0 := by
  have hbot : (ℚ⟮(q : ℂ)⟯ : IntermediateField ℚ ℂ) = ⊥ :=
    IntermediateField.adjoin_simple_eq_bot_iff.2
      (IntermediateField.mem_bot.2 ⟨q, by simp⟩)
  rw [redDeg, hbot, IntermediateField.finrank_bot]

/-- ★★★**`redDeg = 0` は「有理点である」と同値**。

原文 (NCBelyi p.5):
> f and only if m(S) = 0 if and only if S ⊆P1(Q). -/
theorem redDeg_eq_zero_iff {x : ℂ} (hx : IsIntegral ℚ x) :
    redDeg x = 0 ↔ ∃ q : ℚ, x = (q : ℂ) := by
  haveI : FiniteDimensional ℚ ℚ⟮x⟯ := IntermediateField.adjoin.finiteDimensional hx
  constructor
  · intro h
    have h1 : Module.finrank ℚ ℚ⟮x⟯ = 1 := by
      have := one_le_finrank_adjoin x hx
      rw [redDeg] at h
      omega
    have hbot : (ℚ⟮x⟯ : IntermediateField ℚ ℂ) = ⊥ :=
      IntermediateField.finrank_eq_one_iff.1 h1
    have hxm : x ∈ (⊥ : IntermediateField ℚ ℂ) := by
      rw [← hbot]
      exact IntermediateField.mem_adjoin_simple_self ℚ x
    obtain ⟨q, hq⟩ := IntermediateField.mem_bot.1 hxm
    exact ⟨q, by simpa using hq.symm⟩
  · rintro ⟨q, rfl⟩
    exact redDeg_ratCast q

/-- ★★★★**多項式の像は reduced degree を上げない**。 -/
theorem redDeg_aeval_le {x : ℂ} (hx : IsIntegral ℚ x) (g : ℚ[X]) :
    redDeg (Polynomial.aeval x g) ≤ redDeg x := by
  have := finrank_adjoin_aeval_le hx g
  rw [redDeg, redDeg]
  omega

/-- ★★★★★**`f₀(S₀)` の点の reduced degree は `d₀ − 1` 未満**。

★`f₀` が `α₀` の最小多項式なら `d₀ = [ℚ(α₀):ℚ]` なので、
これは `redDeg (f₀(w)) < redDeg α₀` である——★★**最大次数の層に入らない**。 -/
theorem redDeg_eval_root_derivative_lt {f : ℚ[X]} (hdeg : 0 < f.natDegree) {w : ℂ}
    (hw : Polynomial.aeval w (derivative f) = 0) :
    redDeg (Polynomial.aeval w f) < f.natDegree - 1 := by
  have hd0 : derivative f ≠ 0 := by
    intro hc
    have := (Polynomial.derivative_eq_zero (p := f)).1 hc
    omega
  have hint : IsIntegral ℚ w := isIntegral_of_aeval_eq_zero hd0 hw
  have hlt := finrank_adjoin_eval_root_derivative_lt hdeg hw
  have hpos : 1 ≤ Module.finrank ℚ ℚ⟮Polynomial.aeval w f⟯ :=
    one_le_finrank_adjoin _ (isIntegral_aeval hint f)
  rw [redDeg]
  omega

/-! ## ★出典の紐付け(`.src`) -/

def finrank_adjoin_le_natDegree.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(f₀′ の次数が ≤ d₀ − 1 であることから来る次数の限界)",
    sectionId := "ncbelyi-lemma-2-4" }

def finrank_adjoin_aeval_le.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(多項式の像は次数を上げない)",
    sectionId := "ncbelyi-lemma-2-4" }

def finrank_adjoin_eval_root_derivative_lt.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4([ℚ(α′) : ℚ] < d₀ for every α′ ∈ f₀(S₀))",
    sectionId := "ncbelyi-lemma-2-4" }

def redDeg.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(reduced degree ≝ [F : ℚ] − 1)",
    sectionId := "ncbelyi-lemma-2-4" }

def redDeg_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(d(S) = 0 ⟺ m(S) = 0 ⟺ S ⊆ ℙ¹(ℚ))",
    sectionId := "ncbelyi-lemma-2-4" }

end ABC3.Found.NCBelyi
