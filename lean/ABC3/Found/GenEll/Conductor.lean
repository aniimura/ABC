import ABC3.Found.GenEll.ProductFormula

/-!
# [GenEll] Definition 1.5, (ii)(iv) —— 被約化と導手の次数(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.8–p.9。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.8):
> is also an effective Cartier divisor. We shall say that E is reduced if E = Ered.

## ★★何を取るか —— `Proposition 1.6` の非アルキメデス側

原文 (GenEll p.9) の `Proposition 1.6` の証明はこう分かれる:

> the contributions at the **nonarchimedean** primes, from the definition of log-cond_D
> [i.e., involving "(−)_red"] in Definition 1.5, (iv), and, for the contributions at the
> **archimedean** primes, from the fact that the continuous function |s|_L on the compact
> topological space X^arc is bounded.

★**アルキメデス側は `X^arc`(複素解析空間)を要求する**——mathlib に 0 件。
★★**しかし非アルキメデス側は算術因子の上で閉じる**——それが本ファイルである。

> **有効な算術因子 `a` について `deg((a)_red) ≤ deg(a)`**

★これが「`(−)_red` から従う」の中身であり、`X` も `X^arc` も要らない。

## ★`(−)_red` の算術因子版

原文の `(D_x)_red` は「有効 Cartier 因子の被約化」だが、
`𝕍(F)^non` に台を持つ算術因子の水準では **各有限素点の係数を `1` に潰す**ことである
(`ADivRed`)。★アルキメデス側は対象外(原文の `f_x^D` は `𝕍(F)^non` に台を持つ)。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable {F : Type*} [Field F] [NumberField F]

/-! ## ★被約化 -/

/-- **[GenEll] Definition 1.5, (ii)** の `(−)_red`、算術因子の上での形。

原文 (GenEll p.8):
> is also an effective Cartier divisor. We shall say that E is reduced if E = Ered.

★各有限素点の係数を、正なら `1`、そうでなければ `0` に潰す。
`Finsupp.mapRange` が `0 ↦ 0` を要求するので、その条件はここで満たされている。
★アルキメデス側は `0`——原文の `f_x^D` は `𝕍(F)^non` に台を持つからである。 -/
noncomputable def ADivRed (a : ADiv F) : ADiv F :=
  (Finsupp.mapRange (fun n : ℤ => if 0 < n then (1 : ℤ) else 0) (by simp) a.fin, 0)

@[simp] theorem ADivRed_fin (a : ADiv F) (v : FinitePlace F) :
    (ADivRed a).fin v = if 0 < a.fin v then (1 : ℤ) else 0 := rfl

@[simp] theorem ADivRed_arc (a : ADiv F) : (ADivRed a).arc = 0 := rfl

/-- ★`(−)_red` は**冪等**である。原文の `E = E_red` が「被約化の不動点」であることの確認。 -/
theorem adivRed_idem (a : ADiv F) : ADivRed (ADivRed a) = ADivRed a := by
  refine Prod.ext ?_ rfl
  ext v
  simp only [ADivRed, ADiv.fin, Finsupp.mapRange_apply]
  split_ifs <;> omega

/-- ★被約化は**有効**である。 -/
theorem adivRed_isEffective (a : ADiv F) : (ADivRed a).IsEffective := by
  constructor
  · intro v
    rw [ADivRed_fin]
    split_ifs <;> omega
  · intro v
    simp

/-! ## ★★被約化は次数を増やさない —— `Proposition 1.6` の非アルキメデス側 -/

/-- ★各有限素点の剰余体は 2 元以上なので `log q_v > 0`。 -/
theorem log_residueCard_pos (v : FinitePlace F) : 0 < Real.log (residueCard v) := by
  have h1 : 1 < residueCard v := by
    simpa [residueCard] using NumberField.HeightOneSpectrum.one_lt_absNorm v
  exact Real.log_pos (by exact_mod_cast h1)

/-- ★★**`deg((a)_red) ≤ deg(a)`**(`a` が有効なとき)。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

★**これが `Proposition 1.6` の証明の非アルキメデス側**である
(「the contributions at the nonarchimedean primes, from the definition of log-cond_D
[i.e., involving "(−)_red"]」)。

★機構は単純: 有効なら各係数 `n ≥ 0` で、被約化は `n` を `min(n,1)` に潰す。
`log q_v > 0` なので各項が減り、アルキメデス側は被約化で `0` になるが元は非負である。 -/
theorem deg_adivRed_le (a : ADiv F) (ha : a.IsEffective) : deg (ADivRed a) ≤ deg a := by
  classical
  -- 有限側: 台を `a.fin.support` に揃えて項ごとに比べる
  have hsub : (ADivRed a).fin.support ⊆ a.fin.support := by
    intro v hv
    simp only [Finsupp.mem_support_iff, ADivRed_fin] at hv ⊢
    intro hcon
    rw [hcon] at hv
    simp at hv
  have hfin1 : (ADivRed a).fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v))
      = ∑ v ∈ a.fin.support,
          (((ADivRed a).fin v : ℤ) : ℝ) * Real.log (residueCard v) :=
    Finsupp.sum_of_support_subset _ hsub _ (fun v _ => by simp)
  have hfin2 : a.fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v))
      = ∑ v ∈ a.fin.support, ((a.fin v : ℤ) : ℝ) * Real.log (residueCard v) := rfl
  have hle : ∑ v ∈ a.fin.support,
        (((ADivRed a).fin v : ℤ) : ℝ) * Real.log (residueCard v)
      ≤ ∑ v ∈ a.fin.support, ((a.fin v : ℤ) : ℝ) * Real.log (residueCard v) := by
    refine Finset.sum_le_sum fun v _ => ?_
    have hn : 0 ≤ a.fin v := ha.1 v
    have hcoef : (((ADivRed a).fin v : ℤ) : ℝ) ≤ ((a.fin v : ℤ) : ℝ) := by
      rw [ADivRed_fin]
      split_ifs with h
      · exact_mod_cast h
      · exact_mod_cast hn
    exact mul_le_mul_of_nonneg_right hcoef (log_residueCard_pos v).le
  -- アルキメデス側: 被約化は `0`、元は非負
  have harc : (0 : ℝ) ≤ a.arc.sum (fun _ r => r) :=
    Finset.sum_nonneg fun v _ => ha.2 v
  have harc0 : (ADivRed a).arc.sum (fun _ r => r) = 0 := by
    rw [ADivRed_arc]; simp
  rw [deg, deg, harc0, hfin1, hfin2]
  linarith

/-! ## ★出典の紐付け(`.src`) -/

/-- ★**条を明示する**(2026-08-17 夕に修正)。

`ADivRed` は `ADiv F`(= `Spec 𝓞_F` 上の算術因子)の被約化であり、
原文 (ii) が言う**一般の正規ネーター scheme `Z`** の場合ではない。
★`Spec 𝓞_F` は Dedekind なので非零イデアルはすべて可逆であり、
**正則性も Auslander–Buchsbaum も要らない**——だから先に取れた。

★★一般の場合は `RadicalCartier.lean` にあり、そこには
「各茎が UFD」という仮定が残っている(正則 ⟹ UFD が mathlib に無いため)。
★**ラベルが条を書いていないと、一般の場合まで取れたように読める。** -/
def ADivRed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)(Spec 𝓞_F の場合のみ——一般の Z は RadicalCartier.lean)",
    sectionId := "genell-def-1-5" }

/-- ★★**条なしにしてはならない。**

`deg_adivRed_le` は `Proposition 1.6` の**非アルキメデス側だけ**である。
アルキメデス側(「コンパクト空間 `X^arc` 上の連続関数 `|s|_L` が有界」)は
複素解析空間を要求し、mathlib に 0 件である。

★★**2026-08-17 に自分で条なしを付けてしまい、`genell-progress` が 4/24 → 5/24 に
誤って動いた。** 直後の自己監査で発覚し、ここを条つきに直した。
★これは `.src` の 2 値規則が**まさに防ぐために存在する**誤りである。 -/
def deg_adivRed_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.6(証明の非アルキメデス側のみ)",
    sectionId := "genell-prop-1-6" }

end ABC3.Found.GenEll
