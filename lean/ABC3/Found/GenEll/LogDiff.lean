import ABC3.Found.GenEll.ArithDiv
import Mathlib.RingTheory.DedekindDomain.Different
import Mathlib.RingTheory.DedekindDomain.Factorization

/-!
# [GenEll] Definition 1.5, (iii) —— `log-diff` の算術因子(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.8。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.8):
> (iii) Let x ∈ X(F ) ⊆ X(Q), where F is a minimal field of definition of x. Then the diﬀerent ideal of F determines an eﬀective arithmetic divisor

## ★★実装して分かった最も重要なこと —— `δ_x` は `X` にも `x` にも依らない

原文は `δ_x` と書き、`x` の添字を付ける。しかし定義を読むと

> **the different ideal of `F`** determines an effective arithmetic divisor `δ_x ∈ ADiv(F)`

であり、**使っているのは `F`(= `x` の最小定義体)の差積イデアルだけ**である。
★すなわち **`δ_x` は `F` の関数であって、`X` にも `x` にも依らない**。
添字 `x` は「`F` を `x` から取る」という手順を表しているにすぎない。

★★**これは形式化の負荷を大きく変える。**
`log-diff_X` の定義には **scheme `X` が要らない**——数体 `F` の差積イデアルだけでよい。
(いっぽう (iv) の `log-cond_D` は因子 `D ⊆ X` の引き戻しを要求するので、
そちらは scheme 論が要る。★**(iii) と (iv) は要求する道具が違う。**)

## ★mathlib 実測(2026-08-17)

| 要るもの | mathlib |
|---|---|
| 差積イデアル | ★**ある**(`IsDedekindDomain.differentIdeal`) |
| イデアルの素イデアル分解の重複度 | ★**ある**(`Associates.mk v.asIdeal |>.count`) |
| その台の有限性 | ★**ある**(`Associates.finite_factors`) |
| 算術因子 `ADiv(F)` と正規化次数 | ★**我々が実装済み**(`ArithDiv.lean`) |

★ゆえに `Definition 1.5, (iii)` の**算術因子側は今すぐ組める**。
これが `ResearchPaper/genell-goal.md` §9-4 の分解でいう **A1** である。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable (F : Type*) [Field F] [NumberField F]

/-! ## ★イデアルから算術因子へ -/

open scoped Classical in
/-- **イデアルが定める算術因子**。

有限素点 `v` での係数は `v` が `I` の分解に現れる重複度、アルキメデス側は `0`。
★`I = 0` のときは `0` と置く(原文は `I ≠ 0` の場合しか使わない)。 -/
noncomputable def idealADiv (I : Ideal (𝓞 F)) : ADiv F :=
  if h : I = 0 then 0
  else
    (Finsupp.ofSupportFinite
      (fun v : FinitePlace F =>
        ((Associates.mk v.asIdeal).count (Associates.mk I).factors : ℤ))
      (by
        have hfin := Associates.finite_factors (R := 𝓞 F) h
        rw [Filter.eventually_cofinite] at hfin
        exact hfin), 0)

theorem idealADiv_arc (I : Ideal (𝓞 F)) : (idealADiv F I).arc = 0 := by
  classical
  unfold idealADiv
  split_ifs with h
  · rfl
  · rfl

/-- ★**原文「effective arithmetic divisor」** —— 係数は非負である。

重複度は自然数なので、`ℤ` へ落としても非負のままである。 -/
theorem idealADiv_isEffective (I : Ideal (𝓞 F)) : (idealADiv F I).IsEffective := by
  classical
  constructor
  · intro v
    unfold idealADiv
    split_ifs with h
    · exact le_rfl
    · exact Int.natCast_nonneg _
  · intro v
    rw [idealADiv_arc]
    exact le_rfl

/-! ## ★差積イデアルの算術因子 -/

/-- **[GenEll] Definition 1.5, (iii)** の `δ_x`。

原文 (GenEll p.8):
> (iii) Let x ∈ X(F ) ⊆ X(Q), where F is a minimal field of definition of x. Then the diﬀerent ideal of F determines an eﬀective arithmetic divisor

★★**添字は `x` だが、実際に依存しているのは `F` だけである**(上の docstring)。
ゆえに本実装は `F` のみを引数に取る。 -/
noncomputable def differentADiv : ADiv F :=
  idealADiv F (differentIdeal ℤ (𝓞 F))

/-- ★原文「effective」。 -/
theorem differentADiv_isEffective : (differentADiv F).IsEffective :=
  idealADiv_isEffective F _

/-- **[GenEll] Definition 1.5, (iii)** の `log-diff`。

原文 (GenEll p.8):
> determines a well-deﬁned log-diﬀerent function log-diﬀX on X(Q).

`log-diff_X(x) ≝ deg_F(δ_x)`(★`deg` は**下線つき** = 正規化次数)。

★★**`X` は引数に現れない。** 上で見たとおり `δ_x` が `F` だけの関数だからである。
原文が `log-diff_X` と書くのは、`x ∈ X(ℚ̄)` の側から `F` を取る手順を示すためである。

★`Remark 1.5.1` が「`log-diff_X` はスキーム `X_ℚ` だけに依る」と述べるのは、
この観察の**上位互換**である——実際には `X_ℚ` にすら依らず、点の最小定義体だけに依る。 -/
noncomputable def logDiffOfField : ℝ :=
  degNormalized (differentADiv F)

/-! ## ★出典の紐付け(条つき——`Definition 1.5` は (i)(ii)(iv) が残る) -/

def differentADiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5, (iii)",
    sectionId := "genell-def-1-5" }

def logDiffOfField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5, (iii)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
