/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.NCBelyi.BelyiPoly
import ABC3.Found.NCBelyi.Lemma22
import ABC3.Meta.Claim

/-!
# 第 1448 ブロック —— **★★★★★★★★
`[NCBelyi] Theorem 2.5` の証明本体を `Found/` へ**

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5–p.6。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★★これは何か——**新しいゲート G8 が拾った 1 件**である

2026-09-03 に `tools/check.mjs` へ足した G8
「`Skeleton/` に `sorry` 0 かつ条なし `.src` の命題があれば `Found/` へ移す」が、
木全体で**ちょうど 1 件**を拾った。それが本定理である。

☆規約は `Skeleton/` = 原典の主張を写す場所、`Found/` = `sorry` の無い実装であり、
`Lemma 4.1` / `Lemma 4.2` / `Corollary 4.3` / `Corollary 4.4`(第 1444)/
`Theorem 3.8`(第 1445)は既にその形になっている。
★中身は 1 文字も変えていない——`exists_belyi_poly` と `lemma_2_2` を組むだけである。

## ★★取れている範囲は変わらない

原文の証明 4 段のうち**段 2–4 が閉じており、段 1(Riemann–Roch / Serre 双対性で
`ψ : X → ℙ¹` を作り `T = {β}` へ帰着させる段)だけが未了**である。
☆逐語の読みと逸脱の宣言は `Skeleton/NCBelyi/Theorem25.lean` の `theorem_2_5` にある。
-/

namespace ABC3.Found.NCBelyi

open ABC3.Meta

/-- **[NCBelyi] Theorem 2.5**(Noncritical Belyi Maps)—— ★**証明本体**
(第 1448 で `Skeleton` から移した)。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

☆statement の逐語の読みと逸脱の宣言(`ℙ¹` のみ・`T = {β}` など)は
`Skeleton/NCBelyi/Theorem25.lean` の `theorem_2_5` にある。 -/
theorem theorem_2_5 (S : Finset ℂ) (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    -- ★(a)(b) の `ℙ¹` の場合 —— 無条件に取れる(古典的 Belyi)
    (∃ F : Polynomial ℚ, 0 < F.natDegree
        ∧ (∀ x ∈ S, Polynomial.aeval x F = 0 ∨ Polynomial.aeval x F = 1)
        ∧ (∀ w : ℂ, Polynomial.aeval w (Polynomial.derivative F) = 0 →
            Polynomial.aeval w F = 0 ∨ Polynomial.aeval w F = 1))
    -- ★(c) noncritical —— 正規化が済んでいれば 1 点 `β` を `{0,1}` の外へ逃がせる
  ∧ (∀ (S₀ : Finset ℚ) (β : ℚ), (0 : ℚ) ∈ S₀ → (∀ α ∈ S₀, α ≠ 0 → 0 < α) →
        β ∉ S₀ → β ≠ 0 → (∀ α ∈ S₀, α ≠ 0 → 2 * α ≤ β) →
        ∃ f : Polynomial ℚ, 0 < f.natDegree
          ∧ (∀ α ∈ S₀, f.eval α = 0 ∨ f.eval α = 1)
          ∧ f.eval β ≠ 0 ∧ f.eval β ≠ 1
          ∧ (∀ x : ℂ, (Polynomial.derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
              (f.map (algebraMap ℚ ℂ)).eval x = 0
                ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1)) :=
  ⟨ABC3.Found.NCBelyi.exists_belyi_poly S hint,
   fun S₀ β h0 hpos hβ hβ0 hratio =>
     ABC3.Found.NCBelyi.lemma_2_2 S₀ β h0 hpos hβ hβ0 hratio⟩

/-! ## ★出典の紐付け(`.src`) -/

def theorem_2_5.src : Source :=
  { paper := "NCBelyi", pdfPage := 5, item := "Theorem 2.5",
    sectionId := "ncbelyi-thm-2-5" }

def theorem_2_5.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_belyi_poly((a)(b) の ℙ¹ の場合。証明済み)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.exists_belyi_poly") 1,
    .citation "[ABC3]" "lemma_2_2((c) noncritical。証明済み)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.lemma_2_2") 1,
    .implicitStep
      ("★★★★**2026-09-03(第 1448)——場所を規約どおりに直した**。" ++
       "☆同日に `tools/check.mjs` へ足したゲート G8" ++
       "(`Skeleton/` に `sorry` 0 かつ条なし `.src` の命題があれば `Found/` へ移す)が" ++
       "木全体で**ちょうど 1 件**拾ったのが本定理である。" ++
       "★中身は 1 文字も変えていない。" ++
       "★★取れている範囲も変わらない——原文の 4 段のうち段 2-4 が閉じており、" ++
       "段 1(Riemann–Roch / Serre 双対性)だけが未了である。") 6 ]

end ABC3.Found.NCBelyi
