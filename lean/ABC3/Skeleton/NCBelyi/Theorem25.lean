import ABC3.Meta.Claim
import ABC3.Interface.NCBelyi.BelyiSetup
import ABC3.Found.NCBelyi.Separation
import ABC3.Found.NCBelyi.Lemma22
import ABC3.Found.NCBelyi.BelyiPoly
import ABC3.Found.NCBelyi.Lemma24Chain
import ABC3.Found.NCBelyi.Thm25Step3

/-!
# [NCBelyi] Theorem 2.5 —— Belyi Maps Noncritical at Prescribed Points(`Skeleton`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5–p.6。
**260 dpi 目視確認 2026-08-17**。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★この転写が変えること

`[GenEll] Theorem 2.1` の `.needs` にあった

```
.citation "[Mzk1]" "noncritical Belyi maps" (.absent "mathlib に `Belyi` は 0 件") 11
```

は、**原典が `0_Source` にあるのに「不在」と分類していた**点で正確でなかった。
★本ファイルで statement を型に固定したので、
`.otherPaper "[NCBelyi]" "Theorem 2.5"` へ**格上げできる**。

★★**これは「解決」ではない**——`sorry` は残る。
変わったのは**分類**である:「Lean エコシステムに無い」から
「原典はあるが我々がまだ実装していない」へ。**後者は前者より軽い。**

## ★★記法の実測(この論文に固有、2026-08-17)

★**`ℚ̄`(上線つき)と `ℚ` の区別が主張の内容である。**
`Theorem 2.5` は **`ℚ̄` 上の曲線**について `φ : X → ℙ¹_ℚ̄` を出す(`ℙ¹_ℚ` ではない)。
`pdftotext` は上線を落とすので **`.txt` では両者が同じ `Q` になる**。

★**この論文では小文字ギリシャ文字 `φ` が `pdftotext` に残る**
——[GenEll] では `α β γ δ ε λ μ` がすべて空白に潰れたのと**逆**である。
★★**「ギリシャ文字は落ちる」は論文ごとの性質であって、一般則ではない。**
-/

namespace ABC3.Skeleton.NCBelyi

open ABC3.Meta ABC3.Interface.NCBelyi

/-- **[NCBelyi] Theorem 2.5**(Belyi Maps Noncritical at Prescribed Points)。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★★★★★★★★★★★★★ 2026-08-27——**構成に載せ替えた**(第 425 ブロック)

★**旧 statement(`∀ B : BelyiSetup, …`)は偽であった**
——`Check/NCBelyi/Thm25AxiomGap.lean` の `theorem_2_5_false` で機械検証済み。
`Mor : Curve → Type` は**公理を持たない型のフィールド**なので、
`Mor X := Empty` と置けば存在主張は必ず偽になる。
★★**posit した型は空にもなれる**——`Remark 1.5.1` / `Proposition 1.7` と同じ病の別の面である。

## ★★原文の証明の 4 段のうち、どこまで取れたか

原文 p.6–p.7 の証明は 4 段である:

| 段 | 中身 | 状態 |
|---|---|---|
| 1 | Riemann–Roch / Serre 双対性で `ψ : X → ℙ¹` を作り `T = {β}` に帰着 | ❌ `.needs` |
| 2 | `Lemma 2.4` で `S ⊆ ℙ¹(ℚ)` に帰着 | △ 多項式の段は済(第 417–418) |
| 3 | `Lemma 2.3`(`C = 4`)と `x ↦ ν·x + μ` で `Lemma 2.2` の仮定を作る | ❌ `.needs` |
| 4 | `Lemma 2.2` で閉じる | ✅ 済(第 398–404) |

★★★**本 statement が取るのは 2 つ**である:

- **(a)(b) の `ℙ¹` の場合** —— `S` が代数的数の有限集合なら Belyi 多項式が
  **無条件に**取れる(`Found/NCBelyi/BelyiPoly.lean` の `exists_belyi_poly`、第 424)。
  ★これは古典的な Belyi の定理そのものである。
- **(c) の noncritical の部分** —— 正規化が済んでいれば指定した 1 点 `β` を
  `{0,1}` の外へ逃がせる(`Lemma 2.2`、第 398–404)。

## ★★★★逸脱(明示)

| 項 | 原典 | 形式化 | 理由 |
|---|---|---|---|
| 量化する対象 | `∀ B : BelyiSetup` | **`ℚ[x]` の多項式と `ℂ` の点** | 前者では偽だから |
| 曲線 `X` | 一般の滑らかな固有連結曲線 | **`ℙ¹` のみ** | 段 1(RR)が未了 |
| `T` | 有限集合 | **正規化後の 1 点 `β`** | 段 1 が `T = {β}` へ帰着させる |
| 定義体 `F` | 「`F` 上に取れる」 | **含めない** | 段 1–3 が未了 |

★★★★★**落としたものは `.needs` に全部書いてある**。
段 1(Riemann–Roch / Serre 双対性)は `[Stacks]` 53.4 / 53.5 / 48.27 に原典があり、
第 419 で手元にあることを実測した。 -/
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


/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def theorem_2_5.src : Source :=
  { paper := "NCBelyi", pdfPage := 5, item := "Theorem 2.5",
    sectionId := "ncbelyi-thm-2-5" }

/-- ★原文 p.6 の証明を通読して数えた(260 dpi 目視 2026-08-17)。

★骨格: `T` に点を足して `|T| ≥ 2g_X + 1` とし、`D ≝ Σ_{t∈T}[t]` の直線束 `L = O_X(D)` を取る。
`deg(ω_X ⊗ L^{-1}(x)) ≤ −2` から **Serre 双対性**で `H¹(X, L(−x)) = 0`、
ゆえに `Γ(X,L) ↠ L ⊗ k(x)`。基点を持たない線型系が出て有限射 `ψ : X → ℙ¹_ℚ̄` を得る。
あとは `X = ℙ¹` の場合に帰着し、`Lemma 2.3`・`Lemma 2.4` の**入れ子帰納法**で閉じる。 -/
def theorem_2_5.needs : List ProofObligation :=
  [ .implicitStep
      "★★★本 statement は原文の証明の 4 段のうち**段 4(Lemma 2.2)と ℙ¹ の場合の (a)(b)** を取っている。段 1(Riemann–Roch で T = {β} へ帰着)と段 3(Lemma 2.3 + アフィンで正規化)は未了である。★空虚でないことは Check/NCBelyi/Thm25Witness.lean で確かめてある" 6,
    .implicitStep
      "★★★★★★★★★★2026-08-29 の到達点(第 519-526、§9-981〜988): 原文の証明 4 段のうち**3 段が閉じた**。段 2(Lemma 2.4)は Found/NCBelyi/Lemma24Chain.lean の lemma_2_4_chain、段 3(Lemma 2.3(C=4) ＋ x ↦ νx+μ)は Found/NCBelyi/Thm25Step3.lean の exists_lemma22_setup / theorem_2_5_p1_rat、段 4(Lemma 2.2)は第 398-404。★★★**残るのは段 1(Riemann–Roch / Serre 双対性で ψ : X → ℙ¹ を作り T = {β} に帰着させる段)だけ**である——[Stacks] 53.4 / 53.5 / 48.27 に原典があるが、mathlib には曲線の直線束・Serre 双対性・基点自由な線型系がほぼ無い(LineBundle 0 件、2026-08-16 実測)" 9,
    .otherPaper "[Stacks]" "53.5 Riemann-Roch / 53.4 Duality / 48.27 Duality for proper schemes over fields(段 1 の原典側——第 419 で手元にあることを実測した)" 4441,
    .folklore
      "曲線上の直線束の次数と Riemann–Roch(deg(ω_X ⊗ L^{-1}(x)) ≤ (2g_X−2)−(2g_X+1)+1 = −2 の計算)。★これが段 1 であり、ℙ¹ への帰着を与える" 6,
    .otherPaper "[NCBelyi]" "Lemma 2.4(S ⊆ ℙ¹(ℚ) の場合への帰着)——★★★★★2026-08-29 に**閉じた**(Found/NCBelyi/Lemma24Chain.lean の lemma_2_4_chain、第 526)。有理関数は Chain(Möbius と多項式の交互合成)として表し、極は ∞ へ写るので臨界点から除いてある" 5,
    .otherPaper "[NCBelyi]" "Lemma 2.3(有理係数の自己同型による正規化)——★有理版は第 409 で取った(lemma_2_3_rat)。★★x ↦ ν·x + μ との組み合わせ(段 3)も 2026-08-29 に**閉じた**(Found/NCBelyi/Thm25Step3.lean、第 519)" 4,
    .otherPaper "[NCBelyi]" "Lemma 2.2(入れ子帰納法の受け皿)——★第 398-404 で完成し、sorry 無しで本ファイルにある" 3,
    .implicitStep
      "★★旧 statement(∀ B : BelyiSetup, …)は**偽**であった——Check/NCBelyi/Thm25AxiomGap.lean の theorem_2_5_false で機械検証済み。Mor : Curve → Type を Empty に取れば存在主張は必ず偽になる。2026-08-27 に構成へ載せ替えた(第 425 ブロック)" 5,
    .implicitStep
      "★逸脱: 曲線を ℙ¹ に限り、T を正規化後の 1 点 β に限り、定義体 F の条件を含めていない。いずれも段 1 ・段 3 が未了だからである" 6 ]

/-! ## Lemma 2.3 —— 点の集まりの分離(★実装済み) -/

/-- **[NCBelyi] Lemma 2.3**(Separation of Collections of Points)。

原文 (NCBelyi p.4):
> Lemma 2.3. (Separation of Collections of Points) Let

★★**本 statement は `sorry` ではない**——`Found/NCBelyi/Separation.lean` の実装を参照する。

★原文の証明は 1 行(「it suffices to take λ such |λ − β| is sufficiently small」)だが、
実装では **`ε ≤ δ/4` かつ `C·ε ≤ δ/4`** として「sufficiently small」を明示した。 -/
theorem lemma_2_3 (S : Finset ABC3.Found.NCBelyi.P1C) (C : ℝ) (hC : 0 < C) (b : ℂ)
    (hb : (some b) ∉ S) :
    ∃ lam : ℂ, lam ≠ b ∧ (some lam ∉ S)
      ∧ ∀ p ∈ S, C * ABC3.Found.NCBelyi.absInvShift lam p
          ≤ ABC3.Found.NCBelyi.absInvShift lam (some b) :=
  ABC3.Found.NCBelyi.lemma_2_3 S C hC b hb

def lemma_2_3.src : Source :=
  { paper := "NCBelyi", pdfPage := 4, item := "Lemma 2.3",
    sectionId := "ncbelyi-lemma-2-3" }

/-- ★**空リストは省略ではなく主張である** —— 原文の証明は外部依存を持たない
(「`|λ − β|` を十分小さく取ればよい」だけ)。 -/
def lemma_2_3.needs : List ProofObligation := []

/-! ## Lemma 2.2 —— 有理点で非臨界な Belyi 写像(★実装済み) -/

section Lemma22
open Polynomial

/-- **[NCBelyi] Lemma 2.2**(Belyi Maps Noncritical at Prescribed Rational Points)。

原文 (NCBelyi p.3):
> (Belyi Maps Noncritical at Prescribed Rational Points)

★★**本 statement は `sorry` ではない**——`Found/NCBelyi/Lemma22.lean` の実装を参照する。

`S ⊆ ℙ¹(ℚ)` が (i) `0, ∞ ∈ S`、(ii) `α ∈ S∖{0,∞} → α > 0` を満たし、
`β ∈ ℚ∖S` が (iii) `β/α ≥ 2`(`∀ α ∈ S∖{0,∞}`)を満たすとき、
**非定数多項式 `f ∈ ℚ[x]`** が存在して
(a) `φ(S) ⊆ {0,1,∞}`、(b) `φ(β) ∉ {0,1,∞}`、(c) `φ` は `ℙ¹_ℚ∖{0,1,∞}` 上不分岐。

## ★原文からの読み替え(2 点、いずれも後続に影響しない)

1. **`∞` を落とした。** 多項式は常に `∞ ↦ ∞` で、`∞ ∈ {0,1,∞}` だから
   (a)(c) は `∞` で自動的に満たされる。よってアフィン部分だけで書いてある。
2. **`β/α ≥ 2` を `2·α ≤ β` と書いた。** `α > 0` なので同値であり、
   割り算の場合分けを持ち込まずに済む。

## ★★原文が書いていない段を 1 つ足した

原文の証明は『**so long as |S| ≥ 4**』と書くだけで**基底段を書かない**。
`|S∖{0,∞}| ≤ 1` では 1 次式 `x ↦ c·x` で足りる(`exists_base_scale`)。
★★**1 次式は臨界点を持たない**ので条件 (c) は空虚に成り立つ——これが基底段が易しい理由である。 -/
theorem lemma_2_2 (S : Finset ℚ) (β : ℚ)
    (h0S : (0 : ℚ) ∈ S)
    (hpos : ∀ α ∈ S, α ≠ 0 → 0 < α)
    (hβ : β ∉ S) (hβ0 : β ≠ 0)
    (hratio : ∀ α ∈ S, α ≠ 0 → 2 * α ≤ β) :
    ∃ f : Polynomial ℚ, 0 < f.natDegree
      ∧ (∀ α ∈ S, f.eval α = 0 ∨ f.eval α = 1)
      ∧ f.eval β ≠ 0 ∧ f.eval β ≠ 1
      ∧ (∀ x : ℂ, (derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
          (f.map (algebraMap ℚ ℂ)).eval x = 0 ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1) :=
  ABC3.Found.NCBelyi.lemma_2_2 S β h0S hpos hβ hβ0 hratio

def lemma_2_2.src : Source :=
  { paper := "NCBelyi", pdfPage := 3, item := "Lemma 2.2",
    sectionId := "ncbelyi-lemma-2-2" }

/-- ★原文の証明が実際に依拠しているのは `Lemma 2.1` だけである
(『apply Lemma 2.1 [with, say, C = 2] to the set λ · S』)。 -/
def lemma_2_2.needs : List ProofObligation :=
  [ .otherPaper "[NCBelyi]" "Lemma 2.1(Separating Properties of Belyi Maps)" 2 ]

end Lemma22

end ABC3.Skeleton.NCBelyi
