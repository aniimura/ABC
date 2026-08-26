import ABC3.Meta.Claim
import ABC3.Interface.NCBelyi.BelyiSetup
import ABC3.Found.NCBelyi.Separation
import ABC3.Found.NCBelyi.Lemma22

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

原文 (NCBelyi p.6):
> Moreover, if X, S, and T are defined over a number field F, then φ may be taken to be defined over F.

`X` を `ℚ̄` 上の滑らかな固有連結曲線、`S, T ⊆ X(ℚ̄)` を `S ∩ T = ∅` なる有限集合とすると、
**射 `φ : X → ℙ¹_ℚ̄` が存在して**
(a) `φ` は `ℙ¹_ℚ̄ \ {0,1,∞}` の上で不分岐、
(b) `φ(S) ⊆ {0,1,∞}`、
(c) `φ(T) ∩ {0,1,∞} = ∅`。
さらに `X, S, T` が数体 `F` 上定義されているなら **`φ` も `F` 上に取れる**。

★★**(c) が「noncritical」の中身である**——古典的な Belyi の定理は (a)(b) だけを言う。
指定した点集合 `T` を `{0,1,∞}` の**外**に逃がせる、という強化がこの論文の主結果。

★**最後の一文(`F` 上に降りる)が `[GenEll] Theorem 2.1` にとって本質的**である
——そこでは**数体上の曲線**を扱うからである。ゆえに statement に含めた。 -/
theorem theorem_2_5 (B : BelyiSetup) (X : B.Curve) (S T : Set (B.Point X))
    (hS : S.Finite) (hT : T.Finite) (hST : S ∩ T = ∅) :
    (∃ φ : B.Mor X, B.UnramifiedOutsideThree φ
        ∧ (∀ s ∈ S, B.app φ s ∈ B.three)
        ∧ (∀ t ∈ T, B.app φ t ∉ B.three))
  ∧ (∀ F : B.NumField, B.CurveOver X F → B.PointsOver S F → B.PointsOver T F →
        ∃ φ : B.Mor X, B.MorOver φ F ∧ B.UnramifiedOutsideThree φ
          ∧ (∀ s ∈ S, B.app φ s ∈ B.three)
          ∧ (∀ t ∈ T, B.app φ t ∉ B.three)) := by
  sorry

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
  [ .folklore
      "曲線上の直線束の次数と Riemann–Roch(`deg(ω_X ⊗ L^{-1}(x)) ≤ (2g_X−2)−(2g_X+1)+1 = −2` の計算)" 6,
    .citation "[NCBelyi]" "Serre 双対性(`Γ(X, ω_X ⊗ L^{-1}(x))` が `H¹(X, L(−x))` の双対)"
      (.unmeasured) 6,
    .otherPaper "[NCBelyi]" "Lemma 2.4(`S ⊆ ℙ¹(ℚ)` の場合への帰着)" 5,
    .otherPaper "[NCBelyi]" "Lemma 2.3(有理係数の自己同型による正規化)" 4,
    .otherPaper "[NCBelyi]" "Lemma 2.2(入れ子帰納法の受け皿)" 3,
    .implicitStep
      "★原文 p.5 の Lemma 2.4 の証明は `m(S)`・`d(S)`(既約次数の最大値と、その最大値を取る点の既約次数の和)についての**入れ子帰納法**である。『for a suitable choice of C』『bounded by some fixed expression in d₀』が定量的に詰められていない" 5,
    .implicitStep
      "★statement の語彙(ℚ̄ 上の曲線・ℙ¹_ℚ̄ への射・不分岐・定義体)を Interface/NCBelyi/BelyiSetup.lean に posit した。**我々は作っていない**" 5 ]

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
    ∃ f : ℚ[X], 0 < f.natDegree
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
