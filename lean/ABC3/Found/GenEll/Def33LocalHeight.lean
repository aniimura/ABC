/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateCurveWitness

/-!
# [GenEll] Definition 3.3 —— **局所高さ `v_K(q_E)`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★何が揃っていたか

`Found/GaloisRep/TateCurveWitness.lean`（並行セッションの第 309 系列）が
**Tate 母数と一意化を無条件に構成している**。★本ファイルはそれを
`Definition 3.3` の主張の形に束ねるだけである。

| 原文が言うこと | 宣言 |
|---|---|
| `q_E` は `𝔪_K` の元 | `tateParamR_mem` |
| `v_K(q_E)` は**正**の整数 | ★`localHeight_pos_of_split` |
| `v_K(q_E)` は `q_E` の付値そのもの | `localHeight_eq_of_split` |
| `q_E` が Tate 母数であること（一意化） | ★`uniformization_of_split` |

## ★★★原文の仮定と実装の仮定が一致している

原文 (GenEll p.15):
> is proper, while the special fiber Ek of E is isomorphic to (Gm)k, the multiplicative

★「特殊ファイバーが `(𝔾_m)_k` **に同型**」——これは**分裂**乗法還元である
（捻れではない）。実装の `W.HasSplitMultiplicativeReduction R` がそれに当たる。

原文 (GenEll p.15):
> that obtained by extracting an l-th root of the Tate parameter qE ∈ mK . The

★「`q_E ∈ 𝔪_K`」——これが正性の出どころであり、`tateParamR_mem` そのものである。
★★**したがって正性は posit ではなく定理になった**
（`Skeleton/GenEll/Section3.lean` の `localHeight_pos.needs` が
「`vq_pos` を posit したこと自体が仕事の先送りである」と書いていた段）。

## ★★★★逸脱の記録（CLAUDE.md の「逸脱」）

★原文は半アーベル scheme `E → Spec(𝓞_K)` で語るが、実装は
**DVR `R` 上の極小 Weierstrass モデル `W`** で語る（`[W.IsMinimal R]`）。
★★両者の対応（Néron モデル／極小モデルの一意性）は本項目の主張には要らない
——原文が `Definition 3.3` で使うのは `q_E` の付値だけである。

★★★**`𝓞_K/(q_E)` が `M_ell` への分類射による引き戻しに等しい**という
原文の性質（p.15）は**含んでいない**。それは `q_E` の特徴づけであって、
局所高さの定義には要らない。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep

/-- ★★★★★★★★**[GenEll] Definition 3.3**（局所高さ）—— 原文の 4 条をまとめて。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★**「positive integer」が定理になっている**のが本項目の要点である
（`Interface/GenEll/TateLocal.lean` の `vq_pos` は posit だった）。 -/
theorem definition_3_3_bundle
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) :
    tateParamR W h ∈ IsLocalRing.maximalIdeal R
      ∧ 0 < localHeightOf W h
      ∧ (IsDiscreteValuationRing.maximalIdeal R).valuation K ((tateParamK W h : Kˣ) : K)
          = (Multiplicative.ofAdd (-(localHeightOf W h : ℤ)) : Multiplicative ℤ)
      ∧ Nonempty (W.toAffine.Point ≃+ Additive (Kˣ ⧸ Subgroup.zpowers (tateParamK W h))) :=
  ⟨tateParamR_mem W h, localHeight_pos_of_split W h, localHeight_eq_of_split W h,
    uniformization_of_split W h⟩

/-! ### ★★★★★★★★項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
（`tools/genell-progress.mjs` の規則）。 -/

/-- ★★★★★★★★**[GenEll] Definition 3.3**（Local Heights）—— 実装された。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★主張

| 原文 | 宣言 |
|---|---|
| `q_E ∈ 𝔪_K`（p.15） | `tateParamR_mem` |
| `v_K(q_E) ∈ ℤ_{>0}` ——★**正**であること | ★`localHeight_pos_of_split` |
| `v_K(q_E)` が局所高さである（付値そのもの） | `localHeight_eq_of_split` |
| `q_E` が Tate 母数であること（`E(K) ≅ Kˣ/q^ℤ`） | ★`uniformization_of_split` |
| 4 条をまとめたもの | ★`definition_3_3_bundle`（本ファイル） |

## ★★★★逸脱の記録

### 1. 半アーベル scheme ではなく極小 Weierstrass モデルで語る

原文は `E → Spec(𝓞_K)`（特殊ファイバーが `(𝔾_m)_k` に同型な 1 次元半アーベル
scheme）で語るが、実装は DVR `R` 上の極小 Weierstrass モデル `W` で語る。
★両者の対応（Néron モデルの一意性）は本項目の主張には要らない
——原文が `Definition 3.3` で使うのは `q_E` の付値だけである。

★★原文の「特殊ファイバーが `(𝔾_m)_k` **に同型**」は
**分裂**乗法還元（捻れでない）であり、実装の
`W.HasSplitMultiplicativeReduction R` と一致している。

### 2. `𝓞_K/(q_E)` の `M_ell` 的な特徴づけは含まない

原文 (p.15) は `q_E` を「`𝓞_K/(q_E)` が分類射 `Spec(𝓞_K) → M_ell` による
引き戻しに等しい」という性質で特徴づけるが、★それは `q_E` の**一意性**の話であって
局所高さの**定義**には要らない。★★本実装の `q_E` は
`uniformization_of_split`（`E(K) ≅ Kˣ/q^ℤ`）で特徴づけている。 -/
def definition_3_3.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3",
    sectionId := "genell-def-3-3" }

def definition_3_3.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "localHeight_pos_of_split(★原文の「positive integer」——posit ではなく定理)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.localHeight_pos_of_split") 15,
    .citation "[ABC3]" "tateParamR_mem(q_E ∈ 𝔪_K ——正性の出どころ)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateParamR_mem") 15,
    .citation "[ABC3]" "uniformization_of_split(Tate 一意化 E(K) ≅ Kˣ/q^ℤ)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.uniformization_of_split") 15,
    .citation "[ABC3]" "localHeight_eq_of_split(局所高さは q_E の付値そのもの)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.localHeight_eq_of_split") 15,
    .implicitStep
      ("★逸脱 1: 半アーベル scheme E → Spec(𝓞_K) ではなく DVR 上の極小 Weierstrass" ++
       "モデル W で語る。両者の対応(Néron モデルの一意性)は本項目の主張には要らない" ++
       "——原文が Definition 3.3 で使うのは q_E の付値だけである") 15,
    .implicitStep
      ("★逸脱 2: 𝓞_K/(q_E) が分類射 Spec(𝓞_K) → M_ell による引き戻しに等しいという" ++
       "原文 p.15 の特徴づけは含まない。それは q_E の一意性の話であって局所高さの" ++
       "定義には要らない。本実装の q_E は uniformization_of_split で特徴づけている") 15,
    .implicitStep
      ("★★原文の「特殊ファイバーが (𝔾_m)_k に同型」は**分裂**乗法還元(捻れでない)であり、" ++
       "実装の W.HasSplitMultiplicativeReduction R と一致している。非分裂の場合は" ++
       "不分岐 2 次拡大で分裂させると v_K が変わらないが、その段は含んでいない") 15 ]

def definition_3_3_bundle.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(4 条をまとめたもの)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GenEll
