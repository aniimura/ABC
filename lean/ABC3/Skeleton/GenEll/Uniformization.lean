import ABC3.Meta.Claim
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.VariableChange

/-!
# スケルトン —— **複素楕円曲線の一意化**(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★なぜこの節点を立てるのか——G8 の `htFalt` の穴

第 329 で `FaltingsHeightData` の witness を組んだが、その `htFalt` は **`deg∞/12`** であって
Faltings 高さではない(`Check/GaloisRep/HtFaltNotPinned.lean`)。
★界面が `ht^Falt` を固定していないのは、固定する条件が
**アルキメデス素点での計量**に依るからである。

★★真の Faltings 高さは

    12·d·ht^Falt(E) = Σ_p v_p(Δ_min)·log N(p) − Σ_{v|∞} log( (2π)^{12}·|Δ(Λ_v)| )

で、右辺第 2 項の `Λ_v` は `E ⊗_L ℂ_v ≅ ℂ/Λ_v` の**周期束**である。
★★★したがって `ht^Falt` を固定するには、まず**任意の複素楕円曲線が `ℂ/Λ` であること**
——一意化——が要る。**それが本節点である。**

## ★★★★★2026-08-26 の実測(mathlib の在庫)

| 段 | mathlib |
|---|---|
| 周期対 `PeriodPair`・束 `lattice` | ✅ `Analysis/SpecialFunctions/Elliptic/Weierstrass.lean` |
| `℘`・`℘'`・`g₂`・`g₃` | ✅ 同上 |
| **`℘'² = 4℘³ − g₂℘ − g₃`** | ✅ `derivWeierstrassP_sq` |
| Eisenstein 級数・Dedekind η・モジュラー判別式 `Δ = η²⁴` | ✅ `NumberTheory/ModularForms/` |
| **`℘` の加法定理**(`ℂ/Λ → E(ℂ)` が群準同型) | ★**0 件**(`weierstrassP_add` で 0) |
| **一意化**(任意の `E/ℂ` が `ℂ/Λ`) | ★**0 件** |

★★★★**「束 → 曲線」は揃っており、欠けているのは「曲線 → 束」である。**
★古典的には `j` の全射性(モジュラー関数論)で出る。

## ★★下流

    本節点(一意化)
      → 周期束の共体積(アルキメデス norm)
      → `ω_E` を**計量つき**算術直線束にする
      → (D1)(D2)(D3) の `deg` に載せる
      → `ht^Falt = deg(ω_E)` を固定し、G8 の欠陥 #6 を塞ぐ

★★★★★**2026-08-26 の訂正**: (D1)(D2)(D3) は**すでに witness を持つ**
(`Found/Arakelov/APicWitness.lean`・`ADegBase.lean`・`AHeightWitness.lean`)。
★`Interface/Arakelov/APic.lean` の `waiting` は**古い**。
★★したがって **`deg` の機構はもうある**——`APic` の witness は
`Pic X × Multiplicative (arcCM X)` で、アルキメデス側のデータを持っている。
★★★★★★**残る障害はこの節点(一意化)と、`ω_E` のアルキメデス norm
(周期束の共体積)の組み上げ、そして `Proposition 3.4` の解析的内容である。**
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta

/-- ★束 `Λ` に対応する Weierstrass 曲線 `y² = x³ − (g₂/4)x − (g₃/4)`。

★`℘'² = 4℘³ − g₂℘ − g₃` で `x = ℘`、`y = ℘'/2` と置いた形である。 -/
noncomputable def latticeCurve (P : PeriodPair) : WeierstrassCurve ℂ :=
  ⟨0, 0, 0, -P.g₂ / 4, -P.g₃ / 4⟩

/-- ★★★★★★★**複素楕円曲線の一意化**——任意の `E/ℂ` はある周期束の曲線と同型。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★これが `ht^Falt` を固定するための最下流の葉である。 -/
theorem exists_periodPair (W : WeierstrassCurve ℂ) (hell : W.IsElliptic) :
    ∃ (P : PeriodPair) (C : WeierstrassCurve.VariableChange ℂ), C • W = latticeCurve P := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def latticeCurve.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_periodPair.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_periodPair.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★**2026-08-26 の実測**: mathlib は `PeriodPair`・`lattice`・`℘`・`℘'`・`g₂`・`g₃` と **`℘'² = 4℘³ − g₂℘ − g₃`**(`derivWeierstrassP_sq`)を持つ(`Analysis/SpecialFunctions/Elliptic/Weierstrass.lean`)。★★★つまり**「束 → 曲線」は揃っている**。★欠けているのは逆向き——任意の `E/ℂ` に周期束を付けること——であり、古典的には `j` の全射性(モジュラー関数論)で出る。★★`weierstrassP_add`(加法定理)は 0 件で、`ℂ/Λ → E(ℂ)` が群準同型であることも未証明である" 17,
    .implicitStep
      "★★★★★**なぜ要るか**: 真の Faltings 高さは `12·d·ht^Falt(E) = Σ_p v_p(Δ_min)·log N(p) − Σ_{v|∞} log((2π)^12·|Δ(Λ_v)|)` であり、右辺第 2 項の `Λ_v` は `E ⊗_L ℂ_v ≅ ℂ/Λ_v` の周期束である。★★第 329 の witness の `htFalt := deg∞/12` は右辺第 2 項を捨てた形にあたり、界面が `ht^Falt` を固定していない(`Check/GaloisRep/HtFaltNotPinned.lean`)ことの原因はここにある" 17,
    .implicitStep
      "★★★★**下流の鎖**: 本節点 → 周期束の共体積(アルキメデス norm)→ `ω_E` を計量つき算術直線束にする → (D1)(D2)(D3) の `deg` に載せる → `ht^Falt = deg(ω_E)` を固定 → G8 の欠陥 #6 を塞ぐ。★★★★★**2026-08-26 の訂正**: (D1)(D2)(D3) はすでに witness を持つ(`APicWitness`・`ADegBase`・`AHeightWitness`)。`Interface/Arakelov/APic.lean` の `waiting` は古い。★したがって **`deg` の機構はもうある**——`APic` の witness は `Pic X × Multiplicative (arcCM X)` でアルキメデス側を持つ。★★残る障害は本節点と、`ω_E` のアルキメデス normの組み上げ、そして `Proposition 3.4` の解析的内容である" 17,
    .implicitStep
      "★★★★★★**2026-08-26 の測定(本節点の道筋)**: mathlib は **`WeierstrassCurve.exists_variableChange_of_j_eq`**(分離閉体上で `j` が等しければ同型、`AlgebraicGeometry/EllipticCurve/IsomOfJ.lean`)を持つ。★★したがって本節点は **「j の全射性」**に帰着する——与えられた `j₀ : ℂ` に対し `j(latticeCurve P) = j₀` なる周期対 `P` を作ればよい。★★★しかし **mathlib にモジュラー `j` 関数は 0 件**(`Mathlib/NumberTheory/ModularForms/` を `jInvariant|def j` で検索して 0)。★★★★足場はある——Eisenstein 級数、`Δ = η²⁴`(`Discriminant.lean`)、レベル 1 の次元公式と Sturm 境界(`LevelOne/DimensionFormula.lean`)。★★★★★古典的な道は `M₁₂ = ⟨E₄³, Δ⟩` と値の式で、`E₄³ − λΔ` が ℍ に零点を持つことを出す。見積もり **10-30 ブロック**" 17,
    .implicitStep
      "★★★★★★★**2026-08-26 の追加測定(見積もりの上振れ)**: mathlib の `Elliptic/Weierstrass.lean` は **微分方程式で止まっている**。★同ファイルを `27|discrim|halfPeriod|Injective|surjective` で検索して **0 件**——すなわち (i) 束の判別式 `g₂³ − 27g₃² ≠ 0`、(ii) 半周期と `e₁,e₂,e₃`、(iii) `℘` の 2 対 1 性・全射性 はいずれも無い。★★したがって本節点は「j の全射性」だけではなく、**束の判別式の非消失と `ℂ/Λ ≅ E(ℂ)` の群同型**も必要とする。★★★★見積もりを **25-60 ブロック**に上振れさせる。★★★★★これは上流(mathlib)に入るべき仕事であり、**Galois の義務の中では閉じない**" 17,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, VI.5(複素楕円曲線の一意化)"
      (.absent "mathlib に一意化(任意の E/ℂ が ℂ/Λ)は 0 件(2026-08-26 実測)") 17,
    .otherPaper "GenEll" "Proposition 3.4(Faltings 高さと無限遠因子)" 17 ]

end ABC3.Skeleton.GenEll
