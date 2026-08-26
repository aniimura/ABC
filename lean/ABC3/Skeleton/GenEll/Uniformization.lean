import ABC3.Found.GenEll.JSurjective

/-!
# スケルトン —— **複素楕円曲線の一意化**(`Skeleton`)★★★★★★★★**閉じた**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★なぜこの節点を立てたのか——G8 の `htFalt` の穴

第 329 で `FaltingsHeightData` の witness を組んだが、その `htFalt` は **`deg∞/12`** であって
Faltings 高さではない(`Check/GaloisRep/HtFaltNotPinned.lean`)。
★界面が `ht^Falt` を固定していないのは、固定する条件が
**アルキメデス素点での計量**に依るからである。

★★真の Faltings 高さは

    12·d·ht^Falt(E) = Σ_p v_p(Δ_min)·log N(p) − Σ_{v|∞} log( (2π)^{12}·|Δ(Λ_v)| )

で、右辺第 2 項の `Λ_v` は `E ⊗_L ℂ_v ≅ ℂ/Λ_v` の**周期束**である。
★★★したがって `ht^Falt` を固定するには、まず**任意の複素楕円曲線が `ℂ/Λ` であること**
——一意化——が要る。**それが本節点である。**

## ★★★★★★★★2026-08-26 —— **証明が付いた**(第 332-348、17 ブロック)

`sorry` は `Found/GenEll/JSurjective.lean` の `exists_periodPair_of_isElliptic` に置き換わった。

| 段 | ブロック | 内容 |
|---|---|---|
| (i) 判別式の非消失 | 332-346 | `latticeDisc_ne_zero`——`Δ_lat = 4096π¹²·Δ(τ) ≠ 0` |
| `j` の同定 | 346 | `j(ℤ+τℤ) = E₄³/Δ`(古典的モジュラー `j`) |
| (iii) `j` の全射性 | 347-348 | `jFun_surjective`——**開かつ閉** |
| 一意化 | 348 | + mathlib `exists_variableChange_of_j_eq` |

★★★★★**(ii) `ℂ/Λ ≅ E(ℂ)` の群同型は要らなかった**
——本節点は**曲線の同型**(変数変換)であって群同型ではない。
★第 331 の見積もり(25-60 ブロック、「Galois の義務の中では閉じない」)が過大だったのは、
(ii) を要求に数えていたからである。

★★★`j` の全射性は **valence 公式(留数計算)を使わずに**出た:
mathlib の `isCompact_truncatedFundamentalDomain`(切り詰めた基本領域のコンパクト性)と
`AnalyticOnNhd.is_constant_or_isOpen`(開写像定理)で、像が開かつ閉になる。

## ★★下流(ここから先が G8 の残り)

    本節点(一意化)★済
      → 周期束の共体積(アルキメデス norm)
      → `ω_E` を**計量つき**算術直線束にする
      → (D1)(D2)(D3) の `deg` に載せる ★機構は既にある
      → `ht^Falt = deg(ω_E)` を固定し、G8 の欠陥 #6 を塞ぐ

★★(D1)(D2)(D3) は witness を持つ
(`Found/Arakelov/APicWitness.lean`・`ADegBase.lean`・`AHeightWitness.lean`)。
`Interface/Arakelov/APic.lean` の `waiting` は古い。
★★★★★★**残る障害は `ω_E` のアルキメデス norm(周期束の共体積)の組み上げと、
`Proposition 3.4` の解析的内容である。**
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll

/-- ★★★★★★★**複素楕円曲線の一意化**——任意の `E/ℂ` はある周期束の曲線と同型。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★★★★★★2026-08-26: `Found/GenEll/JSurjective.lean` で**証明された**。 -/
theorem exists_periodPair (W : WeierstrassCurve ℂ) (hell : W.IsElliptic) :
    ∃ (P : PeriodPair) (C : WeierstrassCurve.VariableChange ℂ), C • W = latticeCurve P :=
  ABC3.Found.GenEll.exists_periodPair_of_isElliptic W hell

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_periodPair.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_periodPair.needs : List ProofObligation :=
  [ .implicitStep
      "★★★★★★★**2026-08-26 の実測**: mathlib は `PeriodPair`・`lattice`・`℘`・`℘'`・`g₂`・`g₃` と **`℘'² = 4℘³ − g₂℘ − g₃`**(`derivWeierstrassP_sq`)を持つ(`Analysis/SpecialFunctions/Elliptic/Weierstrass.lean`)。★★★つまり**「束 → 曲線」は揃っている**。★欠けていたのは逆向き——任意の `E/ℂ` に周期束を付けること——であり、古典的には `j` の全射性(モジュラー関数論)で出る" 17,
    .implicitStep
      "★★★★★**なぜ要るか**: 真の Faltings 高さは `12·d·ht^Falt(E) = Σ_p v_p(Δ_min)·log N(p) − Σ_{v|∞} log((2π)^12·|Δ(Λ_v)|)` であり、右辺第 2 項の `Λ_v` は `E ⊗_L ℂ_v ≅ ℂ/Λ_v` の周期束である。★★第 329 の witness の `htFalt := deg∞/12` は右辺第 2 項を捨てた形にあたり、界面が `ht^Falt` を固定していない(`Check/GaloisRep/HtFaltNotPinned.lean`)ことの原因はここにある" 17,
    .implicitStep
      "★★★★**下流の鎖**: 本節点 → 周期束の共体積(アルキメデス norm)→ `ω_E` を計量つき算術直線束にする → (D1)(D2)(D3) の `deg` に載せる → `ht^Falt = deg(ω_E)` を固定 → G8 の欠陥 #6 を塞ぐ。★(D1)(D2)(D3) はすでに witness を持つ(`APicWitness`・`ADegBase`・`AHeightWitness`)。`Interface/Arakelov/APic.lean` の `waiting` は古い。★★★★★★**2026-08-26 に本節点が閉じたので、残る障害は `ω_E` のアルキメデス norm の組み上げと、`Proposition 3.4` の解析的内容である**" 17,
    .implicitStep
      "★★★★★★**2026-08-26 の証明(第 347-348)**: mathlib の **`WeierstrassCurve.exists_variableChange_of_j_eq`**(分離閉体上で `j` が等しければ同型、`AlgebraicGeometry/EllipticCurve/IsomOfJ.lean`)に帰着させ、**`j : ℍ → ℂ` の全射性**を建てた。★★mathlib にモジュラー `j` 関数は無かったので `Found/GenEll/JFunction.lean` で `jFun := E₄³/Δ` として定義した。★★★全射性は **valence 公式(留数計算)を使わない**——(a) `j` は上半平面で正則で定数でない(カスプで発散)ので像は**開**(`AnalyticOnNhd.is_constant_or_isOpen`)、(b) `SL(2,ℤ)` 不変性で基本領域に移し、`‖j‖ → ∞` で `im` を抑えると `isCompact_truncatedFundamentalDomain` から像は**閉**、(c) `ℂ` は連結。★★★★★段 (b) のコンパクト性が mathlib にあったのが決め手である" 17,
    .implicitStep
      "★★★★★★★**2026-08-26 の見積もり訂正**: 本節点は当初 **25-60 ブロック**、「上流(mathlib)に入るべき仕事であり、Galois の義務の中では閉じない」と見積もっていた。★★実際には **第 332-348 の 17 ブロック**で閉じた。★★★★差の大半は **(ii) `ℂ/Λ ≅ E(ℂ)` の群同型は要らなかった**ことによる——本節点は**曲線の同型**(変数変換)であって群同型ではない。★見積もりが (ii) を要求に数えていたのが過大の原因である。★★★★★**残る (ii)(`℘` の加法定理、mathlib に 0 件)は本節点の下流ではない**" 17,
    .citation "[Silverman]" "The Arithmetic of Elliptic Curves, VI.5(複素楕円曲線の一意化)"
      (.inProject "ABC3" "Found/GenEll/JSurjective.lean の `exists_periodPair_of_isElliptic`(2026-08-26、sorry なし)") 17,
    .otherPaper "GenEll" "Proposition 3.4(Faltings 高さと無限遠因子)" 17 ]

end ABC3.Skeleton.GenEll
