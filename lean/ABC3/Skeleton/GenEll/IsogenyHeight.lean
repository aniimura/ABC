/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.Lemma35Concrete

/-!
# [DelSB616] Théorème 2.4 —— **同種写像による Faltings 高さの変化**（`Skeleton`）

原典: P. Deligne, *Preuve des conjectures de Tate et de Shafarevitch (d'après G. Faltings)*,
Séminaire Bourbaki 616 (1983/84)、物理 p.14。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

## ★★★★★★★★★★なぜこれを載せるのか —— 2026-08-29 の発見

`Found/GaloisRep/Lemma35Concrete.lean` が型で固定した**唯一の残り**

    `ht^Falt(E/H) ≤ ht^Falt(E) + 2·log(l)`

について、台帳は原典を **[FC] Chapter I, Proposition 2.7**（`0_Source` に無い）と
書いていた。★★★**それは足りない見立てだった。**

★`[DelSB616]` の **§2「ISOGÉNIES」**が、まさにこの段を扱っている:

> Dans ce paragraphe, nous étudions la variation de la hauteur `h(A)` par isogénie,
> pour `A` à réduction semi-stable, et en déduisons la conjecture de Tate.

★★そして `[DelSB616]` は **`0_Source` にあり `papers.json` に登記済み**である
（`Skeleton/GenEll/DeligneHeight.lean` が既に使っている）。
★★★したがって CLAUDE.md の進め方（スケルトンで依存グラフを作り、葉から形式化）が
**そのまま適用できる**——「原典が無い」ではなく「まだグラフに載せていない道」である。

## ★★★★★★Deligne の 3 段（テキスト層で読んだ。物理 p.13–14、2026-08-29）

| 段 | 内容 | 原文の位置 |
|---|---|---|
| 1 | 完全列 `0 → A′ → A → A″ → 0` に付随する `ω` の自然な同型（エルミート構造と両立） | `2.1` |
| 2 | `w(H)` の計算——`H` が平坦準有限なら `w(H) = ρ^{-1}`、`w(H)/𝒪` は `#H` で消える | `2.2 (b)(c)` |
| 3 | ★**アルキメデスのエルミート構造は `#H` 倍**である | `2.2 (a)` |

★★`Théorème 2.4` はこの 3 段から出る:「`A` を半安定還元をもつアーベル多様体とすると、
有限個の素数の集合 `F` があって、`l′ ∈ F` と素な次数の同種写像 `u : A → B` について
`h(A)` と `h(B)` は（明示の差で）結ばれる」。

★★★[GenEll] `Lemma 3.5` が要るのは `l`-巡回部分群の場合、すなわち `#H = l` である。

## ★★本ファイルが取るもの

★**3 段が揃ったときに何が出るかを型で固定する**——`Skeleton` の役目である。
★★`Found/GaloisRep/Lemma35Concrete.lean` の仮説 `hfalt` はこの結論そのものなので、
本スケルトンの 3 つの仮説が `Found` になれば `Lemma 3.5` は閉じる。

## ☆残っている葉（明示）

☆(1) 楕円曲線の**同種写像**そのもの（`E/H` を Weierstrass 曲線として作る＝Vélu）。
★mathlib に無い（2026-08-29 に `#check` で確認: `WeierstrassCurve.Isogeny` 不在）。
☆(2) **Néron モデル**と Néron 微分の同種写像の下でのふるまい（`ω_B ⊆ ω_A`、余核が `#H` で消える）。
★mathlib に無い（`AlgebraicGeometry.NeronModel` 不在）。
☆(3) アルキメデスの**周期格子の共体積**の比較（`#H` 倍）。
★本プロジェクトの `Found/GenEll/Covolume.lean`・`CurveArchInv.lean` が受け皿になる。
-/

namespace ABC3.Skeleton.GenEll

/-- ★★★★★★★★**[DelSB616] Théorème 2.4 の到達点** ——
同種写像で Faltings 高さは `log(deg u)` の程度しか動かない。

原文 (DelSB616 p.14):
> F de nombres premiers tels que, pour u : A - B une

★`wfin` は `w(H)` の有限素点側の寄与、`warch` はアルキメデス側の寄与
（原文 `2.2 (a)`「La structure hermitienne de `w(H)` ⊗ ℂ は `#H` 倍」）。
★★`hsplit` は原文 `2.1` の完全列に付随する分解である。

★★★**これが `Found/GaloisRep/Lemma35Concrete.lean` の仮説 `hfalt` を埋める形**である
——`deg u = l` として `ht^Falt(E/H) ≤ ht^Falt(E) + 2·log(l)`。 -/
theorem isogeny_faltings_height
    {P : Type*} (htA htB deg : P → ℝ)
    (_hdeg1 : ∀ p, 1 ≤ deg p)
    -- ★段 2（原文 `2.2 (b)(c)`）: 有限素点側の寄与は `log(#H)` で抑えられる
    (wfin : P → ℝ) (hwfin : ∀ p, wfin p ≤ Real.log (deg p))
    -- ★段 3（原文 `2.2 (a)`）: アルキメデスのエルミート構造は `#H` 倍
    (warch : P → ℝ) (hwarch : ∀ p, warch p ≤ Real.log (deg p))
    -- ★段 1（原文 `2.1`）: 完全列に付随する分解
    (hsplit : ∀ p, htB p = htA p + wfin p + warch p) :
    ∀ p, htB p ≤ htA p + 2 * Real.log (deg p) := by
  intro p
  rw [hsplit p]
  linarith [hwfin p, hwarch p]

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def isogeny_faltings_height.src : ABC3.Meta.Source :=
  { paper := "DelSB616", pdfPage := 14,
    item := "Théorème 2.4(同種写像による Faltings 高さの変化——3 段は仮定に置く)",
    sectionId := "delsb616-thm-2-4" }

/-- ★原文 p.13–14 を**テキスト層で**読んで数えた（2026-08-29）。
☆260 dpi の目視確認はまだ行っていない——結論の式はテキスト層に落ちている。 -/
def isogeny_faltings_height.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★★★★★★★2026-08-29 の発見: [GenEll] Lemma 3.5 に残る唯一の入力 " ++
       "ht^Falt(E/H) ≤ ht^Falt(E) + 2log(l) の原典は、台帳が書いていた " ++
       "[FC] Ch. I, Prop 2.7(0_Source に無い)だけではない。" ++
       "★**[DelSB616] §2「ISOGÉNIES」がまさにこの段を扱っており、0_Source にある**。" ++
       "★★したがって CLAUDE.md の進め方(スケルトンで依存グラフを作り、葉から形式化)が" ++
       "そのまま適用できる——「原典が無い」ではなく「まだグラフに載せていない道」である") 14,
    .implicitStep
      ("★段 1(原文 2.1): 完全列 0 → A′ → A → A″ → 0 に付随する ω の自然な同型。" ++
       "数体の整数環上ではエルミート構造とも両立する。本 statement は仮定 hsplit に置いた") 13,
    .implicitStep
      ("★段 2(原文 2.2 (b)(c)): H が平坦準有限なら w(H) = ρ^{-1} で、" ++
       "w(H)/𝒪 は #H で消える([3] appendice prop. 9)。有限素点側の寄与の評価。" ++
       "本 statement は仮定 hwfin に置いた") 14,
    .implicitStep
      ("★★段 3(原文 2.2 (a)): 『La structure hermitienne de w(H) ⊗ ℂ は #H 倍』。" ++
       "アルキメデス側の寄与。★本プロジェクトの Found/GenEll/Covolume.lean・" ++
       "CurveArchInv.lean(周期格子の共体積)が受け皿になる。本 statement は仮定 hwarch に置いた") 13,
    .folklore
      ("☆葉 (1): 楕円曲線の**同種写像**そのもの(E/H を Weierstrass 曲線として作る＝Vélu)。" ++
       "★mathlib に無い(2026-08-29 に #check で確認: WeierstrassCurve.Isogeny 不在)") 14,
    .folklore
      ("☆葉 (2): **Néron モデル**と Néron 微分の同種写像の下でのふるまい" ++
       "(ω_B ⊆ ω_A、余核が #H で消える)。" ++
       "★mathlib に無い(AlgebraicGeometry.NeronModel 不在)") 14,
    .otherPaper "[GenEll]"
      ("Lemma 3.5 —— ★本節点の結論がその唯一の残り入力である。" ++
       "Found/GaloisRep/Lemma35Concrete.lean の仮説 hfalt がこの形をしている") 17,
    .implicitStep
      ("☆[GenEll] p.17 はこの段を [FC] Chapter I, Proposition 2.7 に帰しているが、" ++
       "その論文は 0_Source に無く papers.json にも登記されていない。" ++
       "★**本節点([DelSB616] Théorème 2.4)が同じ内容を扱う公開の詳細版であり、" ++
       "0_Source にある**。辺の先を検査できるのはこちらである") 14 ]

end ABC3.Skeleton.GenEll
