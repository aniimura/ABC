/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.FiniteFromNorthcott

/-!
# [GenEll] Lemma 3.5 の組み立て —— **3 つの入力から結論が出る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Thus, Lemma 3.5 follows from Proposition 3.4.

## ★★★原文の最後の 1 文

`Lemma 3.5` の証明は、3 つを並べて `Proposition 3.4` に流すだけである:

| 入力 | 出どころ | 状態 |
|---|---|---|
| `deg(E_H) = l·deg(E)` | ★`Lemma 3.2, (ii)` | ★★**済**（`Lemma32QuotMu.lean` の `vAdd_pow`） |
| `ht_Falt(E_H) = ht_Falt(E) + 2·log(l)` | 原文「integrating a (1,1)-form … differs … by a factor of `l`」 | 未（アルキメデス解析） |
| `κ·deg(E_H) ≤ ht_Falt(E_H) + C` | `Proposition 3.4` を `E_H` に当てる | 未（`[Silv2]`） |

★★**本ファイルはその組み立てだけを取る**——代入 2 回である。

## ★★★★なぜ組み立てを別に取るのか

★原文が「Thus, … follows from …」と書いた段は、**入力が揃ったときに何が出るかを
型で固定する**ことに意味がある。★★入力の 1 つ（`deg(E_H) = l·deg(E)`）は
既に手元にあるので、**残り 2 つが何かがこれで明示される**。

★★★`Lemma 3.7` / `Theorem 3.8` が使うのは結論そのものではなく
**`l` の上界**の形なので、そちらも取る（`l_le_of_lemma_3_5`）。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★実数の水準で述べる——曲線の型は取らない。★★原文の `1/(12(1+ε))` は
`κ` として受ける（正であることだけ使う）。
-/

namespace ABC3.Found.GenEll

/-- ★★★★★★**[GenEll] Lemma 3.5 の組み立て**。

原文 (GenEll p.17):
> Thus, Lemma 3.5 follows from Proposition 3.4.

★3 つの入力（`deg(E_H) = l·deg(E)`、`ht_Falt(E_H) = ht_Falt(E) + 2log l`、
`Proposition 3.4` を `E_H` に当てたもの）から結論が出る。★★代入 2 回である。

★★★1 つめは既に手元にある——`Lemma32QuotMu.lean` の `vAdd_pow`
（`v_K(q^l) = l·v_K(q)`）がそれである。 -/
theorem lemma_3_5_combine (degE degH htFalt htFaltH : ℝ) (l : ℕ) (κ C : ℝ)
    (hdegH : degH = (l : ℝ) * degE)
    (hfaltH : htFaltH = htFalt + 2 * Real.log l)
    (hprop34 : κ * degH ≤ htFaltH + C) :
    κ * ((l : ℝ) * degE) ≤ htFalt + 2 * Real.log l + C := by
  rw [← hdegH]
  calc κ * degH ≤ htFaltH + C := hprop34
    _ = htFalt + 2 * Real.log l + C := by rw [hfaltH]

/-- ★★★★★**`l` の上界の形** —— `Lemma 3.7` / `Theorem 3.8` が使うのはこちら。

原文 (GenEll p.17):
> Thus, Lemma 3.5 follows from Proposition 3.4.

★`Lemma 3.5` の不等式を `l` について解いた形である。
★★`deg(E) > 0`（`Definition 3.3` の局所高さの正性から出る）が要る。 -/
theorem l_le_of_lemma_3_5 (degE htFalt : ℝ) (l : ℕ) (κ C : ℝ)
    (hκ : 0 < κ) (hdeg : 0 < degE)
    (h : κ * ((l : ℝ) * degE) ≤ htFalt + 2 * Real.log l + C) :
    (l : ℝ) ≤ (htFalt + 2 * Real.log l + C) / (κ * degE) := by
  rw [le_div_iff₀ (mul_pos hκ hdeg)]
  calc (l : ℝ) * (κ * degE) = κ * ((l : ℝ) * degE) := by ring
    _ ≤ htFalt + 2 * Real.log l + C := h

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Lemma 3.5` には
`ht_Falt(E_H) = ht_Falt(E) + 2log l`（アルキメデス解析）と
`Proposition 3.4`（`[Silv2]`）が残っている。 -/

def lemma_3_5_combine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(組み立て——3 つの入力から結論。入力自体は含まない)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_combine.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[GenEll]"
      ("Lemma 3.2, (ii)(deg(E_H) = l·deg(E))——★済。Lemma32QuotMu.lean の vAdd_pow " ++
       "(v_K(q^l) = l·v_K(q))がそれである") 15,
    .implicitStep
      ("原文「since integrating a (1,1)-form over Ev differs from integrating over " ++
       "(EH)v by a factor of l」——ht_Falt(E_H) = ht_Falt(E) + 2·log(l)。" ++
       "★アルキメデス解析であり未実装(並行セッションが Covolume.lean 等で構築中)") 17,
    .otherPaper "[GenEll]" "Proposition 3.4 を E_H に当てる——★[Silv2] に依るので未" 17,
    .implicitStep
      ("★原文が「Thus, … follows from …」と書いた段は、入力が揃ったときに何が出るかを" ++
       "型で固定することに意味がある。入力の 1 つは既に手元にあるので、" ++
       "残り 2 つが何かがこれで明示される") 17,
    .implicitStep
      ("★逸脱: 実数の水準で述べる——曲線の型は取らない。原文の 1/(12(1+ε)) は κ として" ++
       "受ける(正であることだけ使う)") 17 ]

def l_le_of_lemma_3_5.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l の上界の形——Lemma 3.7 / Theorem 3.8 が使う形)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
