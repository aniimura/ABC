/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Thm21Extract
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★`e` を大きく取れば次数の比は `1` に近づく（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.12。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★★★★★★★★★★★★★★★★★これは何か

原文 p.12 は、段 A（`D = ∅` への帰着）で

> moreover, by Proposition 1.7, (ii), it follows that by choosing e to be sufficiently
> large, we may assume that deg(ωX(D)|Y ) = deg(ωY (E)) ≤(1+ϵ′)·deg(ωY ).

と言う。★★★**本ファイルはこの段を取る**——`§9-976` が残した幾何の 3 点のうちの
**第 3 点（(d)）**である。

## ★機構 —— Riemann–Hurwitz を認めれば**純粋に算術**である

`n ≔ deg(Y/X)`、`A ≔ 2g_X − 2`、`Dg ≔ #D`、分岐は `D` の上だけで指数はちょうど `e`。
このとき `#E = n·Dg/e` であり、Riemann–Hurwitz は

    `deg(ω_Y) = n·A + (e − 1)·#E = n·(A + (1 − 1/e)·Dg)`

を与える。★したがって

    `deg(ω_Y(E)) = deg(ω_Y) + #E = n·(A + Dg) = deg(ω_X(D)|_Y)`

——★★**原文の等式 `deg(ω_X(D)|_Y) = deg(ω_Y(E))` はこれである**（`degOmega_of_riemannHurwitz`）。

★★★あとは `e` を大きく取れば

    `(A + Dg) / (A + (1 − 1/e)·Dg) = S / (S − Dg/e) → 1`   （`S ≔ A + Dg > 0`）

なので、任意の `ϵ′ > 0` に対し `e ≥ (1+ϵ′)·Dg/(ϵ′·S)` で
`A + Dg ≤ (1+ϵ′)(A + (1 − 1/e)Dg)` になる（`exists_ram_index_deg_ratio`）。
★同時に `A + (1 − 1/e)Dg > 0`——すなわち **`U_Y` も双曲的**である。

## ★★これで `Theorem 2.1` に残るのは 2 点

| 入力 | 状態 |
|---|---|
| (a) 分岐指数がちょうど `e` の連結有限エタール Galois 被覆 | ☆残る（folklore／[Stacks] 58.6） |
| (b) noncritical Belyi 写像（[NCBelyi] `Theorem 2.5`） | ☆残る |
| (c) 局所コンパクト性から `Ξ`・`Ξ_v` を取る段 | ★`§9-977` |
| ★(d) `deg(ω_X(D)|_Y) = deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)` | ★★**本ファイル**（Riemann–Hurwitz を受ける） |

★★★Riemann–Hurwitz 自体（`deg(ω_Y) = n·(2g_X−2) + ∑(e_P−1)`）は
[Stacks] 53.12 にあることを第 419 で実測してある。**本ファイルはそれを受ける**。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★★★★★★Riemann–Hurwitz の帳簿 -/

/-- ★★★★★★★★★★**原文の等式 `deg(ω_X(D)|_Y) = deg(ω_Y(E))` の中身**。

`n ≔ deg(Y/X)`、`A ≔ 2g_X − 2`、`Dg ≔ #D`、`#E = n·Dg/e` として:

* `deg(ω_Y) = n·A + (e−1)·#E = n·(A + (1 − 1/e)·Dg)`（Riemann–Hurwitz）
* `deg(ω_Y(E)) = deg(ω_Y) + #E = n·(A + Dg) = deg(ω_X(D)|_Y)`

★★分岐が `D` の上だけで指数がちょうど `e` のとき、**`e` は消える**
——これが原文の等式である。 -/
theorem degOmega_of_riemannHurwitz (n A Dg e : ℝ) (he : 0 < e) :
    n * A + (e - 1) * (n * Dg / e) = n * (A + (1 - 1 / e) * Dg)
  ∧ (n * A + (e - 1) * (n * Dg / e)) + n * Dg / e = n * (A + Dg) := by
  have hne : e ≠ 0 := ne_of_gt he
  constructor
  · field_simp
  · field_simp
    ring

/-! ## ★★★★★★★★★★★★★★★★★★★`e` を大きく取る -/

/-- ★★★★★★★★★★★★★★★★★★★**`e` を十分大きく取れば
`deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)` になる**（かつ `U_Y` は双曲的のまま）。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

★原文 p.12 の「`by choosing e to be sufficiently large, we may assume that
deg(ωX(D)|Y ) = deg(ωY (E)) ≤(1+ϵ′)·deg(ωY )`」がこれである。

`A ≔ 2g_X − 2`、`Dg ≔ #D`、`S ≔ A + Dg > 0`（`U_X` の双曲性）とすると

    `A + (1 − 1/e)·Dg = S − Dg/e`

なので、`e ≥ (1+ϵ′)·Dg/(ϵ′·S)` かつ `e > 2Dg/S` で両方が出る。 -/
theorem exists_ram_index_deg_ratio (A Dg ep : ℝ) (hDg : 0 ≤ Dg)
    (hhyp : 0 < A + Dg) (hep : 0 < ep) :
    ∃ e₀ : ℕ, 0 < e₀ ∧ ∀ e : ℕ, e₀ ≤ e →
      0 < A + (1 - 1 / (e : ℝ)) * Dg
      ∧ A + Dg ≤ (1 + ep) * (A + (1 - 1 / (e : ℝ)) * Dg) := by
  obtain ⟨N, hN⟩ := exists_nat_gt (max ((1 + ep) * Dg / (ep * (A + Dg))) (2 * Dg / (A + Dg)))
  refine ⟨N + 1, Nat.succ_pos N, fun e he => ?_⟩
  have heN : (N : ℝ) < (e : ℝ) := by
    have h : N < e := lt_of_lt_of_le (Nat.lt_succ_self N) he
    exact_mod_cast h
  have hmax1 : (1 + ep) * Dg / (ep * (A + Dg)) < (e : ℝ) :=
    lt_of_le_of_lt (le_max_left _ _) (lt_trans hN heN)
  have hmax2 : 2 * Dg / (A + Dg) < (e : ℝ) :=
    lt_of_le_of_lt (le_max_right _ _) (lt_trans hN heN)
  have hepos : (0:ℝ) < (e : ℝ) :=
    lt_of_le_of_lt (div_nonneg (by positivity) (le_of_lt hhyp)) hmax2
  have hdiv : Dg / (e : ℝ) ≤ (A + Dg) / 2 := by
    rw [div_le_div_iff₀ hepos (by norm_num)]
    rw [div_lt_iff₀ hhyp] at hmax2
    nlinarith
  have hrw : A + (1 - 1 / (e : ℝ)) * Dg = (A + Dg) - Dg / (e : ℝ) := by
    field_simp
    ring
  refine ⟨by rw [hrw]; linarith, ?_⟩
  rw [hrw]
  have hkey : (1 + ep) * Dg / (e : ℝ) ≤ ep * (A + Dg) := by
    rw [div_le_iff₀ hepos]
    rw [div_lt_iff₀ (by positivity)] at hmax1
    nlinarith
  have hsplit : (1 + ep) * (Dg / (e:ℝ)) = (1 + ep) * Dg / (e:ℝ) := by ring
  nlinarith [hkey]

/-- ★★★★★**(d) の結論の形**——`deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)`。

★`hRH`・`hYE` は `degOmega_of_riemannHurwitz` が与える形、
`hle` は `exists_ram_index_deg_ratio` が与える形である。 -/
theorem deg_omega_le_of_ratio (n A Dg ep e degOmY degOmYE : ℝ) (hn : 0 < n)
    (hRH : degOmY = n * (A + (1 - 1 / e) * Dg))
    (hYE : degOmYE = n * (A + Dg))
    (hpos : 0 < A + (1 - 1 / e) * Dg)
    (hle : A + Dg ≤ (1 + ep) * (A + (1 - 1 / e) * Dg)) :
    0 < degOmY ∧ degOmYE ≤ (1 + ep) * degOmY := by
  refine ⟨by rw [hRH]; positivity, ?_⟩
  rw [hRH, hYE]
  nlinarith [mul_le_mul_of_nonneg_left hle (le_of_lt hn)]

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def degOmega_of_riemannHurwitz.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1(deg(ω_X(D)|_Y) = deg(ω_Y(E))——Riemann–Hurwitz の帳簿)",
    sectionId := "genell-thm-2-1" }

def exists_ram_index_deg_ratio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1(e を十分大きく取れば次数の比は 1 に近づく)",
    sectionId := "genell-thm-2-1" }

def deg_omega_le_of_ratio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1(deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y))",
    sectionId := "genell-thm-2-1" }

def exists_ram_index_deg_ratio.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[Stacks]" "53.12 Riemann-Hurwitz(deg(ω_Y) = n·(2g_X−2) + ∑(e_P−1))" 4441,
    .otherPaper "[GenEll]" "Proposition 1.7, (ii)(原文がここで引く項目)" 9,
    .implicitStep
      ("★★★★★測定(2026-08-29): 原文 p.12 の『by choosing e to be sufficiently large, " ++
       "we may assume that deg(ω_X(D)|_Y) = deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)』は、" ++
       "**Riemann–Hurwitz を認めれば純粋に算術**である。" ++
       "★分岐が D の上だけで指数がちょうど e なら #E = n·Dg/e であり、" ++
       "deg(ω_Y) = n(A + (1−1/e)Dg)、deg(ω_Y(E)) = n(A + Dg)——**e が消える**。" ++
       "★★あとは S − Dg/e → S(S ≔ A + Dg > 0)だけである") 6,
    .implicitStep
      ("★★同時に A + (1−1/e)Dg > 0——すなわち **U_Y も双曲的**である" ++
       "(原文が『Y is a hyperbolic curve』と言う箇所)") 4,
    .implicitStep
      ("★★★これで Theorem 2.1 に残るのは**2 点**である: " ++
       "(a) 分岐指数がちょうど e の連結有限エタール Galois 被覆(folklore／[Stacks] 58.6)、" ++
       "(b) noncritical Belyi 写像([NCBelyi] Theorem 2.5)。" ++
       "★★★算術・解析の段(§9-976 の鎖、§9-977 のコンパクト性、本ファイルの次数)は" ++
       "**すべて取れた**") 7 ]

end ABC3.Found.GenEll
