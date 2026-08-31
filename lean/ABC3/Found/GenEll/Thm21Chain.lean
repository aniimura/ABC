/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Theorem 2.1` の 2 本の鎖（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.11–p.13。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 原文の計算そのもの

原文 p.12–p.13 の `(ii) ⟹ (i)` は**2 段の帰着**でできている。
★どちらの段も、中身は**不等式の鎖の計算**である。本ファイルはその 2 本を取る。

### ★段 A（`D = ∅` への帰着、原文 p.12）

被覆 `U_Y → U_X`（`E ≔ (D ×_X Y)_red` の各点で分岐指数がちょうど `e`）を取り、
`d′ ≔ d·deg(Y/X)` と置く。原文の鎖は

    `ht_{ω_X(D)} ≲ (1+ϵ′)·ht_{ω_Y} ≲ (1+ϵ′)²·log-diff_Y`
                 `≲ (1+ϵ′)²(log-diff_X + log-cond_D) ≲ (1+ϵ)(log-diff_X + log-cond_D)`

（最後は `(1+ϵ′)² ≤ 1+ϵ`）。★★`thm_2_1_stepA` がこれである。

### ★段 B（`D = ∅` の場合の背理法、原文 p.12–p.13）

noncritical Belyi 写像 `φ : X → P` を取り、`E ≔ φ^{-1}(C)_red`
（したがって `φ^*ω_P(C) ≅ ω_X(E)`）と置く。原文の鎖は

    `ht_{ω_X} ≈ ht_{ω_X(E)} − ht_E ≈ ht_{ω_P(C)} − ht_E`
             `≲ (1+ϵ′)(log-diff_P + log-cond_C) − ht_E`
             `≲ (1+ϵ′)(log-diff_X + log-cond_E) − ht_E`
             `≲ (1+ϵ′)(log-diff_X + ht_E) − ht_E ≈ (1+ϵ′)·log-diff_X + ϵ′·ht_E`
             `≈ (1+ϵ′)·log-diff_X + ϵ′·(deg(E)/deg(ω_X))·ht_{ω_X}`

——すなわち `(1 − ϵ′·r)·ht_{ω_X} ≲ (1+ϵ′)·log-diff_X`（`r ≔ deg(E)/deg(ω_X)`）。
`1 + ϵ′ ≤ (1+ϵ)(1 − ϵ′·r)` を満たす `ϵ′` を取れば `ht_{ω_X} ≲ (1+ϵ)·log-diff_X` である。
★★`thm_2_1_stepB` がこれである。

## ★★★★★★★★★★★★★★何が残るか（測定、2026-08-29）

★★**本ファイルは `Theorem 2.1` を閉じない。** 残るのは**幾何の入力**である:

| 入力 | 原文での出どころ | 状態 |
|---|---|---|
| `Prop 1.7, (i)`（`log-diff_P + log-cond_C ≲ log-diff_X + log-cond_E`） | [GenEll] | ★★**本日 `§9-975` で取れた** |
| `Prop 1.6`（`log-cond_E ≲ ht_E`） | [GenEll] | ★**手元にある**（`Skeleton` の `prop_1_6`） |
| `Prop 1.4, (i)(iii)`（高さの加法性・`BD`-類が同型類だけに依ること） | [GenEll] | ★**手元にある**（`§9-962`） |
| `E ≔ (D ×_X Y)_red` の各点で分岐指数 `e` の連結有限エタール Galois 被覆 | folklore／[Stacks] 58.6 | ☆**無い** |
| noncritical Belyi 写像 `φ : X → P` | [NCBelyi] `Theorem 2.5` | ☆**無い**（`Skeleton/NCBelyi/Theorem25.lean` に見出しだけ） |
| 局所コンパクト性から `Ξ`・`Ξ_v` を取る段 | 原文 p.12 | ☆**無い** |
| `deg(ω_X(D)|_Y) = deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)` | `Prop 1.7, (ii)`（Riemann–Hurwitz） | ☆**無い**（`deg` が語彙の外） |

★★★★★したがって `Theorem 2.1` に残るのは**幾何の 4 点**であり、
**算術の鎖は本ファイルで尽きている**。★`.src` は条つき（`Theorem 2.1(…)`）にしてある
——指標には**数えない**。

## ★向きについて

★原文の `≲` を印字どおり（`BDle`）に読むと鎖は逆を向き、`abc` の向きにならない。
★★`Check/GenEll/Prop17Direction.lean` で**印字どおりの読みが `Proposition 1.7` を偽にする**
ことを機械検証した（2026-08-29）。本ファイルは**通常の読み**（`BDge`、`α ≤ β + C`）で取る。
`Gap/GenEll/BDDirection.lean` を参照。
-/

namespace ABC3.Found.GenEll

/-! ## ★`BD` の道具 -/

/-- ★`BDge` の推移律。★鎖を繋ぐのに要る。 -/
theorem bdge_trans {Pt : Type} {α β γ : Pt → ℝ} (h1 : BDge α β) (h2 : BDge β γ) :
    BDge α γ := by
  obtain ⟨C1, hC1⟩ := h1
  obtain ⟨C2, hC2⟩ := h2
  exact ⟨C1 + C2, fun x => by have := hC1 x; have := hC2 x; linarith⟩

/-- ★★**`(i) ⟹ (ii)` の中身** —— `BD`-類の不等式は部分集合へ制限できる。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

★原文 p.12 は「`The fact that (i) ⟹ (ii) is immediate from the definitions.`」と書く。
★原文が「定義から直ちに」と言うのはこれである——**同じ定数 `C` がそのまま使える**。 -/
theorem bdge_restrict {Pt : Type} (α β : Pt → ℝ) (T : Set Pt)
    (h : BDge α β) :
    BDge (fun x : ↥T => α x.1) (fun x : ↥T => β x.1) := by
  obtain ⟨C, hC⟩ := h
  exact ⟨C, fun x => hC x.1⟩

/-! ## ★★★★★★★★★★★★★★段 A —— `D = ∅` への帰着（原文 p.12） -/

/-- ★★★★★★★★★★★★★★**原文 p.12 の鎖**（`D = ∅` への帰着）。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

    `ht_{ω_X(D)} ≲ (1+ϵ′)·ht_{ω_Y} ≲ (1+ϵ′)²·log-diff_Y`
                 `≲ (1+ϵ′)²(log-diff_X + log-cond_D) ≲ (1+ϵ)(log-diff_X + log-cond_D)`

★`g1` は `Prop 1.4, (i)(ii)(iii)` と次数の比較（`Prop 1.7, (ii)`）から、
`g2` は `U_Y` の上で仮定した不等式、`g3` は `Prop 1.7, (i)` から来る。
★★**本補題はその 3 本を繋ぐ計算だけ**である。 -/
theorem thm_2_1_stepA {Pt : Type}
    (htOmXD htOmY logDiffY logDiffX logCondD : Pt → ℝ)
    (eps ep : ℝ) (hep : 0 ≤ ep)
    (hchoice : (1 + ep) ^ 2 ≤ 1 + eps)
    (hXnn : ∀ x, 0 ≤ logDiffX x + logCondD x)
    (g1 : BDge htOmXD (fun x => (1 + ep) * htOmY x))
    (g2 : BDge htOmY (fun x => (1 + ep) * logDiffY x))
    (g3 : BDge logDiffY (fun x => logDiffX x + logCondD x)) :
    BDge htOmXD (fun x => (1 + eps) * (logDiffX x + logCondD x)) := by
  obtain ⟨C1, hC1⟩ := g1
  obtain ⟨C2, hC2⟩ := g2
  obtain ⟨C3, hC3⟩ := g3
  refine ⟨C1 + (1 + ep) * C2 + (1 + ep) ^ 2 * C3, fun x => ?_⟩
  have a1 := hC1 x
  have a2 := hC2 x
  have a3 := hC3 x
  simp only at a1 a2 a3
  have hep1 : (0:ℝ) ≤ 1 + ep := by linarith
  have hstep1 : (1 + ep) * htOmY x ≤ (1 + ep) * ((1 + ep) * logDiffY x + C2) :=
    mul_le_mul_of_nonneg_left (by linarith) hep1
  have hstep2 : (1 + ep) ^ 2 * logDiffY x
      ≤ (1 + ep) ^ 2 * ((logDiffX x + logCondD x) + C3) :=
    mul_le_mul_of_nonneg_left (by linarith) (by positivity)
  have hstep3 : (1 + ep) ^ 2 * (logDiffX x + logCondD x)
      ≤ (1 + eps) * (logDiffX x + logCondD x) :=
    mul_le_mul_of_nonneg_right hchoice (hXnn x)
  nlinarith [hstep1, hstep2, hstep3]

/-! ## ★★★★★★★★★★★★★★★★★★段 B —— `D = ∅` の場合（原文 p.12–p.13） -/

/-- ★★★★★★★★★★★★★★★★★★**原文 p.13 の鎖**（`D = ∅` の場合）。

原文 (GenEll p.13):
> htωX ≈htωX(E) −htE ≈htωP (C) −htE

    `ht_{ω_X} ≈ ht_{ω_X(E)} − ht_E ≈ ht_{ω_P(C)} − ht_E`
             `≲ (1+ϵ′)(log-diff_P + log-cond_C) − ht_E`
             `≲ (1+ϵ′)(log-diff_X + log-cond_E) − ht_E`
             `≲ (1+ϵ′)(log-diff_X + ht_E) − ht_E`
             `≈ (1+ϵ′)·log-diff_X + ϵ′·(deg(E)/deg(ω_X))·ht_{ω_X}`

すなわち `(1 − ϵ′·r)·ht_{ω_X} ≲ (1+ϵ′)·log-diff_X`。
★原文の `ϵ′` の選び方（`1 + ϵ′ ≤ (1+ϵ)(1 − ϵ′·deg(E)/deg(ω_X))`）が `hchoice` である。

★仮定の出どころ:

| 仮定 | 出どころ |
|---|---|
| `h1` | `Prop 1.4, (i)`（高さの加法性） |
| `h2` | `φ^*ω_P(C) ≅ ω_X(E)` と `Prop 1.4, (iii)` |
| `h3` | ★**(ii) そのもの**（`φ(Ξ) ⊆ K_V` の上で） |
| `h4` | `Prop 1.7, (i)`（`e = 1` として）——★`§9-975` |
| `h5` | `Prop 1.6`（`log-cond_E ≲ ht_E`） |
| `h6` | `ht_E ≈ (deg(E)/deg(ω_X))·ht_{ω_X}`（`Prop 1.4, (i)`） | -/
theorem thm_2_1_stepB {Pt : Type}
    (htOm htOmXE htOmPC htE logDiffX logCondE logDiffP logCondC : Pt → ℝ)
    (eps ep r : ℝ) (hep : 0 ≤ ep) (hlt : ep * r < 1)
    (hchoice : 1 + ep ≤ (1 + eps) * (1 - ep * r))
    (hdiffnn : ∀ x, 0 ≤ logDiffX x)
    (h1 : BDeq htOm (fun x => htOmXE x - htE x))
    (h2 : BDeq htOmXE htOmPC)
    (h3 : BDge htOmPC (fun x => (1 + ep) * (logDiffP x + logCondC x)))
    (h4 : BDge (fun x => logDiffP x + logCondC x) (fun x => logDiffX x + logCondE x))
    (h5 : BDge logCondE htE)
    (h6 : BDge htE (fun x => r * htOm x)) :
    BDge htOm (fun x => (1 + eps) * logDiffX x) := by
  obtain ⟨C1, hC1⟩ := h1
  obtain ⟨C2, hC2⟩ := h2
  obtain ⟨C3, hC3⟩ := h3
  obtain ⟨C4, hC4⟩ := h4
  obtain ⟨C5, hC5⟩ := h5
  obtain ⟨C6, hC6⟩ := h6
  have hd : (0:ℝ) < 1 - ep * r := by linarith
  have hep1 : (0:ℝ) ≤ 1 + ep := by linarith
  refine ⟨((1 + ep) * (C4 + C5) + C1 + C2 + C3 + ep * C6) / (1 - ep * r), fun x => ?_⟩
  have a1 : htOm x ≤ htOmXE x - htE x + C1 := by linarith [(abs_le.mp (hC1 x)).2]
  have a2 : htOmXE x ≤ htOmPC x + C2 := by linarith [(abs_le.mp (hC2 x)).2]
  have a3 := hC3 x
  have a4 := hC4 x
  have a5 := hC5 x
  have a6 := hC6 x
  simp only at a3 a4 a6
  have s1 : (1 + ep) * (logDiffP x + logCondC x)
      ≤ (1 + ep) * ((logDiffX x + logCondE x) + C4) :=
    mul_le_mul_of_nonneg_left (by linarith) hep1
  have s2 : (1 + ep) * ((logDiffX x + logCondE x) + C4)
      ≤ (1 + ep) * ((logDiffX x + htE x + C5) + C4) :=
    mul_le_mul_of_nonneg_left (by linarith) hep1
  have s3 : ep * htE x ≤ ep * (r * htOm x + C6) :=
    mul_le_mul_of_nonneg_left (by linarith) hep
  have key : (1 - ep * r) * htOm x
      ≤ (1 + ep) * logDiffX x + ((1 + ep) * (C4 + C5) + C1 + C2 + C3 + ep * C6) := by
    nlinarith [s1, s2, s3]
  have s4 : (1 + ep) * logDiffX x ≤ ((1 + eps) * (1 - ep * r)) * logDiffX x :=
    mul_le_mul_of_nonneg_right hchoice (hdiffnn x)
  show htOm x - (1 + eps) * logDiffX x
      ≤ ((1 + ep) * (C4 + C5) + C1 + C2 + C3 + ep * C6) / (1 - ep * r)
  rw [le_div_iff₀ hd]
  nlinarith [key, s4]

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def bdge_restrict.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1((i) ⟹ (ii)——BD-類の不等式は部分集合へ制限できる)",
    sectionId := "genell-thm-2-1" }

def thm_2_1_stepA.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 12,
    item := "Theorem 2.1(段 A——D = ∅ への帰着の鎖)",
    sectionId := "genell-thm-2-1" }

def thm_2_1_stepB.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 13,
    item := "Theorem 2.1(段 B——D = ∅ の場合の鎖)",
    sectionId := "genell-thm-2-1" }

def thm_2_1_stepB.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "prop_1_7((i) の両側、§9-975)——h4 の出どころ"
      (.inProject "ABC3" "ABC3.Found.GenEll.prop_1_7") 4,
    .otherPaper "[GenEll]" "Proposition 1.6(log-cond_E ≲ ht_E)——h5 の出どころ" 8,
    .otherPaper "[GenEll]" "Proposition 1.4, (i)(iii)(高さの加法性・同型類だけに依ること)" 6,
    .otherPaper "[NCBelyi]"
      ("Theorem 2.5(Belyi Maps Noncritical at Prescribed Points)——★★段 B の幾何の核。" ++
       "★Skeleton/NCBelyi/Theorem25.lean に見出しだけがある。**取っていない**") 5,
    .folklore
      ("任意の正整数 e に対し、E ≙ (D ×_X Y)_red の各点で分岐指数がちょうど e となる" ++
       "連結有限エタール Galois 被覆 U_Y → U_X が存在する(段 A の幾何の核)。" ++
       "★原典側は [Stacks] 58.6 Fundamental groups。**取っていない**") 11,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 本ファイルは Theorem 2.1 を**閉じない**。" ++
       "★原文 p.12-p.13 の (ii) ⟹ (i) は 2 段の帰着であり、" ++
       "**どちらも中身は不等式の鎖の計算**である——それを本ファイルが取った。" ++
       "★★残るのは幾何の 4 点である: (a) 分岐指数 e の連結有限エタール Galois 被覆、" ++
       "(b) noncritical Belyi 写像([NCBelyi] Theorem 2.5)、" ++
       "(c) 局所コンパクト性から Ξ・Ξ_v を取る段、" ++
       "(d) deg(ω_X(D)|_Y) = deg(ω_Y(E)) ≤ (1+ϵ′)·deg(ω_Y)(Prop 1.7, (ii)——deg が語彙の外)。" ++
       "★★★算術の鎖は本ファイルで尽きている") 9,
    .implicitStep
      ("★★★★向きについて: 原文の ≲ を印字どおり(BDle)に読むと鎖は逆を向き abc の向きにならない。" ++
       "★Check/GenEll/Prop17Direction.lean が『印字どおりの読みは Proposition 1.7 を偽にする』" ++
       "ことを機械検証した(2026-08-29)。本ファイルは**通常の読み(BDge)**で取る") 6 ]

end ABC3.Found.GenEll
