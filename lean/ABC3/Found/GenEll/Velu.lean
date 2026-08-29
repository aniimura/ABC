/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.EllipticCurve.Weierstrass
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★Vélu の公式——`E/H` の Weierstrass モデル（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か —— `Lemma 3.5` の葉

`Lemma 3.5` は「`H_F ⊆ E_F` を `l`-巡回部分群スキームとし、`(E_H)_F ≝ E_F/H_F` と書く」
という設定から始まる。★本ファイルはその **`E/H` を Weierstrass 曲線として作る**
——Vélu の公式である。

★★mathlib に楕円曲線の**同種写像も商も無い**（2026-08-29 に `#check` で確認:
`WeierstrassCurve.Isogeny`・`WeierstrassCurve.velu`・`WeierstrassCurve.quotient` いずれも不在）。

## ★なぜこれが要るのか —— アルキメデス項が消えるから

`Found/GaloisRep/IsogenyReduction.lean`（`§9-1024`・`§9-1026`）が示すとおり、
`Lemma 3.5` に残る唯一の入力 `ht^Falt(E/H) ≤ ht^Falt(E) + 2·log(l)` は

    `(l−1)·d·deg∞(E) − (archSum(E′) − archSum(E)) ≤ 24·d·log(l)`

に帰着する。★★Vélu の公式は **`φ^*(ω_{E′}) = ω_E`**（不変微分が引き戻しで一致）と
正規化されているので、ℂ 上では周期格子がそのまま `Λ ⊆ Λ′` になり、
**アルキメデス項が消えて純粋に有限素点の主張に落ちる**（`§9-1027`、第 585 の測定）。

## ★★★本ファイルが取るもの / 取らないもの

★**取るもの**（すべて多項式の恒等式——`ring` で閉じる）:

| | 内容 |
|---|---|
| `veluGx`・`veluGy` | 点 `Q` での 1 階・2 階のデータ |
| `veluV`・`veluU`・`veluW` | Vélu の `v_Q`・`u_Q`・`w_Q` |
| `veluCurve W v w` | `Y² + a₁XY + a₃Y = X³ + a₂X² + (a₄−5v)X + (a₆−b₂v−7w)` |
| `veluCurve_b₂`〜`c₆` | 不変量の変化（`c₄` は `+240v`、`c₆` は `+504b₂v+6048w`） |
| `veluQuotient W S` | 代表点の有限集合 `S` にわたる和で作った商 |

☆**取らないもの**: 「`veluQuotient W S` が本当に `E/H` である」こと
——すなわち次数 `#H` の同種写像 `φ : E → veluQuotient W S` の構成と
`φ^*(ω′) = ω` の証明。★これには群法則・`S` が `(H∖{O})/±` の代表系であること・
Galois 安定性が要る。☆**本ファイルは定義と不変量の帳簿だけである。**

★★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve

variable {R : Type*} [CommRing R]

/-! ## ★★★★★点ごとの量 -/

/-- ★★★★★Vélu の `g^x_Q = 3x² + 2a₂x + a₄ − a₁y`。 -/
def veluGx (W : WeierstrassCurve R) (x y : R) : R :=
  3 * x ^ 2 + 2 * W.a₂ * x + W.a₄ - W.a₁ * y

/-- ★★★★★Vélu の `g^y_Q = −2y − a₁x − a₃`。 -/
def veluGy (W : WeierstrassCurve R) (x y : R) : R := -2 * y - W.a₁ * x - W.a₃

/-- ★★★★★Vélu の `v_Q = 2g^x_Q − a₁·g^y_Q`（`Q` が 2-捩れでないとき）。 -/
def veluV (W : WeierstrassCurve R) (x y : R) : R :=
  2 * veluGx W x y - W.a₁ * veluGy W x y

/-- ★★★★Vélu の `v_Q = g^x_Q`（`Q` が 2-捩れのとき）。

☆原文の場合分けの片側。本ファイルは**両方を定義するだけ**で、
どちらを使うかは `S` の取り方（代表系）に委ねる。 -/
def veluV2 (W : WeierstrassCurve R) (x y : R) : R := veluGx W x y

/-- ★★★★★Vélu の `u_Q = (g^y_Q)²`。 -/
def veluU (W : WeierstrassCurve R) (x y : R) : R := (veluGy W x y) ^ 2

/-- ★★★★★Vélu の `w_Q = u_Q + v_Q·x_Q`。 -/
def veluW (W : WeierstrassCurve R) (x y : R) : R := veluU W x y + veluV W x y * x

/-! ## ★★★★★★★★商曲線 -/

/-- ★★★★★★★★**Vélu の商曲線**

    `Y² + a₁XY + a₃Y = X³ + a₂X² + (a₄ − 5v)X + (a₆ − b₂v − 7w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`v`・`w` は代表点にわたる `v_Q`・`w_Q` の和である。 -/
def veluCurve (W : WeierstrassCurve R) (v w : R) : WeierstrassCurve R where
  a₁ := W.a₁
  a₂ := W.a₂
  a₃ := W.a₃
  a₄ := W.a₄ - 5 * v
  a₆ := W.a₆ - W.b₂ * v - 7 * w

@[simp] theorem veluCurve_a₁ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).a₁ = W.a₁ := rfl
@[simp] theorem veluCurve_a₂ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).a₂ = W.a₂ := rfl
@[simp] theorem veluCurve_a₃ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).a₃ = W.a₃ := rfl
@[simp] theorem veluCurve_a₄ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).a₄ = W.a₄ - 5 * v := rfl
@[simp] theorem veluCurve_a₆ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).a₆ = W.a₆ - W.b₂ * v - 7 * w := rfl

/-- ★★★`v = w = 0` なら曲線は変わらない（`H = {O}` の場合）。 -/
theorem veluCurve_zero (W : WeierstrassCurve R) : veluCurve W 0 0 = W := by
  cases W; simp [veluCurve]

/-! ## ★★★★★★不変量の帳簿 -/

/-- ★★★★★★**`b₂` は変わらない**。 -/
theorem veluCurve_b₂ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).b₂ = W.b₂ := by
  simp [WeierstrassCurve.b₂]

/-- ★★★★★**`b₄ ↦ b₄ − 10v`**。 -/
theorem veluCurve_b₄ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).b₄ = W.b₄ - 10 * v := by
  simp [WeierstrassCurve.b₄]; ring

/-- ★★★★★**`b₆ ↦ b₆ − 4b₂v − 28w`**。 -/
theorem veluCurve_b₆ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).b₆ = W.b₆ - 4 * W.b₂ * v - 28 * w := by
  simp [WeierstrassCurve.b₆, veluCurve]; ring

/-- ★★★★★**`b₈ ↦ b₈ − b₂²v − 7b₂w + 5b₄v − 25v²`**。 -/
theorem veluCurve_b₈ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).b₈
      = W.b₈ - W.b₂ ^ 2 * v - 7 * W.b₂ * w + 5 * W.b₄ * v - 25 * v ^ 2 := by
  simp [WeierstrassCurve.b₈, WeierstrassCurve.b₂, WeierstrassCurve.b₄, veluCurve]; ring

/-- ★★★★★★★**`c₄ ↦ c₄ + 240v`**——Vélu の古典的な形。 -/
theorem veluCurve_c₄ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).c₄ = W.c₄ + 240 * v := by
  simp [WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄, veluCurve]; ring

/-- ★★★★★★★**`c₆ ↦ c₆ + 504b₂v + 6048w`**。 -/
theorem veluCurve_c₆ (W : WeierstrassCurve R) (v w : R) :
    (veluCurve W v w).c₆ = W.c₆ + 504 * W.b₂ * v + 6048 * w := by
  simp [WeierstrassCurve.c₆, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    veluCurve]
  ring

/-! ## ★★★★★代表点の和 -/

/-- ★★★★★代表点の集合 `S` にわたる `v = Σ v_Q`。 -/
noncomputable def veluVSum (W : WeierstrassCurve R) (S : Finset (R × R)) : R :=
  ∑ Q ∈ S, veluV W Q.1 Q.2

/-- ★★★★★代表点の集合 `S` にわたる `w = Σ w_Q`。 -/
noncomputable def veluWSum (W : WeierstrassCurve R) (S : Finset (R × R)) : R :=
  ∑ Q ∈ S, veluW W Q.1 Q.2

/-- ★★★★★★**Vélu の商** `E/H`——`S` は `(H∖{O})/±` の代表系のつもりである。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆「これが本当に `E/H` である」ことは**本ファイルでは示していない**
——同種写像 `φ : E → veluQuotient W S` の構成が残る。 -/
noncomputable def veluQuotient (W : WeierstrassCurve R) (S : Finset (R × R)) :
    WeierstrassCurve R :=
  veluCurve W (veluVSum W S) (veluWSum W S)

/-- ★★★`S = ∅`（`H = {O}`）なら曲線は変わらない。 -/
theorem veluQuotient_empty (W : WeierstrassCurve R) : veluQuotient W ∅ = W := by
  rw [veluQuotient, veluVSum, veluWSum, Finset.sum_empty, Finset.sum_empty, veluCurve_zero]

/-! ## ★出典の紐付け(`.src`)——★★**条つき。定義と帳簿だけである** -/

def veluCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(E/H の Weierstrass モデル——Vélu の公式の定義)",
    sectionId := "genell-lemma-3-5" }

def veluCurve_c₄.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の不変量——c₄ ↦ c₄ + 240v)",
    sectionId := "genell-lemma-3-5" }

def veluQuotient.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商 E/H——代表点の和で作る。同種写像の構成は含まない)",
    sectionId := "genell-lemma-3-5" }

def veluQuotient.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore
      ("☆★★★★★★★★**残っている中核**: 「`veluQuotient W S` が本当に `E/H` である」" ++
       "こと——次数 `#H` の同種写像 `φ : E → veluQuotient W S` の構成と " ++
       "`φ^*(ω′) = ω`(不変微分が引き戻しで一致すること)の証明。" ++
       "★これには群法則(mathlib の WeierstrassCurve.Affine.Point)・" ++
       "`S` が `(H∖{O})/±` の代表系であること・Galois 安定性が要る。" ++
       "★★本ファイルは**定義と不変量の帳簿だけ**である") 9,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 586): mathlib に楕円曲線の商が無い" ++
       "(#check で確認: WeierstrassCurve.velu・quotient・Isogeny いずれも不在)ので、" ++
       "Vélu の公式を本プロジェクトに置いた。" ++
       "★不変量の変化はすべて多項式の恒等式で、ring で閉じる: " ++
       "b₂ は不変、b₄ ↦ b₄−10v、b₆ ↦ b₆−4b₂v−28w、" ++
       "b₈ ↦ b₈−b₂²v−7b₂w+5b₄v−25v²、c₄ ↦ c₄+240v、c₆ ↦ c₆+504b₂v+6048w") 8,
    .implicitStep
      ("★なぜこれが Lemma 3.5 の葉なのか: §9-1024・§9-1026 が示すとおり" ++
       "残る入力は (l−1)d·deg∞(E) − (archSum(E′) − archSum(E)) ≤ 24d·log(l) であり、" ++
       "★★Vélu の正規化 φ^*ω′ = ω により ℂ 上で周期格子が Λ ⊆ Λ′ になるので" ++
       "**アルキメデス項が消えて純粋に有限素点の主張に落ちる**" ++
       "(§9-1027、第 585 の測定)") 9,
    .implicitStep
      ("☆2-捩れ点の場合分け(veluV2 = veluGx)は定義しただけで、" ++
       "どちらを使うかは S の取り方(代表系)に委ねている。" ++
       "★同種写像を構成するときに確定する") 7 ]

end ABC3.Found.GenEll
