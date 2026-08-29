/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.AlgebraicGeometry.EllipticCurve.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.VariableChange
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Formula
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

/-! ## ★★★★★★★★★★モデルの取り替えとの両立 -/

/-- ★スケーリングだけの変数変換 `(u, 0, 0, 0)`。 -/
def scaleChange (u : Rˣ) : VariableChange R := ⟨u, 0, 0, 0⟩

/-- ★★★★★**`g^x_Q` は重さ 4**。 -/
theorem veluGx_scale (W : WeierstrassCurve R) (x y : R) (u : Rˣ) :
    veluGx (scaleChange u • W) (((u⁻¹ : Rˣ) : R) ^ 2 * x) (((u⁻¹ : Rˣ) : R) ^ 3 * y)
      = ((u⁻¹ : Rˣ) : R) ^ 4 * veluGx W x y := by
  simp [veluGx, scaleChange, WeierstrassCurve.variableChange_a₁,
    WeierstrassCurve.variableChange_a₂, WeierstrassCurve.variableChange_a₄]
  ring

/-- ★★★★★**`g^y_Q` は重さ 3**。 -/
theorem veluGy_scale (W : WeierstrassCurve R) (x y : R) (u : Rˣ) :
    veluGy (scaleChange u • W) (((u⁻¹ : Rˣ) : R) ^ 2 * x) (((u⁻¹ : Rˣ) : R) ^ 3 * y)
      = ((u⁻¹ : Rˣ) : R) ^ 3 * veluGy W x y := by
  simp [veluGy, scaleChange, WeierstrassCurve.variableChange_a₁,
    WeierstrassCurve.variableChange_a₃]
  ring

/-- ★★★★★★**`v_Q` は重さ 4**——`a₄` と同じ。 -/
theorem veluV_scale (W : WeierstrassCurve R) (x y : R) (u : Rˣ) :
    veluV (scaleChange u • W) (((u⁻¹ : Rˣ) : R) ^ 2 * x) (((u⁻¹ : Rˣ) : R) ^ 3 * y)
      = ((u⁻¹ : Rˣ) : R) ^ 4 * veluV W x y := by
  rw [veluV, veluGx_scale, veluGy_scale, veluV,
    WeierstrassCurve.variableChange_a₁, scaleChange]
  ring

/-- ★★★★★★**`w_Q` は重さ 6**——`a₆` と同じ。 -/
theorem veluW_scale (W : WeierstrassCurve R) (x y : R) (u : Rˣ) :
    veluW (scaleChange u • W) (((u⁻¹ : Rˣ) : R) ^ 2 * x) (((u⁻¹ : Rˣ) : R) ^ 3 * y)
      = ((u⁻¹ : Rˣ) : R) ^ 6 * veluW W x y := by
  rw [veluW, veluU, veluGy_scale, veluV_scale, veluW, veluU]
  ring

/-- ★★★★★★★★★★**Vélu の商はスケーリングと両立する**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `veluCurve (u • W) (u⁻⁴v) (u⁻⁶w) = u • veluCurve W v w`

★★これが `neronExp` の帳簿に効く——`Found/GaloisRep/IsogenyReduction.lean` の
`hArch` は与えられたモデルに依るので、**モデルを取り替えたときの挙動**が要る。
★`veluV_scale`・`veluW_scale` と重さが揃っている（4 と 6）のが要点である。 -/
theorem veluCurve_scale (W : WeierstrassCurve R) (v w : R) (u : Rˣ) :
    veluCurve (scaleChange u • W) (((u⁻¹ : Rˣ) : R) ^ 4 * v) (((u⁻¹ : Rˣ) : R) ^ 6 * w)
      = scaleChange u • veluCurve W v w := by
  ext
  · simp [veluCurve, scaleChange, WeierstrassCurve.variableChange_a₁]
  · simp [veluCurve, scaleChange, WeierstrassCurve.variableChange_a₂]
  · simp [veluCurve, scaleChange, WeierstrassCurve.variableChange_a₃]
  · simp [veluCurve, scaleChange, WeierstrassCurve.variableChange_a₄]
    ring
  · simp [veluCurve, scaleChange, WeierstrassCurve.variableChange_a₆,
      WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₂,
      WeierstrassCurve.b₂]
    ring

/-! ## ★★★★★★★★★★★★★★★★★★`l = 2` の場合——写像が本当に商へ落ちる -/

section TwoIsogeny

variable {F : Type*} [Field F]

/-- ★★★★★**2-同種写像の Vélu 写像（`X` 座標）**。

★2-捩れ点 `Q = (x₀,y₀)` では `g^y_Q = 0`（したがって `u_Q = 0`）なので、
`X = x + v_Q/(x − x₀)` と簡単になる。 -/
noncomputable def velu2X (W : WeierstrassCurve F) (x₀ y₀ x : F) : F :=
  x + veluV2 W x₀ y₀ / (x - x₀)

/-- ★★★★★**2-同種写像の Vélu 写像（`Y` 座標）**。 -/
noncomputable def velu2Y (W : WeierstrassCurve F) (x₀ y₀ x y : F) : F :=
  y - veluV2 W x₀ y₀ * (W.a₁ * (x - x₀) + y - y₀) / (x - x₀) ^ 2

set_option maxRecDepth 40000 in
set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★★★★★**Vélu の 2-同種写像は本当に商へ落ちる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`Q = (x₀,y₀)` を `E` 上の 2-捩れ点、`P = (x,y)` を `x ≠ x₀` なる点とすると、
**`(velu2X, velu2Y)` は商曲線 `veluCurve W v_Q (v_Q·x₀)` 上にある**。

★★★これが「Vélu の構成が本当に `E/H` を与える」ことの、`#H = 2` での実証である。
★恒等式は

    `Quot-eq(X,Y)·(x−x₀)⁴ = (v_Q − (x−x₀)²)²·E-eq(x,y)`

——★★余因子が `(v_Q − (x−x₀)²)²` と**平方**になるのが要点である
（sympy で係数を求めてから `linear_combination` に渡した）。

☆一般の `l` については、代表系の取り方と群法則が要る。 -/
theorem velu2_equation (W : WeierstrassCurve F) (x₀ y₀ x y : F)
    (h0 : W.toAffine.Equation x₀ y₀) (h2 : 2 * y₀ + W.a₁ * x₀ + W.a₃ = 0)
    (h : W.toAffine.Equation x y) (hx : x ≠ x₀) :
    (veluCurve W (veluV2 W x₀ y₀) (veluV2 W x₀ y₀ * x₀)).toAffine.Equation
      (velu2X W x₀ y₀ x) (velu2Y W x₀ y₀ x y) := by
  rw [WeierstrassCurve.Affine.equation_iff] at h0 h ⊢
  have hne : x - x₀ ≠ 0 := sub_ne_zero.2 hx
  have ha3 : W.a₃ = -(2 * y₀ + W.a₁ * x₀) := by linear_combination h2
  have ha6 : W.a₆ = y₀ ^ 2 + W.a₁ * x₀ * y₀ + W.a₃ * y₀
      - x₀ ^ 3 - W.a₂ * x₀ ^ 2 - W.a₄ * x₀ := by linear_combination -h0
  simp only [velu2X, velu2Y, veluCurve, veluV2, veluGx, WeierstrassCurve.b₂]
  rw [ha6, ha3] at h ⊢
  field_simp
  linear_combination ((3 * x₀ ^ 2 + 2 * W.a₂ * x₀ + W.a₄ - W.a₁ * y₀) - (x - x₀) ^ 2) ^ 2 * h

/-- ★★★★★★★★★★★★★★★★★★★★**`l = 2`: `φ^*(ω′) = ω`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

不変微分は `ω = dx/(2y + a₁x + a₃)` なので、`φ^*(ω′) = ω` は

    `(dX/dx)·(2y + a₁x + a₃) = 2Y + a₁X + a₃`

と同じことである（`dX/dx = 1 − v_Q/(x−x₀)²`）。

★★★**これが Vélu の正規化**であり、`§9-1027`（第 585）の見立て
「アルキメデス項が消える」の中身である
——ℂ 上で `ω = dz` なら周期格子がそのまま `Λ ⊆ Λ′` になる。
★2-捩れ条件がちょうど効く（残差は `(2y₀+a₁x₀+a₃)·(−v_Q)`）。 -/
theorem velu2_omega (W : WeierstrassCurve F) (x₀ y₀ x y : F)
    (h2 : 2 * y₀ + W.a₁ * x₀ + W.a₃ = 0) (hx : x ≠ x₀) :
    (1 - veluV2 W x₀ y₀ / (x - x₀) ^ 2) * (2 * y + W.a₁ * x + W.a₃)
      = 2 * velu2Y W x₀ y₀ x y + W.a₁ * velu2X W x₀ y₀ x + W.a₃ := by
  have hne : x - x₀ ≠ 0 := sub_ne_zero.2 hx
  simp only [velu2X, velu2Y, veluV2, veluGx]
  field_simp
  linear_combination (-(3 * x₀ ^ 2 + 2 * W.a₂ * x₀ + W.a₄ - W.a₁ * y₀)) * h2

end TwoIsogeny

/-! ## ★★★★★★★★★★★★★★★★★★★★`l = 3` の場合 -/

section ThreeIsogeny

variable {F : Type*} [Field F]

/-- ★★★★★**3-同種写像の Vélu 写像（`X` 座標）**。

★位数 3 の点 `Q` では `H = {O, Q, −Q}` で代表系は `S = {Q}` である。 -/
noncomputable def velu3X (W : WeierstrassCurve F) (xQ yQ x : F) : F :=
  x + veluV W xQ yQ / (x - xQ) + veluU W xQ yQ / (x - xQ) ^ 2

/-- ★★★★★**3-同種写像の Vélu 写像（`Y` 座標）**。 -/
noncomputable def velu3Y (W : WeierstrassCurve F) (xQ yQ x y : F) : F :=
  y - veluU W xQ yQ * (2 * y + W.a₁ * x + W.a₃) / (x - xQ) ^ 3
    - veluV W xQ yQ * (W.a₁ * (x - xQ) + y - yQ) / (x - xQ) ^ 2
    - (W.a₁ * veluU W xQ yQ - veluGx W xQ yQ * veluGy W xQ yQ) / (x - xQ) ^ 2

set_option maxRecDepth 100000 in
set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★**Vélu の 3-同種写像は本当に商へ落ちる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`Q = (xQ,yQ)` を**位数 3** の点（3-分点多項式
`ψ₃ = 3x⁴ + b₂x³ + 3b₄x² + 3b₆x + b₈` が `xQ` で消える）、
`P = (x,y)` を `x ≠ xQ` なる点とすると、
**`(velu3X, velu3Y)` は商曲線 `veluCurve W v_Q w_Q` 上にある**。

★★★これで Vélu の構成は `#H = 2` と `#H = 3` の両方で実証された
——★`l = 2` は `u_Q = 0` で退化した場合、`l = 3` は `u_Q ≠ 0` の**一般形**である。

★★係数（`h` と `hpsi` にかける多項式）は sympy の `reduced` で求めた
——`field_simp` 後の多項式は手に負えない大きさである。

☆一般の `l` については代表系 `S` の取り方と群法則が要る。 -/
theorem velu3_equation (W : WeierstrassCurve F) (xQ yQ x y : F)
    (hq : W.toAffine.Equation xQ yQ)
    (hpsi : 3 * xQ ^ 4 + W.b₂ * xQ ^ 3 + 3 * W.b₄ * xQ ^ 2 + 3 * W.b₆ * xQ + W.b₈ = 0)
    (h : W.toAffine.Equation x y) (hx : x ≠ xQ) :
    (veluCurve W (veluV W xQ yQ) (veluW W xQ yQ)).toAffine.Equation
      (velu3X W xQ yQ x) (velu3Y W xQ yQ x y) := by
  rw [WeierstrassCurve.Affine.equation_iff] at hq h ⊢
  have hne : x - xQ ≠ 0 := sub_ne_zero.2 hx
  have ha6 : W.a₆ = yQ ^ 2 + W.a₁ * xQ * yQ + W.a₃ * yQ
      - xQ ^ 3 - W.a₂ * xQ ^ 2 - W.a₄ * xQ := by linear_combination -hq
  simp only [WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈] at hpsi
  simp only [velu3X, velu3Y, veluCurve, veluV, veluU, veluW, veluGx, veluGy,
    WeierstrassCurve.b₂]
  rw [ha6] at hpsi h ⊢
  field_simp
  linear_combination
    (-W.a₁^2*x*xQ - W.a₁^2*xQ^2 - W.a₁*W.a₃*x - 3*W.a₁*W.a₃*xQ - 8*W.a₁*xQ*yQ - 4*W.a₂*x*xQ
      + 4*W.a₂*xQ^2 - 2*W.a₃^2 - 8*W.a₃*yQ - 2*W.a₄*x + 2*W.a₄*xQ + x^3 - 3*x^2*xQ
      - 3*x*xQ^2 + 5*xQ^3 - 8*yQ^2) ^ 2 * h
    + (-(x - xQ) ^ 2 * (-2*W.a₁^2*x*xQ - W.a₁^2*xQ^2 - 2*W.a₁*W.a₃*x - 4*W.a₁*W.a₃*xQ
      - 12*W.a₁*xQ*yQ - 8*W.a₂*x*xQ + 8*W.a₂*xQ^2 - 3*W.a₃^2 - 12*W.a₃*yQ - 4*W.a₄*x
      + 4*W.a₄*xQ + 6*x^3 - 18*x^2*xQ + 6*x*xQ^2 + 6*xQ^3 - 12*yQ^2)) * hpsi

/-- ★★★★★★★★★★★★★★★★★★★★★★**`l = 3`: `φ^*(ω′) = ω`——★仮定なし**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `(dX/dx)·(2y + a₁x + a₃) = 2Y + a₁X + a₃`

（`dX/dx = 1 − v_Q/(x−xQ)² − 2u_Q/(x−xQ)³`）。

★★★**`l = 3` では位数の条件すら要らず、恒等式として成り立つ**
——Vélu の写像は最初から `ω` を保つように書かれている。
★★これが `§9-1027`（第 585）の見立ての中核であり、
ℂ 上で周期格子が `Λ ⊆ Λ′` になる（スケーリングが入らない）ことを意味する。 -/
theorem velu3_omega (W : WeierstrassCurve F) (xQ yQ x y : F) (hx : x ≠ xQ) :
    (1 - veluV W xQ yQ / (x - xQ) ^ 2 - 2 * veluU W xQ yQ / (x - xQ) ^ 3)
        * (2 * y + W.a₁ * x + W.a₃)
      = 2 * velu3Y W xQ yQ x y + W.a₁ * velu3X W xQ yQ x + W.a₃ := by
  have hne : x - xQ ≠ 0 := sub_ne_zero.2 hx
  simp only [velu3X, velu3Y, veluV, veluU, veluGx, veluGy]
  field_simp
  ring

end ThreeIsogeny

/-! ## ★出典の紐付け(`.src`)——★★**条つき。定義と帳簿だけである** -/

def velu2_equation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の 2-同種写像は本当に商へ落ちる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def velu3_equation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の 3-同種写像は本当に商へ落ちる——u_Q ≠ 0 の一般形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def velu3_omega.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の正規化 φ^*(ω′) = ω——l = 3 では仮定なし。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluCurve_scale.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商はスケーリングと両立する——v は重さ 4、w は重さ 6)",
    sectionId := "genell-lemma-3-5" }

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

section GeneralL

variable {F : Type*} [Field F]

/-! ## ★★★★★★★★★★★★★★★★★★`±` の不変性——「代表系」という言葉が正しい理由 -/

/-- ★★★★★点の `±`——`−Q = (x_Q, negY(x_Q, y_Q))`。

★mathlib の `WeierstrassCurve.Affine.negY x y = −y − a₁x − a₃` をそのまま使う。 -/
def negPt (W : WeierstrassCurve R) (Q : R × R) : R × R := (Q.1, W.toAffine.negY Q.1 Q.2)

/-- ★★★`−(−Q) = Q`。 -/
theorem negPt_involutive (W : WeierstrassCurve R) : Function.Involutive (negPt W) := by
  intro Q
  simp only [negPt, WeierstrassCurve.Affine.negY_negY]

/-- ★★★★★**`g^y` は `±` で符号が変わる**——`g^y_{−Q} = −g^y_Q`。 -/
theorem veluGy_negY (W : WeierstrassCurve R) (x y : R) :
    veluGy W x (W.toAffine.negY x y) = -veluGy W x y := by
  simp only [veluGy, WeierstrassCurve.Affine.negY]; ring

/-- ★★★★★**`g^x_{−Q} = g^x_Q − a₁·g^y_Q`**。 -/
theorem veluGx_negY (W : WeierstrassCurve R) (x y : R) :
    veluGx W x (W.toAffine.negY x y) = veluGx W x y - W.a₁ * veluGy W x y := by
  simp only [veluGx, veluGy, WeierstrassCurve.Affine.negY]; ring

/-- ★★★★★★**`u_Q` は `±` で不変**——`(g^y)²` だから。 -/
theorem veluU_negY (W : WeierstrassCurve R) (x y : R) :
    veluU W x (W.toAffine.negY x y) = veluU W x y := by
  simp only [veluU, veluGy_negY]; ring

/-- ★★★★★★★★★★**`v_Q` は `±` で不変**。

`v_{−Q} = 2g^x_{−Q} − a₁g^y_{−Q} = 2(g^x − a₁g^y) − a₁(−g^y) = 2g^x − a₁g^y = v_Q`。

★★★**これが「`(H∖{O})/±` の代表系にわたる和」という言い方が正しい理由である。** -/
theorem veluV_negY (W : WeierstrassCurve R) (x y : R) :
    veluV W x (W.toAffine.negY x y) = veluV W x y := by
  simp only [veluV, veluGx_negY, veluGy_negY]; ring

/-- ★★★★★★★★★★**`w_Q` は `±` で不変**——`x_{−Q} = x_Q` だから。 -/
theorem veluW_negY (W : WeierstrassCurve R) (x y : R) :
    veluW W x (W.toAffine.negY x y) = veluW W x y := by
  simp only [veluW, veluU_negY, veluV_negY]

theorem veluV_negPt (W : WeierstrassCurve R) (Q : R × R) :
    veluV W (negPt W Q).1 (negPt W Q).2 = veluV W Q.1 Q.2 := veluV_negY W Q.1 Q.2

theorem veluW_negPt (W : WeierstrassCurve R) (Q : R × R) :
    veluW W (negPt W Q).1 (negPt W Q).2 = veluW W Q.1 Q.2 := veluW_negY W Q.1 Q.2

/-- ★★★★★★★★★★★★**代表系を取り替えても `v` の和は変わらない**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`e` は「各代表 `Q` を `Q` か `−Q` に送る全単射」である。 -/
theorem veluVSum_congr_negPt (W : WeierstrassCurve R) (S T : Finset (R × R))
    (e : R × R → R × R)
    (hmem : ∀ Q ∈ S, e Q ∈ T)
    (hinj : ∀ Q ∈ S, ∀ Q' ∈ S, e Q = e Q' → Q = Q')
    (hsurj : ∀ P ∈ T, ∃ Q, ∃ _ : Q ∈ S, e Q = P)
    (he : ∀ Q ∈ S, e Q = Q ∨ e Q = negPt W Q) :
    veluVSum W S = veluVSum W T := by
  refine Finset.sum_bij (fun Q _ => e Q) (fun Q hQ => hmem Q hQ)
    (fun Q hQ Q' hQ' h => hinj Q hQ Q' hQ' h) hsurj ?_
  intro Q hQ
  rcases he Q hQ with h | h
  · rw [h]
  · rw [h, veluV_negPt]

/-- ★★★★★★★★★★★★**代表系を取り替えても `w` の和は変わらない**。 -/
theorem veluWSum_congr_negPt (W : WeierstrassCurve R) (S T : Finset (R × R))
    (e : R × R → R × R)
    (hmem : ∀ Q ∈ S, e Q ∈ T)
    (hinj : ∀ Q ∈ S, ∀ Q' ∈ S, e Q = e Q' → Q = Q')
    (hsurj : ∀ P ∈ T, ∃ Q, ∃ _ : Q ∈ S, e Q = P)
    (he : ∀ Q ∈ S, e Q = Q ∨ e Q = negPt W Q) :
    veluWSum W S = veluWSum W T := by
  refine Finset.sum_bij (fun Q _ => e Q) (fun Q hQ => hmem Q hQ)
    (fun Q hQ Q' hQ' h => hinj Q hQ Q' hQ' h) hsurj ?_
  intro Q hQ
  rcases he Q hQ with h | h
  · rw [h]
  · rw [h, veluW_negPt]

/-- ★★★★★★★★★★★★★★★★★★**Vélu の商は代表系の取り方に依らない**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これで `veluQuotient W S` の `S` は本当に `(H∖{O})/±` の**代表系**であってよい
——どの代表を選んでも同じ曲線が出る。☆本ファイルの散文の主張が定理になった。 -/
theorem veluQuotient_congr_negPt (W : WeierstrassCurve R) (S T : Finset (R × R))
    (e : R × R → R × R)
    (hmem : ∀ Q ∈ S, e Q ∈ T)
    (hinj : ∀ Q ∈ S, ∀ Q' ∈ S, e Q = e Q' → Q = Q')
    (hsurj : ∀ P ∈ T, ∃ Q, ∃ _ : Q ∈ S, e Q = P)
    (he : ∀ Q ∈ S, e Q = Q ∨ e Q = negPt W Q) :
    veluQuotient W S = veluQuotient W T := by
  rw [veluQuotient, veluQuotient,
    veluVSum_congr_negPt W S T e hmem hinj hsurj he,
    veluWSum_congr_negPt W S T e hmem hinj hsurj he]

/-- ★★★★★★★★★★**全部の代表を `−` に取り替えても商は同じ**——`e = negPt W` の場合。 -/
theorem veluQuotient_image_negPt (W : WeierstrassCurve R) (S : Finset (R × R))
    [DecidableEq (R × R)] :
    veluQuotient W S = veluQuotient W (S.image (negPt W)) := by
  refine veluQuotient_congr_negPt W S _ (negPt W) (fun Q hQ => Finset.mem_image_of_mem _ hQ)
    (fun Q _ Q' _ h => (negPt_involutive W).injective h) ?_ (fun Q _ => Or.inr rfl)
  intro P hP
  obtain ⟨Q, hQ, rfl⟩ := Finset.mem_image.1 hP
  exact ⟨Q, hQ, rfl⟩

def veluQuotient_congr_negPt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は代表系の取り方に依らない——v_Q・u_Q は ± で不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★一般の `l`——代表系 `S` にわたる和 -/

/-- ★★★★★点 `Q` ごとの `X` 補正 `v_Q/(x−x_Q) + u_Q/(x−x_Q)²`。 -/
noncomputable def veluXterm (W : WeierstrassCurve F) (xQ yQ x : F) : F :=
  veluV W xQ yQ / (x - xQ) + veluU W xQ yQ / (x - xQ) ^ 2

/-- ★★★★★点 `Q` ごとの `Y` 補正。 -/
noncomputable def veluYterm (W : WeierstrassCurve F) (xQ yQ x y : F) : F :=
  veluU W xQ yQ * (2 * y + W.a₁ * x + W.a₃) / (x - xQ) ^ 3
    + veluV W xQ yQ * (W.a₁ * (x - xQ) + y - yQ) / (x - xQ) ^ 2
    + (W.a₁ * veluU W xQ yQ - veluGx W xQ yQ * veluGy W xQ yQ) / (x - xQ) ^ 2

/-- ★★★★★点 `Q` ごとの `dX/dx` 補正 `v_Q/(x−x_Q)² + 2u_Q/(x−x_Q)³`。 -/
noncomputable def veluDterm (W : WeierstrassCurve F) (xQ yQ x : F) : F :=
  veluV W xQ yQ / (x - xQ) ^ 2 + 2 * veluU W xQ yQ / (x - xQ) ^ 3

/-- ★★★★★★★★**正規化は点ごとに成り立つ**——仮定なし。

    `(dX/dx)_Q · (2y + a₁x + a₃) = 2·Y_Q − a₁·X_Q`

★これが `velu3_omega`（第 590）の中身であり、★★**和を取ればそのまま一般の `l` になる**。 -/
theorem velu_omega_term (W : WeierstrassCurve F) (xQ yQ x y : F) (hx : x ≠ xQ) :
    veluDterm W xQ yQ x * (2 * y + W.a₁ * x + W.a₃)
      = 2 * veluYterm W xQ yQ x y - W.a₁ * veluXterm W xQ yQ x := by
  have hne : x - xQ ≠ 0 := sub_ne_zero.2 hx
  simp only [veluDterm, veluYterm, veluXterm, veluV, veluU, veluGx, veluGy]
  field_simp
  ring

/-- ★★★★★★**Vélu の `X`**——代表系 `S` にわたる和。 -/
noncomputable def veluXGen (W : WeierstrassCurve F) (S : Finset (F × F)) (x : F) : F :=
  x + ∑ Q ∈ S, veluXterm W Q.1 Q.2 x

/-- ★★★★★★**Vélu の `Y`**——代表系 `S` にわたる和。 -/
noncomputable def veluYGen (W : WeierstrassCurve F) (S : Finset (F × F)) (x y : F) : F :=
  y - ∑ Q ∈ S, veluYterm W Q.1 Q.2 x y

/-- ★★★★★★**Vélu の `dX/dx`**——代表系 `S` にわたる和。 -/
noncomputable def veluDGen (W : WeierstrassCurve F) (S : Finset (F × F)) (x : F) : F :=
  1 - ∑ Q ∈ S, veluDterm W Q.1 Q.2 x

/-- ★★★★★★★★★★★★★★**一般の `l` での正規化 `φ^*(ω′) = ω`**。

不変微分は `ω = dx/(2y + a₁x + a₃)` なので、

    `(dX/dx)·(2y + a₁x + a₃) = 2Y + a₁X + a₃`

は `φ^*(ω_{E′}) = ω_E` と同じことである。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★**仮定は `x ∉ S の x 座標` だけ**——群法則も `S` が部分群の代表系であることも要らない。
★★Vélu の写像は**最初から `ω` を保つように書かれている**。

★★★★★これが `§9-1027`（第 585 の測定）の見立て
「アルキメデス項が消える」の中身である——ℂ 上で `ω = dz` なら
周期格子がそのまま `Λ ⊆ Λ′` になり、スケーリングが入らない。 -/
theorem velu_omega_gen (W : WeierstrassCurve F) (S : Finset (F × F)) (x y : F)
    (hS : ∀ Q ∈ S, x ≠ Q.1) :
    veluDGen W S x * (2 * y + W.a₁ * x + W.a₃)
      = 2 * veluYGen W S x y + W.a₁ * veluXGen W S x + W.a₃ := by
  have key : (∑ Q ∈ S, veluDterm W Q.1 Q.2 x) * (2 * y + W.a₁ * x + W.a₃)
      = 2 * (∑ Q ∈ S, veluYterm W Q.1 Q.2 x y)
        - W.a₁ * (∑ Q ∈ S, veluXterm W Q.1 Q.2 x) := by
    rw [Finset.sum_mul, Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun Q hQ => velu_omega_term W Q.1 Q.2 x y (hS Q hQ)
  simp only [veluDGen, veluXGen, veluYGen]
  rw [sub_mul, one_mul, key]
  ring

def velu_omega_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の正規化を一般の l へ——代表系 S にわたる和。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★`S = ∅` なら写像は恒等。 -/
@[simp] theorem veluXGen_empty (W : WeierstrassCurve F) (x : F) : veluXGen W ∅ x = x := by
  simp [veluXGen]

/-- ★★★`S = ∅` なら写像は恒等。 -/
@[simp] theorem veluYGen_empty (W : WeierstrassCurve F) (x y : F) : veluYGen W ∅ x y = y := by
  simp [veluYGen]

/-- ★★★`S = ∅` なら `dX/dx = 1`。 -/
@[simp] theorem veluDGen_empty (W : WeierstrassCurve F) (x : F) : veluDGen W ∅ x = 1 := by
  simp [veluDGen]

/-- ★★★★★**`S` が 1 点なら `velu3X`・`velu3Y` に一致**。

★これで第 589（`l = 3`）が一般形の特別な場合であることが確かめられる。 -/
theorem veluXGen_singleton (W : WeierstrassCurve F) (xQ yQ x : F) :
    veluXGen W {(xQ, yQ)} x = velu3X W xQ yQ x := by
  simp only [veluXGen, veluXterm, velu3X, Finset.sum_singleton]
  ring

/-- ★★★★★**`S` が 1 点なら `velu3Y` に一致**。 -/
theorem veluYGen_singleton (W : WeierstrassCurve F) (xQ yQ x y : F) :
    veluYGen W {(xQ, yQ)} x y = velu3Y W xQ yQ x y := by
  simp only [veluYGen, veluYterm, velu3Y, Finset.sum_singleton]
  ring

end GeneralL

end ABC3.Found.GenEll
