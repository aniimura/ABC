import ABC3.Found.GaloisRep.InfinityKer
import ABC3.Found.GaloisRep.RedHom

/-!
# Galois (G5) 第 163 ブロック —— **★★★★★★★★点の還元は群準同型である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★3 つの場合が合流した

第 160 は**両方の座標が付値環に入る**場合の加法性であった。第 162 で残りが埋まり、
本ブロックで 3 つの場合が合流する:

| 場合 | 判定 | 根拠 |
|---|---|---|
| 整 + 整 | `w x₁ ≤ 1`、`w x₂ ≤ 1` | 第 160(`redPoint_add_some`) |
| 無限遠 + 無限遠 | `1 < w x₁`、`1 < w x₂` | 第 162(`one_lt_val_addX_of_infinity`) |
| 無限遠 + 整 | 混合 | 第 162 + **第 160 の再利用** |

### ★★★★★★混合の場合が第 160 から出る仕組み

`S₁ ∈ E₁`、`S₂ ∉ E₁` を直接評価しようとすると `w(ℓ)² = w(x₁)` の打ち消しに突き当たる。
★しかし第 162 の `val_addX_le_of_mixed` で `T := S₁ + S₂` が付値環に入ることが分かる。
★★そこで

    T + (−S₂) = S₁

に**第 160 の整な加法性**を当てると `red T + red(−S₂) = red S₁ = 0`、すなわち
`red T = red S₂`。★★★**混合の場合を直接評価する必要はない**。

## ★★★★これで `n • ` との可換性が出る

`redHom` を `AddMonoidHom` として束ねれば `map_nsmul`・`map_zsmul` が mathlib から来る。
★D2 の残りは「`P' = XYIdeal(n·Q_v)` への組み立て」と「類の計算」だけになった。

## ★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `redPoint_eq_zero_iff` | ★★★★`red S = 0 ⟺ w(x) > 1` |
| `redPoint_neg` | ★★★還元は符号を保つ |
| `redPoint_add_of_infinity` | ★★★★★無限遠 + 無限遠 |
| `redPoint_add_of_mixed` | ★★★★★★無限遠 + 整 |
| `redPoint_add` | ★★★★★★★★**還元は加法を保つ(全場合)** |
| `redHom` | ★★★★★★★★**還元写像(`AddMonoidHom`)** |
| `redPoint_nsmul` / `redPoint_zsmul` | ★★★★★★★**`n • ` と可換** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)
  {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

variable [W.IsElliptic]

/-! ## ★★★★還元の核の判定 -/

/-- ★★★★**還元が 0 になるのは `x` の付値が 1 を超えるときに限る**。

★`w x ≤ 1` なら第 162 の `val_y_le_of_val_x_le` で `w y ≤ 1` も従うので、
`redPoint` の定義の `then` 枝に入り、`Point.some` は `0` ではない。 -/
theorem redPoint_eq_zero_iff {x y : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular x y) :
    redPoint W v h hv (Point.some x y hns) = 0 ↔ 1 < v.valuation W.FunctionField x := by
  constructor
  · intro hz
    by_contra hcon
    rw [not_lt] at hcon
    rw [redPoint_some W v h hv hns hcon (val_y_le_of_val_x_le W v hns.1 hcon)] at hz
    exact absurd hz (by simp)
  · intro hx
    exact redPoint_eq_zero_of_not_le W v h hv hns (fun hcon => absurd hcon.1 (not_le.2 hx))

/-- ★★★**還元は符号を保つ**。 -/
theorem redPoint_neg (S : (W.map (algebraMap F W.FunctionField)).Point) :
    redPoint W v h hv (-S) = -redPoint W v h hv S := by
  match S with
  | 0 => rw [neg_zero, redPoint_zero]; simp
  | Point.some x y hns =>
    rw [Point.neg_some]
    by_cases hx : v.valuation W.FunctionField x ≤ 1
    · have hy := val_y_le_of_val_x_le W v hns.1 hx
      rw [redPoint_some W v h hv _ hx (val_negY_le W v hx hy),
        redPoint_some W v h hv hns hx hy, Point.neg_some]
      simp only [redConst_negY W v h hv hx hy]
    · rw [not_le] at hx
      rw [redPoint_eq_zero_of_not_le W v h hv _ (fun hcon => absurd hcon.1 (not_le.2 hx)),
        redPoint_eq_zero_of_not_le W v h hv hns (fun hcon => absurd hcon.1 (not_le.2 hx))]
      simp

variable [DecidableEq F]

/-! ## ★★★★★残る 2 つの場合 -/

/-- ★★★★★**両方が無限遠に落ちるなら、和も無限遠に落ちる**(第 162)。 -/
theorem redPoint_add_of_infinity {x₁ y₁ x₂ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂ y₂)
    (hx₁ : 1 < v.valuation W.FunctionField x₁) (hx₂ : 1 < v.valuation W.FunctionField x₂) :
    redPoint W v h hv (Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂) = 0 := by
  by_cases hzero : x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂
  · rw [Point.add_of_Y_eq hzero.1 hzero.2, redPoint_zero]
  · rw [Point.add_some hzero]
    exact redPoint_eq_zero_of_not_le W v h hv _ (fun hcon => absurd hcon.1
      (not_le.2 (one_lt_val_addX_of_infinity W v h₁.1 h₂.1 hzero hx₁ hx₂)))

/-- ★★★★★★**片方だけ無限遠なら、和の還元は残りの点の還元に等しい**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 162 の `val_addX_le_of_mixed` で `T := S₁+S₂` が付値環に入ることを言い、
`T + (−S₂) = S₁` に**第 160 の整な加法性**を当てる。
★★直接評価すると `w(ℓ)² = w(x₁)` の打ち消しに突き当たるが、この回り道で完全に避けられる。 -/
theorem redPoint_add_of_mixed (h2 : IsUnit (2 : F)) {x₁ y₁ x₂ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂ y₂)
    (hx₁ : 1 < v.valuation W.FunctionField x₁)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1) :
    redPoint W v h hv (Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂)
      = redPoint W v h hv (Point.some x₂ y₂ h₂) := by
  have hy₂ : v.valuation W.FunctionField y₂ ≤ 1 := val_y_le_of_val_x_le W v h₂.1 hx₂
  have hne : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂) := by
    rintro ⟨rfl, -⟩; exact absurd hx₂ (not_le.2 hx₁)
  have hx₃ := val_addX_le_of_mixed W v h2 h₁.1 h₂.1 hx₁ hx₂ hy₂ hne
  have hns₃ := nonsingular_add h₁ h₂ hne
  have hy₃ := val_y_le_of_val_x_le W v hns₃.1 hx₃
  have hnsN : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂
      ((W.map (algebraMap F W.FunctionField)).negY x₂ y₂) := (nonsingular_neg x₂ y₂).mpr h₂
  have hstep := redPoint_add_some W v h hv hns₃ hnsN hx₃ hy₃ hx₂ (val_negY_le W v hx₂ hy₂)
  rw [← Point.add_some hne] at hstep
  rw [show Point.some x₂ ((W.map (algebraMap F W.FunctionField)).negY x₂ y₂) hnsN
      = -Point.some x₂ y₂ h₂ from (Point.neg_some h₂).symm,
    add_neg_cancel_right, redPoint_neg W v h hv,
    (redPoint_eq_zero_iff W v h hv h₁).2 hx₁] at hstep
  have hsym := hstep.symm
  rwa [← sub_eq_add_neg, sub_eq_zero] at hsym

/-! ## ★★★★★★★★合流 -/

/-- ★★★★★★★★**点の還元は加法を保つ**(全場合)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★整 + 整(第 160)、無限遠 + 無限遠(第 162)、混合(第 162 + 第 160)の 3 つ。 -/
theorem redPoint_add (h2 : IsUnit (2 : F))
    (S₁ S₂ : (W.map (algebraMap F W.FunctionField)).Point) :
    redPoint W v h hv (S₁ + S₂) = redPoint W v h hv S₁ + redPoint W v h hv S₂ := by
  match S₁, S₂ with
  | 0, S₂ => rw [zero_add, redPoint_zero, zero_add]
  | Point.some x₁ y₁ h₁, 0 => rw [add_zero, redPoint_zero, add_zero]
  | Point.some x₁ y₁ h₁, Point.some x₂ y₂ h₂ =>
    by_cases hx₁ : v.valuation W.FunctionField x₁ ≤ 1
    · by_cases hx₂ : v.valuation W.FunctionField x₂ ≤ 1
      · exact redPoint_add_some W v h hv h₁ h₂ hx₁ (val_y_le_of_val_x_le W v h₁.1 hx₁)
          hx₂ (val_y_le_of_val_x_le W v h₂.1 hx₂)
      · rw [not_le] at hx₂
        rw [add_comm, redPoint_add_of_mixed W v h hv h2 h₂ h₁ hx₂ hx₁,
          (redPoint_eq_zero_iff W v h hv h₂).2 hx₂, add_zero]
    · rw [not_le] at hx₁
      by_cases hx₂ : v.valuation W.FunctionField x₂ ≤ 1
      · rw [redPoint_add_of_mixed W v h hv h2 h₁ h₂ hx₁ hx₂,
          (redPoint_eq_zero_iff W v h hv h₁).2 hx₁, zero_add]
      · rw [not_le] at hx₂
        rw [redPoint_add_of_infinity W v h hv h₁ h₂ hx₁ hx₂,
          (redPoint_eq_zero_iff W v h hv h₁).2 hx₁,
          (redPoint_eq_zero_iff W v h hv h₂).2 hx₂, add_zero]

/-- ★★★★★★★★**還元写像**——加法群の準同型として束ねたもの。 -/
noncomputable def redHom (h2 : IsUnit (2 : F)) :
    (W.map (algebraMap F W.FunctionField)).Point →+ W.Point where
  toFun := redPoint W v h hv
  map_zero' := redPoint_zero W v h hv
  map_add' := redPoint_add W v h hv h2

@[simp] theorem redHom_apply (h2 : IsUnit (2 : F))
    (S : (W.map (algebraMap F W.FunctionField)).Point) :
    redHom W v h hv h2 S = redPoint W v h hv S := rfl

/-- ★★★★★★★**還元は `n • ` と可換**——D2 に必要だった段。 -/
theorem redPoint_nsmul (h2 : IsUnit (2 : F)) (n : ℕ)
    (S : (W.map (algebraMap F W.FunctionField)).Point) :
    redPoint W v h hv (n • S) = n • redPoint W v h hv S := by
  simpa using map_nsmul (redHom W v h hv h2) n S

/-- ★★★★★★整数倍でも同じ。 -/
theorem redPoint_zsmul (h2 : IsUnit (2 : F)) (n : ℤ)
    (S : (W.map (algebraMap F W.FunctionField)).Point) :
    redPoint W v h hv (n • S) = n • redPoint W v h hv S := by
  simpa using map_zsmul (redHom W v h hv h2) n S

/-! ## ★出典の紐付け(`.src`) -/

def redPoint_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——点の還元が群準同型であること)",
    sectionId := "genell-thm-3-8" }

def redPoint_nsmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——点の還元が n 倍と可換であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
