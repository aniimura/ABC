/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.HeightOneBridge
import ABC3.Found.Divisor.CartierMonoid

/-!
# Cartier 因子(局所主因子)—— 鎖 `cartier` の `cartier-def` / `cartier-subgroup` / `qcartier`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> then we shall say that DK is K-Q-Cartier. In the following, we shall assume that

## ★★Cartier 因子を Weil 因子の言葉で書く

原文は `D_K` に **`K`-`Q`-Cartier** 性を仮定する。その土台になるのが
「**局所主因子**」——各点のまわりで `div(f)` になっている Weil 因子である:

  `IsCartierDiv D ⟺ ∀ x, ∃ 開近傍 U と f ∈ K(X)^×, U 上で D = div(f)`

★★**部分群になることが要点**である。`Found/Divisor/CartierMonoid.lean` の
単系論はすべて `Γ : AddSubgroup (S →₀ ℤ)` からしか出発しないので、
**ここが `Γ` の供給源**になる。

## ★★★中身は `ord` の乗法性だけ

| フィールド | 根拠 |
|---|---|
| `0 ∈ Γ` | `f = 1`、`ordPt_one` |
| `D₁ + D₂ ∈ Γ` | `U₁ ⊓ U₂` の上で `f₁·f₂`、`ordPt_mul` |
| `-D ∈ Γ` | 同じ `U` の上で `f⁻¹`、`ordPt_inv` |

★**新しい幾何は 1 つも要らない** —— `SchemeWeilOrd.lean` の `ordPt_mul` / `ordPt_one` だけ。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `ordPt_inv` | `ord_x(f⁻¹) = -ord_x(f)` |
| `IsCartierDiv` | ★局所主因子の述語(`cartier-def`) |
| `cartierSubgroup` | ★★**Cartier 因子は部分群**(`cartier-subgroup`) |
| `IsQCartierDiv` / `isQCartierSubgroup_cartierSubgroup` | ★`Q`-Cartier(`qcartier`) |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry

universe u

variable {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X]

/-! ## ★1. `ord` の逆元則 -/

/-- ★**`ord_x(f⁻¹) = -ord_x(f)`**。 -/
theorem ordPt_inv (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    {f : X.functionField} (hf : f ≠ 0) :
    ordPt X hnorm x f⁻¹ = -ordPt X hnorm x f := by
  have hinv : f⁻¹ ≠ 0 := inv_ne_zero hf
  have h := ordPt_mul hnorm x hf hinv
  rw [mul_inv_cancel₀ hf, ordPt_one] at h
  omega

/-! ## ★2. 局所主因子(`cartier-def`) -/

/-- ★★★**Cartier 因子(局所主因子)** —— 各点のまわりで `div(f)` になっている Weil 因子。

原文 (FrdI p.109):
> identified with the group of Cartier divisors on V [L], and that -/
def IsCartierDiv (hnorm : IsNormalScheme X) (D : WeilDiv X) : Prop :=
  ∀ x : X, ∃ (U : X.Opens) (_ : x ∈ U) (f : X.functionField) (_ : f ≠ 0),
    ∀ v : PrimeDivisorPt X, (v : X) ∈ U → D v = ordPt X hnorm v f

/-! ## ★3. 部分群であること(`cartier-subgroup`) -/

/-- ★★★★**Cartier 因子は `WeilDiv X` の部分群**。

★中身は `ordPt_one` / `ordPt_mul` / `ordPt_inv` だけである。 -/
def cartierSubgroup (hnorm : IsNormalScheme X) : AddSubgroup (WeilDiv X) where
  carrier := {D | IsCartierDiv hnorm D}
  zero_mem' := by
    intro x
    exact ⟨⊤, trivial, 1, one_ne_zero, fun v _ => by simp [ordPt_one]⟩
  add_mem' := by
    rintro D₁ D₂ h₁ h₂ x
    obtain ⟨U₁, hx₁, f₁, hf₁, hD₁⟩ := h₁ x
    obtain ⟨U₂, hx₂, f₂, hf₂, hD₂⟩ := h₂ x
    refine ⟨U₁ ⊓ U₂, ⟨hx₁, hx₂⟩, f₁ * f₂, mul_ne_zero hf₁ hf₂, fun v hv => ?_⟩
    rw [Finsupp.add_apply, hD₁ v hv.1, hD₂ v hv.2, ordPt_mul hnorm v hf₁ hf₂]
  neg_mem' := by
    rintro D h x
    obtain ⟨U, hx, f, hf, hD⟩ := h x
    refine ⟨U, hx, f⁻¹, inv_ne_zero hf, fun v hv => ?_⟩
    rw [Finsupp.neg_apply, hD v hv, ordPt_inv hnorm v hf]

@[simp] theorem mem_cartierSubgroup (hnorm : IsNormalScheme X) (D : WeilDiv X) :
    D ∈ cartierSubgroup hnorm ↔ IsCartierDiv hnorm D := Iff.rfl

/-! ## ★4. `Q`-Cartier(`qcartier`) -/

/-- ★★**`Q`-Cartier** —— ある正の整数倍が Cartier。

原文の「some positive integer multiple of each of the prime divisors」そのもの。 -/
def IsQCartierDiv (hnorm : IsNormalScheme X) (D : WeilDiv X) : Prop :=
  ∃ n : ℕ, 0 < n ∧ IsCartierDiv hnorm (n • D)

/-- ★★★**各素因子が `Q`-Cartier なら、Cartier 因子の群は `IsQCartierSubgroup`**。

★★これが `Found/Divisor/CartierMonoid.lean` の単系論への**入口**である
(単系側の `IsQCartierSubgroup Γ := ∀ s, ∃ n > 0, single s n ∈ Γ` にそのまま入る)。 -/
theorem isQCartierSubgroup_cartierSubgroup (hnorm : IsNormalScheme X)
    (h : ∀ v : PrimeDivisorPt X, IsQCartierDiv hnorm (Finsupp.single v 1)) :
    ABC3.Found.FrdI.IsQCartierSubgroup (cartierSubgroup hnorm) := by
  intro v
  obtain ⟨n, hn, hc⟩ := h v
  refine ⟨n, hn, ?_⟩
  have hsm : (n : ℤ) • (Finsupp.single v (1 : ℤ)) = Finsupp.single v (n : ℤ) := by
    rw [Finsupp.smul_single, smul_eq_mul, mul_one]
  have hnat : (n • (Finsupp.single v (1 : ℤ)) : WeilDiv X)
      = Finsupp.single v (n : ℤ) := by
    rw [← hsm, ← Nat.cast_smul_eq_nsmul ℤ]
  rw [← hnat]
  exact hc

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の Cartier 因子(局所主因子)と `Q`-Cartier 性。 -/
def IsCartierDiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Cartier 因子(局所主因子)と Q-Cartier 性",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
