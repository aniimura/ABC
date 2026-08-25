import ABC3.Found.GaloisRep.LatticeInvariant

/-!
# Galois (G6) 第 216 ブロック —— **★★★★★★Lipschitz 公式を格子の各段に当てる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★道 β の段 5(の前半)

第 215 で段 4(`g₂, g₃` の q 展開)を取った。残るのは:

| 段 | 内容 | 状態 |
|---|---|---|
| 4 | 格子 `ℤ + τℤ` の `g₂, g₃` を Eisenstein の q 展開に繋ぐ | 済(第 215) |
| 5 | `℘` の q 展開(Lipschitz 公式を `n ∈ ℤ` の各段に当てる) | ★**本ブロックはその前半** |
| 6 | 「`ℤ` 係数の形式級数が関数として 0 なら形式的に 0」 | 未着手 |

## ★★★★★格子を τ の段に分ける

格子 `ℤτ + ℤ` は `n ∈ ℤ` ごとの**横一列**(`nτ + ℤ`)に分かれる。
`℘` の二重和を段ごとに切り、各段に Lipschitz 公式(第 214 の `lipschitz_tateXterm`)を当てると:

    ∑_{m ∈ ℤ} 1/(z + nτ + m)²  =  (2πi)² · tateXterm(qⁿ u)        (n ≥ 0)
    ∑_{m ∈ ℤ} 1/(z − (n+1)τ + m)² = (2πi)² · tateXterm(q^{n+1} u⁻¹)  (n ≥ 0)

ここで `q = exp(2πiτ)`、`u = exp(2πiz)`。
★★**これがちょうど形式側 `tateXpair a w` の `(a, w)` の対の形である**——
`a` が `u` に、`w` が `q` に対応し、上段が `tateXterm (qⁿ a)`、下段が `tateXterm (qⁿ⁺¹ a⁻¹)`。

## ★★★負の段は「折り返し」で取る

`n < 0` の段は `z − (n+1)τ` の虚部が負になりうるので、上半平面の点として直接は扱えない。
★`∑' m, 1/(w + m)² = ∑' m, 1/(−w + m)²`(`m ↦ −m` の置換)で折り返してから
`(n+1)τ − z`(これは `im z < im τ` なら上半平面にある)に当てる。
★★★**平方だから折り返しで符号が消える**——これが `tsum_inv_sq_neg` の要点である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `shiftUp` | `z + nτ` を上半平面の点として |
| `exp_add_nat_mul` | `exp(2πi(z+nτ)) = qⁿ·u` |
| `lipschitz_shift` | ★★★★★**上向きの段の Lipschitz 公式** |
| `tsum_inv_sq_neg` | ★★★★**平方和は折り返しで不変** |
| `shiftDown` / `exp_sub_nat_mul` | `(n+1)τ − z` の側 |
| `lipschitz_shift_neg` | ★★★★★★**下向きの段の Lipschitz 公式** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-- ★`z + nτ` は上半平面にある(`n ≥ 0`)。 -/
noncomputable def shiftUp (z τ : UpperHalfPlane) (n : ℕ) : UpperHalfPlane :=
  ⟨(z : ℂ) + n * (τ : ℂ), by
    have hz : 0 < (z : ℂ).im := z.2
    have hτ : 0 < (τ : ℂ).im := τ.2
    simp only [Complex.add_im, Complex.mul_im, Complex.natCast_re, Complex.natCast_im]
    nlinarith [Nat.cast_nonneg (α := ℝ) n]⟩

theorem shiftUp_coe (z τ : UpperHalfPlane) (n : ℕ) :
    ((shiftUp z τ n : UpperHalfPlane) : ℂ) = (z : ℂ) + n * (τ : ℂ) := rfl

/-- ★`exp(2πi(z + nτ)) = qⁿ · u`。 -/
theorem exp_add_nat_mul (z τ : ℂ) (n : ℕ) :
    Complex.exp (2 * ↑π * I * (z + n * τ))
      = Complex.exp (2 * ↑π * I * τ) ^ n * Complex.exp (2 * ↑π * I * z) := by
  rw [mul_add, Complex.exp_add, ← Complex.exp_nat_mul]
  ring_nf

/-- ★★★★★**上向きの段の Lipschitz 公式**——`n` 段目は `tateXterm (qⁿ u)` になる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem lipschitz_shift (z τ : UpperHalfPlane) (n : ℕ) :
    ∑' m : ℤ, 1 / ((z : ℂ) + n * (τ : ℂ) + m) ^ 2
      = (2 * ↑π * I) ^ 2 * tateXterm
        (Complex.exp (2 * ↑π * I * τ) ^ n * Complex.exp (2 * ↑π * I * z)) := by
  have h := lipschitz_tateXterm (shiftUp z τ n)
  rw [shiftUp_coe] at h
  rw [h, exp_add_nat_mul]

/-- ★★★★**平方和は折り返しで不変**——`m ↦ −m` の置換で `w` を `−w` に替えられる。

★平方だから符号が消える。負の段を上半平面に持ち上げるのに使う。 -/
theorem tsum_inv_sq_neg (w : ℂ) :
    ∑' m : ℤ, 1 / (w + m) ^ 2 = ∑' m : ℤ, 1 / (-w + m) ^ 2 := by
  rw [← (Equiv.neg ℤ).tsum_eq (fun m : ℤ => 1 / (-w + m) ^ 2)]
  refine tsum_congr fun m => ?_
  simp only [Equiv.neg_apply, Int.cast_neg]
  rw [show (-w + -(m : ℂ)) ^ 2 = (w + m) ^ 2 by ring]

/-- ★`(n+1)τ − z` は上半平面にある(`im z < im τ` のとき)。 -/
noncomputable def shiftDown (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (n : ℕ) :
    UpperHalfPlane :=
  ⟨(n + 1 : ℕ) * (τ : ℂ) - (z : ℂ), by
    have hτ : 0 < (τ : ℂ).im := τ.2
    simp only [Complex.sub_im, Complex.mul_im, Complex.natCast_re, Complex.natCast_im]
    push_cast
    nlinarith [Nat.cast_nonneg (α := ℝ) n]⟩

theorem shiftDown_coe (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (n : ℕ) :
    ((shiftDown z τ him n : UpperHalfPlane) : ℂ) = (n + 1 : ℕ) * (τ : ℂ) - (z : ℂ) := rfl

/-- ★`exp(2πi((n+1)τ − z)) = q^{n+1} · u⁻¹`。 -/
theorem exp_sub_nat_mul (z τ : ℂ) (n : ℕ) :
    Complex.exp (2 * ↑π * I * ((n + 1 : ℕ) * τ - z))
      = Complex.exp (2 * ↑π * I * τ) ^ (n + 1) * (Complex.exp (2 * ↑π * I * z))⁻¹ := by
  rw [mul_sub, Complex.exp_sub, ← Complex.exp_nat_mul, div_eq_mul_inv]
  congr 2
  push_cast
  ring

/-- ★★★★★★**下向きの段の Lipschitz 公式**——`−(n+1)` 段目は `tateXterm (q^{n+1} u⁻¹)` になる。

★折り返し(`tsum_inv_sq_neg`)で上半平面に持ち上げてから当てる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem lipschitz_shift_neg (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (n : ℕ) :
    ∑' m : ℤ, 1 / ((z : ℂ) - (n + 1 : ℕ) * (τ : ℂ) + m) ^ 2
      = (2 * ↑π * I) ^ 2 * tateXterm
        (Complex.exp (2 * ↑π * I * τ) ^ (n + 1) * (Complex.exp (2 * ↑π * I * z))⁻¹) := by
  have hrefl := tsum_inv_sq_neg ((z : ℂ) - (n + 1 : ℕ) * (τ : ℂ))
  have h := lipschitz_tateXterm (shiftDown z τ him n)
  rw [shiftDown_coe] at h
  rw [hrefl, show -((z : ℂ) - (n + 1 : ℕ) * (τ : ℂ)) = (n + 1 : ℕ) * (τ : ℂ) - (z : ℂ) by ring,
    h, exp_sub_nat_mul]

/-! ## ★出典の紐付け(`.src`) -/

def lipschitz_shift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子の上向きの段の Lipschitz 公式)",
    sectionId := "genell-def-3-3" }

def tsum_inv_sq_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——平方和の折り返し)",
    sectionId := "genell-def-3-3" }

def lipschitz_shift_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子の下向きの段の Lipschitz 公式)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
