import ABC3.Found.GaloisRep.TateClassMap
import Mathlib.NumberTheory.ModularForms.EisensteinSeries.QExpansion
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass

/-!
# Galois (G6) 第 214 ブロック —— **★★★★★葉 (b) の道を測って選んだ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★葉 (b) には道が 2 本ある

第 212 で Weierstrass 方程式は `I` を法として取れた。★厳密な等式には道が 2 本ある。

### 道 α —— 形式的・代数的

`ℤ[u, u⁻¹, (1−u)⁻¹][[q]]` の中で係数ごとに確かめる。★★`q` の各係数は `u` の有理式の
恒等式になるが、**有限の帰着が見えない**——係数ごとに別の議論が要る。

### 道 β —— 複素解析から移す

`℘` の理論で `(X, Y)` が方程式を満たすことは古典的に分かる。★係数が `ℤ` なので、
関数としての等式から**形式級数の等式**が従う。

## ★★★★★★測ったら道 β の在庫が厚かった(2026-08-20)

| 段 | 内容 | mathlib 在庫 |
|---|---|---|
| 1 | `℘` の構成と **`℘'² = 4℘³ − g₂℘ − g₃`** | ★★★★`PeriodPair.derivWeierstrassP_sq` ✓ |
| 2 | **Lipschitz 公式** `Σ_{n∈ℤ} 1/(z+n)^{k+1} = ((−2πi)^{k+1}/k!)·Σ n^k qⁿ` | ★★★★`EisensteinSeries.qExpansion_identity` ✓ |
| 3 | Eisenstein の q 展開 `E_k = 1 − (2k/B_k)Σσ_{k−1}(n)qⁿ` | ★★★`EisensteinSeries.q_expansion_bernoulli` ✓ |
| 4 | 格子 `ℤ + τℤ` の `g₂, g₃` を `E₄, E₆` に繋ぐ | ★無い |
| 5 | `℘` の q 展開(2 を `n ∈ ℤ` の各段に当てる) | ★無い |
| 6 | 「`ℤ` 係数の形式級数が関数として 0 なら形式的に 0」 | ★無い |

★★★**解析の難所(`℘` の構成・微分方程式・Lipschitz 公式)はすべて mathlib に在る**。
残る 4-6 は機械的な段であり、**道 β を採る**。

## ★★★★★★本ブロックで橋を 1 本架けた

Lipschitz 公式の `k = 1` は、右辺がちょうど **`f(t) = t/(1−t)²`**——
すなわち我々の `tateXterm`——の形になる:

    Σ_{n∈ℤ} 1/(z+n)² = (2πi)² · f(u),      u = e^{2πiz}

★★`w := z + kτ` に当てれば `f(qᵏu)` が出る——これが `℘` の q 展開の各段である。
★★★**解析側の `f` と形式側の `tateXterm` が同じものであることを型で確かめた**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `lipschitz_two` | ★★★★★**Lipschitz 公式(`k = 1`)** |
| `lipschitz_tateXterm` | ★★★★★★**解析側の `f` は `tateXterm` である** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-- ★★★★★**Lipschitz の公式**(`k = 1`)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`℘` の q 展開の各段はこれである(`z` を `z + kτ` に取り替える)。 -/
theorem lipschitz_two (z : UpperHalfPlane) :
    ∑' n : ℤ, 1 / ((z : ℂ) + n) ^ 2
      = (2 * ↑π * I) ^ 2
        * (Complex.exp (2 * ↑π * I * z) / (1 - Complex.exp (2 * ↑π * I * z)) ^ 2) := by
  have h := EisensteinSeries.qExpansion_identity (k := 1) le_rfl z
  simp only [Nat.factorial_one, Nat.cast_one, div_one, pow_one] at h
  have hnorm : ‖Complex.exp (2 * ↑π * I * (z : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one z
  rw [h, tsum_coe_mul_geometric_of_norm_lt_one hnorm]
  congr 1
  rw [show ((-2 : ℂ) * ↑π * I) ^ 2 = (2 * ↑π * I) ^ 2 by ring]

/-- ★★★★★★**解析側の `f` は形式側の `tateXterm` である**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★これが道 β の入口である——`℘` の各段が `tateXterm (qᵏu)` に化ける。 -/
theorem lipschitz_tateXterm (z : UpperHalfPlane) :
    ∑' n : ℤ, 1 / ((z : ℂ) + n) ^ 2
      = (2 * ↑π * I) ^ 2 * tateXterm (Complex.exp (2 * ↑π * I * z)) := by
  have hnorm : ‖Complex.exp (2 * ↑π * I * (z : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one z
  have hne : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 1 := by
    intro h
    rw [h, norm_one] at hnorm
    exact lt_irrefl 1 hnorm
  rw [lipschitz_two z, tateXterm_eq_div hne]

/-! ## ★出典の紐付け(`.src`) -/

def lipschitz_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Lipschitz の公式)",
    sectionId := "genell-def-3-3" }

def lipschitz_tateXterm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——解析側の項が tateXterm に一致すること)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
