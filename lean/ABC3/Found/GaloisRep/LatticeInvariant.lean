import ABC3.Found.GaloisRep.LipschitzBridge

/-!
# Galois (G6) 第 215 ブロック —— **★★★★★★格子不変量 `g₂, g₃` の q 展開**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★道 β の段 4

第 214 で道 β(複素解析から移す)を選び、残る段を 3 つに絞った:

| 段 | 内容 | 状態 |
|---|---|---|
| 4 | 格子 `ℤ + τℤ` の `g₂, g₃` を Eisenstein の q 展開に繋ぐ | ★**本ブロック** |
| 5 | `℘` の q 展開(Lipschitz 公式を `n ∈ ℤ` の各段に当てる) | 未着手 |
| 6 | 「`ℤ` 係数の形式級数が関数として 0 なら形式的に 0」 | 未着手 |

## ★★★★★格子和と Eisenstein 和は同じもの

mathlib は 2 つの世界を別々に持っている:

| 世界 | 和のとり方 |
|---|---|
| `PeriodPair.G n` | `∑' l : L.lattice, (lⁿ)⁻¹`(**格子の元**にわたる) |
| `EisensteinSeries.eisSummand` | `(v₀·z + v₁)^{−k}`(**整数の対**にわたる) |

★`τ` から周期対 `⟨τ, 1⟩` を作り、`latticeEquivProd : L.lattice ≃ₗ[ℤ] ℤ × ℤ` と
`finTwoArrowEquiv : (Fin 2 → ℤ) ≃ ℤ × ℤ` を噛ませると**同じ和になる**。
★★`l = 0` の項は `(0ⁿ)⁻¹ = 0` なので勝手に落ちる——除外の細工は要らなかった。

### ★★★★★★そこから q 展開まで一直線

    G_k = ∑' v, eisSummand k v τ            (本ブロック)
        = ζ(k) · eisensteinSeries 0 k τ      (mathlib)
        = 2ζ(k) · E_k(τ)                     (E は 1/2 倍で定義されている)
        = 2ζ(k) · (1 − (2k/B_k)·∑ σ_{k−1}(n)qⁿ)   (mathlib)

★★`g₂ = 60·G₄`・`g₃ = 140·G₆` なので、そのまま `σ₃`・`σ₅` の級数で書ける
——**これは我々の `tateA4 = −5·sigmaSeries 3`・`tateA6` と同じ形**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tauPair` | `τ` から周期対 `⟨τ, 1⟩` |
| `G_eq_tsum_eisSummand` | ★★★★**格子和は Eisenstein 和である** |
| `G_eq_two_zeta_mul_E` | ★★★★★`G_k = 2ζ(k)·E_k` |
| `G_qExpansion` | ★★★★★★**`G_k` の q 展開** |
| `g2_qExpansion` / `g3_qExpansion` | ★★★★★★**`g₂, g₃` の q 展開** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair EisensteinSeries ModularForm

/-- ★上半平面の点から周期対 `⟨τ, 1⟩` を作る。 -/
noncomputable def tauPair (τ : UpperHalfPlane) : PeriodPair where
  ω₁ := (τ : ℂ)
  ω₂ := 1
  indep := by
    rw [linearIndependent_fin2]
    refine ⟨by simp, ?_⟩
    intro a h
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one] at h
    have him : ((τ : ℂ)).im = 0 := by
      rw [← h]
      simp
    exact absurd him (ne_of_gt τ.2)

/-- ★★★★**格子和 `G_n` は Eisenstein 和である**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`l = 0` の項は `(0ⁿ)⁻¹ = 0` なので勝手に落ちる。 -/
theorem G_eq_tsum_eisSummand (τ : UpperHalfPlane) (n : ℕ) :
    (tauPair τ).G n = ∑' v : Fin 2 → ℤ, eisSummand n v τ := by
  have h1 := ((tauPair τ).latticeEquivProd.toEquiv.symm).tsum_eq
    (fun l : (tauPair τ).lattice => ((l : ℂ) ^ n)⁻¹)
  have h2 := (finTwoArrowEquiv ℤ).tsum_eq
    (fun x : ℤ × ℤ => (((x.1 : ℂ) * (τ : ℂ) + (x.2 : ℂ)) ^ n)⁻¹)
  rw [PeriodPair.G, ← h1]
  have hpt : ∀ x : ℤ × ℤ,
      ((((tauPair τ).latticeEquivProd.toEquiv.symm x : (tauPair τ).lattice) : ℂ) ^ n)⁻¹
        = (((x.1 : ℂ) * (τ : ℂ) + (x.2 : ℂ)) ^ n)⁻¹ := by
    intro x
    congr 2
    have heq : (tauPair τ).latticeEquivProd.toEquiv.symm x
        = (tauPair τ).latticeEquivProd.symm x := rfl
    rw [heq]
    simpa [tauPair] using (tauPair τ).latticeEquiv_symm_apply x
  rw [tsum_congr hpt, ← h2]
  refine tsum_congr fun v => ?_
  rw [eisSummand]
  simp only [finTwoArrowEquiv_apply, piFinTwoEquiv]
  rw [zpow_neg, zpow_natCast]

/-- ★★★★★**`G_k` は `2ζ(k)·E_k` である**。

★`ModularForm.E` は `(1/2) • eisensteinSeriesSIF` として定義されているので 2 が出る。 -/
theorem G_eq_two_zeta_mul_E (τ : UpperHalfPlane) {k : ℕ} (hk : 3 ≤ k) :
    (tauPair τ).G k = 2 * riemannZeta k * ModularForm.E hk τ := by
  rw [G_eq_tsum_eisSummand τ k, tsum_eisSummand_eq_riemannZeta_mul_eisensteinSeries hk τ]
  have hE : ModularForm.E hk τ = (1 / 2 : ℂ) • eisensteinSeriesSIF (N := 1) 0 k τ := rfl
  rw [hE, eisensteinSeriesSIF_apply]
  simp only [smul_eq_mul]
  ring

/-- ★★★★★★**`G_k` の q 展開**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem G_qExpansion (τ : UpperHalfPlane) {k : ℕ} (hk : 3 ≤ k) (hk2 : Even k) :
    (tauPair τ).G k = 2 * riemannZeta k *
      (1 - (2 * k / bernoulli k) *
        ∑' n : ℕ+, ((ArithmeticFunction.sigma (k - 1) (n : ℕ) : ℕ) : ℂ)
          * Complex.exp (2 * ↑π * I * τ) ^ (n : ℤ)) := by
  rw [G_eq_two_zeta_mul_E τ hk, EisensteinSeries.q_expansion_bernoulli hk hk2 τ]

/-- ★★★★★★**`g₂` の q 展開**——`σ₃` の級数で書ける。 -/
theorem g2_qExpansion (τ : UpperHalfPlane) :
    (tauPair τ).g₂ = 60 * (2 * riemannZeta 4 *
      (1 - (8 / bernoulli 4) *
        ∑' n : ℕ+, ((ArithmeticFunction.sigma 3 (n : ℕ) : ℕ) : ℂ)
          * Complex.exp (2 * ↑π * I * τ) ^ (n : ℤ))) := by
  rw [PeriodPair.g₂, G_qExpansion τ (k := 4) (by norm_num) (by decide)]
  norm_num

/-- ★★★★★★**`g₃` の q 展開**——`σ₅` の級数で書ける。 -/
theorem g3_qExpansion (τ : UpperHalfPlane) :
    (tauPair τ).g₃ = 140 * (2 * riemannZeta 6 *
      (1 - (12 / bernoulli 6) *
        ∑' n : ℕ+, ((ArithmeticFunction.sigma 5 (n : ℕ) : ℕ) : ℂ)
          * Complex.exp (2 * ↑π * I * τ) ^ (n : ℤ))) := by
  rw [PeriodPair.g₃, G_qExpansion τ (k := 6) (by norm_num) (by decide)]
  norm_num

/-! ## ★出典の紐付け(`.src`) -/

def G_eq_tsum_eisSummand.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子和が Eisenstein 和であること)",
    sectionId := "genell-def-3-3" }

def g2_qExpansion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子不変量 g2 の q 展開)",
    sectionId := "genell-def-3-3" }

def g3_qExpansion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子不変量 g3 の q 展開)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
