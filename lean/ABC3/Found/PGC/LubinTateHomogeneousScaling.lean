import ABC3.Found.PGC.LubinTateLinearization

/-!
# Lubin-Tate 形式群の存在補題: g 側の準備(斉次式のスケーリング、`sorry` 無し)

f 側の線形化(`Found/PGC/LubinTateLinearization.lean`)に続き、g 側
(`Φ.subst(g,g)` の変化分 `φ.subst(g,g)`)の計算に必要な最初の部品——
**斉次式を定数倍した変数へ代入すると、その定数のべきがそのまま前に出る**
——を確立する。

## ★見つけた近道

自分で「有限積からスカラーを汲み出す」補題を手で組もうとしたところ、
mathlib に**すでにそのものずばりがあった**: `MvPowerSeries.
rescale_homogeneous_eq_smul`(`RingTheory/MvPowerSeries/Substitution.lean`)——
`rescale (const σ r) f = r^n • f`(`f` が次数 `n` の斉次式のとき)。
`rescale_eq_subst`(`rescale a f = subst (a•X) f`)と合わせるだけで
`homogeneous_subst_const_smul'` が数行で出た——config手で `Finsupp.prod` の
スカラー汲み出しを組もうとして苦労した後だったので、「名前を知らないと
たどり着けない」の実例がまた1つ増えた(`expand`/`homogeneousComponent` に
続く3例目、`memory` に記録する価値がある)。

## 続報(2026-09-04): g 側の線形化は完成した

`g ≡ πX (mod deg 2)` のとき `φ.subst(g,g)`(`φ` が次数 `n+1` の斉次式)を
`π^{n+1}•φ` で近似する誤差を次数 `n+2` 以上に押し出す段は、
`Found/PGC/LubinTateGLinearization.lean::coeff_subst_g_linearize` として
sorry 無しで完成した。当初の3段階計画通り: (1) 一般化した
`order_pow_sub_pow_ge'`(`k≥1` すべてで `order(a-a')` を変数のまま持ち回る
版、`LubinTateLinearization.lean::order_pow_sub_pow_ge` の `k≥2` 限定版の
一般化)、(2) 2変数の telescoping `order_prod_pow_sub_prod_pow_ge`
(`a_0^{d_0}a_1^{d_1}−a_0'^{d_0}a_1'^{d_1}` の分解)、(3) それを
`coeff_subst` の有限和に適用する `coeff_ad_eq_of_degree` を経て、
`coeff_subst_g_linearize` 本体が数行で出た——`coeff_subst_linearize`(f 側)
と並行の構造。これで帰納法の1ステップに要る両側の線形化(f・g)が
どちらも sorry 無しで揃った。
-/

namespace ABC3.Found.PGC

/-- ★★**斉次式を定数倍した変数へ代入すると、定数のべきが前に出る**。`φ` が
`Found/PGC/LubinTateLinearization.lean` と同じ「係数の消失」による次数 `p`
の斉次性を持つとき、`subst (c•X) φ = c^p • φ`。`MvPowerSeries.
rescale_homogeneous_eq_smul` + `rescale_eq_subst` から数行で出る。 -/
theorem homogeneous_subst_const_smul {A : Type*} [CommRing A] {φ : MvPowerSeries (Fin 2) A} {p : ℕ}
    (hφ : ∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ p → MvPowerSeries.coeff d φ = 0) (c : A) :
    MvPowerSeries.subst (fun i => c • (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)) φ = c ^ p • φ := by
  have hsupp : ∀ d ∈ (φ : MvPowerSeries (Fin 2) A).support, Finsupp.degree d = p := by
    intro d hd
    by_contra hne
    exact hd (hφ d hne)
  have h1 : MvPowerSeries.rescale (Function.const (Fin 2) c) φ = c ^ p • φ :=
    MvPowerSeries.rescale_homogeneous_eq_smul hsupp
  rw [← h1, MvPowerSeries.rescale_eq_subst]
  congr 1

end ABC3.Found.PGC
