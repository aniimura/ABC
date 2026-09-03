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

## まだ無いもの(次に戻るときの具体的な計画)

`g ≡ πX (mod deg 2)` のとき `φ.subst(g,g)`(`φ` が次数 `n+1` の斉次式)を
`π^{n+1}•φ` で近似する誤差を次数 `n+2` 以上に押し出す段は、まだ残る。
手計算で以下まで詰めた:

1. **一般化した `order_pow_sub_pow_ge`**: `order(a),order(a')≥1` かつ
   `order(a-a')≥2` のとき、`order(a^k-a'^k) ≥ order(a-a')+(k-1)`
   (`k≥0` すべてで成立、`k=0` は自明、`k≥1` は `geom_sum₂_mul` の
   因数分解 `a^k-a'^k=(a-a')·Σaⁱa'^{k-1-i}` + 各項の次数 `≥k-1` から)。
   これは `LubinTateLinearization.lean::order_pow_sub_pow_ge`(`k≥2` 限定・
   `φ` 側の次数だけ使う版)の一般化にあたる。
2. **2変数の telescoping**: `a_0^{d_0}a_1^{d_1} − a_0'^{d_0}a_1'^{d_1}
   = (a_0^{d_0}−a_0'^{d_0})·a_1^{d_1} + a_0'^{d_0}·(a_1^{d_1}−a_1'^{d_1})`
   の各項に部品1を当てると、`d_0+d_1=n+1` のとき両方とも次数 `≥n+2`。
3. これを `MvPowerSeries.coeff_subst` の有限和(`φ` が次数 `n+1` の斉次式
   なので `d` は `degree d=n+1` の項しか効かない)に適用すれば、
   `coeff e (φ.subst(g,g)) = π^{n+1}·coeff e φ`(次数 `e≤n+1` の範囲)が出る
   ——`coeff_subst_linearize`(f 側)と並行の構造。
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
