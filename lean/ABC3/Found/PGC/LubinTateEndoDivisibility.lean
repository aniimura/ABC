import ABC3.Found.PGC.LubinTateDivisibility
import Mathlib.RingTheory.PowerSeries.Expand

/-!
# Lubin-Tate 形式群の`𝒪_K`作用への拡張: 1変数版の可除性(`sorry` 無し)

`ResearchPaper/blocked-leaves.json` の `progress2026_09_04q`/`r` で分解した
節目(1)——各 `a∈𝒪_K` について `[a]_f(X)≡aX(mod deg2)` かつ `f∘[a]_f=[a]_f∘f`
(合成として可換)を満たす冪級数 `[a]_f` の存在——に向けた最初の一歩。

この存在補題は `Found/PGC/LubinTateF` 系列(2変数の `Φ.subst(g,g)=f.subst(Φ)`)
と同じ次数ごとの `Nat.rec` 近似列の型だが、対象が **1変数** の冪級数 `φ`
(基底段が `X` ではなく `aX`)であり、満たすべき等式も `f.subst(φ)=φ.subst(g)`
(`f(φ(X))=φ(g(X))`)という**1変数どうしの合成**である点が違う。本ファイルは
その**可除性**(次数ごとの再帰で現れる剰余項 `R_n:=f.subst(φ_n)−φ_n.subst(g)`
が一意化元 `π` で割り切れること、剰余体への還元が0であること)を確立する
——2変数版(`Found/PGC/LubinTateDivisibility.lean::residue_divides_R`)と
全く同じ議論(有限体上の Frobenius 恒等式)が、`PowerSeries R` が
`MvPowerSeries Unit R` の別名(`abbrev`)であることを利用してほぼそのまま
移植できた。

## 鍵になった発見

- `PowerSeries.expand`(mathlib、`RingTheory/PowerSeries/Expand.lean`)は
  **定義そのものが** `MvPowerSeries.expand`(`fun p hp => MvPowerSeries.expand
  p hp`)——`PowerSeries.map`・`PowerSeries.constantCoeff` も同様に
  `MvPowerSeries` の対応物として直接定義されている。これにより
  `Found/PGC/FrobeniusExpand.lean::mvPowerSeries_pow_card_eq_expand`
  (`σ` について一般)を `σ:=Unit` として**そのまま**適用でき、1変数版の
  Frobenius 恒等式 `powerSeries_pow_card_eq_expand` は実質1行で済んだ。
- `Found/PGC/LubinTateDivisibility.lean` の補助補題(`powerSeries_map_subst`・
  `expand_eq_subst_pow`・`hasSubst_X_i`・`substXpow_eq_pow`)は元から一般
  (`τ` について一般、`Fin 2` に固定されていない)だったので、そのまま
  再利用できた——`residue_divides_R` 自身だけが `MvPowerSeries (Fin 2) A`
  に固定されていた。

## まだ無いもの

可除性(1変数版)のみを確立する。両側の線形化・1ステップの組み立て・
近似列・極限という、2変数版と同型だが1変数へ作り直す必要がある残りの
工程は別途の課題として残る(`ResearchPaper/blocked-leaves.json` の
`progress2026_09_04r` 参照)。
-/

namespace ABC3.Found.PGC

/-- **有限体上の1変数冪級数の Frobenius 恒等式**——`mvPowerSeries_pow_card_
eq_expand` を `σ:=Unit`(`PowerSeries F = MvPowerSeries Unit F` という
`abbrev`)としてそのまま適用しただけ。 -/
theorem powerSeries_pow_card_eq_expand {F : Type*} [Field F] [Fintype F] {p : ℕ}
    [ExpChar F p] {f : ℕ} (hq : Fintype.card F = p ^ f) (φ : PowerSeries F) :
    φ ^ (p ^ f) = PowerSeries.expand (p ^ f) (pow_ne_zero f (expChar_ne_zero F p)) φ :=
  mvPowerSeries_pow_card_eq_expand hq φ

/-- ★★★**`𝒪_K` 作用への拡張(1変数版)の可除性(剰余体版)**。`f, g` が
どちらも剰余体上で `X^q` に還元されるとき、任意の(定数項0の)1変数冪級数
`φ` について `f.subst(φ) − φ.subst(g)`(`f(φ(X))−φ(g(X))`)の剰余体への
還元は恒等的に0——2変数版(`residue_divides_R`)と全く同じ議論
(`φ̄^q = φ̄.subst(X^q) = expand(q)(φ̄)` という Frobenius 恒等式)がそのまま
使える。`φ` の構成の仕方に依存する帰納的な不変量は不要。 -/
theorem residue_divides_R_endo {A κ : Type*} [CommRing A] [Field κ] [Fintype κ] [Algebra A κ]
    {pp : ℕ} [ExpChar κ pp] {ff : ℕ} (hq : Fintype.card κ = pp ^ ff)
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0) (f : PowerSeries A)
    (hf : PowerSeries.map (algebraMap A κ) f = PowerSeries.X ^ (pp ^ ff))
    (hg : PowerSeries.map (algebraMap A κ) g = PowerSeries.X ^ (pp ^ ff))
    (φ : PowerSeries A) (hφ0 : PowerSeries.constantCoeff φ = 0) :
    PowerSeries.map (algebraMap A κ)
      (PowerSeries.subst φ f - PowerSeries.subst g φ) = 0 := by
  set residue : A →+* κ := algebraMap A κ with hres
  set q := pp ^ ff with hqdef
  set φbar := PowerSeries.map residue φ with hφbar
  have hφ0' : MvPowerSeries.constantCoeff φ = 0 := hφ0
  have hg0' : MvPowerSeries.constantCoeff g = 0 := hg0
  have hφHasSubst : PowerSeries.HasSubst φ := by
    show IsNilpotent (MvPowerSeries.constantCoeff φ); rw [hφ0']; exact IsNilpotent.zero
  have hgHasSubst : PowerSeries.HasSubst g := by
    show IsNilpotent (MvPowerSeries.constantCoeff g); rw [hg0']; exact IsNilpotent.zero
  have hφbarHasSubst : PowerSeries.HasSubst φbar := by
    show IsNilpotent (MvPowerSeries.constantCoeff (MvPowerSeries.map residue φ))
    rw [MvPowerSeries.constantCoeff_map, hφ0', map_zero]
    exact IsNilpotent.zero
  show MvPowerSeries.map residue (PowerSeries.subst φ f - PowerSeries.subst g φ) = 0
  rw [map_sub]
  have hterm1 : MvPowerSeries.map residue (PowerSeries.subst φ f) = φbar ^ q := by
    rw [powerSeries_map_subst hφHasSubst residue f]
    show PowerSeries.subst φbar (PowerSeries.map residue f) = _
    rw [hf]
    exact substXpow_eq_pow hφbarHasSubst
  have hterm2 : MvPowerSeries.map residue (PowerSeries.subst g φ) = φbar ^ q := by
    rw [powerSeries_map_subst hgHasSubst residue φ]
    show PowerSeries.subst (PowerSeries.map residue g) φbar = _
    rw [hg]
    rw [powerSeries_pow_card_eq_expand hq φbar]
    exact (expand_eq_subst_pow (pow_ne_zero ff (expChar_ne_zero κ pp)) φbar).symm
  rw [hterm1, hterm2, sub_self]

end ABC3.Found.PGC
