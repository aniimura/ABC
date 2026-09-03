import ABC3.Found.PGC.FrobeniusExpand
import Mathlib.RingTheory.MvPowerSeries.Substitution
import Mathlib.RingTheory.PowerSeries.Substitution

/-!
# Lubin-Tate 形式群の存在補題の可除性(`sorry` 無し、一般の有限体で)

Lubin-Tate 形式群の構成(`ResearchPaper/blocked-leaves.json` の
`[pGC] Proposition 1.2 / ...` で局所類体論の相互律への迂回路として検討している)
の**存在補題**の核心にある可除性——次数ごとの再帰構成で現れる剰余項
`R_n := Φ_n(g,…,g) − f(Φ_n)` が一意化元 `π` で割り切れること——を、
剰余体への還元という形で一般に確立する。

## 数学的な内容

`f, g ∈ 𝔉_π`(`f ≡ g ≡ X^q (mod π)`、`q` は剰余体の元の個数)とすると、
剰余体 `κ = A/π` 上で `f̄(T) = T^q`・`ḡ(X) = X^q` となる。任意の(次数ごとの
構成の途中経過である)`Φ`(定数項 0)について:

```
R̄_n = Φ̄(ḡ,ḡ) − f̄(Φ̄) = Φ̄(X^q, X^q) − Φ̄^q
```

ここで `Φ̄(X^q, X^q)` は「`Φ̄` の各変数を `q` 乗で置き換えたもの」——これは
`Found/PGC/FrobeniusExpand.lean::mvPowerSeries_pow_card_eq_expand`
(有限体上の冪級数の Frobenius 恒等式、`φ^q = expand(q)(φ)`)が
**`Φ̄` が何であっても**そのまま与える `Φ̄^q` に一致する。したがって
`R̄_n = Φ̄^q − Φ̄^q = 0` ——**`Φ` の構成の仕方に依存する帰納的な不変量は
一切不要**で、可除性は既存の1個の補題から即座に出る
(`ResearchPaper/blocked-leaves.json` の `progress2026_09_04e` で最初に
気づいた簡略化を、実際に形式化したもの)。

## 使った mathlib の道具(実測、2026-09-04)

- `MvPowerSeries.map_subst` / 本ファイルの `powerSeries_map_subst`
  (1変数版、`MvPowerSeries.map_subst` から `PowerSeries.subst_def` 経由で導出)——
  還元写像 `MvPowerSeries.map` が `subst` と可換であること。
- `MvPowerSeries.subst_self`・`expand_subst`・`expand_X` から
  本ファイルの `expand_eq_subst_pow`(`expand p hp φ = subst (X^p) φ`)を導出——
  mathlib に直接の形では無かった。
- `Found/PGC/FrobeniusExpand.lean::mvPowerSeries_pow_card_eq_expand`
  (本セッション既存)。

## まだ無いもの

本ファイルは可除性(`R_n ≡ 0 (mod π)` の剰余体版)のみを確立する。次数ごとの
係数の**再帰的構成そのもの**(`φ_{n+1} := R_n/(π−π^{n+1})` を実際に定義し、
`Nat.rec` で組み立て、最終的な `F` が関数等式を exact に満たすことを示す)は
別途の課題として残る。
-/

namespace ABC3.Found.PGC

open MvPowerSeries

/-! ### 部品1: `PowerSeries.subst` の還元可換性(1変数版、`MvPowerSeries.map_subst` から) -/

/-- `MvPowerSeries.map` は `PowerSeries.subst`(1変数版)とも可換である。
`PowerSeries.subst a f = MvPowerSeries.subst (fun _ : Unit => a) f`
(`PowerSeries.subst_def`)を経由して `MvPowerSeries.map_subst` から導く。 -/
theorem powerSeries_map_subst {R S τ : Type*} [CommRing R] [CommRing S] {a : MvPowerSeries τ R}
    (ha : PowerSeries.HasSubst a) (h : R →+* S) (f : PowerSeries R) :
    (MvPowerSeries.map h) (PowerSeries.subst a f) =
      PowerSeries.subst ((MvPowerSeries.map h) a) ((MvPowerSeries.map h) f) := by
  have ha' : MvPowerSeries.HasSubst (fun _ : Unit => a) := by
    constructor
    · intro _; exact ha
    · intro d; exact Set.toFinite _
  rw [PowerSeries.subst_def, PowerSeries.subst_def]
  exact MvPowerSeries.map_subst ha' f

/-! ### 部品2: `expand` を `subst` で書き直す -/

/-- `expand p hp φ = subst (fun i ↦ X i ^ p) φ`。mathlib に直接の形では無かった
(`expand_subst` を `f := X`(恒等代入、`subst_self`)に特殊化して得る)。 -/
theorem expand_eq_subst_pow {σ R : Type*} [CommRing R] {p : ℕ} (hp : p ≠ 0)
    (φ : MvPowerSeries σ R) :
    expand p hp φ = subst (fun i : σ => (X i : MvPowerSeries σ R) ^ p) φ := by
  have h := expand_subst (τ := σ) (f := (X : σ → MvPowerSeries σ R)) HasSubst.X (p := p) (hp := hp)
    (φ := φ)
  rw [subst_self] at h
  simp only [id] at h
  rw [h]
  congr 1
  funext i
  exact expand_X p hp i

/-! ### 部品3: 補助事実(`X i` の `HasSubst`・`subst` による `X^q ↦ (⋅)^q`) -/

theorem hasSubst_X_i {S τ : Type*} [CommRing S] (i : τ) :
    PowerSeries.HasSubst (X i : MvPowerSeries τ S) := by
  show IsNilpotent (constantCoeff (X i : MvPowerSeries τ S))
  rw [MvPowerSeries.constantCoeff_X]
  exact IsNilpotent.zero

theorem substXpow_eq_pow {R S τ : Type*} [CommRing R] [CommRing S] [Algebra R S] {q : ℕ}
    {b : MvPowerSeries τ S} (hb : PowerSeries.HasSubst b) :
    PowerSeries.subst b (PowerSeries.X ^ q : PowerSeries R) = b ^ q := by
  rw [PowerSeries.subst_pow hb PowerSeries.X q, PowerSeries.subst_X hb]

/-! ### ★★★可除性の核心 -/

/-- ★★★**Lubin-Tate 存在補題の可除性(剰余体版)**。`f, g` がどちらも剰余体上で
`X^q`(`q` = 剰余体の元の個数)に還元されるとき、`Φ`(定数項 0、任意)について
`Φ(g,g) − f(Φ)` の剰余体への還元は**恒等的に 0** ——`Φ` の構成の仕方に依存する
帰納的な不変量は一切不要で、`Found/PGC/FrobeniusExpand.lean` の
Frobenius 恒等式1個から即座に出る。

Lubin-Tate 形式群 `F` を次数ごとに再帰構成する際、`R_n := Φ_n(g,g) − f(Φ_n)` が
`π` で割り切れること(⟺ 剰余体への還元が 0)の核心はこれ——`Φ_n` が
どのように構成された途中経過であっても成り立つ。 -/
theorem residue_divides_R {A κ : Type*} [CommRing A] [Field κ] [Fintype κ] [Algebra A κ]
    {pp : ℕ} [ExpChar κ pp] {ff : ℕ} (hq : Fintype.card κ = pp ^ ff)
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0) (f : PowerSeries A)
    (hf : PowerSeries.map (algebraMap A κ) f = PowerSeries.X ^ (pp ^ ff))
    (hg : PowerSeries.map (algebraMap A κ) g = PowerSeries.X ^ (pp ^ ff))
    (Φ : MvPowerSeries (Fin 2) A) (hΦ0 : MvPowerSeries.constantCoeff Φ = 0) :
    MvPowerSeries.map (algebraMap A κ)
      (MvPowerSeries.subst (fun i => PowerSeries.subst (X i) g) Φ - PowerSeries.subst Φ f) = 0 := by
  set residue : A →+* κ := algebraMap A κ with hres
  set q := pp ^ ff with hqdef
  set a : Fin 2 → MvPowerSeries (Fin 2) A := fun i => PowerSeries.subst (X i) g with ha_def
  have haM : MvPowerSeries.HasSubst a := by
    constructor
    · intro i
      show IsNilpotent (constantCoeff (a i))
      have heq0 : constantCoeff (a i) = 0 :=
        PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i) g hg0
      rw [heq0]; exact IsNilpotent.zero
    · intro d; exact Set.toFinite _
  set Φbar := MvPowerSeries.map residue Φ with hΦbar
  have hΦHasSubst : PowerSeries.HasSubst Φ := by
    show IsNilpotent (constantCoeff Φ); rw [hΦ0]; exact IsNilpotent.zero
  have hΦbarHasSubst : PowerSeries.HasSubst Φbar := by
    show IsNilpotent (constantCoeff Φbar)
    rw [hΦbar, MvPowerSeries.constantCoeff_map, hΦ0, map_zero]
    exact IsNilpotent.zero
  rw [map_sub]
  have hai : ∀ i : Fin 2,
      MvPowerSeries.map residue (a i) = (X i : MvPowerSeries (Fin 2) κ) ^ q := by
    intro i
    show MvPowerSeries.map residue (PowerSeries.subst (X i) g) = _
    rw [powerSeries_map_subst (hasSubst_X_i i) residue g]
    have hgi : MvPowerSeries.map residue (X i : MvPowerSeries (Fin 2) A)
        = (X i : MvPowerSeries (Fin 2) κ) := MvPowerSeries.map_X residue i
    rw [hgi]
    show PowerSeries.subst (X i) (PowerSeries.map residue g) = _
    rw [hg]
    exact substXpow_eq_pow (hasSubst_X_i i)
  have hterm1 : MvPowerSeries.map residue (MvPowerSeries.subst a Φ) = Φbar ^ q := by
    rw [MvPowerSeries.map_subst haM]
    have hfun : (fun i => MvPowerSeries.map residue (a i))
        = (fun i : Fin 2 => (X i : MvPowerSeries (Fin 2) κ) ^ q) := funext hai
    rw [hfun, ← expand_eq_subst_pow (pow_ne_zero ff (expChar_ne_zero κ pp))]
    exact (mvPowerSeries_pow_card_eq_expand hq Φbar).symm
  have hterm2 : MvPowerSeries.map residue (PowerSeries.subst Φ f) = Φbar ^ q := by
    rw [powerSeries_map_subst hΦHasSubst residue f]
    show PowerSeries.subst Φbar (PowerSeries.map residue f) = _
    rw [hf]
    exact substXpow_eq_pow hΦbarHasSubst
  rw [hterm1, hterm2, sub_self]

end ABC3.Found.PGC
