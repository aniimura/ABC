import ABC3.Found.PGC.LubinTateActionAdd

/-!
# `𝒪_K` 作用への拡張: `[π^n]_f` は `f` の `n` 回自己合成(`sorry` 無し)

節目(3)(捩れ点の構成 `Λ_n:=ker[π^n]_f`)へ向けた最初の一歩。古典的な
Lubin-Tate 理論の基本事実——`[π]_f = f` 自身、`[π^n]_f` は `f` の
`n` 回自己合成——を、既に確立済みの乗法側の環準同型性
(`Found/PGC/LubinTateActionMul.lean::LubinTateAction_comp`)から
形式的に導く。

## 証明の筋

`[π]_f = f` は、`f` 自身が(a) 次数1の係数 `π`、(b) 自分自身との
関数等式(自明)を満たすことから、1変数の一意性補題で `LubinTateAction
(π) = f` と結論する——`LubinTateAction` の定義そのものが「この2条件を
満たす唯一の冪級数」だから、ほとんど同語反復に近い。

`[1]_f = X`(恒等冪級数)も同様——`X` は次数1の係数1・自分自身との
自明な関数等式(`X.subst(f)=f=f.subst(X)`)を満たす。

`[π^(n+1)]_f = f([π^n]_f)` は `LubinTateAction_comp`(`[ab]_f=[a]_f∘
[b]_f`)を `a:=π^n,b:=π` に適用し、`[π]_f=f` で置き換え、`[π^n]_f` が
`f` と可換であること(`LubinTateAction_functional_equation`)で
「`f`を右から合成」から「`f`を左から合成」へ向きを変えるだけで出る。
これは `n` 回自己合成の関数 `iteratedLubinTate`(`Nat.rec` による定義、
`f^{(0)}:=X`・`f^{(n+1)}:=f^{(n)}.subst(f)`)とちょうど一致する漸化式
なので、帰納法で `LubinTateAction(π^n) = iteratedLubinTate f n` が
出る——新しい次数ごとの構成は一切不要だった。
-/

namespace ABC3.Found.PGC

/-! ### 部品0: `f` の `n` 回自己合成 -/

/-- `f` の `n` 回自己合成: `f^{(0)}:=X`・`f^{(n+1)}:=f^{(n)}.subst(f)`
(すなわち `f(f^{(n)})`)。 -/
noncomputable def iteratedLubinTate {A : Type*} [CommRing A] (f : PowerSeries A) : ℕ → PowerSeries A
  | 0 => PowerSeries.X
  | n + 1 => PowerSeries.subst (iteratedLubinTate f n) f

/-! ### 部品1: `[π]_f = f`・`[1]_f = X` -/

/-- ★★★★★**`[π]_f = f`**——`f` 自身が「次数1の係数 `π`・自分自身との
(自明な)関数等式」を満たすので、`LubinTateAction` の定義する唯一性
から直ちに従う。 -/
theorem LubinTateAction_pi_eq_f {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) :
    LubinTateAction hq hπmax f hf0 hf1 hf π = f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  apply powerSeries_uniqueness hπmax hπne0 hf0c hf1
    (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf π) hf0c
  · rw [coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf π, hf1]
  · exact (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf π).symm
  · rfl

/-- **`[1]_f = X`**(恒等冪級数)——`X` が「次数1の係数1・自分自身との
(自明な)関数等式 `X.subst(f)=f=f.subst(X)`」を満たすことから。 -/
theorem LubinTateAction_one_eq_X {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) :
    LubinTateAction hq hπmax f hf0 hf1 hf 1 = PowerSeries.X := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hfHS : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0c]; exact IsNilpotent.zero
  apply powerSeries_uniqueness hπmax hπne0 hf0c hf1
    (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf 1) PowerSeries.constantCoeff_X
  · rw [coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf 1]
    simp
  · exact (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf 1).symm
  · rw [PowerSeries.subst_X hfHS, PowerSeries.X_subst]

/-! ### 部品2: `[π^n]_f = f` の `n` 回自己合成 -/

/-- ★★★★★★★**`[π^n]_f` は `f` の `n` 回自己合成**——古典的な Lubin-Tate
理論の基本事実(捩れ点 `Λ_n:=ker[π^n]_f` の定義の出発点)。
`LubinTateAction_comp`・`LubinTateAction_pi_eq_f`・`LubinTateAction_
functional_equation`(可換性)だけから帰納法で出る。 -/
theorem LubinTateAction_pi_pow {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n) = iteratedLubinTate f n := by
  induction n with
  | zero =>
      show LubinTateAction hq hπmax f hf0 hf1 hf (π ^ 0) = PowerSeries.X
      rw [pow_zero]
      exact LubinTateAction_one_eq_X hq hπmax hπne0 f hf0 hf1 hf
  | succ n ih =>
    rw [pow_succ, LubinTateAction_comp hq hπmax hπne0 f hf0 hf1 hf (π ^ n) π,
      LubinTateAction_pi_eq_f hq hπmax hπne0 f hf0 hf1 hf,
      (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf (π ^ n)).symm, ih]
    rfl

end ABC3.Found.PGC
