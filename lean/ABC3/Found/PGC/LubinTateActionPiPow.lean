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

/-! ### 部品3: `iteratedLubinTate f n` は mod `π` で `X^(q^n)` -/

/-- ★★★★★★★**`[π^n]_f` は mod `π` で `X^(q^n)`**——古典的な Lubin-Tate
理論の基本事実(`n=1` の場合が仮定 `hf` そのもの)。捩れ点 `Λ_n` の
議論で `[π^n]_f` を distinguished 多項式として扱う際の出発点になる
見込み(`RingTheory/PowerSeries/WeierstrassPreparation.lean` の
`PowerSeries.IsWeierstrassFactorization` 等)。`PowerSeries.map_subst`
(写像は代入と可換)を軸に、`f` の mod `π` での姿(仮定 `hf`)を
`n` 回適用するだけの帰納法——新しい構成は一切不要だった。 -/
theorem iteratedLubinTate_map_residue {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    PowerSeries.map (IsLocalRing.residue A) (iteratedLubinTate f n) = PowerSeries.X ^ (pp ^ ff) ^ n := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hppffpos : 0 < pp ^ ff := hq ▸ Fintype.card_pos
  induction n with
  | zero =>
      show PowerSeries.map (IsLocalRing.residue A) PowerSeries.X = PowerSeries.X ^ (pp ^ ff) ^ 0
      rw [pow_zero, pow_one, PowerSeries.map_X]
  | succ n ih =>
      have hiter0 : PowerSeries.constantCoeff (iteratedLubinTate f n) = 0 := by
        rw [← LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf n]
        exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)
      have hvHS : PowerSeries.HasSubst (iteratedLubinTate f n) := by
        show IsNilpotent (PowerSeries.constantCoeff (iteratedLubinTate f n))
        rw [hiter0]; exact IsNilpotent.zero
      have hstep : PowerSeries.map (IsLocalRing.residue A) (PowerSeries.subst (iteratedLubinTate f n) f) =
          PowerSeries.subst (PowerSeries.map (IsLocalRing.residue A) (iteratedLubinTate f n))
            (PowerSeries.map (IsLocalRing.residue A) f) :=
        PowerSeries.map_subst hvHS f
      show PowerSeries.map (IsLocalRing.residue A) (PowerSeries.subst (iteratedLubinTate f n) f) =
        PowerSeries.X ^ (pp ^ ff) ^ (n + 1)
      rw [hstep, hf, ih]
      have hXHS : PowerSeries.HasSubst
          (PowerSeries.X ^ (pp ^ ff) ^ n : PowerSeries (IsLocalRing.ResidueField A)) := by
        show IsNilpotent (PowerSeries.constantCoeff
          (PowerSeries.X ^ (pp ^ ff) ^ n : PowerSeries (IsLocalRing.ResidueField A)))
        rw [map_pow, PowerSeries.constantCoeff_X, zero_pow (pow_ne_zero n hppffpos.ne')]
        exact IsNilpotent.zero
      rw [substXpow_eq_pow hXHS, ← pow_mul, pow_succ]

/-! ### 部品4: `[π^n]_f ∣ [π^m]_f`(`n≤m`)——`Λ_n⊆Λ_m` の由来 -/

/-- `[π^n]_f` 自身の定数項は0(`LubinTateAction_pi_pow` 経由で
`constantCoeff_LubinTateAction` に帰着)。 -/
theorem constantCoeff_iteratedLubinTate {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    PowerSeries.constantCoeff (iteratedLubinTate f n) = 0 := by
  rw [← LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf n]
  exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)

/-- **`[π^{a+b}]_f = [π^a]_f∘[π^b]_f`**(合成として)——`LubinTateAction_
comp` を `π^a・π^b=π^{a+b}` に適用するだけ。 -/
theorem iteratedLubinTate_add {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : ℕ) :
    iteratedLubinTate f (a + b) =
      PowerSeries.subst (iteratedLubinTate f b) (iteratedLubinTate f a) := by
  rw [← LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf (a + b),
    ← LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf a,
    ← LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf b,
    pow_add, LubinTateAction_comp hq hπmax hπne0 f hf0 hf1 hf (π ^ a) (π ^ b)]

/-- ★★★★★★★**`[π^b]_f ∣ [π^{a+b}]_f`**——古典的な Lubin-Tate 理論で
`Λ_b⊆Λ_{a+b}`(捩れ点の包含)の由来になる事実。`[π^a]_f` の定数項が0
(`X∣[π^a]_f`)であることと `iteratedLubinTate_add` を組み合わせ、
`subst` が積を保つ(`PowerSeries.subst_mul`)ことから従う——新しい
次数ごとの構成は一切不要だった。 -/
theorem iteratedLubinTate_dvd_iteratedLubinTate_add {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : ℕ) :
    iteratedLubinTate f b ∣ iteratedLubinTate f (a + b) := by
  obtain ⟨ra, hra⟩ := PowerSeries.X_dvd_iff.mpr (constantCoeff_iteratedLubinTate hq hπmax hπne0 f hf0 hf1 hf a)
  have hbHS : PowerSeries.HasSubst (iteratedLubinTate f b) := by
    show IsNilpotent (PowerSeries.constantCoeff (iteratedLubinTate f b))
    rw [constantCoeff_iteratedLubinTate hq hπmax hπne0 f hf0 hf1 hf b]
    exact IsNilpotent.zero
  refine ⟨PowerSeries.subst (iteratedLubinTate f b) ra, ?_⟩
  rw [iteratedLubinTate_add hq hπmax hπne0 f hf0 hf1 hf a b, hra,
    PowerSeries.subst_mul hbHS, PowerSeries.subst_X hbHS]

/-- `n≤m` ならば `[π^n]_f ∣ [π^m]_f`——`iteratedLubinTate_dvd_
iteratedLubinTate_add` を `m=k+n` の形に書き直すだけ。 -/
theorem iteratedLubinTate_dvd_of_le {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) {n m : ℕ} (hnm : n ≤ m) :
    iteratedLubinTate f n ∣ iteratedLubinTate f m := by
  obtain ⟨k, hk⟩ := Nat.le.dest hnm
  rw [← hk, add_comm]
  exact iteratedLubinTate_dvd_iteratedLubinTate_add hq hπmax hπne0 f hf0 hf1 hf k n

end ABC3.Found.PGC
