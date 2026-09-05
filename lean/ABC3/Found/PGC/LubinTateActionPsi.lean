import ABC3.Found.PGC.LubinTateActionWeierstrass

/-!
# `ψ_n` は既約(古典的 Lubin-Tate 理論の核心定理、`sorry` 無し)

`Found/PGC/LubinTateActionWeierstrass.lean::irreducible_iteratedLubinTatePrimitive_one`
(`n=1` 限定)を**任意の `n≥1`** へ一般化する。`ψ_n:=D_n/D_{n-1}`
(`n=1` の場合は `ψ_1=φ_1`、`D_0=X` なので)——原始 `π^n`-捩れ点
(`Λ_n\Λ_{n-1}`)を統べる、次数 `q^n-q^{n-1}` の多項式——が**既約**で
あることを示す。

## 証明の筋

`[π^n]_f = [π]_f∘[π^{n-1}]_f = f([π^{n-1}]_f(X))` と `f(T)=T・h(T)`
(`h:=f/X`、`h(0)=f'(0)=π` ちょうど、`f` の次数1の係数が `π` である
ことの言い換え)を組み合わせ、`[π^n]_f = [π^{n-1}]_f・r_n`
(`r_n:=h([π^{n-1}]_f)`、`r_n` の定数項は**ちょうど `π`**——`n=1` の
場合の `φ_1` の議論(`π^n・単元` という形)より単純)という分解を得る。
`r_n` を Weierstrass 分解し(`r_n=ψ_n・U'_n`、`ψ_n` は distinguished)、
`(D_{n-1}・ψ_n)・(U_{n-1}・U'_n)` が `[π^n]_f` のもう1つの Weierstrass
分解であることから一意性で `D_n=D_{n-1}・ψ_n` を得る。`ψ_n` は
distinguished(弱Eisenstein込み)なので、残るは `ψ_n` の定数項が
`𝔪²` に属さないことだけ——`ψ_n.coeff0・U'_n定数項=π`(`r_n` の定数項が
ちょうど `π` であることから)から、`n=1` の場合と全く同じ背理法
(属すと仮定すると `π` 自身が単元になる)で閉じる。
-/

namespace ABC3.Found.PGC

/-! ### 部品0: 代入と定数項の一般公式(`coeff_one_mul_eq` は既出、
`LubinTateActionWeierstrass.lean` から使う) -/

/-- **代入は定数項を保つ**(代入する値 `v` の定数項が0のとき):
`constantCoeff(subst v p) = constantCoeff p`——次数 `≥1` の項は
代入後の定数項に効かない、という標準的な次数勘定
(`coeff_one_subst_1var` の次数0版)。 -/
theorem constantCoeff_subst_1var {A : Type*} [CommRing A] {v p : PowerSeries A}
    (hv0 : PowerSeries.constantCoeff v = 0) :
    PowerSeries.constantCoeff (PowerSeries.subst v p) = PowerSeries.constantCoeff p := by
  have hvHS : PowerSeries.HasSubst v := by
    show IsNilpotent (PowerSeries.constantCoeff v); rw [hv0]; exact IsNilpotent.zero
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, ← PowerSeries.coeff_zero_eq_constantCoeff_apply]
  show MvPowerSeries.coeff (Finsupp.single () 0) (PowerSeries.subst v p) = _
  rw [PowerSeries.coeff_subst hvHS]
  rw [finsum_eq_sum_of_support_subset _ (s := ({0} : Finset ℕ)) (fun d hd => by
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_singleton, Set.mem_singleton_iff]
    by_contra hcon
    apply hd
    have hdpos : 0 < d := Nat.pos_of_ne_zero hcon
    have hz : MvPowerSeries.coeff (Finsupp.single () (0 : ℕ)) (v ^ d : PowerSeries A) = 0 := by
      have hvorder : (1 : ℕ∞) ≤ MvPowerSeries.order (v : PowerSeries A) :=
        MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hv0
      have hvdorder : ((d : ℕ) : ℕ∞) ≤ MvPowerSeries.order (v ^ d : PowerSeries A) := by
        calc ((d : ℕ) : ℕ∞) = d • (1 : ℕ∞) := by simp
          _ ≤ d • MvPowerSeries.order (v : PowerSeries A) := by gcongr
          _ ≤ _ := MvPowerSeries.le_order_pow d
      have hvdorder' : ((d : ℕ) : ℕ∞) ≤ PowerSeries.order (v ^ d : PowerSeries A) := by
        rw [PowerSeries.order_eq_order]; exact hvdorder
      have hpos : (0 : ℕ∞) < d := by exact_mod_cast hdpos
      exact PowerSeries.coeff_of_lt_order 0 (lt_of_lt_of_le hpos hvdorder')
    rw [hz, smul_zero])]
  rw [Finset.sum_singleton, pow_zero]
  show PowerSeries.coeff 0 p • MvPowerSeries.coeff (Finsupp.single () (0 : ℕ)) (1 : PowerSeries A) = _
  simp

/-- `f=X・r` のとき `r` の定数項は `f` の次数1の係数——`X` の低次係数
(`coeff0=0`・`coeff1=1`)から `coeff_one_mul_eq` を評価するだけ。 -/
theorem ra_constantCoeff_eq_coeff_one {A : Type*} [CommRing A]
    (g ra : PowerSeries A) (c : A) (hg1 : PowerSeries.coeff 1 g = c) (hra : g = PowerSeries.X * ra) :
    PowerSeries.constantCoeff ra = c := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
  have h1 := coeff_one_mul_eq (PowerSeries.X : PowerSeries A) ra
  rw [← hra, hg1] at h1
  have hX0 : PowerSeries.coeff 0 (PowerSeries.X : PowerSeries A) = 0 := by simp
  have hX1 : PowerSeries.coeff 1 (PowerSeries.X : PowerSeries A) = 1 := by simp
  rw [hX0, hX1, zero_mul, one_mul, zero_add] at h1
  exact h1.symm

/-! ### 部品1: `[π^n]_f = [π^{n-1}]_f・r_n`(`r_n` の定数項はちょうど `π`) -/

/-- `n≥1` のとき、`[π^n]_f = [π^{n-1}]_f・r_n`(`r_n` の定数項はちょうど
`π`)という分解が存在する——`f=X・r_a`(`r_a` の定数項は `f` の次数1
係数=`π`)を `[π^{n-1}]_f` へ代入して得る。 -/
theorem iteratedLubinTate_step_decomp {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    ∃ r : PowerSeries A, iteratedLubinTate f n = iteratedLubinTate f (n - 1) * r ∧
      PowerSeries.constantCoeff r = π := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  obtain ⟨ra, hra⟩ := PowerSeries.X_dvd_iff.mpr hf0c
  have hraconst : PowerSeries.constantCoeff ra = π := ra_constantCoeff_eq_coeff_one f ra π hf1 hra
  set gprev := iteratedLubinTate f (n - 1) with hgprev_def
  have hgprev0 : PowerSeries.constantCoeff gprev = 0 :=
    constantCoeff_iteratedLubinTate hq hπmax hπne0 f hf0 hf1 hf (n - 1)
  have hgprevHS : PowerSeries.HasSubst gprev := by
    show IsNilpotent (PowerSeries.constantCoeff gprev); rw [hgprev0]; exact IsNilpotent.zero
  refine ⟨PowerSeries.subst gprev ra, ?_, ?_⟩
  · have hstep : iteratedLubinTate f (1 + (n - 1)) =
        PowerSeries.subst gprev (iteratedLubinTate f 1) :=
      iteratedLubinTate_add hq hπmax hπne0 f hf0 hf1 hf 1 (n - 1)
    rw [show 1 + (n - 1) = n by omega, iteratedLubinTate_one_eq_f hq hπmax hπne0 f hf0 hf1 hf] at hstep
    rw [hstep, hra, PowerSeries.subst_mul hgprevHS, PowerSeries.subst_X hgprevHS]
  · rw [constantCoeff_subst_1var hgprev0, hraconst]

/-! ### 部品2: `r_n` の Weierstrass 分解と `D_n=D_{n-1}・ψ_n` -/

/-- `[π^n]_f=[π^{n-1}]_f・r` の商 `r`(いずれか)の mod `π` の像は非零
——`iteratedLubinTate_dvd_map_residue` を `n-1≤n` に適用するだけ。 -/
theorem iteratedLubinTate_step_map_residue_ne_zero {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ)
    (r : PowerSeries A) (hr : iteratedLubinTate f n = iteratedLubinTate f (n - 1) * r) :
    PowerSeries.map (IsLocalRing.residue A) r ≠ 0 := by
  have := iteratedLubinTate_dvd_map_residue hq hπmax hπne0 f hf0 hf1 hf (Nat.sub_le n 1) r hr
  rw [this]; exact pow_ne_zero _ PowerSeries.X_ne_zero

/-- `iteratedLubinTate_step_decomp` の商 `r`(選択関数による代表)。 -/
noncomputable def iteratedLubinTateStepR {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    PowerSeries A :=
  (iteratedLubinTate_step_decomp hq hπmax hπne0 f hf0 hf1 hf n hn).choose

theorem iteratedLubinTateStepR_eq {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    iteratedLubinTate f n = iteratedLubinTate f (n - 1) * iteratedLubinTateStepR hq hπmax hπne0 f hf0 hf1 hf n hn :=
  (iteratedLubinTate_step_decomp hq hπmax hπne0 f hf0 hf1 hf n hn).choose_spec.1

theorem iteratedLubinTateStepR_constantCoeff {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    PowerSeries.constantCoeff (iteratedLubinTateStepR hq hπmax hπne0 f hf0 hf1 hf n hn) = π :=
  (iteratedLubinTate_step_decomp hq hπmax hπne0 f hf0 hf1 hf n hn).choose_spec.2

theorem iteratedLubinTateStepR_map_residue_ne_zero {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    PowerSeries.map (IsLocalRing.residue A) (iteratedLubinTateStepR hq hπmax hπne0 f hf0 hf1 hf n hn) ≠ 0 :=
  iteratedLubinTate_step_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n _
    (iteratedLubinTateStepR_eq hq hπmax hπne0 f hf0 hf1 hf n hn)

/-- `r_n` の Weierstrass 標準分解の distinguished 多項式部分——
`ψ_n:=D_n/D_{n-1}`、原始 `π^n`-捩れ点を統べる多項式。 -/
noncomputable def iteratedLubinTatePsi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    Polynomial A :=
  (iteratedLubinTateStepR hq hπmax hπne0 f hf0 hf1 hf n hn).weierstrassDistinguished
    (iteratedLubinTateStepR_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n hn)

/-- `r_n` の Weierstrass 標準分解の単元部分。 -/
noncomputable def iteratedLubinTateStepU {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    PowerSeries A :=
  (iteratedLubinTateStepR hq hπmax hπne0 f hf0 hf1 hf n hn).weierstrassUnit
    (iteratedLubinTateStepR_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n hn)

theorem iteratedLubinTateStepR_eq_psi_mul_U {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    iteratedLubinTateStepR hq hπmax hπne0 f hf0 hf1 hf n hn =
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn : PowerSeries A) *
        iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn :=
  PowerSeries.eq_weierstrassDistinguished_mul_weierstrassUnit
    (iteratedLubinTateStepR_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n hn)

theorem isDistinguishedAt_iteratedLubinTatePsi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).IsDistinguishedAt (IsLocalRing.maximalIdeal A) :=
  PowerSeries.isDistinguishedAt_weierstrassDistinguished
    (iteratedLubinTateStepR_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n hn)

theorem isUnit_iteratedLubinTateStepU {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    IsUnit (iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn) :=
  PowerSeries.isUnit_weierstrassUnit
    (iteratedLubinTateStepR_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n hn)

/-- `ψ_n` の定数項と `U'_n` の定数項の積は `π`——`r_n=ψ_n・U'_n` の定数項
を比較し、`r_n` の定数項がちょうど `π` であること(`iteratedLubinTateStepR_
constantCoeff`)と組み合わせるだけ。 -/
theorem iteratedLubinTatePsi_coeff_zero_mul {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0 *
      PowerSeries.constantCoeff (iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn) = π := by
  have h1 := iteratedLubinTateStepR_constantCoeff hq hπmax hπne0 f hf0 hf1 hf n hn
  rw [iteratedLubinTateStepR_eq_psi_mul_U hq hπmax hπne0 f hf0 hf1 hf n hn, map_mul,
    Polynomial.constantCoeff_coe] at h1
  exact h1

/-- ★★★★★★★★★**`D_n=D_{n-1}・ψ_n`**——`(D_{n-1}・ψ_n)・(U_{n-1}・U'_n)` が
`[π^n]_f` のもう1つの Weierstrass 分解であることを確認し、一意性
から結論する(`iteratedLubinTateDistinguished_dvd_of_le` と同じ手法、
`a=1` の場合に特化)。 -/
theorem iteratedLubinTateDistinguished_eq_mul_psi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n =
      iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1) *
        iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn := by
  have hcombined : iteratedLubinTate f n =
      ((iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1) *
        iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn : Polynomial A) : PowerSeries A) *
        (iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf (n - 1) *
          iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn) := by
    rw [iteratedLubinTateStepR_eq hq hπmax hπne0 f hf0 hf1 hf n hn,
      iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf (n - 1),
      iteratedLubinTateStepR_eq_psi_mul_U hq hπmax hπne0 f hf0 hf1 hf n hn]
    push_cast
    ring
  have hfactcombined : (iteratedLubinTate f n).IsWeierstrassFactorization
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1) *
        iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)
      (iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf (n - 1) *
        iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn) :=
    { isDistinguishedAt :=
        (isDistinguishedAt_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1)).mul
          (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)
      isUnit :=
        (isUnit_iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf (n - 1)).mul
          (isUnit_iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn)
      eq_mul := hcombined }
  have huniq := hfactcombined.unique (iteratedLubinTate_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n)
  exact huniq.1.symm

/-- `ψ_n` の次数は `q^n-q^{n-1}`——`D_n=D_{n-1}・ψ_n` と両者の次数公式
(`natDegree_iteratedLubinTateDistinguished`)から従う。 -/
theorem natDegree_iteratedLubinTatePsi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).natDegree = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  have hdegsplit : (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).natDegree =
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1)).natDegree +
        (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).natDegree := by
    rw [iteratedLubinTateDistinguished_eq_mul_psi hq hπmax hπne0 f hf0 hf1 hf n hn]
    exact Polynomial.natDegree_mul
      (monic_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1)).ne_zero
      (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic.ne_zero
  rw [natDegree_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n,
    natDegree_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1)] at hdegsplit
  omega

/-! ### 部品3: `ψ_n` は既約(Eisenstein の判定法) -/

/-- ★★★★★★★★★**`ψ_n` の定数項は `maximalIdeal A ^ 2` に属さない**——
`n=1`(`φ_1`)の場合と全く同じ背理法: 属すと仮定すると `π` 自身が
単元になってしまい、`π∈𝔪` と矛盾する。 -/
theorem iteratedLubinTatePsi_coeff_zero_notMem_sq {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0 ∉
      IsLocalRing.maximalIdeal A ^ 2 := by
  intro hmem
  rw [hπmax, Ideal.span_singleton_pow] at hmem
  obtain ⟨w, hw⟩ := Ideal.mem_span_singleton'.mp hmem
  set U := PowerSeries.constantCoeff (iteratedLubinTateStepU hq hπmax hπne0 f hf0 hf1 hf n hn) with hU_def
  have hone : (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).coeff 0 * U = π :=
    iteratedLubinTatePsi_coeff_zero_mul hq hπmax hπne0 f hf0 hf1 hf n hn
  rw [← hw] at hone
  have hcancel : π * (π * (w * U)) = π * 1 := by rw [mul_one]; linear_combination hone
  have hπv : π * (w * U) = 1 := mul_left_cancel₀ hπne0 hcancel
  have hπunit : IsUnit π := IsUnit.of_mul_eq_one (w * U) hπv
  have hπmem : π ∈ IsLocalRing.maximalIdeal A := by
    rw [hπmax]; exact Ideal.subset_span rfl
  exact (IsLocalRing.maximalIdeal.isMaximal A).ne_top
    (Ideal.eq_top_of_isUnit_mem _ hπmem hπunit)

/-- ★★★★★★★★★★**`ψ_n` は Eisenstein**——distinguished(モニック・弱 Eisenstein)
であることと `coeff 0 ∉ 𝔪²` を束ねただけ。

★2026-09-05: `irreducible_iteratedLubinTatePsi` の証明の中に埋まっていたものを
**独立した定理として取り出した**。`Found/PGC/TotallyRamified.lean::
isTotallyRamifiedAdjoin_of_eisenstein`(Eisenstein の根が生成する拡大は完全分岐)
に渡すため。 -/
theorem isEisensteinAt_iteratedLubinTatePsi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).IsEisensteinAt
      (IsLocalRing.maximalIdeal A) := by
  have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
  constructor
  · rw [hmonic.leadingCoeff]
    exact fun h => (IsLocalRing.maximalIdeal.isMaximal A).ne_top
      (Ideal.eq_top_of_isUnit_mem _ h isUnit_one)
  · exact fun {k} hk =>
      (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).mem hk
  · exact iteratedLubinTatePsi_coeff_zero_notMem_sq hq hπmax hπne0 f hf0 hf1 hf n hn

/-- ★★★★★★★★★★★**`ψ_n` は既約(任意の `n≥1`)**——古典的な Lubin-Tate
理論の核心定理。原始 `π^n`-捩れ点(`Λ_n\Λ_{n-1}`)を1つ添加した拡大
`K(α)/K(Λ_{n-1})` は次数 `q^n-q^{n-1}` の**完全分岐拡大**である、
という基本定理の多項式版。Eisenstein の判定法
(`Polynomial.IsEisensteinAt.irreducible`)を、`ψ_n` が distinguished
であること(モニック・弱Eisenstein)と `coeff0∉𝔪²` に適用するだけ。
`n=1`(`irreducible_iteratedLubinTatePrimitive_one`)の一般化。 -/
theorem irreducible_iteratedLubinTatePsi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    Irreducible (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) := by
  have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
  have heis := isEisensteinAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn
  have hdegpos : 0 < (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).natDegree := by
    rw [natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]
    have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
    have hlt : (pp ^ ff) ^ (n - 1) < (pp ^ ff) ^ n := by
      apply Nat.pow_lt_pow_right h2; omega
    omega
  exact heis.irreducible (IsLocalRing.maximalIdeal.isMaximal A).isPrime hmonic.isPrimitive hdegpos

end ABC3.Found.PGC
