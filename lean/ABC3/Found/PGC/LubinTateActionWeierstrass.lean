import ABC3.Found.PGC.LubinTateActionPiPow
import Mathlib.RingTheory.PowerSeries.WeierstrassPreparation

/-!
# `[π^n]_f` の Weierstrass 標準分解(`sorry` 無し)

節目(3)(捩れ点 `Λ_n:=ker[π^n]_f` の構成)へ向けた第三歩。`A` が
`π`-進完備(`[IsAdicComplete (maximalIdeal A) A]`——古典的な Lubin-Tate
理論の標準的な設定、`A=𝒪_K`)であるという、この節目で新たに必要になる
仮定のもとで、`[π^n]_f`(`=iteratedLubinTate f n`、`f`のn回自己合成)
に mathlib の Weierstrass 標準分解定理
(`RingTheory/PowerSeries/WeierstrassPreparation.lean`、前回発見済み)
を適用し、次数 `q^n`(`q:=pp^ff`)の distinguished 多項式を取り出す。

## 証明の筋

`iteratedLubinTate_map_residue`(前回確立、`[π^n]_f` は mod `π` で
`X^(q^n)`)から `[π^n]_f` の mod `π` の像が非零であることが従い、これが
`PowerSeries.exists_isWeierstrassFactorization` の唯一の前提。あとは
mathlib の標準分解定理(`eq_weierstrassDistinguished_mul_weierstrassUnit`・
`isDistinguishedAt_weierstrassDistinguished`・`isUnit_weierstrassUnit`)を
そのまま適用するだけ——Weierstrass 標準分解そのものを自前で構築する
必要は一切無かった。

次数が `q^n` に一致することは、mathlib の `Polynomial.IsDistinguishedAt.
coe_natDegree_eq_order_map`(distinguished多項式の次数はmod単元での
元の冪級数の位数に一致する)を、`iteratedLubinTate_map_residue` が
与える具体的な位数(`X^(q^n)` の位数 `=q^n`)と組み合わせて得る。

古典的な Lubin-Tate 理論では、この distinguished 多項式の根(が生成する
extension)が捩れ点 `Λ_n` そのものであり、その次数 `q^n` が
`|Λ_n|=q^n`(`Λ_n≅(𝒪_K/π^n)`、加法群として)の由来になる。
-/

namespace ABC3.Found.PGC

/-! ### 部品0: `[π^n]_f` の mod `π` の像は非零 -/

/-- `iteratedLubinTate_map_residue` と `X^k≠0` から直ちに従う。 -/
theorem iteratedLubinTate_map_residue_ne_zero {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    PowerSeries.map (IsLocalRing.residue A) (iteratedLubinTate f n) ≠ 0 := by
  rw [iteratedLubinTate_map_residue hq hπmax hπne0 f hf0 hf1 hf n]
  exact pow_ne_zero _ PowerSeries.X_ne_zero

/-! ### 部品1: distinguished 多項式・単元への分解 -/

/-- `[π^n]_f` の Weierstrass 標準分解の distinguished 多項式の部分。 -/
noncomputable def iteratedLubinTateDistinguished {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    Polynomial A :=
  (iteratedLubinTate f n).weierstrassDistinguished
    (iteratedLubinTate_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n)

/-- `[π^n]_f` の Weierstrass 標準分解の単元の部分。 -/
noncomputable def iteratedLubinTateUnit {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    PowerSeries A :=
  (iteratedLubinTate f n).weierstrassUnit
    (iteratedLubinTate_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n)

/-! ### 部品2: 分解の性質 -/

/-- ★★★★★★★**`[π^n]_f` の Weierstrass 標準分解**:
`[π^n]_f = D_n · U_n`(`D_n:=`distinguished 多項式、`U_n:=`単元)。 -/
theorem iteratedLubinTate_eq_distinguished_mul_unit {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    iteratedLubinTate f n =
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n : PowerSeries A) *
        iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n :=
  PowerSeries.eq_weierstrassDistinguished_mul_weierstrassUnit
    (iteratedLubinTate_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n)

/-- `D_n` は distinguished 多項式(`IsDistinguishedAt`、`maximalIdeal A`
に関して)——係数が最高次以外すべて `maximalIdeal A` に属し、
定数項は `maximalIdeal A ^ 2` に属さない(Eisenstein 型の条件)。 -/
theorem isDistinguishedAt_iteratedLubinTateDistinguished {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).IsDistinguishedAt
      (IsLocalRing.maximalIdeal A) :=
  PowerSeries.isDistinguishedAt_weierstrassDistinguished
    (iteratedLubinTate_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n)

/-- `U_n` は単元(冪級数として)。 -/
theorem isUnit_iteratedLubinTateUnit {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    IsUnit (iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n) :=
  PowerSeries.isUnit_weierstrassUnit
    (iteratedLubinTate_map_residue_ne_zero hq hπmax hπne0 f hf0 hf1 hf n)

/-! ### 部品3: `D_n` の次数は `q^n` -/

/-- ★★★★★★★★**`D_n` の次数は `q^n`**(`q:=pp^ff`)——古典的な Lubin-Tate
理論で `|Λ_n|=q^n` の由来になる事実。distinguished 多項式の次数は
mod `maximalIdeal A` での元の冪級数の位数に一致する
(`Polynomial.IsDistinguishedAt.coe_natDegree_eq_order_map`)ことと、
`iteratedLubinTate_map_residue`(`[π^n]_f` は mod `π` で `X^(q^n)`、
位数 `q^n`)を組み合わせて得る。 -/
theorem natDegree_iteratedLubinTateDistinguished {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).natDegree = (pp ^ ff) ^ n := by
  have hunotmem : PowerSeries.constantCoeff (iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n) ∉
      IsLocalRing.maximalIdeal A := by
    rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff, not_not]
    exact PowerSeries.isUnit_constantCoeff _ (isUnit_iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n)
  have hkey := Polynomial.IsDistinguishedAt.coe_natDegree_eq_order_map
    (iteratedLubinTate f n) (iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n)
    (isDistinguishedAt_iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)
    hunotmem
    (iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf n)
  have hkey' : ((iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).natDegree : ℕ∞) =
      (((pp ^ ff) ^ n : ℕ) : ℕ∞) := by
    rw [hkey]
    show ((PowerSeries.map (IsLocalRing.residue A)) (iteratedLubinTate f n)).order = ((pp ^ ff) ^ n : ℕ)
    rw [iteratedLubinTate_map_residue hq hπmax hπne0 f hf0 hf1 hf n, PowerSeries.order_X_pow]
  exact_mod_cast hkey'

/-! ### 部品4: `D_n` の定数項は0(`X ∣ D_n`) -/

/-- `[π^n]_f` 自身の定数項が0であること(`LubinTateAction_pi_pow` 経由で
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

/-- ★★★★★★★**`D_n` の定数項は正確に0**——`[π^n]_f` 自身の定数項が0で、
`U_n` の定数項が単元(`≠0`、`A` は整域)であることから、
`[π^n]_f = D_n·U_n` の定数項を比較して直ちに従う。これは
`D_n(X) = X·φ_n(X)`(`φ_n` が原始 `π^n`-捩れ点を統べる多項式)という
古典的な構造の出発点。 -/
theorem iteratedLubinTateDistinguished_coeff_zero {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).coeff 0 = 0 := by
  have hconst0 : PowerSeries.constantCoeff (iteratedLubinTate f n) = 0 :=
    constantCoeff_iteratedLubinTate hq hπmax hπne0 f hf0 hf1 hf n
  have heq := iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf n
  have hprod : PowerSeries.constantCoeff (iteratedLubinTate f n) =
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n).coeff 0 *
        PowerSeries.constantCoeff (iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n) := by
    rw [heq, map_mul, Polynomial.constantCoeff_coe]
  rw [hconst0] at hprod
  have hune0 : PowerSeries.constantCoeff (iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n) ≠ 0 :=
    (PowerSeries.isUnit_constantCoeff _ (isUnit_iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n)).ne_zero
  rcases mul_eq_zero.mp hprod.symm with h | h
  · exact h
  · exact absurd h hune0

/-- `X ∣ D_n`——`D_n` の定数項が0であることの言い換え。 -/
theorem X_dvd_iteratedLubinTateDistinguished {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) :
    Polynomial.X ∣ iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n :=
  Polynomial.X_dvd_iff.mpr (iteratedLubinTateDistinguished_coeff_zero hq hπmax hπne0 f hf0 hf1 hf n)

end ABC3.Found.PGC
