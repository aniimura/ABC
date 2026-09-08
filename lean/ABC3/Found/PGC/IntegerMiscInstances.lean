import ABC3.Found.PGC.IntegerResidueBase

/-!
# [pGC] `HasseArfStrongInduction` の残り **9 項目を全部**載せる

直前の波で `hresA`(`exists_sub_mem_maximalIdeal`)と `hadj`(`adjoin_pi_eq_top`)が載った。
残っていたのは
`HasseArfStrongInduction.lean:447-459`
`exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top` の**9 項目**である。
本ファイルはその **9 項目すべて**を閉じる。

| 項目 | 本ファイルの宣言 | 見積もり→実測 |
| --- | --- | --- |
| 7 `hAinj` | `injective_baseRingHom` | 安 → 3 行 |
| 4 `[SMulCommClass G A B]` | `smulCommClass_base` | 安 → 5 行 |
| 6 `[CharP (ResidueField B) p]` | `charP_residueField` | 安 → 10 行 |
| 5 `[FaithfulSMul G B]` | `faithfulSMul_integer` | 中 → 25 行 |
| 8 `habel` | `mul_comm_of_forall_zpow`（抽象核） | 安 → 4 行 |
| 9 `h1 : G_1 = ⊤` | `lowerRamificationGroup_one_eq_top` | 中 → 7 行 |
| 1 `[IsDiscreteValuationRing A]` | `isDiscreteValuationRing_baseIntegerSubring` | 高いと見た→**7 行** |
| 2 `[IsNoetherian A B]` | `isNoetherian_baseIntegerSubring` | 高いと見た→**4 行**（前提は `Module.Finite`）|
| 3 `[Fintype G]` | `fintype_algEquiv`（測定のみ） | 無償だった |

## ★★前波の自分の断定を 2 つ覚す（訂正）

★他ファイルの docstring は書き換えない。ここに訂正を書く。

1. **「`K` にノルムが無いから `𝒪_K` を DVR にできない」は偽。**
   `BaseIntegerAlgebra.lean` では「像のノルムで部分環を作ればよい」まで進めたが、
   実は `NormedField.induced`(`Analysis/Normed/Field/Basic.lean:335`)で
   ★**`K` 自体をノルム体にできる**ので、DVR の証明は 1 行も複写しなくてよい。
   測定した事実(これが鍵):
   * `letI := NormedField.induced K M (algebraMap K M) _` の下で
     `‖a‖ = ‖algebraMap K M a‖` は ★**`rfl`**（`fast_instance%` でも射影は還元される）。
   * ★★`integerSubring K = baseIntegerSubring K M` が ★**`Subring.ext (fun _ => Iff.rfl)`**。
   ⇒ `IntegerDVR.isDiscreteValuationRing_integerSubring` を `M := K` でそのまま使える。

2. **「項目 1 ・2 は高い」も偽。** 項目 1 は 7 行、1 は 2 の前提だが
   `isNoetherian_of_isNoetherianRing_of_finite`(`RingTheory/Noetherian/Basic.lean:324`)が
   ★**instance** なので `letI` 2 つと `infer_instance` で閉じた。

## ★在庫の測定（コマンドを残す）

```
grep -n "NormedField.induced" .cache/mathlib-index.txt
  → Analysis/Normed/Field/Basic.lean:335  ★`R S` は明示引数（索引の行には出ない、#297）
grep -n "isUltrametricDist_of_forall_norm_natCast_le_one" .cache/mathlib-index.txt
  → Analysis/Normed/Field/Ultra.lean:103  ★`‖(n : K)‖ ≤ 1` だけで IsUltrametricDist が出る
grep -n "AlgEquiv.fintype" .cache/mathlib-index.txt
  → FieldTheory/Fixed.lean:318  ★**instance**。項目 3 は無償
grep -n "isNoetherian_of_isNoetherianRing_of_finite" .cache/mathlib-index.txt
  → RingTheory/Noetherian/Basic.lean:324  ★**instance**
grep -n "CharP.charP_iff_prime_eq_zero" .cache/mathlib-index.txt   → Algebra/CharP/Basic.lean:103
grep -n "IsLocalRing.residue_eq_zero_iff" .cache/mathlib-index.txt
  → RingTheory/LocalRing/ResidueField/Basic.lean:38
```

★「無いと思ったが在った」: `smul_pow` ではなく ★**`smul_pow'`**
(`Algebra/Group/Action/Basic.lean:188`, `r • x ^ n = (r • x) ^ n`)。
`smul_pow`(`Defs.lean:487`)は `(r • x) ^ n = r ^ n • x ^ n` で**向きも形も違う**。

## ★抽象核と具体層を別の宣言にしたところ

* `mul_comm_of_forall_zpow` / `eq_top_of_generator_mem` —— ★**純群論**。
  分岐・付値・Galois の語が 1 語も出ない。
  `h1` の中身は「部分群なので生成元だけ見ればよい」という群論である。
* `inducedNormedField` / `integerSubring_eq_baseIntegerSubring` ——
  ★体の拡大 `M/K` に依らず、単射な環準同型 1 つで済む。

## 逸脱の記録

1. `MulSemiringAction` を仮説 `hiso` から作るので `instance` にせず
   `def`(`@[implicit_reducible]`)＋`letI`。前波までと同じ流儀。
2. `smulCommClass_base` は `hfix : g • algebraMap K M a = algebraMap K M a` を
   ★**仮説で受ける**。`G` が抽象な `MulSemiringAction` なので
   `K` を固定することは型からは出ない（`M ≃ₐ[K] M` なら `AlgEquiv.commutes` で自動）。
3. `lowerRamificationGroup_one_eq_top` は `G := M ≃ₐ[K] M` に固定している。
   `RamNormBridge.mem_ramification_iff` が `M ≃ₐ[K] M` で書かれているため。
4. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace IntegerNorm

section Misc

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★`hAinj` —— `𝒪_K → 𝒪_M` は単射。 -/
theorem injective_baseRingHom : Function.Injective (baseRingHom K M) := by
  intro a b hab
  apply Subtype.ext
  have h := congrArg (fun y : ↥(integerSubring M) => (y : M)) hab
  simp only [coe_baseRingHom] at h
  exact (algebraMap K M).injective h

/-- ★`[SMulCommClass G 𝒪_K 𝒪_M]` —— `G` が底の像を固定するなら成り立つ。 -/
theorem smulCommClass_base {G : Type*} [Monoid G] [MulSemiringAction G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖)
    (hfix : ∀ (g : G) (a : K), g • algebraMap K M a = algebraMap K M a) :
    letI := integerMulSemiringAction hiso
    letI := baseAlgebra K M
    SMulCommClass G ↥(baseIntegerSubring K M) ↥(integerSubring M) := by
  letI := integerMulSemiringAction hiso
  letI := baseAlgebra K M
  refine ⟨fun g a b => ?_⟩
  apply Subtype.ext
  show g • ((algebraMap K M (a : K)) * (b : M))
    = (algebraMap K M (a : K)) * (g • (b : M))
  rw [smul_mul', hfix]

/-- ★`[CharP (ResidueField 𝒪_M) p]` —— `‖p‖ < 1` から出る。 -/
theorem charP_residueField {p : ℕ} (hp : p.Prime) (hpM : ‖((p : ℕ) : M)‖ < 1) :
    letI := isLocalRing_integerSubring (M := M)
    CharP (IsLocalRing.ResidueField ↥(integerSubring M)) p := by
  letI := isLocalRing_integerSubring (M := M)
  refine (CharP.charP_iff_prime_eq_zero hp).mpr ?_
  have hmem : ((p : ℕ) : ↥(integerSubring M)) ∈ IsLocalRing.maximalIdeal ↥(integerSubring M) := by
    rw [IsLocalRing.mem_maximalIdeal]
    intro hu
    have h1 := isUnit_iff_norm_eq_one.mp hu
    have hcoe : (((p : ℕ) : ↥(integerSubring M)) : M) = ((p : ℕ) : M) := by push_cast; ring
    rw [hcoe] at h1
    linarith
  have hres := (IsLocalRing.residue_eq_zero_iff
    (R := ↥(integerSubring M)) ((p : ℕ) : ↥(integerSubring M))).mpr hmem
  rw [map_natCast] at hres
  exact hres


/-- ★`[FaithfulSMul G 𝒪_M]` —— `M` 上で忠実なら `𝒪_M` 上でも忠実。

★`z = (z·π^n)/π^n` で `M` の元を `𝒪_M` の元 2 つの商に書く。 -/
theorem faithfulSMul_integer {G : Type*} [Monoid G] [MulSemiringAction G M] [FaithfulSMul G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖) {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hπmem : π ∈ integerSubring M) :
    letI := integerMulSemiringAction hiso
    FaithfulSMul G ↥(integerSubring M) := by
  letI := integerMulSemiringAction hiso
  refine ⟨fun {g h} hgh => ?_⟩
  have hpi : g • π = h • π := by
    have := hgh ⟨π, hπmem⟩
    exact congrArg (fun y : ↥(integerSubring M) => (y : M)) this
  refine eq_of_smul_eq_smul (α := M) (fun z => ?_)
  obtain ⟨n, hn⟩ := exists_pow_lt_of_lt_one (show (0:ℝ) < 1 / (‖z‖ + 1) by positivity) hπ1
  have hx1 : ‖z * π ^ n‖ ≤ 1 := by
    rw [norm_mul, norm_pow]
    have hz0 : (0:ℝ) ≤ ‖z‖ := norm_nonneg z
    have h2 : ‖z‖ * ‖π‖ ^ n ≤ ‖z‖ * (1 / (‖z‖ + 1)) := by
      refine mul_le_mul_of_nonneg_left (le_of_lt hn) hz0
    have h3 : ‖z‖ * (1 / (‖z‖ + 1)) ≤ 1 := by
      rw [mul_one_div, div_le_one (by linarith)]
      linarith
    linarith
  have hx := hgh ⟨z * π ^ n, hx1⟩
  have hxc : g • (z * π ^ n) = h • (z * π ^ n) :=
    congrArg (fun y : ↥(integerSubring M) => (y : M)) hx
  rw [smul_mul', smul_mul', smul_pow', smul_pow', hpi] at hxc
  have hne : (h • π) ^ n ≠ 0 := by
    refine pow_ne_zero n (norm_pos_iff.mp ?_)
    rw [hiso]
    exact hπ0
  exact mul_right_cancel₀ hne hxc

end Misc

/-! ## §2 ★★★`𝒪_K` を DVR にする —— ノルムを `K` に引き戻す -/

section BaseDVR

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★**`K` に `M` からノルムを引き戻す**。

★在庫の測定: `NormedField.induced`(`Analysis/Normed/Field/Basic.lean:335`)。
★★`R S` は**明示引数**(索引の行には出ない `variable (R S)`、#297)。 -/
@[implicit_reducible] def inducedNormedField (K M : Type*) [Field K] [NormedField M]
    [Algebra K M] : NormedField K :=
  NormedField.induced K M (algebraMap K M) (algebraMap K M).injective

omit [IsUltrametricDist M] in
/-- 引き戻したノルムは `rfl` で `‖algebraMap K M a‖`。 -/
theorem norm_induced (a : K) :
    letI := inducedNormedField K M
    ‖a‖ = ‖algebraMap K M a‖ := rfl

/-- 引き戻したノルムも非アルキメデス。

★`isUltrametricDist_of_forall_norm_natCast_le_one`(`Analysis/Normed/Field/Ultra.lean:103`)を
使うと `‖(n : K)‖ ≤ 1` だけで済む。 -/
@[implicit_reducible] def inducedIsUltrametricDist (K M : Type*) [Field K] [NormedField M]
    [IsUltrametricDist M] [Algebra K M] :
    letI := inducedNormedField K M
    IsUltrametricDist K := by
  letI := inducedNormedField K M
  refine IsUltrametricDist.isUltrametricDist_of_forall_norm_natCast_le_one (fun n => ?_)
  show ‖algebraMap K M (n : K)‖ ≤ 1
  rw [map_natCast]
  exact IsUltrametricDist.norm_natCast_le_one M n

/-- ★★**引き戻したノルムの整数環は `baseIntegerSubring` そのもの**。

`Iff.rfl` で通る（台集合が defeq、`fast_instance%` でも射影は還元される）。 -/
theorem integerSubring_eq_baseIntegerSubring :
    letI := inducedNormedField K M
    letI := inducedIsUltrametricDist K M
    integerSubring K = baseIntegerSubring K M := Subring.ext (fun _ => Iff.rfl)

/-- ★★★**`[IsDiscreteValuationRing 𝒪_K]`**（項目 1）。

★ノルムを `K` に引き戻して `IntegerDVR.isDiscreteValuationRing_integerSubring` を
`M := K` でそのまま使う。★**証明を複写しない**。

★前波の私の断定「`K` にノルムが無いから `𝒪_K` が作れない」は**偽**であり、
さらに「像のノルムで定義すればよい」と書いた後も、実は
**`NormedField.induced` で `K` 自体をノルム体にできる**ので写しは一切不要だった。 -/
theorem isDiscreteValuationRing_baseIntegerSubring {πK : K}
    (hπ0 : 0 < ‖algebraMap K M πK‖) (hπ1 : ‖algebraMap K M πK‖ < 1)
    (hval : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖algebraMap K M πK‖ ^ m)
    (hπmem : πK ∈ baseIntegerSubring K M) :
    IsDiscreteValuationRing ↥(baseIntegerSubring K M) := by
  letI := inducedNormedField K M
  letI := inducedIsUltrametricDist K M
  have hmem : πK ∈ integerSubring K := hπmem
  have h := isDiscreteValuationRing_integerSubring (M := K) (π := πK) hπ0 hπ1 hval hmem
  rw [integerSubring_eq_baseIntegerSubring] at h
  exact h

end BaseDVR


/-! ## §3 `habel` と `h1` -/

section Gen

/-- ★**抽象核(純群論)** —— 巡回群は可換。

`hgen` は「`g` が生成する」を `∃ n : ℤ, h = g ^ n` の形で受ける
(`IsCyclic` を経由しないのは、木の側で `zpowers g = ⊤` の形で出るため)。 -/
theorem mul_comm_of_forall_zpow {G : Type*} [Group G] {g : G}
    (hgen : ∀ h : G, ∃ n : ℤ, h = g ^ n) (x y : G) : x * y = y * x := by
  obtain ⟨n, rfl⟩ := hgen x
  obtain ⟨m, rfl⟩ := hgen y
  exact (Commute.refl g).zpow_zpow n m

/-- ★**抽象核(純群論)** —— 生成元が入れば部分群は `⊤`。 -/
theorem eq_top_of_generator_mem {G : Type*} [Group G] {g : G} {H : Subgroup G}
    (hgen : ∀ h : G, ∃ n : ℤ, h = g ^ n) (hg : g ∈ H) : H = ⊤ := by
  rw [Subgroup.eq_top_iff']
  intro h
  obtain ⟨n, rfl⟩ := hgen h
  exact zpow_mem hg n

end Gen

section H1

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★**`h1 : G_1 = ⊤`**(項目 9)。

`G_1` は部分群なので、生成元 `g` について測れば済む(`eq_top_of_generator_mem`)。
生成元の側は `RamNormBridge.mem_ramification_iff` の逆向きで
`1 ≤ t` から `∀ z, ‖z‖ ≤ 1 → ‖g z − z‖ ≤ ‖π‖²` が出る。
`t = u 0` なので `1 ≤ t` はちょうど `harith` (1) そのもの。 -/
theorem lowerRamificationGroup_one_eq_top [FiniteDimensional K M]
    {π : M} {n t : ℕ} {g : M ≃ₐ[K] M}
    (hiso : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ • z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hπmem : π ∈ integerSubring M)
    (hbr : ‖g π - π‖ = ‖π‖ ^ (t + 1)) (ht : 1 ≤ t)
    (hgen : ∀ h : M ≃ₐ[K] M, ∃ j : ℤ, h = g ^ j) :
    letI := isLocalRing_integerSubring (M := M)
    letI := integerMulSemiringAction hiso
    lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) 1 = ⊤ := by
  letI := isLocalRing_integerSubring (M := M)
  letI := integerMulSemiringAction hiso
  refine eq_top_of_generator_mem hgen ?_
  rw [mem_lowerRamificationGroup_iff_norm hiso hπ0 hπ1 hval hπmem]
  exact (RamNormBridge.mem_ramification_iff (i := 1) g (hiso g) hπ0 hπ1 hn hvalK hbr).mpr ht

end H1

/-! ## §4 ★★`Module.Finite 𝒪_K 𝒪_M` と `[IsNoetherian 𝒪_K 𝒪_M]` -/

section Noeth

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

open Finset JumpMono in
/-- ★★**`𝒪_M` は `𝒪_K` 上 `π^0, …, π^{n−1}` で生成される**。

★中身は既存の 2 本の合成であり、新しい不等式は 1 つも使わない:
* `JumpMono.exists_coeff_norm_le`(`JumpStrictMono.lean:116`) —— `z = Σ c_l π^l`、`‖c_l π^l‖ ≤ ‖z‖`
* `RamNormBridge.norm_algebraMap_le_one_of_le_one`(`RamificationGroupNormBridge.lean:202`)
  —— `‖c_l π^l‖ ≤ 1` ⇒ `‖c_l‖ ≤ 1`、すなわち `c_l ∈ 𝒪_K`。 -/
theorem module_finite_baseIntegerSubring [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hπmem : π ∈ integerSubring M) :
    letI := baseAlgebra K M
    Module.Finite ↥(baseIntegerSubring K M) ↥(integerSubring M) := by
  letI := baseAlgebra K M
  classical
  refine ⟨⟨image (fun l : Fin n => (⟨π, hπmem⟩ : ↥(integerSubring M)) ^ (l : ℕ)) univ, ?_⟩⟩
  rw [eq_top_iff]
  rintro x -
  obtain ⟨c, hc, hnorm⟩ := exists_coeff_norm_le hπ0 hπ1 hn hvalK (x : M)
  have hcmem : ∀ l : Fin n, c l ∈ baseIntegerSubring K M := by
    intro l
    rw [mem_baseIntegerSubring]
    refine RamNormBridge.norm_algebraMap_le_one_of_le_one hπ0 hπ1 hvalK l.isLt ?_
    have hl := hnorm l
    rw [Algebra.smul_def] at hl
    exact le_trans hl x.2
  set f : Fin n → ↥(integerSubring M) := fun l =>
    (⟨c l, hcmem l⟩ : ↥(baseIntegerSubring K M)) •
      ((⟨π, hπmem⟩ : ↥(integerSubring M)) ^ (l : ℕ)) with hf
  have hcoe : ((∑ l, f l : ↥(integerSubring M)) : M) = ∑ l, ((f l : ↥(integerSubring M)) : M) :=
    map_sum (integerSubring M).subtype f univ
  have hxeq : x = ∑ l, f l := by
    apply Subtype.ext
    rw [hcoe, ← hc]
    refine sum_congr rfl (fun l _ => ?_)
    rw [Algebra.smul_def]
    rfl
  rw [hxeq]
  refine Submodule.sum_mem _ (fun l _ => Submodule.smul_mem _ _ (Submodule.subset_span ?_))
  exact mem_coe.mpr (mem_image_of_mem _ (mem_univ l))

/-- ★★★**`[IsNoetherian 𝒪_K 𝒪_M]`**（項目 2）。

★在庫の測定: `isNoetherian_of_isNoetherianRing_of_finite`
(`RingTheory/Noetherian/Basic.lean:324`)は **instance** なので、
`IsNoetherianRing 𝒪_K`（DVR → PID → Noether）と `Module.Finite` を `letI` で置けば
`inferInstance` で出る。 -/
theorem isNoetherian_baseIntegerSubring [FiniteDimensional K M] {π : M} {πK : K} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hπmem : π ∈ integerSubring M)
    (hK0 : 0 < ‖algebraMap K M πK‖) (hK1 : ‖algebraMap K M πK‖ < 1)
    (hvalKK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖algebraMap K M πK‖ ^ m)
    (hKmem : πK ∈ baseIntegerSubring K M) :
    letI := baseAlgebra K M
    IsNoetherian ↥(baseIntegerSubring K M) ↥(integerSubring M) := by
  letI := baseAlgebra K M
  letI := isDiscreteValuationRing_baseIntegerSubring hK0 hK1 hvalKK hKmem
  letI := module_finite_baseIntegerSubring hπ0 hπ1 hn hvalK hπmem
  infer_instance

end Noeth

/-! ## §5 ★`[Fintype G]` は無償（項目 3）の測定 -/

section FintypeMeasure

/-- ★**測定** —— `AlgEquiv.fintype`(`FieldTheory/Fixed.lean:318`)は **instance** なので
`[FiniteDimensional K M]` から `Fintype (M ≃ₐ[K] M)` はそのまま出る。
★前波の表では「未測定」としていた項目 3 はこれで閉じた。 -/
theorem fintype_algEquiv (K M : Type*) [Field K] [Field M] [Algebra K M]
    [FiniteDimensional K M] : Nonempty (Fintype (M ≃ₐ[K] M)) := ⟨inferInstance⟩

end FintypeMeasure

/-! ## `.src` と 公理 -/

def isDiscreteValuationRing_baseIntegerSubring.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def lowerRamificationGroup_one_eq_top.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms injective_baseRingHom
#print axioms smulCommClass_base
#print axioms charP_residueField
#print axioms faithfulSMul_integer
#print axioms mul_comm_of_forall_zpow
#print axioms eq_top_of_generator_mem
#print axioms inducedNormedField
#print axioms integerSubring_eq_baseIntegerSubring
#print axioms isDiscreteValuationRing_baseIntegerSubring
#print axioms lowerRamificationGroup_one_eq_top
#print axioms module_finite_baseIntegerSubring
#print axioms isNoetherian_baseIntegerSubring
#print axioms fintype_algEquiv

end IntegerNorm

end ABC3.Found.PGC
