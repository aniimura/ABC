import ABC3.Found.PGC.HarithAssembly
import ABC3.Found.PGC.LocalFieldNorm
import ABC3.Found.PGC.PureStepSetup

/-!
# [pGC] ★★★★★★★★`ℚ_p` を下に敷くと 5 本のうち 3 本が落ちる

前波の `HarithAssembly.exists_norm_sub_algebraMap_le_prod_axDecay_of_norm` は
ノルムの仮説を 5 本要求していた:
`hiso` / `hnK` / `hvalK` / `hnormp` / `heM`（と `he : 0 < e`）。

★★本ファイルで **`hiso` / `hnormp` / `heM` の 3 本と `he` が消えた**。
★★★さらに **`e`（絶対分岐指数）は入力でなくなり、値群から内部で作られる。**

## ★★安い順に測った結果（持ち場の問い）

| 本 | 結果 | 供給元 |
|---|---|---|
| `hnormp` | ★**3 行で落ちた** | `norm_algebraMap'`（mathlib）＋`Padic.norm_p` |
| `hiso` | ★**1 行で落ちた** | `PureStepSetup.norm_algHom_eq`（`PureStepSetup.lean:291`） |
| `heM` と `he` | ★**落ちた**（`e` は出力） | 本ファイル `exists_absRamIndex` |
| `hnK` | ★**残る** | —— |
| `hvalK` | ★**残る** | —— |

## ★★台帳の前回の測定を 1 つ覚す（訂正）

`decisions-pending.md` の

    VERDICT[GT-c]: 半分 — heM は定理になったが hnormp は落ちない

は ★**`[NormedAlgebra ℚ_[p] M]` がある設定では偽**である。
`LocalFieldNorm.lean:78 normedAlgebra` が `PAdicLocalField` にこの instance を与えており、
`norm_algebraMap' M x : ‖algebraMap ℚ_[p] M x‖ = ‖x‖`（mathlib、`[NormOneClass]` だけ）と
`Padic.norm_p` で 3 行である。

同じく前回の「② `harith` の下界には剰余体と等長性が要る」は、
★剰余体は本日 `one_le_jump_zero_pow` で不要になり、
★等長性も本ファイルで落ちた。⇒ ★★**② は完全に死んだ。**

## ★★前回の①と③が今も生きているか（持ち場の問い）

* ① **不分岐側（`f = p`）の 1 段** —— ★**生きている。本波は触っていない。**
  本ファイルの `hvalK` は「全分岐」そのものであり、不分岐層は対象外である。
* ③ **`IntermediateField` の層（#59 の危険区間）** ——
  ★★**出口の形を書くのには不要になった**。
  §4 は `M` を `IntermediateField.adjoin` ではなく
  ★**もう一つの `PAdicLocalField` の `carrier`** として取るので、
  `AdjoinPAdicLocalField.lean:44-49` が警告している instance diamond
  （`closureNormedField` 対 `spectralNorm.normedField ℚ_[p] _`）に入らない。
  ★ただし**一般の `K` と `x` からその `L` を作る**側には依然として `adjoin` が要る。
  ⇒ ★③ は**構成側にのみ生きている**。

## ★★残る 2 本の意味（次の 1 点）

`hnK : Module.finrank K M = p^{k+1}` と
`hvalK : ∀ a ≠ 0, ∃ m, ‖algebraMap K M a‖ = ‖π‖^{p^{k+1}·m}` は合わせて
★**「`M/K` は素元 `π` を持つ次数 `p^{k+1}` の全分岐拡大」**ということであり、
★これは「一般の `K` と `x` から塔を作る」という**内容そのもの**である。
★すなわちこれ以上は「配管」では落ちない。

## 逸脱の記録

1. §4 は `IsScalarTower ℚ_[p] K.carrier L.carrier` を★**仮説で受ける**。
   `PAdicLocalField` は `Algebra ℚ_[p] carrier` を各々独立に持つので、
   二つの体の間の塔の整合性は型からは出ない（D13/D30–D32）。
2. `Algebra.IsAlgebraic.of_finite` を local instance にしている（`LocalFieldNorm` と同じ流儀）。
3. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace HarithPAdic

open ABC3.Skeleton.PGC

/-! ## §1 ★★`hnormp` と `hiso` は `ℚ_p` を下に敷けば落ちる -/

section Supply

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**`hnormp` が落ちた** —— `‖(p : M)‖ = p⁻¹`。

★★**前回の測定の訂正**: 台帳 `decisions-pending.md` の
`VERDICT[GT-c]: 半分 — heM は定理になったが hnormp は落ちない` は
★**`[NormedAlgebra ℚ_[p] M]` がある設定では偽**である。
`norm_algebraMap'`（mathlib）と `Padic.norm_p` で 3 行で出る。
（`LocalFieldNorm.normedAlgebra` が `PAdicLocalField` にこの instance を与えている。） -/
theorem norm_natCast_p {M : Type*} [NormedField M] [NormedAlgebra ℚ_[p] M] :
    ‖((p : ℕ) : M)‖ = ((p : ℝ))⁻¹ := by
  have h : ((p : ℕ) : M) = algebraMap ℚ_[p] M ((p : ℕ) : ℚ_[p]) := by rw [map_natCast]
  rw [h, norm_algebraMap' M, Padic.norm_p]

/-- ★★★**`hiso` が落ちた** —— `K`-同型は等長。

★`PureStepSetup.norm_algHom_eq`（`PureStepSetup.lean:291`）を
`k := ℚ_[p]`、`F := K` で使うだけ。
★前波は抽象層（素の `[Field K]`）でこれが出ないと測ったが、
`ℚ_p` を下に敷けば完備性が手に入るので出る。 -/
theorem norm_algEquiv {K M : Type*} [Field K] [Algebra ℚ_[p] K] [NormedField M]
    [NormedAlgebra ℚ_[p] M] [Algebra K M] [IsScalarTower ℚ_[p] K M]
    [FiniteDimensional K M] [Algebra.IsAlgebraic ℚ_[p] M]
    (σ : M ≃ₐ[K] M) (z : M) : ‖σ z‖ = ‖z‖ :=
  PureStepSetup.norm_algHom_eq (k := ℚ_[p]) (F := K) (σ : M →ₐ[K] M) z

end Supply

/-! ## §2 ★★`e` と `heM` は**出力**になる -/

section RamIndex

variable {K M : Type*} [Field K] [NormedField M] [Algebra K M]

/-- ★★★**`he` と `heM` が落ちた** —— 絶対分岐指数 `e` は値群から**出てくる**。

`p ∈ K` なので `hvalK` で `‖(p:M)‖ = ‖π‖^{q·m}`。
`‖(p:M)‖ = p⁻¹ < 1` と `‖π‖ < 1` から `m ≥ 1`。★`e := m.toNat` とすればよい。 -/
theorem exists_absRamIndex {π : M} {q p : ℕ} (hq : 0 < q) (hp1 : 1 < p)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((q : ℤ) * m))
    (hnormp : ‖((p : ℕ) : M)‖ = ((p : ℝ))⁻¹) :
    ∃ e : ℕ, 0 < e ∧ ‖((p : ℕ) : M)‖ = ‖π‖ ^ (q * e) := by
  have hp0R : (0 : ℝ) < (p : ℝ) := by
    have : (0 : ℕ) < p := by omega
    exact_mod_cast this
  have hp1R : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp1
  have hpM0 : ((p : ℕ) : M) ≠ 0 := by
    intro hc
    rw [hc, norm_zero] at hnormp
    have : (0 : ℝ) < ((p : ℝ))⁻¹ := inv_pos.mpr hp0R
    linarith
  have hpK : (p : K) ≠ 0 := by
    intro hc
    apply hpM0
    have := congrArg (algebraMap K M) hc
    rwa [map_natCast, map_zero] at this
  obtain ⟨m, hm⟩ := hvalK (p : K) hpK
  rw [map_natCast] at hm
  have hlt1 : ‖π‖ ^ ((q : ℤ) * m) < 1 := by
    rw [← hm, hnormp]
    rw [inv_lt_one_iff₀]
    right; exact hp1R
  have hpos : 0 < (q : ℤ) * m := (zpow_lt_one_iff_right_of_lt_one₀ hπ0 hπ1).mp hlt1
  have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq
  have hm1 : 1 ≤ m := by nlinarith
  have htn : ((m.toNat : ℕ) : ℤ) = m := Int.toNat_of_nonneg (by omega)
  refine ⟨m.toNat, by omega, ?_⟩
  rw [hm, ← zpow_natCast ‖π‖ (q * m.toNat)]
  congr 1
  push_cast [htn]
  ring

end RamIndex

/-! ## §3 ★★★★★★★★出口から `hiso` / `hnormp` / `he` / `heM` が消えた -/

section Exit

variable {p : ℕ} [Fact p.Prime]

/-- ★★★★★★★★**`ℚ_p` の上の塔での出口**。

前波の `HarithAssembly.exists_norm_sub_algebraMap_le_prod_axDecay_of_norm` は
ノルムの仮説を 5 本（`hiso` / `hnK` / `hvalK` / `hnormp` / `heM`）要求していた。
★★**うち 3 本（`hiso` / `hnormp` / `heM`）と `he` が本定理で消えた。**

★★★★**`e`（絶対分岐指数）はもはや入力ではない**——値群から内部で作られる。
（結論に `e` が現れないので存在量化で消せる。）

残るノルムの仮説は `hnK`（次数）と `hvalK`（全分岐）の **2 本**、
それに `π` が素元であること（`hπ0` / `hπ1`）。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_padic
    {K M : Type*} [Field K] [Algebra ℚ_[p] K]
    [NormedField M] [IsUltrametricDist M] [NormedAlgebra ℚ_[p] M]
    [Algebra K M] [IsScalarTower ℚ_[p] K M] [FiniteDimensional K M]
    [Algebra.IsAlgebraic ℚ_[p] M]
    {k : ℕ} {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (τ : M →+ M) (hτ : ∀ z : M, τ z = g z)
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (x : M) :
    ∃ y : GainedTowerModel.GaloisTower.twr g p k,
      ‖x - algebraMap (GainedTowerModel.GaloisTower.twr g p k) M y‖
        ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp : p.Prime := Fact.out
  have hq0 : 0 < p ^ (k + 1) := pow_pos hp.pos _
  have hnormp : ‖((p : ℕ) : M)‖ = ((p : ℝ))⁻¹ := norm_natCast_p
  obtain ⟨e, he, heM⟩ :=
    exists_absRamIndex (q := p ^ (k + 1)) hq0 hp.one_lt hπ0 hπ1 hvalK hnormp
  exact HarithAssembly.exists_norm_sub_algebraMap_le_prod_axDecay_of_norm
    (p := p) (e := e) (k := k) (π := π) g hg (fun z => norm_algEquiv (p := p) g z) τ hτ s hsg he
    hnK hvalK hnormp heM x

end Exit

/-! ## §4 ★`PAdicLocalField` の 2 つ組での形 -/

section LocalField

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] Algebra.IsAlgebraic.of_finite

/-- ★★**`PAdicLocalField` の塔 `K ⊆ L` での出口**。

★内容は §3 への代入だけである。★本定理の価値は
**`LocalFieldNorm` の scoped instance がそのまま使えることを測った**ことにある:
`NormedField` / `IsUltrametricDist` / `NormedAlgebra ℚ_[p] _` は `LocalFieldNorm` が与え、
`Algebra.IsAlgebraic ℚ_[p] _` は `Algebra.IsAlgebraic.of_finite` を
local instance にすれば出る。

★★`IsScalarTower ℚ_[p] K.carrier L.carrier` は★**仮説で受ける**。
`PAdicLocalField` は `Algebra ℚ_[p] carrier` を各々独立に持つので、
★二つの体の間の塔の整合性は**型からは出ない**（D13/D30–D32 の同型不変性の問題）。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_localField
    (K L : PAdicLocalField p) [Algebra K.carrier L.carrier]
    [IsScalarTower ℚ_[p] K.carrier L.carrier] [FiniteDimensional K.carrier L.carrier]
    {k : ℕ} {π : L.carrier}
    (g : L.carrier ≃ₐ[K.carrier] L.carrier) (hg : orderOf g = p ^ (k + 1))
    (τ : L.carrier →+ L.carrier) (hτ : ∀ z : L.carrier, τ z = g z)
    (s : ℕ → (L.carrier →+* L.carrier))
    (hsg : ∀ j, ∀ z : L.carrier, s j z = (g ^ p ^ j) z)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hnK : Module.finrank K.carrier L.carrier = p ^ (k + 1))
    (hvalK : ∀ a : K.carrier, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K.carrier L.carrier a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (x : L.carrier) :
    ∃ y : GainedTowerModel.GaloisTower.twr g p k,
      ‖x - algebraMap (GainedTowerModel.GaloisTower.twr g p k) L.carrier y‖
        ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ :=
  exists_norm_sub_algebraMap_le_prod_axDecay_of_padic
    (K := K.carrier) (M := L.carrier) g hg τ hτ s hsg hπ0 hπ1 hnK hvalK x

end LocalField




/-! ## `.src` と 公理 -/

def exists_norm_sub_algebraMap_le_prod_axDecay_of_padic.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_natCast_p
#print axioms norm_algEquiv
#print axioms exists_absRamIndex
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_padic
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_localField

end HarithPAdic

end ABC3.Found.PGC
