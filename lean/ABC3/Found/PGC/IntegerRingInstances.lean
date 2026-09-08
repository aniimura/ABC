import ABC3.Found.PGC.IntegerSubringNorm

/-!
# [pGC] `𝒪_M` に型クラスを載せる —— ★**4 つのうち 3 つが載り、橋が繋がった**

前波(第 1120)で残りは「型クラスを載せる作業 4 つ」と特定した。★本ファイルはそのうち
**2・3・4 を閉じ、1(`IsDiscreteValuationRing`)は要らないことを測った**。

## ★★★測定 —— 1 は 4 に要らない(前波の順序づけの訂正)

`LowerRamificationGroup.lean:265` の定義を読むと

```
def lowerRamificationGroup (B : Type*) [CommRing B] [IsLocalRing B] (G : Type*) [Group G]
    [MulSemiringAction G B] (n : ℕ) : Subgroup G := ((maximalIdeal B) ^ (n + 1)).inertia G
```

★**`IsDiscreteValuationRing` を要求していない。**要るのは `IsLocalRing` と
`MulSemiringAction` だけである。⇒ 前波が「1→2→3→4」と並べたのは順序として誤りで、
★**4 は 1 に依らない**。`IsDiscreteValuationRing` が要るのは `herbrandPhiGroup`
(`ramIndex` 経由)だけである。

## 出したもの

| 宣言 | 内容 |
|---|---|
| `isUnit_iff_norm_eq_one` | `𝒪_M` の単元 ⟺ `‖x‖ = 1` |
| ★`isLocalRing_integerSubring` | **`IsLocalRing 𝒪_M`**(`of_isUnit_or_isUnit_one_sub_self` ＋ 超距離) |
| ★★`integerMulSemiringAction` | **`MulSemiringAction G 𝒪_M`**(等長な環作用の制限、項目 2) |
| `norm_le_norm_pi_of_lt_one` | 離散性: `‖x‖ < 1` ⇒ `‖x‖ ≤ ‖π‖` |
| ★★`maximalIdeal_eq_span` | **`𝔪_{𝒪_M} = (π)`**(項目 3) |
| `coe_smul_integer` | 制限した作用の座標は `M` の作用そのもの(`rfl`) |
| ★★★`mem_lowerRamificationGroup_iff_norm` | **`σ ∈ G_i ⟺ ∀ z, ‖z‖ ≤ 1 → ‖σz − z‖ ≤ ‖π‖^{i+1}`**(項目 4) |

★★最後の右辺は `RamificationGroupNormBridge.mem_ramification_iff` の**左辺そのもの**である。
⇒ ★**環の言葉(`lowerRamificationGroup`)と体・ノルムの言葉が繋がった。**

## ★配管(実測)

* `hiso` を仮説で受けるので `MulSemiringAction` は `instance` にできず
  ★**`def` ＋ `letI`**(木の流儀、`lean-idioms.md` #165)。
* ★**class 型の `def` には `@[implicit_reducible]` が要る**:
  `Definition ... of class type must be marked with @[reducible] or @[implicit_reducible]`。
* `letI` を**statement に書く**形(`letI := …` を `:` の前に置く)は
  `RamCard`/`JumpMono` と同じ流儀で通る。

## ★残り(正確に、`file:line` つき)

★**`IsDiscreteValuationRing ↥(integerSubring M)` 1 本だけ**である。
中身は前波の `IntegerSubringNorm.exists_pow_mul_unit`
(`z ≠ 0`, `‖z‖ ≤ 1` ⇒ `z = π^n · w`, `‖w‖ = 1`)で尽きている。
入口の候補(前波で測った):

* `Valuation.valuationSubring_isDiscreteValuationRing`
  (`RingTheory/Valuation/Discrete/Basic.lean:453`、値群が巡回かつ非自明だけで DVR)
* `IsNonarchimedeanLocalField` の `instance`(`NumberTheory/LocalField/Basic.lean:108`)
* `NormedField.toValued`(`Topology/Algebra/Valued/NormedValued.lean:67`)

★これが載ると `HerbrandComposition.lean:455 herbrandPhiGroup_natCast` が使え、
本ファイルの項目 4 ＋ `RamCard.card_eq_pow_of_mem_iff` で `|G_i| = p^{k+1−m}` が入り、
`hrec` から `HasseArfCongruence.dvd_sub_of_phi_intCast` で `harith` の (3) が閉じる。

## 逸脱の記録

1. `MulSemiringAction` / `IsLocalRing` は**大域 `instance` にしていない**(`def`/`theorem` ＋ `letI`)。
   ★`hiso` を仮説で受けるため、および全体 import への影響を避けるため。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

namespace IntegerNorm

section Units

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★`𝒪_M` の単元 ⟺ ノルムが `1`。 -/
theorem isUnit_iff_norm_eq_one {x : ↥(integerSubring M)} : IsUnit x ↔ ‖(x : M)‖ = 1 := by
  constructor
  · rintro ⟨u, hu⟩
    have h1 : ((u : ↥(integerSubring M)) : M) * ((u⁻¹ : (integerSubring M)ˣ) : M) = 1 := by
      have := u.mul_inv
      exact congrArg (fun y : ↥(integerSubring M) => (y : M)) this
    have hx : ((x : M)) = ((u : ↥(integerSubring M)) : M) := by rw [hu]
    have hnorm : ‖((u : ↥(integerSubring M)) : M)‖ * ‖((u⁻¹ : (integerSubring M)ˣ) : M)‖ = 1 := by
      rw [← norm_mul, h1, norm_one]
    have hle1 : ‖((u : ↥(integerSubring M)) : M)‖ ≤ 1 := (u : ↥(integerSubring M)).2
    have hle2 : ‖((u⁻¹ : (integerSubring M)ˣ) : M)‖ ≤ 1 := (u⁻¹ : (integerSubring M)ˣ).1.2
    rw [hx]
    nlinarith [norm_nonneg ((u : ↥(integerSubring M)) : M),
      norm_nonneg ((u⁻¹ : (integerSubring M)ˣ) : M)]
  · intro hx
    have hx0 : (x : M) ≠ 0 := by
      intro h
      rw [h, norm_zero] at hx
      norm_num at hx
    have hmem : ((x : M))⁻¹ ∈ integerSubring M := by
      rw [mem_integerSubring, norm_inv, hx, inv_one]
    refine isUnit_iff_exists_inv.mpr ⟨⟨((x : M))⁻¹, hmem⟩, ?_⟩
    ext
    simpa using mul_inv_cancel₀ hx0

/-- ★★`𝒪_M` は局所環(超距離だけ)。 -/
theorem isLocalRing_integerSubring : IsLocalRing ↥(integerSubring M) := by
  haveI : Nontrivial ↥(integerSubring M) := inferInstance
  refine IsLocalRing.of_isUnit_or_isUnit_one_sub_self (fun a => ?_)
  have ha1 : ‖(a : M)‖ ≤ 1 := a.2
  rcases eq_or_lt_of_le ha1 with h1 | h1
  · exact Or.inl (isUnit_iff_norm_eq_one.mpr h1)
  · refine Or.inr (isUnit_iff_norm_eq_one.mpr ?_)
    have hco : ((1 - a : ↥(integerSubring M)) : M) = 1 - (a : M) := by push_cast; ring
    rw [hco]
    have hne : ‖(1 : M)‖ ≠ ‖-(a : M)‖ := by
      rw [norm_one, norm_neg]
      exact ne_of_gt h1
    have hmax := IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (x := (1 : M))
      (y := -(a : M)) hne
    rw [← sub_eq_add_neg] at hmax
    rw [hmax, norm_one, norm_neg]
    exact max_eq_left (le_of_lt h1)


end Units

/-! ## §2 作用の制限 -/

section Action

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★**等長な環作用は `𝒪_M` に制限できる** —— `MulSemiringAction G ↥(integerSubring M)`。

★木の流儀(`lean-idioms.md` #165)に合わせ **`instance` にせず `def`** にして `letI` で貼る
(`hiso` を仮説で受けるので `instance` にはできない)。 -/
@[implicit_reducible] def integerMulSemiringAction {G : Type*} [Monoid G] [MulSemiringAction G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖) :
    MulSemiringAction G ↥(integerSubring M) where
  smul g x := ⟨g • (x : M), by rw [mem_integerSubring, hiso]; exact x.2⟩
  one_smul x := Subtype.ext (one_smul G (x : M))
  mul_smul g h x := Subtype.ext (mul_smul g h (x : M))
  smul_zero g := Subtype.ext (smul_zero g)
  smul_add g x y := Subtype.ext (smul_add g (x : M) (y : M))
  smul_one g := Subtype.ext (smul_one g)
  smul_mul g x y := Subtype.ext (smul_mul' g (x : M) (y : M))

end Action

/-! ## §3 極大イデアル -/

section Maximal

variable {M : Type*} [NormedField M] [IsUltrametricDist M]
omit [IsUltrametricDist M] in

/-- ★離散性 —— `‖x‖ ≤ 1` かつ `‖x‖ ≠ 1` なら `‖x‖ ≤ ‖π‖`。 -/
theorem norm_le_norm_pi_of_lt_one {π x : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hx : ‖x‖ < 1) :
    ‖x‖ ≤ ‖π‖ := by
  rcases eq_or_ne x 0 with hx0 | hx0
  · rw [hx0, norm_zero]; exact le_of_lt hπ0
  · obtain ⟨m, hm⟩ := hval x hx0
    have hm1 : 1 ≤ m := by
      by_contra hcon
      rw [not_le] at hcon
      have : (1:ℝ) ≤ ‖π‖ ^ m := by
        calc (1:ℝ) = ‖π‖ ^ (0 : ℤ) := (zpow_zero _).symm
          _ ≤ ‖π‖ ^ m := zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (by omega)
      rw [← hm] at this
      linarith
    rw [hm]
    calc ‖π‖ ^ m ≤ ‖π‖ ^ (1 : ℤ) := zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) hm1
      _ = ‖π‖ := zpow_one _

/-- ★★★★**`𝔪_{𝒪_M} = (π)`** —— 極大イデアルは素元 `π` が生成する。 -/
theorem maximalIdeal_eq_span {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hπmem : π ∈ integerSubring M) :
    letI := isLocalRing_integerSubring (M := M)
    IsLocalRing.maximalIdeal ↥(integerSubring M)
      = Ideal.span {(⟨π, hπmem⟩ : ↥(integerSubring M))} := by
  letI := isLocalRing_integerSubring (M := M)
  ext x
  rw [IsLocalRing.mem_maximalIdeal, Ideal.mem_span_singleton]
  constructor
  · intro hx
    have hxn : ‖(x : M)‖ ≠ 1 := fun h => hx (isUnit_iff_norm_eq_one.mpr h)
    have hxlt : ‖(x : M)‖ < 1 := lt_of_le_of_ne x.2 hxn
    have hle : ‖(x : M)‖ ≤ ‖π‖ := norm_le_norm_pi_of_lt_one hπ0 hπ1 hval hxlt
    obtain ⟨w, hw, hxw⟩ := (dvd_iff_norm_le (π := π) (z := (x : M)) hπ0).mpr hle
    exact ⟨⟨w, hw⟩, Subtype.ext hxw⟩
  · rintro ⟨c, hc⟩ hunit
    have hxc : ‖(x : M)‖ ≤ ‖π‖ := by
      refine (dvd_iff_norm_le (π := π) (z := (x : M)) hπ0).mp ⟨(c : M), c.2, ?_⟩
      exact congrArg (fun y : ↥(integerSubring M) => (y : M)) hc
    have h1 : ‖(x : M)‖ = 1 := isUnit_iff_norm_eq_one.mp hunit
    rw [h1] at hxc
    linarith

end Maximal

/-! ## §4 下付き分岐群のノルム版 -/

section Ram

variable {M : Type*} [NormedField M] [IsUltrametricDist M]


/-- 制限した作用の座標は `M` の作用そのもの(`rfl`)。 -/
theorem coe_smul_integer {G : Type*} [Monoid G] [MulSemiringAction G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖) (g : G) (x : ↥(integerSubring M)) :
    letI := integerMulSemiringAction hiso
    ((g • x : ↥(integerSubring M)) : M) = g • (x : M) := rfl

open IsLocalRing in
/-- ★★★★★★**`lowerRamificationGroup` のノルム版** ——
`σ ∈ G_i ⟺ ∀ z, ‖z‖ ≤ 1 → ‖σ z − z‖ ≤ ‖π‖^{i+1}`。

★右辺は `RamNormBridge.mem_ramification_iff` の左辺そのものであり、
これで `harith` (3) の橋の**環側と体側が繋がる**。 -/
theorem mem_lowerRamificationGroup_iff_norm {G : Type*} [Group G] [MulSemiringAction G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖) {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hπmem : π ∈ integerSubring M)
    {i : ℕ} {σ : G} :
    letI := isLocalRing_integerSubring (M := M)
    letI := integerMulSemiringAction hiso
    σ ∈ lowerRamificationGroup ↥(integerSubring M) G i ↔
      ∀ z : M, ‖z‖ ≤ 1 → ‖σ • z - z‖ ≤ ‖π‖ ^ (i + 1) := by
  letI := isLocalRing_integerSubring (M := M)
  letI := integerMulSemiringAction hiso
  have hp0 : (0:ℝ) < ‖π ^ (i + 1)‖ := by rw [norm_pow]; exact pow_pos hπ0 _
  have hpn : ‖π ^ (i + 1)‖ = ‖π‖ ^ (i + 1) := norm_pow π (i + 1)
  rw [mem_lowerRamificationGroup_iff_forall,
    maximalIdeal_eq_span hπ0 hπ1 hval hπmem, Ideal.span_singleton_pow]
  constructor
  · intro hall z hz
    have h := hall ⟨z, hz⟩
    rw [Ideal.mem_span_singleton] at h
    obtain ⟨c, hc⟩ := h
    have hcoe : σ • z - z = π ^ (i + 1) * (c : M) := by
      have := congrArg (fun y : ↥(integerSubring M) => (y : M)) hc
      simpa [coe_smul_integer hiso] using this
    rw [← hpn]
    exact (dvd_iff_norm_le (π := π ^ (i + 1)) (z := σ • z - z) hp0).mp ⟨(c : M), c.2, hcoe⟩
  · intro hall x
    rw [Ideal.mem_span_singleton]
    have hx : ‖σ • (x : M) - (x : M)‖ ≤ ‖π ^ (i + 1)‖ := by
      rw [hpn]; exact hall (x : M) x.2
    obtain ⟨w, hw, hxw⟩ := (dvd_iff_norm_le (π := π ^ (i + 1))
      (z := σ • (x : M) - (x : M)) hp0).mpr hx
    refine ⟨⟨w, hw⟩, ?_⟩
    apply Subtype.ext
    simpa [coe_smul_integer hiso] using hxw

end Ram

/-! ## `.src` と 公理 -/

def maximalIdeal_eq_span.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def mem_lowerRamificationGroup_iff_norm.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms isUnit_iff_norm_eq_one
#print axioms isLocalRing_integerSubring
#print axioms integerMulSemiringAction
#print axioms norm_le_norm_pi_of_lt_one
#print axioms maximalIdeal_eq_span
#print axioms mem_lowerRamificationGroup_iff_norm

end IntegerNorm

end ABC3.Found.PGC
