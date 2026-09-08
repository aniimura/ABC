import ABC3.Found.PGC.TotallyRamifiedAdjoinNorm

/-!
# [pGC] ★★★★★★★代入が通った —— `IsTotallyRamifiedAdjoin ⇒ hvalK`

前波の自分の依頼「部品は全部揃っているは見込みだ。次の波はそこを疑え」に従い、
★**部品表の 4 つを先に 1 本ずつ測った**。

## ★★測定の結果（§1）

| 測定 | 結果 |
|---|---|
| (A) `K(x)` の素元 | ★通る（`letI := nontriviallyNormedField_adjoin K x` を置いて） |
| (B) `K` の素元 | ★通る |
| (C) 包含写像の等長性 | ★`spectralNorm_extends a` で 1 行 |
| (D) 素元のノルム | ★`Valued.integer.norm_irreducible_pos` が体の側でそのまま使える |
| (E) `Irreducible` の移植 | ★**型を書き換えるだけ**（defeq） |
| (F) DVR instance の移植 | ★同じく defeq |
| (G) `ramificationIndex` の形 | ★`rfl` |
| (H) `hcompat` | ★`rfl` |

★★**`UnramifiedExtension.lean:690-695` の「`whnf` でタイムアウトする」という警告は、
`letI := nontriviallyNormedField_adjoin K x` を置けば起きない**（実測 8.6–10.2 秒）。
★警告自体がその対処を書いており、そのとおりにしただけである。

## ★唯一の詰まりとその直し方（実測したエラー）

```
error: failed to synthesize instance of type class
  Algebra ↥(IntegerNorm.integerSubring K.carrier) ↥(IntegerNorm.integerSubring ↥K.carrier⟮x⟯)
```

★TC 探索は `adjoinIntegersAlgebra`（`AdjoinIntegers.lean:319`）を見つけられない。
型が defeq でも**字面が違う**からである。
★直し方は `integerAlgebra`（本ファイル）のように
**欲しい型を明示して defeq で受け直すこと**。
同じ形が `IsDiscreteValuationRing.exists_irreducible` にも出た
（`(adjoinIntegers K x)` ではなく `↥(IntegerNorm.integerSubring …)` を渡す）。

## 出したもの

| 宣言 | 内容 |
|---|---|
| `integerAlgebra` | ★TC が見つけられない代数構造に名前を付ける |
| ★★★★★★★`valK_of_isTotallyRamifiedAdjoin` | **`IsTotallyRamifiedAdjoin ⇒ π と `hvalK`** |

## ★★これで残るのは `hnK` だけか（★未測定の部分を明示する）

`HarithPAdicSupply.exists_norm_sub_algebraMap_le_prod_axDecay_of_padic` の仮説のうち
`hvalK` は本ファイルが供給する。★しかし同定理は他に
`[NormedAlgebra ℚ_[p] M]`・`[IsScalarTower ℚ_[p] K.carrier M]`・
`[Algebra.IsAlgebraic ℚ_[p] M]`・`[IsUltrametricDist M]` を instance で要求する。
★★**`M = K.carrier⟮x⟯`（`K.closure` の中の中間体）にこれらが揃うかは本ファイル §3 で測る**。
★ノルムの出所が違うことに注意: 本ファイルの `M` のノルムは
`closureNormedField K`（= `spectralNorm K.carrier K.closure`）の制限であり、
`HarithPAdicSupply` §4 が使った `LocalFieldNorm.normedField L`（= `spectralNorm ℚ_[p] L.carrier`）
とは**別の道**である（`AdjoinPAdicLocalField.lean:44-49` が警告している diamond）。

## 逸脱の記録

1. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace TotRamSupply

open ABC3.Skeleton.PGC
open scoped NormedField Valued

section Probe

variable {p : ℕ} [Fact p.Prime]

/-- 測定 (A): `K(x)` の素元は取れるか（★`whnf` タイムアウトの警告は生きているか）。 -/
example (K : PAdicLocalField p) (x : K.closure) :
    ∃ ϖ : adjoinIntegers K x, Irreducible ϖ := by
  letI := nontriviallyNormedField_adjoin K x
  haveI : IsDiscreteValuationRing (adjoinIntegers K x) :=
    isDiscreteValuationRing_adjoinIntegers K x
  exact IsDiscreteValuationRing.exists_irreducible _

/-- 測定 (B): `K` の素元は取れるか。 -/
example (K : PAdicLocalField p) : ∃ ϖ : 𝒪[K.carrier], Irreducible ϖ := by
  haveI := isDiscreteValuationRing_carrierIntegers K
  exact IsDiscreteValuationRing.exists_irreducible _

/-- 測定 (C): 包含写像は等長か（`spectralNorm_extends`）。 -/
example (K : PAdicLocalField p) (x : K.closure) (a : K.carrier) :
    ‖algebraMap K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) a‖ = ‖a‖ :=
  spectralNorm_extends a

/-- 測定 (D): 素元のノルムは体の側で書けるか。 -/
example (K : PAdicLocalField p) (x : K.closure)
    (ϖ : adjoinIntegers K x) (hϖ : Irreducible ϖ) :
    letI := nontriviallyNormedField_adjoin K x
    0 < ‖(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ := by
  letI := nontriviallyNormedField_adjoin K x
  exact Valued.integer.norm_irreducible_pos hϖ

/-- 測定 (E): `Irreducible` は `integerSubring` の側に移るか。 -/
example (K : PAdicLocalField p) (x : K.closure)
    (ϖ : adjoinIntegers K x) (hϖ : Irreducible ϖ) :
    Irreducible (⟨(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)), ϖ.2⟩ :
      ↥(IntegerNorm.integerSubring
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) := hϖ

/-- 測定 (F): DVR の instance は `integerSubring` の側に移るか。 -/
example (K : PAdicLocalField p) (x : K.closure) :
    IsDiscreteValuationRing ↥(IntegerNorm.integerSubring
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
  isDiscreteValuationRing_adjoinIntegers K x

/-- ★**整数環の間の代数構造**を `integerSubring` の型で名前を付ける。

★★TC 探索は `adjoinIntegersAlgebra` を見つけられない（型が字面で違う）ので、
★**型を明示して defeq で受け直す**。
実測したエラー:
`failed to synthesize instance of type class
  Algebra ↥(IntegerNorm.integerSubring K.carrier) ↥(IntegerNorm.integerSubring ↥K.carrier⟮ x⟯)` -/
@[implicit_reducible] noncomputable def integerAlgebra (K : PAdicLocalField p) (x : K.closure) :
    Algebra ↥(IntegerNorm.integerSubring K.carrier)
      ↥(IntegerNorm.integerSubring
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
  adjoinIntegersAlgebra K x

/-- 測定 (G): `ramificationIndex` は `integerSubring` の形で書けるか。 -/
example (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    letI := IntegerNorm.isLocalRing_integerSubring (M := K.carrier)
    letI := IntegerNorm.isLocalRing_integerSubring
      (M := IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    letI := integerAlgebra K x
    ramificationIndex K x
      = Ideal.ramificationIdx
          (IsLocalRing.maximalIdeal ↥(IntegerNorm.integerSubring K.carrier))
          (IsLocalRing.maximalIdeal ↥(IntegerNorm.integerSubring
            (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) := rfl

/-- 測定 (H): `hcompat` を `integerSubring` の形で。 -/
example (K : PAdicLocalField p) (x : K.closure)
    (y : ↥(IntegerNorm.integerSubring K.carrier)) :
    letI := integerAlgebra K x
    ((algebraMap ↥(IntegerNorm.integerSubring K.carrier)
        ↥(IntegerNorm.integerSubring
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) y :
        ↥(IntegerNorm.integerSubring
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = algebraMap K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (y : K.carrier) := rfl


end Probe

/-! ## §2 ★★★★★★★代入 —— `IsTotallyRamifiedAdjoin ⇒ hvalK` -/

section Main

variable {p : ℕ} [Fact p.Prime]

/-- ★★★★★★★**完全分岐ならばノルムの側の `hvalK` が出る**。

`HarithPAdicSupply.exists_norm_sub_algebraMap_le_prod_axDecay_of_padic` に残っていた
2 本（`hnK` / `hvalK`）のうち ★**`hvalK` が `PAdicLocalField` の上で出た**。

★素元 `π` も同時に出す（出口が `π` を要求するため）。 -/
theorem valK_of_isTotallyRamifiedAdjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    ∃ π : IntermediateField.adjoin K.carrier ({x} : Set K.closure),
      0 < ‖π‖ ∧ ‖π‖ < 1 ∧
      ∀ a : K.carrier, a ≠ 0 → ∃ m : ℤ,
        ‖algebraMap K.carrier
            (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) a‖
          = ‖π‖ ^ ((Module.finrank K.carrier
              (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : ℤ) * m) := by
  classical
  letI := nontriviallyNormedField_adjoin K x
  haveI hDVRM : IsDiscreteValuationRing ↥(IntegerNorm.integerSubring
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
    isDiscreteValuationRing_adjoinIntegers K x
  haveI hDVRK : IsDiscreteValuationRing ↥(IntegerNorm.integerSubring K.carrier) :=
    isDiscreteValuationRing_carrierIntegers K
  letI := integerAlgebra K x
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible
    ↥(IntegerNorm.integerSubring
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
  obtain ⟨ϖK, hϖK⟩ :=
    IsDiscreteValuationRing.exists_irreducible ↥(IntegerNorm.integerSubring K.carrier)
  have hπ0 : 0 < ‖(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ :=
    Valued.integer.norm_irreducible_pos hϖ
  have hπ1 : ‖(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1 :=
    Valued.integer.norm_irreducible_lt_one hϖ
  have hK0 : 0 < ‖(ϖK : K.carrier)‖ := Valued.integer.norm_irreducible_pos hϖK
  have hK1 : ‖(ϖK : K.carrier)‖ < 1 := Valued.integer.norm_irreducible_lt_one hϖK
  have hπmem : (ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ∈ IntegerNorm.integerSubring
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ϖ.2
  have hKmem : (ϖK : K.carrier) ∈ IntegerNorm.integerSubring K.carrier := ϖK.2
  have hirrM : Irreducible
      (⟨(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)), hπmem⟩ :
        ↥(IntegerNorm.integerSubring
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) := hϖ
  have hirrK : Irreducible (⟨(ϖK : K.carrier), hKmem⟩ :
      ↥(IntegerNorm.integerSubring K.carrier)) := hϖK
  have hvalM : ∀ z : IntermediateField.adjoin K.carrier ({x} : Set K.closure), z ≠ 0 →
      ∃ m : ℤ, ‖z‖ = ‖(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ^ m :=
    fun z hz => TotRamAdjoin.exists_zpow_norm_of_dvr hπ0 hπmem hirrM hz
  have hvalKK : ∀ y : K.carrier, y ≠ 0 → ∃ m : ℤ, ‖y‖ = ‖(ϖK : K.carrier)‖ ^ m :=
    fun y hy => TotRamAdjoin.exists_zpow_norm_of_dvr hK0 hKmem hirrK hy
  have hiso : ∀ a : K.carrier,
      ‖algebraMap K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) a‖ = ‖a‖ :=
    fun a => spectralNorm_extends a
  have hpKne : algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (ϖK : K.carrier) ≠ 0 := by
    refine (map_ne_zero _).mpr ?_
    intro hc
    rw [hc, norm_zero] at hK0
    linarith
  obtain ⟨mc, hmc⟩ := hvalM _ hpKne
  have hmclt : ‖(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ^ mc < 1 := by
    rw [← hmc, hiso]
    exact hK1
  have hmcpos : 0 < mc := (zpow_lt_one_iff_right_of_lt_one₀ hπ0 hπ1).mp hmclt
  have htn : ((mc.toNat : ℕ) : ℤ) = mc := Int.toNat_of_nonneg (by omega)
  have hval : ‖algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (ϖK : K.carrier)‖
        = ‖(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ^ mc.toNat := by
    rw [hmc, ← zpow_natCast _ mc.toNat, htn]
  have he : ramificationIndex K x
      = Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    TotRamAdjoin.ramificationIndex_eq_finrank_of_isTotallyRamified K x ht
  refine ⟨(ϖ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)), hπ0, hπ1, ?_⟩
  exact TotRamAdjoin.valK_of_ramificationIdx_eq (fun _ => rfl) hπ0 hπ1 hπmem hvalM
    hK0 hK1 hKmem hvalKK hval hiso he

end Main

/-! ## §3 ★出口が要求する instance が `K.carrier⟮x⟯` に揃うかの測定 -/

section ExitInstances

variable {p : ℕ} [Fact p.Prime]

example (K : PAdicLocalField p) (x : K.closure) :
    IsUltrametricDist (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := inferInstance

example (K : PAdicLocalField p) (x : K.closure) :
    IsScalarTower ℚ_[p] K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := inferInstance

/-! ### ★★ここで止まった（実測したエラーを逐語で残す）

残り 2 つの instance は ★**TC で見つからない**:

```
error: failed to synthesize instance of type class
  NormedAlgebra ℚ_[p] ↥K.carrier⟮x⟯
error: failed to synthesize instance of type class
  Algebra.IsAlgebraic ℚ_[p] ↥K.carrier⟮x⟯
```

★`Algebra.IsAlgebraic` は `Algebra.IsAlgebraic.of_finite` を local instance にすれば出るはずである
（`LocalFieldNorm.lean:64` が同じ手を使っている）が、★`FiniteDimensional ℚ_[p] K.carrier⟮x⟯`
を先に作る必要があり、本波はそこに降りていない。

★★`NormedAlgebra ℚ_[p] ↥K.carrier⟮x⟯` は★**ノルムの出所の違いそのもの**である。
本ファイルの `K.carrier⟮x⟯` のノルムは `closureNormedField K`
（= `spectralNorm K.carrier K.closure` の制限）であり、
`spectralNorm.normedAlgebra ℚ_[p] _` が与えるのは `spectralNorm ℚ_[p] _` の方である。
★両者は数学的には一致するが、★**definitionally 一致するとは限らない**
—— `AdjoinPAdicLocalField.lean:44-49` がまさにこの diamond を警告している。
⇒ ★★★**次の 1 点は「2 つのスペクトルノルムの一致を橋渡す補題」である。** -/

end ExitInstances



/-! ## `.src` と 公理 -/

def valK_of_isTotallyRamifiedAdjoin.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms integerAlgebra
#print axioms valK_of_isTotallyRamifiedAdjoin

end TotRamSupply

end ABC3.Found.PGC
