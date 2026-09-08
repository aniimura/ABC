import ABC3.Found.PGC.IntegerDVR

/-!
# [pGC] 底の整数環 `𝒪_K` と `Algebra 𝒪_K 𝒪_M` —— 次のノードの第 1 歩

前波(第 1123)で `IsDiscreteValuationRing 𝒪_M` が載り、残りは
`HerbrandComposition.lean:455 herbrandPhiGroup_natCast` の
`[Algebra A B]` / `[SMulCommClass G A B]` / `hadj` / `[Fintype G]` になった。

## ★★測定 —— 「`K` にノルムが無いから `𝒪_K` が作れない」は**偽**

我々の設定の `K` は素の `Field`(ノルムを持たない)。★しかし `𝒪_K` は
**像のノルム**で定義できる:

    baseIntegerSubring K M := {a : K | ‖algebraMap K M a‖ ≤ 1} : Subring K

★超距離だけで `Subring` になる(`K` 上のノルムは要らない)。
★これで `[Algebra A B]`(`A = 𝒪_K`, `B = 𝒪_M`)が `codRestrict` 1 本で載る。

## 出したもの

| 宣言 | 内容 |
|---|---|
| `baseIntegerSubring` | `{a : K \| ‖algebraMap K M a‖ ≤ 1}` が `Subring K` |
| `baseRingHom` | `𝒪_K →+* 𝒪_M` |
| ★`baseAlgebra` | **`Algebra ↥(baseIntegerSubring K M) ↥(integerSubring M)`** |
| `coe_baseAlgebraMap` | 代数射の座標は `algebraMap K M`(`rfl`) |

## ★残り(正確に)

* `[SMulCommClass G 𝒪_K 𝒪_M]` —— `G` が底を固定することから出るはず。★**未**。
* `hadj : Algebra.adjoin 𝒪_K {π'} = ⊤`(`𝒪_M = 𝒪_K[π]`)——
  ★中身は `IntegerSubringNorm.norm_algebraMap_le_one_of_le_one`(係数が `𝒪_K` に入る)と
  `JumpMono.exists_coeff_norm_le`(素元冪が基底)に在るが、★**部分環の statement にしていない**。
* `[Fintype G]`、および Hasse–Arf 本体の `hresA` / `hAinj` / `habel` / `h1`。
  ★`h1 : lowerRamificationGroup B G 1 = ⊤` は第 1122 の
  `mem_lowerRamificationGroup_iff_norm` ＋ 跳びの下界(`1 ≤ u 0`)から出るはずだが★**測っていない**。

## 逸脱の記録

1. `instance` にせず `def`(`@[implicit_reducible]` ＋ `letI`)。前波までと同じ流儀。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

namespace IntegerNorm

section Base

variable (K M : Type*) [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★**底の整数環** —— `K` にノルムが無くても像のノルムで定義できる。 -/
def baseIntegerSubring : Subring K where
  carrier := {a : K | ‖algebraMap K M a‖ ≤ 1}
  zero_mem' := by simp
  one_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq, map_add] at *
    exact le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ha hb)
  neg_mem' := by intro a ha; simpa using ha
  mul_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq, map_mul, norm_mul] at *
    calc ‖algebraMap K M a‖ * ‖algebraMap K M b‖ ≤ 1 * 1 :=
          mul_le_mul ha hb (norm_nonneg _) zero_le_one
      _ = 1 := one_mul 1

@[simp] theorem mem_baseIntegerSubring {a : K} :
    a ∈ baseIntegerSubring K M ↔ ‖algebraMap K M a‖ ≤ 1 := Iff.rfl

/-- `𝒪_K →+* 𝒪_M`。 -/
def baseRingHom : ↥(baseIntegerSubring K M) →+* ↥(integerSubring M) :=
  RingHom.codRestrict ((algebraMap K M).comp (baseIntegerSubring K M).subtype)
    (integerSubring M) (fun a => a.2)

/-- ★★**`Algebra 𝒪_K 𝒪_M`**。 -/
@[implicit_reducible] def baseAlgebra :
    Algebra ↥(baseIntegerSubring K M) ↥(integerSubring M) :=
  (baseRingHom K M).toAlgebra

theorem coe_baseRingHom (a : ↥(baseIntegerSubring K M)) :
    ((baseRingHom K M a : ↥(integerSubring M)) : M) = algebraMap K M (a : K) := rfl

end Base

/-! ## `.src` と 公理 -/

def baseAlgebra.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms baseIntegerSubring
#print axioms baseRingHom
#print axioms baseAlgebra
#print axioms coe_baseRingHom

end IntegerNorm

end ABC3.Found.PGC
