import ABC3.Found.PGC.BaseIntegerAlgebra

/-!
# [pGC] `hresA`(剰余体が底に一致)と `hadj`(`𝒪_M = 𝒪_K[π]`)を**部分環の statement にする**

前波(第 1124)で残った `herbrandPhiGroup_natCast` / Hasse–Arf 本体の仮説のうち、
★**`hresA` と `hadj` の 2 本を閉じた**。

## ★#323 の測定(着手前) —— 部分環の上では**焼き切れない**

`lean-idioms.md #323` は「`NormedField M` のまま `Algebra.adjoin` の所属を書くと
`isDefEq` timeout(200000 heartbeat)で 3 回焼き切れ、`[Field M]` だけの層に切り出すと通った」
と記録している。★今回は `Algebra.adjoin` を**部分環 `↥(integerSubring M)` の上**で書くので
条件が違う。★実測した(probe):

```lean
example {π : M} (hπmem : π ∈ integerSubring M) (l : ℕ) :
    letI := baseAlgebra K M
    (⟨π, hπmem⟩ : ↥(integerSubring M)) ^ l ∈ Algebra.adjoin ↥(baseIntegerSubring K M) {⟨π, hπmem⟩} := by
  letI := baseAlgebra K M
  exact pow_mem (Algebra.self_mem_adjoin_singleton _ _) l
```
★**8.2 秒(baseline)で通った。**⇒ #323 は **`NormedField M` を直接扱う層に固有**であり、
部分環の層には効かない(#332 と同じ形の測定)。

## 出したもの

| 宣言 | 内容 |
|---|---|
| ★★`exists_sub_mem_maximalIdeal` | **`hresA`**: `∀ b : 𝒪_M, ∃ a : 𝒪_K, b − algebraMap a ∈ 𝔪_{𝒪_M}` |
| ★★★`adjoin_pi_eq_top` | **`hadj`**: `Algebra.adjoin 𝒪_K {π'} = ⊤`(＝ `𝒪_M = 𝒪_K[π]`) |

中身はどちらも既存:
* `hresA` = `WildBreak.exists_sub_algebraMap_norm_le`(全分岐なら剰余体が下りる)
  ＋「`a` も整」(超距離)＋「`𝔪 = (π)`」(第 1122)。
* `hadj` = `JumpMono.exists_coeff_norm_le`(素元冪が基底、各項 ≤ 全体)
  ＋ ★`RamNormBridge.norm_algebraMap_le_one_of_le_one`(係数も整)。

★**在庫の測定の訂正(自分の前波の報告)**: 前波は係数の補題を
「`IntegerSubringNorm.norm_algebraMap_le_one_of_le_one`」と書いたが、★**誤り**である。
実際は `RamificationGroupNormBridge.lean:202`(名前空間 `RamNormBridge`)に在る。
`Unknown identifier` で気づいた(#330 の手順どおり `grep -rn "theorem <名前>" lean/ABC3/Found/PGC/*.lean`)。

## ★★残り(`HasseArfStrongInduction.lean:447-459` の全仮説に対する過不足、実測)

| 要求 | 現状 |
|---|---|
| `[CommRing A]` `[IsDomain A]` | ★済(`baseIntegerSubring` は体の部分環) |
| `[IsDiscreteValuationRing A]`(★`𝒪_K` 側) | ★**未**。第 1123 の証明が `𝒪_K` にも効くかは未測定(値群の生成元は `TotallyRamifiedValueGroup.exists_valSub_gen` が出すはず) |
| `[CommRing B]` `[Algebra A B]` `[IsDomain B]` `[IsDiscreteValuationRing B]` | ★済(第 1123・1124) |
| `[IsNoetherian A B]` | ★**未** |
| `[MulSemiringAction G B]` | ★済(第 1122) |
| `[Fintype G]` | ★**未**(容易) |
| `[SMulCommClass G A B]` | ★**未** |
| `[FaithfulSMul G B]` | ★**未** |
| `[CharP (ResidueField B) p]` | ★**未**(`‖p‖ < 1` ⇒ `p ∈ 𝔪` から出るはず) |
| `hπ' : Irreducible π'` | ★済(第 1123 `irreducible_pi`) |
| `hresA` | ★★**本ファイルで済** |
| `hadj` | ★★**本ファイルで済** |
| `hAinj : Function.Injective (algebraMap A B)` | ★**未**(容易: 体の射) |
| `habel` | ★**未**(巡回なので容易) |
| `h1 : lowerRamificationGroup B G 1 = ⊤` | ★**未**。第 1122 `mem_lowerRamificationGroup_iff_norm` ＋ 跳びの下界 `1 ≤ u 0` から出るはずだが**測っていない** |

## 逸脱の記録

1. `letI` を statement に書く形(`baseAlgebra` / `isLocalRing_integerSubring`)。前波までと同じ流儀。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

namespace IntegerNorm

section Res

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★**`hresA`** —— `𝒪_M` の任意の元は `𝒪_K` の元と `𝔪_{𝒪_M}` を法として合同。

★中身は `WildBreak.exists_sub_algebraMap_norm_le`(全分岐なら剰余体が下りる)。
★ここでは「係数が `𝒪_K` に入る」ことと「`𝔪 = (π)`」を足して**部分環の形**にした。 -/
theorem exists_sub_mem_maximalIdeal [FiniteDimensional K M] {π : M} {n : ℕ} (hn1 : 1 < n)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hval : ∀ z : M, z ≠ 0 → ∃ m : ℤ, ‖z‖ = ‖π‖ ^ m) (hπmem : π ∈ integerSubring M)
    (b : ↥(integerSubring M)) :
    letI := isLocalRing_integerSubring (M := M)
    letI := baseAlgebra K M
    ∃ a : ↥(baseIntegerSubring K M),
      b - algebraMap ↥(baseIntegerSubring K M) ↥(integerSubring M) a
        ∈ IsLocalRing.maximalIdeal ↥(integerSubring M) := by
  letI := isLocalRing_integerSubring (M := M)
  letI := baseAlgebra K M
  obtain ⟨a, ha⟩ := WildBreak.exists_sub_algebraMap_norm_le hn1 hπ0 hπ1 hn hvalK b.2
  have hamem : a ∈ baseIntegerSubring K M := by
    rw [mem_baseIntegerSubring]
    have hsplit : algebraMap K M a = (b : M) + (-((b : M) - algebraMap K M a)) := by ring
    rw [hsplit]
    refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le b.2 ?_)
    rw [norm_neg]
    exact le_trans ha (le_of_lt hπ1)
  refine ⟨⟨a, hamem⟩, ?_⟩
  rw [maximalIdeal_eq_span hπ0 hπ1 hval hπmem, Ideal.mem_span_singleton]
  have hcoe : ((b - algebraMap ↥(baseIntegerSubring K M) ↥(integerSubring M)
      ⟨a, hamem⟩ : ↥(integerSubring M)) : M) = (b : M) - algebraMap K M a := rfl
  obtain ⟨w, hw, hbw⟩ := (dvd_iff_norm_le (π := π) (z := (b : M) - algebraMap K M a) hπ0).mpr ha
  exact ⟨⟨w, hw⟩, Subtype.ext (by rw [hcoe]; simpa using hbw)⟩


/-- ★★★★★**`hadj`** —— `𝒪_M = 𝒪_K[π]`(部分環の形)。

★中身は `JumpMono.exists_coeff_norm_le`(素元冪が基底、各項 ≤ 全体)と
`IntegerSubringNorm.norm_algebraMap_le_one_of_le_one`(係数も整)の 2 本。

★#323 の測定: 部分環 `↥(integerSubring M)` の上で `Algebra.adjoin` の所属を書いても
**焼き切れない**(実測 8.2 秒)。#323 は `NormedField M` を**直接**扱う層に固有である。 -/
theorem adjoin_pi_eq_top [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hπmem : π ∈ integerSubring M) :
    letI := baseAlgebra K M
    Algebra.adjoin ↥(baseIntegerSubring K M)
      ({(⟨π, hπmem⟩ : ↥(integerSubring M))} : Set ↥(integerSubring M)) = ⊤ := by
  letI := baseAlgebra K M
  refine eq_top_iff.mpr (fun x _ => ?_)
  obtain ⟨c, hc, hnorm⟩ := JumpMono.exists_coeff_norm_le hπ0 hπ1 hn hvalK (x : M)
  have hcmem : ∀ l : Fin n, c l ∈ baseIntegerSubring K M := by
    intro l
    rw [mem_baseIntegerSubring]
    refine RamNormBridge.norm_algebraMap_le_one_of_le_one hπ0 hπ1 hvalK l.isLt ?_
    have h0 := hnorm l
    rw [Algebra.smul_def] at h0
    exact le_trans h0 x.2
  have hx : x = ∑ l : Fin n, (algebraMap ↥(baseIntegerSubring K M) ↥(integerSubring M)
      ⟨c l, hcmem l⟩) * (⟨π, hπmem⟩ : ↥(integerSubring M)) ^ (l : ℕ) := by
    apply Subtype.ext
    rw [← hc]
    push_cast
    refine Finset.sum_congr rfl (fun l _ => ?_)
    rw [Algebra.smul_def]
    rfl
  rw [hx]
  refine sum_mem (fun l _ => mul_mem (Subalgebra.algebraMap_mem _ _)
    (pow_mem (Algebra.self_mem_adjoin_singleton _ _) _))

end Res

/-! ## `.src` と 公理 -/

def exists_sub_mem_maximalIdeal.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def adjoin_pi_eq_top.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms exists_sub_mem_maximalIdeal
#print axioms adjoin_pi_eq_top

end IntegerNorm

end ABC3.Found.PGC
