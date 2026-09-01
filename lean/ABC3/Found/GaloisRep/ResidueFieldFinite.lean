/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.IntegerCongr

/-!
# 第 989 ブロック —— **★★★★★★★★★★★★★★★★★★★★完備化の剰余体は有限**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か —— mathlib の穴を埋める

第 982（分裂性の二者択一）は剰余体が**有限**であることを要求する。
☆だが `Finite (IsLocalRing.ResidueField (p.adicCompletionIntegers L))` は
mathlib のインスタンスになっていない（第 983・985 で実測）。
★`CompactSpace`・`DenseRange (𝓞 L → 整数環)`・剰余体の同一視のいずれも無い。

★本ブロックはそれを**内製**する。道は 3 段:

| 段 | 中身 | 出どころ |
|---|---|---|
| 1 | `L` は `Lv` で稠密 | `HeightOneSpectrum.denseRange_algebraMap`（mathlib） |
| 2 | `v_p(y) ≤ 1` なら `y ≡ a (mod p)` なる `a ∈ 𝓞 L` がある | 第 988 |
| 3 | よって `𝓞 L / p → ResidueField R` が全射 | `Ideal.finiteQuotientOfFreeOfNeBot`（mathlib） |

☆これは第 897（`IsAdicComplete` の内製）と同じ形の作業である——
真だが mathlib に無いので自分で建てる。

## ★配管の注意（実測）

`Valued.mem_nhds` は `Valued.v.restrict` と `ValueGroup₀` の言葉で述べられている。
★`γ := 1` を渡し、`Valuation.restrict_lt_iff_lt_embedding` と
`simp [Units.val_one, map_one]` で `Valued.v _ < 1` に直すこと。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField

/-- ★**`Valued.v` と `p.valuation` の橋**（`L` の元について）。 -/
theorem valued_algebraMap_eq (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (y : L) :
    Valued.v (algebraMap L (p.adicCompletion L) y) = HeightOneSpectrum.valuation L p y := by
  have h := HeightOneSpectrum.valuedAdicCompletion_eq_valuation' (R := 𝓞 L) (K := L) p y
  convert h using 2
  rfl

/-- ★**付値が `< 1` なら `𝔪` の元**（完備化の整数環で）。 -/
theorem mem_max_of_valued_lt_one (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (z : p.adicCompletionIntegers L)
    (hz : Valued.v (z : p.adicCompletion L) < 1) :
    z ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) := by
  rw [IsLocalRing.mem_maximalIdeal]
  intro hu
  rw [HeightOneSpectrum.adicCompletionIntegers.isUnit_iff_valued_eq_one] at hu
  rw [hu] at hz
  exact lt_irrefl _ hz

/-- ★★★★★★★★★★★★**完備化の整数環の元は `𝓞 L` の元と `𝔪` を法として合同**。

☆稠密性で `y : L` に寄せ、第 988 で `y` を `𝓞 L` の元に寄せる。 -/
theorem exists_integer_congr_completion (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (x : p.adicCompletionIntegers L) :
    ∃ a : 𝓞 L, Valued.v ((x : p.adicCompletion L)
      - algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)) < 1 := by
  have hd := HeightOneSpectrum.denseRange_algebraMap L p
  have hs : {z : p.adicCompletion L | Valued.v ((x : p.adicCompletion L) - z) < 1}
      ∈ nhds (x : p.adicCompletion L) := by
    rw [Valued.mem_nhds]
    refine ⟨1, ?_⟩
    intro z hz
    simp only [Set.mem_setOf_eq] at hz ⊢
    rw [Valuation.restrict_lt_iff_lt_embedding] at hz
    simp only [Units.val_one, map_one] at hz
    have heq : Valued.v ((x : p.adicCompletion L) - z)
        = Valued.v (z - (x : p.adicCompletion L)) := by
      rw [← neg_sub]; exact Valuation.map_neg _ _
    rw [heq]; exact hz
  obtain ⟨y, hy⟩ := hd.mem_nhds hs
  simp only [Set.mem_setOf_eq] at hy
  have hx1 : Valued.v (x : p.adicCompletion L) ≤ 1 :=
    (HeightOneSpectrum.mem_adicCompletionIntegers (𝓞 L) L p).1 x.2
  have hy1 : Valued.v (algebraMap L (p.adicCompletion L) y) ≤ 1 := by
    have heq2 : algebraMap L (p.adicCompletion L) y
        = (x : p.adicCompletion L) - ((x : p.adicCompletion L)
          - algebraMap L (p.adicCompletion L) y) := by ring
    rw [heq2]
    exact le_trans (Valuation.map_sub _ _ _) (max_le hx1 (le_of_lt hy))
  have hyv : HeightOneSpectrum.valuation L p y ≤ 1 := by
    rw [← valued_algebraMap_eq L p y]; exact hy1
  obtain ⟨a, ha⟩ := exists_integer_congr p y hyv
  refine ⟨a, ?_⟩
  have hya : Valued.v (algebraMap L (p.adicCompletion L)
      (y - algebraMap (𝓞 L) L a)) < 1 := by
    rw [valued_algebraMap_eq L p]; exact ha
  have hsplit : (x : p.adicCompletion L)
      - algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)
      = ((x : p.adicCompletion L) - algebraMap L (p.adicCompletion L) y)
        + algebraMap L (p.adicCompletion L) (y - algebraMap (𝓞 L) L a) := by
    rw [map_sub]; ring
  rw [hsplit]
  exact lt_of_le_of_lt (Valuation.map_add _ _ _) (max_lt hy hya)

/-- ★★★★★★★★★★★★★★★★★★★★**完備化の剰余体は有限である**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 989）**——mathlib に無かったので内製した。
☆`𝓞 L / p → ResidueField R` が全射（第 988 ＋ 稠密性）で、
`𝓞 L / p` は有限（`Ideal.finiteQuotientOfFreeOfNeBot`）だから。 -/
instance finite_residueField_adicCompletionIntegers (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) :
    Finite (IsLocalRing.ResidueField (p.adicCompletionIntegers L)) := by
  haveI : Finite (𝓞 L ⧸ p.asIdeal) := Ideal.finiteQuotientOfFreeOfNeBot p.asIdeal p.ne_bot
  have hker : ∀ a ∈ p.asIdeal,
      ((IsLocalRing.residue (p.adicCompletionIntegers L)).comp
        (algebraMap (𝓞 L) (p.adicCompletionIntegers L))) a = 0 := by
    intro a ha
    show IsLocalRing.residue _ (algebraMap (𝓞 L) (p.adicCompletionIntegers L) a) = 0
    rw [IsLocalRing.residue_eq_zero_iff]
    refine mem_max_of_valued_lt_one L p _ ?_
    have hcoe : ((algebraMap (𝓞 L) (p.adicCompletionIntegers L) a
        : p.adicCompletionIntegers L) : p.adicCompletion L)
        = algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a) := rfl
    rw [hcoe, valued_algebraMap_eq L p, HeightOneSpectrum.valuation_of_algebraMap,
      HeightOneSpectrum.intValuation_lt_one_iff_dvd]
    exact Ideal.dvd_span_singleton.2 ha
  refine Finite.of_surjective
    (Ideal.Quotient.lift p.asIdeal
      ((IsLocalRing.residue (p.adicCompletionIntegers L)).comp
        (algebraMap (𝓞 L) (p.adicCompletionIntegers L))) hker) ?_
  intro w
  obtain ⟨x, rfl⟩ := Ideal.Quotient.mk_surjective w
  obtain ⟨a, ha⟩ := exists_integer_congr_completion L p x
  refine ⟨Ideal.Quotient.mk p.asIdeal a, ?_⟩
  show IsLocalRing.residue _ (algebraMap (𝓞 L) (p.adicCompletionIntegers L) a)
    = IsLocalRing.residue _ x
  have hz : IsLocalRing.residue (p.adicCompletionIntegers L)
      (algebraMap (𝓞 L) (p.adicCompletionIntegers L) a - x) = 0 := by
    rw [IsLocalRing.residue_eq_zero_iff]
    refine mem_max_of_valued_lt_one L p _ ?_
    have hcoe : ((algebraMap (𝓞 L) (p.adicCompletionIntegers L) a - x
        : p.adicCompletionIntegers L) : p.adicCompletion L)
        = algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)
          - (x : p.adicCompletion L) := rfl
    rw [hcoe]
    have heq : algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)
        - (x : p.adicCompletion L)
        = -((x : p.adicCompletion L)
          - algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)) := by ring
    rw [heq, Valuation.map_neg]
    exact ha
  rw [map_sub, sub_eq_zero] at hz
  exact hz

/-! ## ★★★★★★★★★★★★第 990 —— 剰余の代表は `𝓞 L` から取れる

★第 982 の二者択一が出す捻り `d` は**完備化の整数環 `R` の元**である。
☆一方 `minDeltaExp_eq_mul_of_nonsplit`（第 929）が受ける捻りは **`L` の元**である。

★第 989 の全射性を取り出しておけば、`d` と同じ剰余をもつ `𝓞 L` の元に取り替えられる。
☆2 次式の係数は `d` に `φ d` を通してしか依存しないので、
剰余が同じなら `Splits` はそのまま移る。 -/

/-- ★★★★★★★★★★★★**完備化の整数環の剰余の代表は `𝓞 L` から取れる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 989 の全射性の中身を取り出したものである。 -/
theorem exists_integer_residue_eq (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (x : p.adicCompletionIntegers L) :
    ∃ a : 𝓞 L, IsLocalRing.residue (p.adicCompletionIntegers L)
        (algebraMap (𝓞 L) (p.adicCompletionIntegers L) a)
      = IsLocalRing.residue (p.adicCompletionIntegers L) x := by
  obtain ⟨a, ha⟩ := exists_integer_congr_completion L p x
  refine ⟨a, ?_⟩
  have hz : IsLocalRing.residue (p.adicCompletionIntegers L)
      (algebraMap (𝓞 L) (p.adicCompletionIntegers L) a - x) = 0 := by
    rw [IsLocalRing.residue_eq_zero_iff]
    refine mem_max_of_valued_lt_one L p _ ?_
    have hcoe : ((algebraMap (𝓞 L) (p.adicCompletionIntegers L) a - x
        : p.adicCompletionIntegers L) : p.adicCompletion L)
        = algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)
          - (x : p.adicCompletion L) := rfl
    rw [hcoe]
    have heq : algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)
        - (x : p.adicCompletion L)
        = -((x : p.adicCompletion L)
          - algebraMap L (p.adicCompletion L) (algebraMap (𝓞 L) L a)) := by ring
    rw [heq, Valuation.map_neg]
    exact ha
  rw [map_sub, sub_eq_zero] at hz
  exact hz

def exists_integer_residue_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(完備化の整数環の剰余の代表は 𝓞 L から取れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★出典の紐付け(`.src`) -/

def exists_integer_congr_completion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(完備化の整数環の元は 𝓞 L の元と 𝔪 を法として合同。★無条件)",
    sectionId := "genell-lemma-3-5" }

def finite_residueField_adicCompletionIntegers.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(完備化の剰余体は有限である。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
