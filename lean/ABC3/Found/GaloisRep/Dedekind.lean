import ABC3.Found.GaloisRep.LocalPrincipal
import Mathlib.RingTheory.DedekindDomain.Dvr
import Mathlib.RingTheory.DiscreteValuationRing.TFAE

/-!
# Galois (G5) 第 137 ブロック —— **★★★★★★★座標環は Dedekind 環**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★これが層 3 の関門であった

§9-442 で「`IsIntegrallyClosed W.CoordinateRing` **だけ**が欠けている(10-25 ブロック)」
と測った。★本ブロックでそれが出る。

    IsDedekindDomainDvr F[W]                     (mathlib の class)
      = IsNoetherian F[W] F[W]                   ✅ mathlib
        + 「0 でない素イデアルでの局所化がすべて DVR」

★★DVR であることは mathlib の `IsDiscreteValuationRing.TFAE` の
**項 0 ⟺ 項 4**(`DVR ⟺ 極大イデアルが単項`)に流す。単項性は第 136 で出ている。

★★★`IsDedekindDomainDvr.isDedekindDomain` は **instance** なので、
そこから `IsDedekindDomain` が出て、`IsIntegrallyClosed` は**副産物**として付いてくる。

## ★★★★★★ここまでの一本道(第 127 → 第 137)

| ブロック | 内容 |
|---|---|
| 127 | `Ring.DimensionLEOne F[W]`(整拡大から) |
| 129 | 平方完成 `z² = Ψ₂Sq(x)` |
| 130 | `Ψ₂Sq` の判別式は `16Δ`(mathlib にあった) |
| 131 | `Ψ₂Sq` は squarefree |
| 132 | 分岐点の局所構造 `z² = (x−c)·e(x)` |
| 133 | 不分岐点の局所構造 `(z−s)(z+s) = (x−c)·g(x)` |
| 134 | 素イデアルの `x` 座標(零点定理なし) |
| 135 | `P = (x − c, y − y₀)` |
| 136 | 局所化の極大イデアルが単項 |
| **137** | **`IsDedekindDomain F[W]`** |

★★★★★★★**11 ブロックで抜けた。**§9-442 の見積もり(10-25 + 5-15 + 10-25)より短い。

## ★★★これで因子機構が開く

mathlib は Dedekind 環の因子機構を完備している——
`IsDedekindDomain.HeightOneSpectrum`・`count`・`finprod_heightOneSpectrum_factorization`・
adic 付値。★`div(f)` は**主分数イデアルの素分解**としてそのまま出る。

★★次は `div(f_P ∘ [n])` が `n` で割れることの因子計算である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_map_eq_span` | ★★★★局所化の極大イデアルは単項 |
| `isDedekindDomain_coordinateRing` | ★★★★★★★**`F[W]` は Dedekind 環** |
| `isIntegrallyClosed_coordinateRing` | ★★★★★**`F[W]` は整閉**(§9-442 で欠けていた 1 件) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F]

set_option maxHeartbeats 1000000 in
/-- ★★★★**0 でない素イデアルの局所化における極大イデアルは単項である**。

★`Ψ₂Sq(c)` が 0 かどうかで第 136 の 2 つの補題に振り分ける。 -/
theorem exists_map_eq_span [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    (P : Ideal W.CoordinateRing) (hP : P ≠ ⊥) [hPp : P.IsPrime]
    (S : Type) [CommRing S] [Algebra W.CoordinateRing S] [IsLocalization P.primeCompl S] :
    ∃ π : W.CoordinateRing,
      Ideal.map (algebraMap W.CoordinateRing S) P
        = Ideal.span {algebraMap W.CoordinateRing S π} := by
  obtain ⟨c, s, hs, hz, hspan⟩ := exists_prime_gens h2 W P hP
  by_cases hs0 : s = 0
  · subst hs0
    refine ⟨genZ W, map_eq_span_genZ h2 W P S ?_ ?_ ?_ hspan⟩
    · rw [← hs]; ring
    · exact hspan ▸ Ideal.subset_span (by simp)
    · rwa [map_zero, sub_zero] at hz
  · refine ⟨genX W - algebraMap F W.CoordinateRing c,
      map_eq_span_genX h2 W P S hs hs0 ?_ hz hspan⟩
    exact hspan ▸ Ideal.subset_span (by simp)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**座標環は Dedekind 環である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`IsDedekindDomainDvr`(Noether + 局所化がすべて DVR)を経由する。
★★DVR 判定は `IsDiscreteValuationRing.TFAE` の項 0 ⟺ 項 4(極大イデアルが単項)。
★★★これで mathlib の**因子機構が丸ごと使えるようになる**。 -/
theorem isDedekindDomain_coordinateRing [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] : IsDedekindDomain W.CoordinateRing := by
  haveI hdvr : IsDedekindDomainDvr W.CoordinateRing := by
    constructor
    intro P hP hPp
    haveI := hPp
    have hmapne : Ideal.map (algebraMap W.CoordinateRing (Localization.AtPrime P)) P ≠ ⊥ := by
      obtain ⟨a, haP, ha0⟩ := (Submodule.ne_bot_iff P).1 hP
      intro hbot
      have hinj : Function.Injective
          (algebraMap W.CoordinateRing (Localization.AtPrime P)) :=
        IsLocalization.injective _ (Ideal.primeCompl_le_nonZeroDivisors P)
      have hz : algebraMap W.CoordinateRing (Localization.AtPrime P) a = 0 := by
        have hmem := Ideal.mem_map_of_mem
          (algebraMap W.CoordinateRing (Localization.AtPrime P)) haP
        rwa [hbot, Ideal.mem_bot] at hmem
      exact ha0 (hinj (by rw [hz, map_zero]))
    have hnf : ¬ IsField (Localization.AtPrime P) := by
      rw [IsLocalRing.isField_iff_maximalIdeal_eq, ← Localization.AtPrime.map_eq_maximalIdeal]
      exact hmapne
    refine ((IsDiscreteValuationRing.TFAE (Localization.AtPrime P) hnf).out 4 0).1 ?_
    rw [← Localization.AtPrime.map_eq_maximalIdeal]
    obtain ⟨π, hπ⟩ := exists_map_eq_span h2 W P hP (Localization.AtPrime P)
    exact ⟨⟨_, hπ⟩⟩
  exact IsDedekindDomainDvr.isDedekindDomain _

set_option maxHeartbeats 1000000 in
/-- ★★★★★**座標環は整閉である**——§9-442 で「これだけが欠けている」と測った 1 件。

★Dedekind 性の副産物として出る。 -/
theorem isIntegrallyClosed_coordinateRing [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] : IsIntegrallyClosed W.CoordinateRing := by
  haveI := isDedekindDomain_coordinateRing h2 W
  infer_instance

/-! ## ★出典の紐付け(`.src`) -/

def isDedekindDomain_coordinateRing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——座標環が Dedekind 環であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
