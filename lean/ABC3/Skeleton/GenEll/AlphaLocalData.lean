/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.ImageSL2FromH2H1
import ABC3.Found.GenEll.BadPrimeFromMultRed
import ABC3.Meta.Claim

/-!
# 第 1352 ブロック —— **悪い素点の惰性は幂単かつ非自明**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——`α` の側の**残り 1 本を節点にする**

`EllModuliWitness` の `imageContainsSL2J_torsionExt`（#4）は
`imageContainsSL2J_of_h2_h1`（第 1321）で **`∃ σ, h2 ∧ h1` 1 本**に落ちている:

* `h2`——`σ` は `T_l E` の上で `mod l` で**幂単**（`(σ−1)² ≡ 0`）
* `h1`——`σ` は `mod l` で**非自明**（`σ ≠ 1`）

☆これは古典的には「乗法還元の素点の**惰性**が Tate 曲線の一意化で
`(1 1; 0 1)` の形に作用する」ことである。

## ★★★★★★★★何が残っているのか——測定（2026-09-02、第 1352）

★道は `exists_h2_h1_of_bad_prime`（第 1320、**証明済み**）で**一本道**である。
そこに渡す局所のデータのうち

| 枚 | 内容 | 状態 |
|---|---|---|
| 1 | 悪い素点 `p` と `v_p(j) < 0` | ★第 1350（証明済み） |
| 2 | `hp`（付値の両立） | ★`valuation_algebraMap_adicCompletion`（在庫） |
| 3 | `E ⊗ L_p` の極小性 | ★第 1351（証明済み） |
| 4 | `E ⊗ L_p` の**乗法**還元 | ★第 1351（証明済み） |
| 5 | **分裂**乗法還元 | ☆非分岐 2 次拡大の段（第 1025-1043 に材料） |
| 6 | `ζ_l ∈ L_p` | ☆円分拡大の段 |

★★★**残るのは 5 と 6 の 2 枚だけ**である。

★★☆**訂正（2026-09-02、第 1354）**——この 2 枚を一度「局所の配管だけ」と書いたが、
**それは誤りである**。どちらも `L_p` の**有限拡大**を取る段であり、
☆**mathlib に「完備離散付値環の有限拡大は再び完備離散付値環」が無い**
（2026-09-02 確認: `Mathlib/RingTheory/DiscreteValuationRing/` は `Basic`・`TFAE` の 2 本だけ、
`Mathlib/NumberTheory/LocalField/Basic.lean` の `IsNonarchimedeanLocalField` は
**他のどのファイルからも使われていない**）。
★したがって `AdjoinRoot f` や `L_p(ζ_l)` に
`IsDiscreteValuationRing`・`IsFractionRing`・`IsAdicComplete` を与える段が丸ごと無い。
☆`ResearchPaper/blocked-leaves.json` に記録した。

★★☆**精査（2026-09-02、第 1355）**——「何も無い」ではなく
**解析の側はある／代数の側の橋が無い**が正確である。

☆mathlib に**ある**もの: `Analysis/Normed/Unbundled/SpectralNorm.lean`
（`spectralNorm`・`spectralNorm_unique_field_norm_ext`・`spectralMulAlgNorm`・
完備体の有限拡大の `CompleteSpace` インスタンス）、`Field/Krasner.lean`、
`Unbundled/FiniteExtension.lean`。★つまり**norm の一意な延長は取れている**。

★**無い**のはそこから `IsDiscreteValuationRing`・`IsFractionRing`・`IsAdicComplete`
へ渡す橋（値群の離散性＝分岐指数の有限性）である。
☆`Topology/Algebra/Valued/LocallyCompact.lean` に
「整数環がコンパクト ↔ 完備かつ離散付値環かつ剰余体が有限」があるので、
**局所コンパクト性を経由する道**が最短かもしれない。

★★★**欠けている primitive を 1 本に特定した（2026-09-02、第 1356）**

> **完備なイデアルに沿って幂等元が持ち上がる**
> （`IsAdicComplete I R` かつ `I ≤ Jacobson R` なら `R/I` の幂等元は `R` へ持ち上がる。
> 同値に「完備局所環の上の加群有限代数は局所環の有限直積」）

☆mathlib の `RingTheory/Idempotents.lean` は**幂零核に沿った持ち上げ**しか持たない
（`lift_of_isNilpotent_ker`・`existsUnique_isIdempotentElem_eq_of_ker_isNilpotent`）。

★これさえ入れば道は繋がる:
`R′ ≔ integralClosure R L` は `IsIntegralClosure.finite`（加群有限）と
`IsIntegralClosure.isDedekindDomain`（Dedekind）を持ち、`R′/m_R R′` は剰余体上有限次元
＝Artin 環。☆`R′` は整域なので非自明な幂等元を持たず、よって `R′/m_R R′` は局所、
`m_R R′ ⊆ Jacobson` なので `R′` は局所。
★最後に `AdicCompletion/LocalRing.lean` の `isLocalRing_of_isAdicComplete_maximal` と
合わせて DVR 構造が出る。

☆5 の材料（`Found/GaloisRep/UnramQuad.lean`）は
`valuation_algebraMap_ext`・`isMinimal_baseChange_ext`・`hasMultiplicativeReduction_ext`・
`hasSplitMultiplicativeReduction_ext` まで揃っているが、
`AdjoinRoot f` が離散付値環であること（`IsDomain`・`IsLocalRing`・`IsDiscreteValuationRing`）
はまだ仮説である。★6 は `L_p(ζ_l)/L_p` が `p ∤ l` なら非分岐であることを使う。

## ★★★★★★★★★★★★★★★★★★★★mathlib の欠落は塞がった（2026-09-02、第 1365-1372）

★上の「分岐する拡大の理論が mathlib に無い」は**自前で作って塞いだ**。

| 第 | 内容 |
|---|---|
| 1357-1364 | 完備なイデアルに沿う幂等元の持ち上げから 8 段 |
| 1365 | 整拡大では `m_R·S ≤ jacobson ⊥` |
| **1366** | ★**完備 DVR の有限分離整閉包は DVR**（`isDiscreteValuationRing_integralClosure`） |
| 1367-1368 | `m_A C`-進完備 → `m_C`-進完備、`IsFractionRing C L` |
| 1369 | 分岐指数 `e`（`m_A C = m_C^e`、`1 ≤ e ≤ [L:K]`） |
| 1370 | `hp` の分岐版 `hpe`（付値指数は `e` 倍） |
| 1371 | `hcop` の出どころを分岐版に（`l ∤ e` なら通る） |
| **1372** | ★**`exists_h2_h1_of_bad_prime_ram`**——`L_v` は `L_p` の任意の有限拡大でよい |

★★★**測定の鍵（第 1370）**——`exists_h2_h1_of_bad_prime` の証明で
不分岐の仮定 `hp` が使われるのは `vAdd_algebraMap_eq_valAdd` の**ただ 1 箇所**である。
☆他（`IsMinimal`・`HasSplitMultiplicativeReduction`・`IsElliptic`）は
インスタンス引数なので呼び出し側が別途与える。

★★残るのは**組み立ての配管**であって mathlib の欠落ではない:

1. `L_v′/L_v` の付値の延長公式 `v_{L_v′}(y) = v_{L_v}(y)^e`（`y ∈ L_v`）
   ——第 1369 の `m_A C = m_C^e` から出るはず。
2. 拡大の上での `IsMinimal`・`IsElliptic`・`HasSplitMultiplicativeReduction`
   ——不分岐版は `Found/GaloisRep/UnramQuad.lean`（第 1029-1033）にある。
3. `ζ_l ∈ L_v′` の構成（`L_v(ζ_l)` は `≤ l−1` 次なので `l ∤ e`）。
4. `SSCurve` の言葉（`E.alg ≃ₐ[E.fld] E.alg`）への翻訳。

## ★★★★★★★★★★★★6（`ζ_l ∈ L_p`）の道の測定（2026-09-02、第 1353）

☆道は 2 つある。

**(A) 局所を伸ばす**——`L_p(ζ_l)` を取る。
★`hp`（付値の両立）が保たれるのは**非分岐のとき**、すなわち `p ∤ l` のときである。
☆`ℚ(ζ_l)/ℚ` は `l` の上でしか分岐しないので、`p` が `l` の上に無ければよい。

**(B) 大域を伸ばす**——`L′ ≔ L(ζ_l)` を取る。
★`Gal(L̄/L′) ⊆ Gal(L̄/L)` なので、`L′` の上で見つけた `σ` は `L` の上の `σ` でもある
（`h2`・`h1` は `T_l E` の上の条件なのでそのまま伝わる）。
☆この場合も `L′` の素点 `p′` が `p` の上で非分岐であること（＝`p ∤ l`）が要る。

★★★**どちらの道でも要るのは「`l` の上に無い悪い素点が 1 つある」ことである**。
☆`Theorem 3.8` の条件 (a) は `l ≥ 100·d·(ht^Falt + C·d^ε)` を与えるので
`l` は十分大きいが、「悪い素点の剰余標数より大きい」ことを出すには
導手と `ht^Falt` の関係を経由する段が要る——**それが次の節点である**。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Interface.GaloisRep ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点の惰性は `T_l E` の上で `mod l` で幂単かつ非自明**（第 1352）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆Tate 曲線の一意化で `σ` が `(1 1; 0 1)` の形に作用することの言い換えである。

★★★これが `imageContainsSL2J_torsionExt`（#4）に残る**ただ 1 つの節点**である。 -/
theorem exists_h2_h1_unipotent (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (hm : E.HasMultRed) (hpr : E.PrimeToLocalHeights l) :
    ∃ σ : E.alg ≃ₐ[E.fld] E.alg,
      (∀ x : E.tate l, ∃ u : E.tate l,
          galTate E.W l σ (galTate E.W l σ x) + x
            = galTate E.W l σ x + galTate E.W l σ x + l • u) ∧
        (∃ x : E.tate l, ∀ u : E.tate l, galTate E.W l σ x ≠ x + l • u) := by
  sorry

/-- ★★★★★★★★★★★★★★★★★★★★
**`ImageContainsSL2J` はこの 1 本から出る**——★（第 1352）。

☆`imageContainsSL2J_of_h2_h1`（第 1321、無条件）に渡すだけである。 -/
theorem imageContainsSL2J_of_multRed (E : SSCurve) (l : ℕ) (hl : Nat.Prime l)
    (hl5 : 5 ≤ l) (hm : E.HasMultRed) (hpr : E.PrimeToLocalHeights l)
    (hno : ¬ HasLCyclicJ E l) : ImageContainsSL2J E l := by
  haveI : Fact l.Prime := ⟨hl⟩
  exact imageContainsSL2J_of_h2_h1 E l hl5 hno (exists_h2_h1_unipotent E l hm hpr)

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_h2_h1_unipotent.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(悪い素点の惰性は T_l E の上で mod l で幂単かつ非自明)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_unipotent.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_h2_h1_of_bad_prime(第 1320、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_bad_prime") 1,
    .citation "[ABC3]" "SSCurve.exists_local_multRed(第 1351、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.SSCurve.exists_local_multRed") 1,
    .citation "[ABC3]" "hasSplitMultiplicativeReduction_ext(第 1033、証明済み。DVR 構造は仮説)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hasSplitMultiplicativeReduction_ext") 1,
    .citation "[ABC3]" "isDiscreteValuationRing_integralClosure(第 1366、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isDiscreteValuationRing_integralClosure") 1,
    .citation "[ABC3]" "exists_h2_h1_of_bad_prime_ram(第 1372、証明済み。分岐版)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_bad_prime_ram") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1365-1372）**——**mathlib の欠落は塞がった**。" ++
       "★完備 DVR の有限分離整閉包は完備 DVR（第 1366-1368）、" ++
       "分岐指数は `e ≤ [L_v′:L_v]`（第 1369）、" ++
       "`hp` の分岐版で `h2`・`h1` が出る（第 1370-1372）。" ++
       "☆残るのは付値の延長公式・拡大の上の局所データ・`ζ_l` の構成・" ++
       "`SSCurve` の言葉への翻訳の**配管**だけである。") 19 ]

def imageContainsSL2J_of_multRed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ImageContainsSL2J はこの 1 本から出る)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_of_multRed.needs : List ProofObligation :=
  [ .citation "[ABC3]" "imageContainsSL2J_of_h2_h1(第 1321、無条件)"
      (.inProject "ABC3" "ABC3.Found.GenEll.imageContainsSL2J_of_h2_h1") 1,
    .citation "[ABC3]" "exists_h2_h1_unipotent(本ファイルの節点)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_h2_h1_unipotent") 1 ]

end ABC3.Skeleton.GenEll
