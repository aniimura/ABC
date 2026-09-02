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

★★★**残るのは 5 と 6 の 2 枚だけ**であり、どちらも**局所の配管**である
——未証明の外部引用は 1 本もない。

☆5 の材料（`Found/GaloisRep/UnramQuad.lean`）は
`valuation_algebraMap_ext`・`isMinimal_baseChange_ext`・`hasMultiplicativeReduction_ext`・
`hasSplitMultiplicativeReduction_ext` まで揃っているが、
`AdjoinRoot f` が離散付値環であること（`IsDomain`・`IsLocalRing`・`IsDiscreteValuationRing`）
はまだ仮説である。★6 は `L_p(ζ_l)/L_p` が `p ∤ l` なら非分岐であることを使う。

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
    .implicitStep
      ("★★★★**2026-09-02（第 1352）**——道は第 1320 で一本道であり、" ++
       "残るのは**分裂乗法還元**（非分岐 2 次拡大）と **`ζ_l ∈ L_p`**（円分拡大）の" ++
       "2 枚の局所の配管だけである。☆未証明の外部引用は 1 本もない。") 19 ]

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
