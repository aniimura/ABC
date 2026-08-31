/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.DetCycloChar
import ABC3.Meta.Claim

/-!
# `Theorem 3.8` に残る Galois 側の葉（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19–p.22。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★これは何か

`Found/GenEll/EllModuliGalois.lean` の

* `imageContainsSL2J_of_alpha`（`§9-1188`、第 762）
* `imageSurjectiveJ_of_containsSL2`（`§9-1187`、第 761）

は、`EllModuliData` の `imageContainsSL2_of_torsionExt` 欄と
`imageSurjective_of_containsSL2` 欄を**それぞれ 1 本の仮説に帰着**した。
本ファイルはその仮説を**依存グラフの節点として固定する**。

## ★★★残る葉

| # | 命題 | 原文の言い方 |
|---|---|---|
| 3 | `α = (1 1 / 0 1)` が mod `l` 像に入る | 『by the local theory (cf. the discussion preceding Lemma 3.2)』 |
| 4 | 円分指標が **`mod l^n` で**全射 | 『`ℚ(ζ_{l^∞})/ℚ` は `l` で完全分岐するので `L/ℚ` と線型無関連』（★逆極限の段は第 773 で済んだ） |
| 5 | `galRep` は連続 | （原文は位相を明示しない。★像が閉であることは連続性に帰着済み——第 765） |

☆どれも**群論ではない**——群論の核（`Lemma 3.1, (iv)`）は
`Found/GenEll/Sl2Padic.lean`・`Thm38Bridge.lean` に**実装済み**である。

## ★★材料

* 葉 3: `Found/GaloisRep/Lemma32Tate.lean`（Tate 一意化と `Lemma 3.2, (i)`）
* 葉 4: mathlib の `IsCyclotomicExtension.autEquivPow`
  （`Irreducible (cyclotomic n K)` から `Gal(L/K) ≃* (ZMod n)ˣ`）。
  ★残るのは `L` の上での `cyclotomic (l^n)` の既約性——`l` が `L` で不分岐であることから
  `L ∩ ℚ(ζ_{l^n}) = ℚ` が出る、という古典的な議論である。
* 葉 5: mathlib の `Mathlib/FieldTheory/KrullTopology.lean`
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Interface.GaloisRep
open Matrix Matrix.SpecialLinearGroup
open scoped MatrixGroups Classical

/-- **[GenEll] `Theorem 3.8` の局所理論の側**——`α` が mod `l` 像に入る。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★★★`Found/GenEll/EllModuliGalois.lean` の `imageContainsSL2J_of_alpha` が
**そのまま消費する形**である。

☆機構: 乗法還元の素点で `E` は Tate 曲線 `E_q` であり、完全列
`0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0` に合わせた基底で mod `l` 像は上三角になる。
`l ∤ v(q)` なら `q` の `l` 乗根の抽出が非自明で、`α = (1 1 / 0 1)` が像に入る。 -/
theorem alpha_in_modl_image (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (hmult : E.HasMultRed) (hpr : E.PrimeToLocalHeights l) :
    ∃ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e).range).map (glRedPadic l) := by
  sorry

def alpha_in_modl_image.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所理論の行列表示——α = (1 1 / 0 1) が mod l 像に入る)",
    sectionId := "genell-thm-3-8" }

def alpha_in_modl_image.needs : List ProofObligation :=
  [ .citation "[FC]" "Degenerations of Abelian Varieties, Ch. III, Cor. 7.3(完全列 0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0)"
      (.absent "mathlib に Tate 曲線・Tate twist・M_l(E) はいずれも無い(2026-08-16 実測)") 10,
    .implicitStep
      ("★完全列に合わせた基底で mod l 像は上三角になり、l ∤ v(q) なら " ++
       "q の l 乗根の抽出が非自明なので α が像に入る。" ++
       "☆材料は Found/GaloisRep/Lemma32Tate.lean(Tate 一意化と Lemma 3.2, (i))") 10,
    .citation "[ABC3]" "sl2_range_basis_transfer（SL₂ ⊆ 像 は基底に依らない、第 825）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sl2_range_basis_transfer") 1,
    .implicitStep
      ("★★★訂正（2026-08-31、第 826）: 以前は基底 `e` を**任意に与えられた**形であったが、"
       ++ "それでは**循環する**——α がすべての共役に入ることは SL₂ ⊆ 像 を既に"
       ++ "知っていないと言えない。★Tate 一意化に適合した基底 (ζ_l, q^{1/l}) を取るのだから"
       ++ "**`∃ e`** が正しい形である。SL₂ は正規部分群なのでこれで十分である") 1 ]

/-! ## ★★★★★★★★★★葉 5 は閉じた（2026-08-31、第 766-774）

★`galRep` の連続性と像の閉性は `Found/GenEll/GalRepContinuity.lean`・
`Found/GenEll/GalRepClosed.lean` で**無条件に証明された**ので、
本ファイルからは削除した。
-/

/-! ## ★★★★★★★★★★★★★★★★★★★★葉 4 も閉じた（2026-08-31、第 780-784）

★円分指標の全射性は `Found/GenEll/CycloDisjoint.lean`・
`Found/GenEll/DetCycloChar.lean` の `cyclotomic_det_surjective` で
**代数的な仮説 `¬ (l : ℤ) ∣ disc L`（= `l` が `L` で不分岐）だけから
証明された**ので、本ファイルからは削除した。

☆道筋は「完全分岐」を一度も使わない:
`disc ℚ(ζ_{l^k}) = ±l^m` と `l ∤ disc L` から判別式が互いに素で、
`NumberField.linearDisjoint_of_isGalois_isCoprime_discr` が線型無関連を出す。
-/

end ABC3.Skeleton.GenEll
