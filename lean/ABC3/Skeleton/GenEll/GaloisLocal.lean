/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Meta.Claim

/-!
# `Theorem 3.8`・`Corollary 4.3` に残る Galois 側の 3 つの葉（`Skeleton`）

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

## ★★★3 つの葉

| # | 命題 | 原文の言い方 |
|---|---|---|
| 3 | `α = (1 1 / 0 1)` が mod `l` 像に入る | 『by the local theory (cf. the discussion preceding Lemma 3.2)』 |
| 4 | 円分指標 `Gal(L̄/L) → ℤ_l^×` が全射 | 『`ℚ(ζ_{l^∞})/ℚ` は `l` で完全分岐するので `L/ℚ` と線型無関連』 |
| 5 | 像は閉部分群 | （原文は位相を明示しない。profinite 群の連続像である） |

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
    (hmult : E.HasMultRed) (hpr : E.PrimeToLocalHeights l)
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) :
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
       "☆材料は Found/GaloisRep/Lemma32Tate.lean(Tate 一意化と Lemma 3.2, (i))") 10 ]

/-- **[GenEll] `Corollary 4.3` の円分指標の側**——`l` が `L` で不分岐なら全射。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★★★★★`Found/GenEll/EllModuliGalois.lean` の `imageSurjectiveJ_of_containsSL2` が
**そのまま消費する形**である。

☆`det ρ(σ)` が円分指標であること自体は `det_cyclotomic_full`（★無条件）で済んでいる。 -/
theorem cyclotomic_det_surjective (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (hunram : True)
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (u : ℤ_[l]ˣ) :
    ∃ σ : E.alg ≃ₐ[E.fld] E.alg,
      Matrix.GeneralLinearGroup.det (galRep E.W l e σ) = u := by
  sorry

def cyclotomic_det_surjective.src : Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(円分指標 Gal(L̄/L) → ℤ_l^× の全射性——l が L で不分岐のとき)",
    sectionId := "genell-cor-4-3" }

def cyclotomic_det_surjective.needs : List ProofObligation :=
  [ .citation "[mathlib]" "IsCyclotomicExtension.autEquivPow(Irreducible (cyclotomic n K) から Gal ≃* (ZMod n)ˣ)"
      (.inMathlib "IsCyclotomicExtension.autEquivPow") 2,
    .implicitStep
      ("★★L の上で cyclotomic (l^n) が既約であること。" ++
       "l が L で不分岐かつ ℚ(ζ_{l^n})/ℚ が l で完全分岐なので L ∩ ℚ(ζ_{l^n}) = ℚ" ++
       "——古典的だが mathlib に直接の形は無い") 8,
    .implicitStep "★各 n での全射性から ℤ_l^× への全射性を出す逆極限の段" 4 ]

/-- **[GenEll] `Theorem 3.8` の位相の側**——`galRep` の像は閉部分群。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`Lemma 3.1, (iv)`（`Found/GenEll/Sl2Padic.lean`）が閉部分群を要求するので要る。
★`Gal(L̄/L)` は Krull 位相で profinite（コンパクト）であり、`galRep` は連続なので像は閉。 -/
theorem galRep_range_isClosed (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) :
    IsClosed (((galRep E.W l e).range : Subgroup (GL (Fin 2) ℤ_[l])) :
      Set (GL (Fin 2) ℤ_[l])) := by
  sorry

def galRep_range_isClosed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galRep の像は閉部分群——profinite 群の連続像)",
    sectionId := "genell-thm-3-8" }

def galRep_range_isClosed.needs : List ProofObligation :=
  [ .citation "[mathlib]" "KrullTopology(Gal(L̄/L) の位相)"
      (.inMathlib "krullTopology") 2,
    .implicitStep
      ("★Gal(L̄/L) は Krull 位相で profinite(コンパクト)であり、galRep は連続。" ++
       "☆連続性は galRep が有限レベル(E[l^n] への作用)の逆極限であることから出る") 6 ]

end ABC3.Skeleton.GenEll
