/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.DetCycloChar
import ABC3.Found.GenEll.AlphaUnipotent
import ABC3.Found.GenEll.AlphaMemImage
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

## ★★★葉の一覧（★2026-09-06 に葉 3 も閉じたので、残る葉は 0）

| # | 命題 | 原文の言い方 |
|---|---|---|
| 3 | `α = (1 1 / 0 1)` が mod `l` 像に入る（★2026-09-06 に閉じた） | 『by the local theory (cf. the discussion preceding Lemma 3.2)』 |
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
`l ∤ v(q)` なら `q` の `l` 乗根の抽出が非自明で、`α = (1 1 / 0 1)` が像に入る。

## ★★★★★★2026-09-06: `sorry` は消えた —— `Found/GenEll/` からの配線 2 行

* `SSCurve.exists_h2_h1_unipotent_of_multRed`(`Found/GenEll/AlphaUnipotent.lean`、`sorry` 0)
  ——仮定 `(hmult) (hpr)` が本定理と完全に一致する。
* `alpha_mem_map_of_galTate`(`Found/GenEll/AlphaMemImage.lean`、`sorry` 0)
  ——結論が本定理と逐語同一である。

## ★★★★逸脱の記録(2026-09-06) —— 仮引数 `e` を足した

`alpha_mem_map_of_galTate` は基底 `e : E.tate l ≃+ (Fin 2 → ℤ_[l])` を
**入力として受け取る**。これを与えないと、本定理は
「`T_l E` が階数 2 の自由 `ℤ_l` 加群である」ことを**暗に主張する**ことになり、
その事実を与える宣言は本プロジェクトに 0 件である(2026-09-06 測定)。

★そこで本定理に `e` を**仮引数として足した**。CLAUDE.md の「逸脱」の条
(「後続の証明に影響が出ないならば前提の追加を許容する」)に従う。

☆**影響が無いことの測定(2026-09-06)**: 消費側 `ImageContainsSL2J` の定義
(`Found/GenEll/EllModuliGalois/Theorem38.lean:52`)は `∀ e, …` の形なので、
基底 `e` は**消費側がゴールから供給する**。生成側では要求しなくてよい。
また本定理を Lean の項として消費している宣言は木の中に 0 件である
(`EllModuliWitness.lean:623` は `.needs` の文字列参照のみで、名前は変えていない)。

☆**第 826 の訂正との関係**: 「`∀ e` → `∃ e`」という訂正は
**消費側の仮定としては正しい**(Tate 一意化に適合した基底を取るのだから、
任意の基底で `α` が出るとは言えない)。しかし `e` を**入力に持たない
生成側の主張**に写すと、上の通り「階数 2 の自由性」まで主張してしまい過剰になる。
★結論の `∃ e₀` はそのまま残してあるので、訂正の内容自体は保たれている。 -/
theorem alpha_in_modl_image (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (hmult : E.HasMultRed) (hpr : E.PrimeToLocalHeights l)
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) :
    ∃ e₀ : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e₀).range).map (glRedPadic l) := by
  obtain ⟨σ, h2, h1⟩ := E.exists_h2_h1_unipotent_of_multRed l hmult hpr
  exact alpha_mem_map_of_galTate E l e σ h2 h1

def alpha_in_modl_image.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所理論の行列表示——α = (1 1 / 0 1) が mod l 像に入る)",
    sectionId := "genell-thm-3-8" }

def alpha_in_modl_image.needs : List ProofObligation :=
  [ .citation "[FC]" "Degenerations of Abelian Varieties, Ch. III, Cor. 7.3（完全列 0 → ᴽ_l(1) → M_l(E) → ᴽ_l → 0）"
      (.absent "mathlib に Tate 曲線・Tate twist・M_l(E) はいずれも無い(2026-08-16 実測)。★2026-09-06 に再測: re:`TateCurve|tateCurve|Tate curve|tateTwist|TateTwist|Tate twist|TateModule`→0") 8,
    .citation "[ABC3]" "tatePhi_map（Tate 一意化 Φ は σ で同変、証明済み）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhi_map") 1,
    .citation "[ABC3]" "sl2_range_basis_transfer（SL₂ ⊆ 像 は基底に依らない、第 825）"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sl2_range_basis_transfer") 1,
    .citation "[ABC3]" "SSCurve.exists_h2_h1_unipotent_of_multRed（悪い素点の惰性は mod l で冪単かつ非自明、第 1384、sorry 0）"
      (.inProject "ABC3" "ABC3.Found.GenEll.SSCurve.exists_h2_h1_unipotent_of_multRed") 19,
    .citation "[ABC3]" "alpha_mem_map_of_galTate（冪単・非自明から α が mod l 像に入る、第 1237、sorry 0）"
      (.inProject "ABC3" "ABC3.Found.GenEll.alpha_mem_map_of_galTate") 19,
    .implicitStep
      ("★★★訂正（2026-08-31、第 826）: 以前は基底 `e` を**任意に与えられた**形であったが、"
       ++ "それでは**循環する**。★Tate 一意化に適合した基底 (ζ_l, q^{1/l}) を取るのだから"
       ++ "**`∃ e`** が正しい形である。SL₂ は正規部分群なのでこれで十分である。"
       ++ "☆★逸脱（2026-09-06）: 結論の `∃ e₀` はそのままに、"
       ++ "基底 `e` を**仮引数として足した**。足さないと T_l E の階数 2 の自由性を"
       ++ "暗に主張してしまい、それを与える宣言は木に 0 件だからである。"
       ++ "消費側 ImageContainsSL2J は `∀ e` なので影響が無いことを測定した") 1,
    .implicitStep
      ("★★機構: 乗法還元の素点で E は Tate 曲線 E_q であり、"
       ++ "E[l] = ⟨ζ_l, q^{1/l}⟩。σ(ζ) = ζ^χ, σ(q^{1/l}) = ζ^a q^{1/l} なので"
       ++ "この基底で mod l 像は上三角。★l ∤ v(q) なら q^{1/l} ∉ K(ζ_l) なので"
       ++ "χ(σ) = 1 かつ a(σ) ≠ 0 なる σ があり、その冪を取れば α = (1 1 / 0 1) が出る") 6,
    .implicitStep
      ("★★★測定（2026-08-31、第 828）——**局所→大域の橋が要る**。"
       ++ "galRep は大域体 E.fld の上の表現であり、上の機構は完備体の上の話である。"
       ++ "★分解群 Gal(L̄_p/L_p) → Gal(L̄/L) と Tate 加群の互換性が要る。"
       ++ "☆本プロジェクトにも mathlib にも decompositionGroup の機構は無い（grep で 0 件）。"
       ++ "★★★2026-09-06: この橋は Found/GenEll/AlphaUnipotent.lean（第 1365-1384）で"
       ++ "渡り終えており、本節点の `sorry` は消えた") 8 ]

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
