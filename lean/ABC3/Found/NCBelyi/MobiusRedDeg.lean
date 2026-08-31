/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.DescendData
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★`ℚ`-Möbius は reduced degree を変えない（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> applying an automorphism (with rational coeﬃcients!) as in Lemma 2.3 and then

## ★★★★★★★★★★★★★★★★★★★★これは何か —— (b) を通すための鍵

`Lemma 2.4` に残っているのは **(b)（分離）** だけであり
（`Found/NCBelyi/Lemma24Package.lean` の測定）、
それを作る原文の段は

> applying an automorphism (with rational coefficients!) as in Lemma 2.3 and then
> multiplying by some positive rational number, we may assume that |α| ≤1, for all
> α ∈S\{∞}, while |β| ≥C

である。★★★**ここで `ℚ` 係数の自己同型を挟むので、入れ子帰納法の測度
（`redDeg`／`maxRedDeg`）が動かないことを先に確かめておく必要がある。**

★★本ファイルはそれを取る:

    **`ℚ⟮ν/(x−λ) + μ⟯ = ℚ⟮x⟯`**   （`λ, ν, μ ∈ ℚ`、`ν ≠ 0`、`x ≠ λ`）

したがって `redDeg` も `maxRedDeg` も**不変**である。

## ★機構 —— 逆 Möbius も `ℚ` 係数

* `ν/(x−λ) + μ ∈ ℚ⟮x⟯` —— 体の演算で閉じているから
* `x = ν/(y−μ) + λ ∈ ℚ⟮y⟯`（`y ≔ ν/(x−λ) + μ`）—— **逆写像も同じ形**

★★★原文が『**with rational coefficients!**』と感嘆符つきで書くのは、
まさにこの「測度を壊さない」ことのためである。

## ★これで `Lemma 2.4` の (b) に必要な部品は

| 部品 | 場所 |
|---|---|
| 正規化（`|α| ≤ 1`、`|β| ≥ C`） | `Found/NCBelyi/RatSeparation.lean` の `exists_rat_normalization` |
| 係数の限界 `≤ d₀^{d₀}` | `Found/NCBelyi/CoeffBound.lean` |
| 入れ子帰納法 | `Found/NCBelyi/NestedInduction.lean`／`DescendData.lean` |
| ★**測度が Möbius で不変** | ★★**本ファイル** |

★★残る仕事は `nested_induction_descend'` を `β` を運ぶ形に書き直すことだけである。
-/

namespace ABC3.Found.NCBelyi

open IntermediateField

/-! ## ★★★★★★★★★★★★★★★★★★★★生成する体は変わらない -/

/-- ★★★★★★★★★★★★★★★★★★★★**`ℚ⟮ν/(x−λ) + μ⟯ = ℚ⟮x⟯`**。

原文 (NCBelyi p.5):
> applying an automorphism (with rational coeﬃcients!) as in Lemma 2.3 and then

★`⊆` は体の演算で閉じていること、`⊇` は**逆写像も同じ形**（`x = ν/(y−μ) + λ`）
であることから出る。
★★原文の『**with rational coefficients!**』の感嘆符はこのためである。 -/
theorem adjoin_mobius_eq (lam nu mu : ℚ) (hnu : nu ≠ 0) (x : ℂ) (hx : x ≠ (lam : ℂ)) :
    ℚ⟮(nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ)⟯ = ℚ⟮x⟯ := by
  have hxl : x - (lam : ℂ) ≠ 0 := sub_ne_zero.mpr hx
  have hnuC : ((nu : ℂ)) ≠ 0 := by exact_mod_cast hnu
  have hratmem : ∀ (K : IntermediateField ℚ ℂ) (q : ℚ), ((q : ℂ)) ∈ K := by
    intro K q
    have h := K.algebraMap_mem q
    simpa using h
  refine le_antisymm ?_ ?_
  · rw [adjoin_simple_le_iff]
    have hxm : x ∈ ℚ⟮x⟯ := mem_adjoin_simple_self ℚ x
    exact add_mem (mul_mem (hratmem _ nu) (inv_mem (sub_mem hxm (hratmem _ lam))))
      (hratmem _ mu)
  · rw [adjoin_simple_le_iff]
    set y : ℂ := (nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ) with hy
    have hym : y ∈ ℚ⟮y⟯ := mem_adjoin_simple_self ℚ y
    have hne : y - (mu : ℂ) ≠ 0 := by
      rw [hy]
      simp only [add_sub_cancel_right]
      exact mul_ne_zero hnuC (inv_ne_zero hxl)
    have hxeq : x = (nu : ℂ) * (y - (mu : ℂ))⁻¹ + (lam : ℂ) := by
      rw [hy]
      simp only [add_sub_cancel_right]
      rw [mul_inv, inv_inv, ← mul_assoc, mul_inv_cancel₀ hnuC, one_mul]
      ring
    rw [hxeq]
    exact add_mem (mul_mem (hratmem _ nu) (inv_mem (sub_mem hym (hratmem _ mu))))
      (hratmem _ lam)

/-! ## ★★★★★★★★★★★★測度は不変 -/

/-- ★★★★★★★★★★★★**`ℚ`-Möbius は reduced degree を変えない**。 -/
theorem redDeg_mobius (lam nu mu : ℚ) (hnu : nu ≠ 0) (x : ℂ) (hx : x ≠ (lam : ℂ)) :
    redDeg ((nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ)) = redDeg x := by
  rw [redDeg, redDeg, adjoin_mobius_eq lam nu mu hnu x hx]

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`ℚ`-Möbius は `maxRedDeg` を変えない**
——入れ子帰納法の測度が正規化で動かない。 -/
theorem maxRedDeg_mobius (lam nu mu : ℚ) (hnu : nu ≠ 0) (S : Finset ℂ)
    (hlam : ((lam : ℂ)) ∉ S) :
    maxRedDeg (S.image (fun x => (nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ)))
      = maxRedDeg S := by
  classical
  rw [maxRedDeg, maxRedDeg, Finset.sup_image]
  refine Finset.sup_congr rfl (fun x hx => ?_)
  exact redDeg_mobius lam nu mu hnu x (fun h => hlam (h ▸ hx))

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def adjoin_mobius_eq.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(ℚ⟮ν/(x−λ)+μ⟯ = ℚ⟮x⟯——with rational coefficients! の中身)",
    sectionId := "ncbelyi-lemma-2-4" }

def redDeg_mobius.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(ℚ-Möbius は reduced degree を変えない)",
    sectionId := "ncbelyi-lemma-2-4" }

def maxRedDeg_mobius.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(ℚ-Möbius は maxRedDeg を変えない——測度が正規化で動かない)",
    sectionId := "ncbelyi-lemma-2-4" }

def maxRedDeg_mobius.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★★測定(2026-08-29): 原文が『applying an automorphism " ++
       "(with rational coefficients!) as in Lemma 2.3』と**感嘆符つきで**書くのは、" ++
       "その自己同型が**入れ子帰納法の測度を壊さない**ためである。" ++
       "★中身は ℚ⟮ν/(x−λ)+μ⟯ = ℚ⟮x⟯ であり、⊇ は逆写像 x = ν/(y−μ)+λ が" ++
       "**同じ形**であることから出る") 6,
    .implicitStep
      ("★★これで Lemma 2.4 の (b)(分離)に必要な部品はそろった: " ++
       "正規化(exists_rat_normalization)・係数の限界(CoeffBound.lean)・" ++
       "入れ子帰納法(NestedInduction.lean／DescendData.lean)・" ++
       "★測度の Möbius 不変性(本ファイル)。" ++
       "★★★残る仕事は nested_induction_descend′ を β を運ぶ形に書き直すことだけである") 6 ]

end ABC3.Found.NCBelyi
