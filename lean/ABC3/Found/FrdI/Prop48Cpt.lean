/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop33UnTr

/-!
# [FrdI] Proposition 4.8, (iii) —— Frobenius-compact 対象の移送(下ごしらえ)

原文 (FrdI p.88):
> (iii) If C is of rationally standard type, then (Cistr)birat is of standard

★★`Frobenius-compact` は `𝒪^×` についての 3 条である。
★`Remark 4.5.1` が与えるのは `((𝒞^istr)^un-tr)^birat` の compact 対象で、
`Proposition 4.8, (iii)` が要求するのは `(𝒞^istr)^birat` のもの ——
**向きが逆**なのが `iii-b-compact` の障害だった。

## ★★★測って分かったこと(2026-08-19)

第 2 条(無限位数の単元が存在する)は**準同型の逆向きに落ちる**
(`infiniteOrder_of_map`、`Prop48Sec.lean`)。
したがって **`𝒪^×` の写像が全射でありさえすれば第 2 条は渡る**。

★その全射性は `Proposition 4.4, (ii)` の記述
(`𝒪^×(A^birat)` は `𝒪^▷(A)^gp` の像)を通し、
**`𝒪^▷(𝒞^istr) ↠ 𝒪^▷((𝒞^istr)^un-tr)`** に帰着する。
★本ファイルはその**最初の 1 本**を取る。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★★★★**`𝒪^▷` は `𝒞^istr → 𝒞^un-tr` で全射**。

★★`𝒞^un-tr` は射を `𝔽_Φ` の像で同一視した商なので、
`Base` も `degFr` も**定義的に一致する** ——
したがって `istrToUnTr` の充満性(在庫)で代表元を取れば、
それは自動的に `𝒪^▷` に入る。 -/
theorem otri_unTr_surjective (Fc : FrobenioidCore P) (A : UnTr P) (β : End A)
    (hβ : β ∈ OTri (unTrPre P Fc) A) :
    ∃ β₀ : End (show Istr P from A),
      β₀ ∈ OTri (istrPre P Fc) (show Istr P from A) ∧ (istrToUnTr P).map β₀ = β := by
  obtain ⟨β₀, hβ₀⟩ := (istrToUnTr P).map_surjective β
  refine ⟨β₀, ⟨?_, ?_⟩, hβ₀⟩
  · have h : (unTrPre P Fc).Base ((istrToUnTr P).map β₀)
        = (unTrPre P Fc).Base (𝟙 A) := by
      rw [hβ₀]
      exact hβ.1
    exact h
  · have h : (unTrPre P Fc).degFr ((istrToUnTr P).map β₀) = 1 := by
      rw [hβ₀]
      exact hβ.2
    exact h

def otri_unTr_surjective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — 𝒪^▷ は 𝒞^istr → 𝒞^un-tr で全射",
    sectionId := "frdi-prop-4-8" }

/-! ## ★★★`Remark 4.5.1` —— **在庫にあった**(在庫の見落とし 11 件目)

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is

★★★`istr_standardType`(`Rmk451.lean:227`)が
**`IsOfStandardType` の 6 フィールドすべてを `𝒞^istr` へ移して完成している**。

★私は `frobIsotropic` / `frobNormalized` / group-like の降下を
**手で書き下してしまった** —— どれも
`istr_frobIsotropicType` / `istr_frobNormalizedType` / `istr_groupLikeType_down`
として在庫にあった。★重複は撤去した。

## ★★対策(記憶に追記済み)

`Remark 4.5.1` を実装しようとするとき、**ファイル名が `Rmk451.lean`** なのだから
そこを最初に開くべきだった。
★★**「原典の番号でファイルを引く」**を検索の第 3 の手として加える
(第 1: `prop_<番号>_.*_of`、第 2: `<構成の接頭辞>_<性質>`)。

★測ったこと自体は無駄ではなく、
`frobIsotropic` が **`𝒞^istr` では自明**(`Dd := A`、`φ := 𝟙 A`)であることは
在庫の `istr_frobIsotropicType` の中身と一致していた。 -/

end ABC3.Found.FrdI
