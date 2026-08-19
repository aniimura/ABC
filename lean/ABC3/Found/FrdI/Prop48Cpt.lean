/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop33UnTr
import ABC3.Found.FrdI.Prop44Otri

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

/-! ## ★`Frobenius-compact` の第 1 条は `𝒞^birat` では**ただで出る**

原文 (FrdI p.23):
> A Frobenius-compact object of C is de

★★**原文の 3 条**(物理 p.23、実測):
1. `𝒪^×(C)` が可換
2. `𝒪^×(C)^pf ≠ 0` —— ★**抽出では `̸=` が `=` に化ける**(記憶の事故型そのもの)。
   第 3 条が `𝒪^×(C)^pf` の上の作用を語るので、`≠` でなければ意味を成さない
3. `Aut_𝒞(C)` の元が `𝒪^×(C)^pf` に `ℚ>0` の元の掛け算として作用するなら、実は自明に作用する

★我々の `IsFrobeniusCompact` の第 2 条(無限位数の単元が存在する)は
`𝒪^×(C)^pf ≠ 0` の言い換えであり、**読みは合っている**。

★★`IsFrobeniusCompact` の 3 条のうち**第 1 条(`𝒪^×` の可換性)**は、
`birationally Frobenius-normalized` 型のとき `otri_mul_comm` から直ちに出る ——
`𝒪^▷(A^birat)` が可換なら、その部分単系 `𝒪^×(A^birat)` も可換である。

★★★`Proposition 4.8, (iii)` の仮定は `rationally standard` 型であり、
そこに **`biratFrobNormalized` が入っている**(`Definition 4.5, (iii), (a)`)ので、
**この条は仮定から直接出る**。 -/

variable {P} {G : Frobenioid P}

/-- ★★★★**`𝒪^×(A^birat)` は可換**(birationally Frobenius-normalized 型のとき)。 -/
theorem birat_otimes_comm
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (A : BiratCat P G) (x y : End A)
    (hx : x ∈ OTimes (biratPre P G) A) (hy : y ∈ OTimes (biratPre P G) A) :
    x * y = y * x :=
  congrArg Subtype.val
    (otri_mul_comm (biratPre P G) (hfn A) ⟨x, hx.1⟩ ⟨y, hy.1⟩)

/-- ★★`Frobenius-compact` の第 1 条は仮定から出る。 -/
theorem birat_frobeniusCompact_cond1
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (A : BiratCat P G) :
    ∀ x y : End A, x ∈ OTimes (biratPre P G) A → y ∈ OTimes (biratPre P G) A →
      x * y = y * x :=
  fun x y hx hy => birat_otimes_comm hfn A x y hx hy

def birat_otimes_comm.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — Frobenius-compact の第 1 条は 𝒞^birat で自動",
    sectionId := "frdi-def-1-2-iv" }


/-! ## ★★★★★★第 3 条は **`c = d` さえ出れば終わる** —— 核の torsion 性は要らない

★★以前の台帳は「`σ` が恒等 ⟹ `Div^gp` の核に入る ⟹ **核が torsion であること**が要る」
と書いていたが、これは**遠回り**であった。

★★★仮定は
```
∀ u ∈ 𝒪^×(A), ∃ k, (endConj θ u)^(d*k) = u^(c*k)
```
であり、結論は `∀ u ∈ 𝒪^×(A), ∃ k, (endConj θ u)^k = u^k` である。
★**`c = d` なら仮定の式そのものが結論の式である**(`k' := d*k` と取ればよい)。
★★したがって残る問いは **`c = d` を出すこと**だけであり、
`𝒪^×(A)` が torsion かどうかは**まったく関係が無い**。
-/

section Cond3

variable {Dd : Type u} [Category.{v} Dd] {Cc : Type u2} [Category.{v2} Cc]
  {Φ₀ : MonoidOn.{v, u, w} Dd} (Q : PreFrobenioid Cc Φ₀)

/-- ★★★★★★**第 3 条は `c = d` から直ちに出る**。

★仮定の `k` に対し `k' := d * k` を取るだけ。 -/
theorem frobeniusCompact_cond3_of_eq {A : Cc} (θ : A ≅ A) (c d : ℕ+) (hcd : c = d)
    (hyp : ∀ u : End A, u ∈ OTimes Q A → ∃ k : ℕ+,
      ((endConj θ u) ^ (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)
        = (u ^ (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)) :
    ∀ u : End A, u ∈ OTimes Q A → ∃ k : ℕ+,
      ((endConj θ u) ^ ((k : ℕ+) : ℕ) : End A) = (u ^ ((k : ℕ+) : ℕ) : End A) := by
  intro u hu
  obtain ⟨k, hk⟩ := hyp u hu
  refine ⟨d * k, ?_⟩
  show ((endConj θ u) ^ (((d * k : ℕ+) : ℕ)) : End A) = (u ^ (((d * k : ℕ+) : ℕ)) : End A)
  rw [show ((d * k : ℕ+) : ℕ) = ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) from rfl]
  rw [hk, hcd]

def frobeniusCompact_cond3_of_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — Frobenius-compact の第 3 条は c = d から出る",
    sectionId := "frdi-def-1-2-iv" }

end Cond3

end ABC3.Found.FrdI
