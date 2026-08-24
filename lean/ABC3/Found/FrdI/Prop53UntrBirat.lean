/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop53Birat
import ABC3.Found.FrdI.Thm51Span
import ABC3.Found.FrdI.Prop44Core

/-!
# [FrdI] `Φ^birat` は `𝒞` と `𝒞^un-tr` で一致する

原文 (FrdI p.103):
> tively, Φpf) and the rational function monoid Φbirat (respectively, Q · Φbirat =

★★`Proposition 5.3` の第 2 文は `𝒞^un-tr` と `(𝒞^un-tr)^pf` の有理関数の単系を
**`Φ^birat`(＝ `𝒞` のもの)** と書く。★我々の ` phiBiratOn` は
**Frobenioid ごとに定義されている**ので、`𝒞^un-tr` で計算したものと
`𝒞` で計算したものが一致することを別に言う必要がある。

## ★短い道 —— birationalization を比較しない

`Theorem 5.1, (i)` のために作った**辞書**

* `exists_preStepPair_of_mem_phiBiratAt` —— `Φ^birat` の元 ⇒ base-equivalent な pre-step の対
* `mem_phiBiratAt_of_preStepPair` —— その逆

は `Φ^birat(A)` を **`𝒞` の言葉だけ**(pre-step の対と `sliceDivGpOf`)で書いている。
★★`𝒞^un-tr` の射は `𝒞` の射の商であり、`Base`・`Div`・`deg_Fr` は **`rfl` で一致する**
(`unTrPre_Base` / `unTrPre_Div` / `unTrPre_degFr`)。
★したがって pre-step の対は両向きに移り、`sliceDivGpOf` の値も変わらない。

★★★**birationalization どうしの比較関手(`Rmk451Birat.lean` の 4 段)は要らない。**
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section UnTrBirat

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★`𝒞^un-tr` の射の pre-step 性は代表元のそれと一致する(`rfl`)。 -/
theorem unTr_isPreStep_iff_rep (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    IsPreStep (unTrPre P Fc) (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      ↔ IsPreStep P α :=
  Iff.rfl

/-- ★`sliceDivGpOf` も代表元のそれと一致する(`rfl`)。 -/
theorem sliceDivGpOf_unTr (Fc : FrobenioidCore P) {A B W : Istr P}
    (a : W.obj ⟶ A.obj) (ha : IsIso (P.Base a)) (φ : W.obj ⟶ B.obj) :
    sliceDivGpOf (P := unTrPre P Fc)
        (toHomUnTr P a : (show UnTr P from W) ⟶ (show UnTr P from A)) ha
        (toHomUnTr P φ : (show UnTr P from W) ⟶ (show UnTr P from B))
      = sliceDivGpOf (P := P) a ha φ :=
  rfl

set_option maxHeartbeats 800000 in
/-- ★★★★★★**`Φ^birat` は `𝒞` と `𝒞^un-tr` で一致する**(点ごと)。

★両側とも `Theorem 5.1, (i)` の辞書で
**base-equivalent な pre-step の対**に翻訳し、
`𝒞^un-tr` の射が `𝒞` の射の商であることで往復する。 -/
theorem phiBiratAt_unTr_eq (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (A : Istr P) :
    phiBiratAt (unTrPre P Fc) (unTr_frobenioid P Fc G) (show UnTr P from A)
      = phiBiratAt P G A.obj := by
  ext x
  constructor
  · -- `𝒞^un-tr` ⟹ `𝒞` : 代表元を取る
    intro hx
    obtain ⟨W, δ₁, δ₂, hs₁, hc₁, hs₂, hbase, hval⟩ :=
      exists_preStepPair_of_mem_phiBiratAt (unTr_frobenioid P Fc G) hx
    obtain ⟨δ₁', hδ₁'⟩ := Quotient.exists_rep δ₁
    obtain ⟨δ₂', hδ₂'⟩ := Quotient.exists_rep δ₂
    subst hδ₁'
    subst hδ₂'
    have hs₁' : IsPreStep P δ₁' := hs₁
    have hs₂' : IsPreStep P δ₂' := hs₂
    have hbase' : P.Base δ₁' = P.Base δ₂' := hbase
    have hc₁' : IsCoAngular P δ₁' :=
      isCoAngular_of_isotropic P G (hiso W.obj) δ₁' hs₁'
    have hc₂' : IsCoAngular P δ₂' :=
      isCoAngular_of_isotropic P G (hiso W.obj) δ₂' hs₂'
    have hmem := mem_phiBiratAt_of_preStepPair G δ₁' δ₂' hc₁' hs₁' hc₂' hs₂' hbase'
    rwa [show sliceDivGpOf (P := P) δ₁' hs₁'.2 δ₂' = x from hval] at hmem
  · -- `𝒞` ⟹ `𝒞^un-tr` : 射を商へ送る
    intro hx
    obtain ⟨W, δ₁, δ₂, hs₁, hc₁, hs₂, hbase, hval⟩ :=
      exists_preStepPair_of_mem_phiBiratAt G hx
    have hW : IsIsotropic P W := hiso W
    have hs₁' : IsPreStep (unTrPre P Fc)
        (toHomUnTr P δ₁ : (show UnTr P from (⟨W, hW⟩ : Istr P)) ⟶ (show UnTr P from A)) := hs₁
    have hs₂' : IsPreStep (unTrPre P Fc)
        (toHomUnTr P δ₂ : (show UnTr P from (⟨W, hW⟩ : Istr P)) ⟶ (show UnTr P from A)) := hs₂
    have hbase' : (unTrPre P Fc).Base
        (toHomUnTr P δ₁ : (show UnTr P from (⟨W, hW⟩ : Istr P)) ⟶ (show UnTr P from A))
      = (unTrPre P Fc).Base
        (toHomUnTr P δ₂ : (show UnTr P from (⟨W, hW⟩ : Istr P)) ⟶ (show UnTr P from A)) := hbase
    have hc₁' : IsCoAngular (unTrPre P Fc) _ :=
      isCoAngular_of_isotropic (unTrPre P Fc) (unTr_frobenioid P Fc G)
        (unTr_isotropic P Fc (show UnTr P from (⟨W, hW⟩ : Istr P))) _ hs₁'
    have hc₂' : IsCoAngular (unTrPre P Fc) _ :=
      isCoAngular_of_isotropic (unTrPre P Fc) (unTr_frobenioid P Fc G)
        (unTr_isotropic P Fc (show UnTr P from (⟨W, hW⟩ : Istr P))) _ hs₂'
    have hmem := mem_phiBiratAt_of_preStepPair (unTr_frobenioid P Fc G)
      (toHomUnTr P δ₁ : (show UnTr P from (⟨W, hW⟩ : Istr P)) ⟶ (show UnTr P from A))
      (toHomUnTr P δ₂) hc₁' hs₁' hc₂' hs₂' hbase'
    rwa [show sliceDivGpOf (P := unTrPre P Fc)
        (toHomUnTr P δ₁ : (show UnTr P from (⟨W, hW⟩ : Istr P)) ⟶ (show UnTr P from A))
        hs₁'.2 (toHomUnTr P δ₂) = x from hval] at hmem

/-- ★★★★★★★**`Φ^birat` は `𝒞` と `𝒞^un-tr` で一致する**(`𝒟` 上の部分関手として)。

★これが `Proposition 5.3` の第 2 文が `Φ^birat` を **`𝒞` のもの**として書くことの根拠である。 -/
theorem phiBiratOn_unTr_eq (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (d : D) :
    phiBiratOn (unTr_frobenioid P Fc G) d = phiBiratOn G d := by
  ext y
  have hA : IsIsotropic P (biratBaseObj G d) := hiso _
  have e : (P.toElem.obj (biratBaseObj G d)).base ≅ d := biratBaseIso G d
  rw [mem_phiBiratOn_iff (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (A := (show UnTr P from (⟨biratBaseObj G d, hA⟩ : Istr P))) e y,
    mem_phiBiratOn_iff G hiso (A := biratBaseObj G d) e y,
    phiBiratAt_unTr_eq Fc G hiso (⟨biratBaseObj G d, hA⟩ : Istr P)]
  rfl

end UnTrBirat

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.3` の第 2 文が `Φ^birat` を
`𝒞` のものとして書くことの根拠。 -/
def phiBiratOn_unTr_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — Φ^birat は 𝒞 と 𝒞^un-tr で一致する",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
