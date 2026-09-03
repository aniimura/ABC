/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ComapOnAffineOpen
import ABC3.Found.GenEll.DivIdealChart
import ABC3.Meta.Claim

/-!
# DivisorOfSectionComap —— `[GenEll] Definition 1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
variable {X : Scheme.{0}}

/-! ## ★★★★★★`divIdeal` を任意の小さいアフィン開で読む -/

/-- ★★★★★★**`divIdeal M s₀ W = span {(s₀/t) の制限}`**
（`W` が自明化つき `V` の中の、`X_t` に含まれるアフィン開なら）。

★`§9-923` の一般化（そこでは `W = X_t ⊓ V` に固定していた）。
★★機構は同じ——`trivValue(s₀)` と `sectionRatio(s₀,t)` は単元倍しか違わない。 -/
theorem divIdeal_eq_span_globalRatio_of_le (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (V : X.Opens) (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s₀ t : (M.obj (op ⊤) : Type))
    {W : X.Opens} (hW : IsAffineOpen W) (hWV : W ≤ V) (hWt : W ≤ nonVanishing M t) :
    divIdeal M s₀ ⟨W, hW⟩
      = Ideal.span {(X.presheaf.map (homOfLE hWt).op).hom (globalRatio M hM s₀ t)} := by
  have hle : W ≤ nonVanishing M t ⊓ V := le_inf hWt hWV
  have h1 : (X.presheaf.map (homOfLE hWt).op).hom (globalRatio M hM s₀ t)
      = (X.presheaf.map (homOfLE hle).op).hom (sectionRatio M V e s₀ t) := by
    rw [← globalRatio_res M hM s₀ t ⟨V, e⟩, res_trans]
  rw [divIdeal_eq M s₀ ⟨W, hW⟩ (trivialOfLe M hWV e), trivValue_restrict M hWV e s₀,
    h1, sectionRatio, map_mul, res_trans]
  refine (Ideal.span_singleton_mul_right_unit ?_ _).symm
  exact ((isUnit_trivValue_res M V e t).unit⁻¹.isUnit).map
    (X.presheaf.map (homOfLE hle).op).hom

/-! ## ★出典の紐付け(`.src`) -/

def divIdeal_eq_span_globalRatio_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(divIdeal は任意の小さいアフィン開で (s₀/t) の制限が生成する)",
    sectionId := "genell-def-1-2" }

end ABC3.Found.GenEll
