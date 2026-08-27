import ABC3.Found.GenEll.NorthcottCoord
import ABC3.Found.GenEll.NorthcottClassical

/-!
# `northcott_of_projModel` の非空虚検査 —— **`ℙ¹`（`j`-線）で実際に効く**

`Found/GenEll/NorthcottCoord.lean` の `northcott_of_projModel` は
「射影埋め込みのデータを与えられたものとして」`Proposition 1.4, (iv)`（Northcott 性）を出す。

★**その仮定の束が空虚でないことは、散文のままでは検査されていない。**
本ファイルは `ℙ¹` の場合に実際に仮定を満たして結論を得る。

**これは原典の主張ではない**（我々のモデルについての事実）ので `.src` を持たない。

## ★★★なぜ `ℙ¹` か —— 消費者の場合だから

`ResearchPaper/mathlib-gap.json` の `ample-and-projective-embedding` に書いたとおり、
`Proposition 1.4, (iv)` の消費者は `Proposition 3.4` であり、そこでの適用先は
`M_ell` である。★その粗モジュライは `j`-線、すなわち `ℙ¹` であり、
**`ℤ` 上の射影モデルを明示的に持つ**。

★★だから「原文の `L_ℚ` は ample」を「`L` は `ℤ` 上の射影埋め込みを与える」と
読み替える逸脱は、**消費者に影響しない**——本ファイルはその読み替えが
実際に機能することを機械にかける。

## ★★★★何を確かめたか

| `northcott_of_projModel` の仮定 | `ℙ¹` での埋め方 |
|---|---|
| 各点の定義体（次数 `≤ d`） | `ℚ⟮α⟯` |
| 同次座標 `crd` と割る成分 `idx` | `![α, 1]` と `idx = 1` |
| 正規化座標が**単射** | `0` 成分を見れば `α` が戻る |
| `H(crd p) ≤ exp(ht p + const)` | ★`mulHeight₁_eq_mulHeight`（`mulHeight₁ x = mulHeight ![x,1]`）で**等号** |

★★★結論は古典的 Northcott（`finite_of_finrank_le_of_mulHeight₁_le` と同じ主張）である
——★**別ルートで出る**ことが、抽象形が空虚でないことの証拠になる。
-/

namespace ABC3.Check.GenEll

open ABC3.Found.GenEll NumberField IntermediateField Height

/-- ★★★★★★**`northcott_of_projModel` は空虚でない** —— `ℙ¹`（`j`-線）で実際に効く。

★仮定の束（定義体・同次座標・単射性・高さの比較）を `ℙ¹` の場合に実際に満たし、
古典的 Northcott の主張を**抽象形から**導く。

★★これが `Proposition 1.4, (iv)` の逸脱（「ample」を「射影埋め込みを与える」と読み替える）が
消費者（`Proposition 3.4` の `M_ell`）で機能することの検査である。 -/
theorem northcott_projLine_nonvacuous (d : ℕ) (C : ℝ) :
    {α : ℂ | ∃ h : IsIntegral ℚ α,
      Module.finrank ℚ ℚ⟮α⟯ ≤ d ∧
      (haveI := IntermediateField.adjoin.finiteDimensional h
       haveI : NumberField ℚ⟮α⟯ := ⟨⟩
       Real.log (Height.mulHeight₁
         (⟨α, IntermediateField.mem_adjoin_simple_self ℚ α⟩ : ℚ⟮α⟯)) ≤ C)}.Finite := by
  classical
  set P : Type := {α : ℂ // ∃ h : IsIntegral ℚ α, Module.finrank ℚ ℚ⟮α⟯ ≤ d} with hPdef
  have hnf : ∀ p : P, NumberField ℚ⟮p.1⟯ := by
    intro p
    obtain ⟨h, _⟩ := p.2
    haveI := IntermediateField.adjoin.finiteDimensional h
    exact ⟨⟩
  set fld : P → IntermediateField ℚ ℂ := fun p => ℚ⟮p.1⟯ with hfld
  set elt : ∀ p : P, fld p := fun p =>
    ⟨p.1, IntermediateField.mem_adjoin_simple_self ℚ p.1⟩ with helt
  set ht : P → ℝ := fun p => haveI := hnf p; Real.log (Height.mulHeight₁ (elt p)) with hht
  set crd : ∀ p : P, Fin 2 → fld p := fun p => ![elt p, 1] with hcrd
  have hdeg : ∀ p : P, Module.finrank ℚ (fld p) ≤ d := fun p => p.2.choose_spec
  have hcmp : ∀ p : P, haveI := hnf p; Height.mulHeight (crd p) ≤ Real.exp (ht p + 0) := by
    intro p
    haveI := hnf p
    rw [add_zero, hht]
    simp only
    rw [Real.exp_log (lt_of_lt_of_le zero_lt_one (one_le_mulHeight₁ _))]
    exact le_of_eq (Height.mulHeight₁_eq_mulHeight (elt p)).symm
  have hinj : Function.Injective
      (fun (p : P) (i : Fin 2) => ((crd p i / crd p 1 : fld p) : ℂ)) := by
    intro p q hpq
    have h0 := congrFun hpq 0
    simp only [hcrd, Matrix.cons_val_zero, Matrix.cons_val_one, div_one] at h0
    exact Subtype.ext h0
  have hfin := northcott_of_projModel (P := P) (ι := Fin 2) d ht fld hnf hdeg crd 1 0 hcmp hinj C
  have himg : {α : ℂ | ∃ h : IsIntegral ℚ α,
      Module.finrank ℚ ℚ⟮α⟯ ≤ d ∧
      (haveI := IntermediateField.adjoin.finiteDimensional h
       haveI : NumberField ℚ⟮α⟯ := ⟨⟩
       Real.log (Height.mulHeight₁
         (⟨α, IntermediateField.mem_adjoin_simple_self ℚ α⟩ : ℚ⟮α⟯)) ≤ C)}
      ⊆ Subtype.val '' {p : P | ht p ≤ C} := by
    rintro α ⟨h, hd, hC⟩
    exact ⟨⟨α, ⟨h, hd⟩⟩, hC, rfl⟩
  exact Set.Finite.subset (hfin.image _) himg

end ABC3.Check.GenEll
