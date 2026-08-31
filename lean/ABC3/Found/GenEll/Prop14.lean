/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottCoordsImage
import ABC3.Found.GenEll.Prop16
import ABC3.Found.GenEll.HeightMetric
import ABC3.Found.GenEll.HeightAdditive
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★[GenEll] Proposition 1.4 —— 4 条を 1 本に（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★4 条の内訳

| 条 | 原文 | 実装 |
|---|---|---|
| (i) 加法性 | `ht_{L̄⊗M̄} = ht_L̄ + ht_M̄` | `htArith_tensor_unconditional`（`HeightAdditive.lean`） |
| (ii) 下に有界 | `ht_L̄` は下に一様に有界 | `prop_1_4_ii`（`Prop16.lean`） |
| (iii) BD-類 | 計量を替えても BD-類は同じ | `htArith_sub_abs_le`（`HeightMetric.lean`） |
| (iv) Northcott | `L_ℚ` が豊富なら有界な点は有限 | ★`northcott_of_isAmple_coords_image`（`§9-961`、本日） |

★★★★★★**(iv) が本日（2026-08-29）閉じたので、4 条が揃った**。

## ★★★★★逸脱の記録（`CLAUDE.md` の「逸脱」）

### 1. 量化する対象

原文は `∀ D : HeightTheoryData` の形で読めるが、**それは偽である**
——`HeightTheoryData` は公理を 1 つも持たないデータであり、
`Check/GenEll/HeightAxiomGap.lean` の `prop_1_4_statement_false` で機械検証済み。
★公理を足すのは `check.mjs` の **B5** そのもの（`Proposition 1.4` の (i)-(iv) は
足すべき公理そのもの）なので、**構成に置き換えた**——`ArcModel` ＋ `ArithCartier`。

### 2. (iii) の形

原文は『生成ファイバーが同じなら』。因子表示では『因子が同じ（計量だけ違う）』である。
★**垂直因子の差の分は含めていない**（`Remark 1.4.1` の側で扱う）。

### 3. (iv) の形 —— ★本日変わったところ

* 以前は**射影モデルをデータとして受けていた**（`northcott_of_projModel`）。
* ★★★本 statement は**豊富性から射影モデルを作る**（`§9-914`〜`§9-920`）ので、
  受けるのは原文どおり `IsAmple` である。
* ★★点の側は「点 `xF p` とその**同次座標** `x p`」だけ
  （同次座標は `§9-941` で**必ず取れる**）。
* ★★★★結論は**正規化座標の像の有限性**である
  ——`hinj`（座標が点を分ける）は `§9-960` で測ったとおり
  結論を添字集合へ移すためだけの条件であり、Northcott の内容ではない。
  ★点の族としての有限性が要るときは `§9-952` の形に `hinj` を足せばよい。

### 4. `Green` 関数の連続性

原文の `|−|_L` が連続なのはエルミート計量の定義から。
★本実装の `GreenFn` は任意の実数値関数なので、連続性は利用者が与える
（`Prop16Proper.lean` と同じ扱い）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★4 条 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] Proposition 1.4**
（Basic Properties of Heights）—— (i)〜(iv) の 4 条。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

* **(i)** 加法性 `ht_{D̄⊗Ē} = ht_D̄ + ht_Ē`
* **(ii)** 下に一様に有界（定数は `F` にも点にも依らない）
* **(iii)** 因子が同じで計量が連続なら高さの差は一様に有界（BD-類）
* **(iv)** ★**豊富なら Northcott**——点の側は点とその同次座標だけ、単射性なし

★★★★★逸脱はファイル冒頭に記録した。 -/
theorem prop_1_4 {V : Type}
    [NormedAddCommGroup V] [NormedSpace ℂ V] [FiniteDimensional ℂ V]
    (M : ArcModel X V) [Nonempty (complexPoints X)]
    (Dv : ArithCartier X) (hg : @Continuous _ _ M.topology _ Dv.green) :
    -- ★(i) 加法性
    (∀ (Ev : ArithCartier X) (F : Type) [Field F] [NumberField F]
        (xF : specRingOfIntegers F ⟶ X),
        pullbackIdeal F Dv.divisor xF ≠ 0 → pullbackIdeal F Ev.divisor xF ≠ 0 →
        htArith F (Dv.tensor Ev) xF = htArith F Dv xF + htArith F Ev xF)
    -- ★(ii) 下に一様に有界
  ∧ (∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
        (xF : specRingOfIntegers F ⟶ X), -C ≤ htArith F Dv xF)
    -- ★(iii) 計量だけ違えば高さの差は一様に有界
  ∧ (∀ Ev : ArithCartier X, Ev.divisor = Dv.divisor →
        @Continuous _ _ M.topology _ Ev.green →
        ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
          (xF : specRingOfIntegers F ⟶ X),
          |htArith F Dv xF - htArith F Ev xF| ≤ C)
    -- ★★★(iv) 豊富なら Northcott（点の側は同次座標だけ、単射性なし）
  ∧ (∀ (Lb : X.PresheafOfModules) (hLb : IsLocallyTrivial X Lb) (hample : IsAmple Lb)
        (x₀ : (X : Type)) (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type)) (d : ℕ),
        ∃ (L : ℕ) (_ : 0 < L) (N : ℕ)
          (s : Fin (N + 1) → ((presheafTensorPow Lb L).obj (op ⊤) : Type))
          (hcov : (⨆ k, nonVanishing (presheafTensorPow Lb L) (s k)) = ⊤),
          (∀ i, IsAffineOpen (nonVanishing (presheafTensorPow Lb L) (s i))) ∧
          ∀ {P : Type}
            (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
            (_ : ∀ p, Module.finrank ℚ (fld p) ≤ d)
            (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
            (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
            (_ : ∀ p, x p ≠ 0) (_ : ∀ p, x p 0 ≠ 0)
            (idx0 : P → Fin (N + 1))
            (hx : ∀ p, haveI := hnf p; Set.range
              (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 (fld p)) (fld p))) ≫ xF p
                ≫ globalToProj (presheafTensorPow Lb L)
                  (isLocallyTrivial_presheafTensorPow hLb L) φ s hcov).base
              ⊆ Set.range (chartA N ℤ (idx0 p)).base)
            (_ : ∀ p, x p (idx0 p) ≠ 0)
            (_ : ∀ p, haveI := hnf p; ∀ k, ((x p k : fld p))
              = projPointCoord N ℤ (fld p)
                  (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 (fld p)) (fld p))) ≫ xF p
                    ≫ globalToProj (presheafTensorPow Lb L)
                      (isLocallyTrivial_presheafTensorPow hLb L) φ s hcov)
                  (idx0 p) (hx p) k * ((x p (idx0 p) : fld p)))
            (idx : Fin (N + 1)) (C : ℝ),
            ((fun (p : P) (k : Fin (N + 1)) =>
              ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)) ''
              {p | haveI := hnf p;
                htArith (fld p) ((hyperplaneArith N).comap
                  (globalToProj (presheafTensorPow Lb L)
                    (isLocallyTrivial_presheafTensorPow hLb L) φ s hcov)) (xF p) ≤ C}).Finite) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · intro Ev F _ _ xF hD hE
    exact htArith_tensor_unconditional F Dv Ev xF hD hE
  · exact prop_1_4_ii M Dv hg
  · intro Ev hdiv hcont
    exact htArith_sub_abs_le M Dv Ev hdiv.symm hg hcont
  · intro Lb hLb hample x₀ φ d
    exact northcott_of_isAmple_coords_image Lb hLb hample x₀ φ d

/-! ## ★出典の紐付け(`.src`) -/

def prop_1_4.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6, item := "Proposition 1.4",
    sectionId := "genell-prop-1-4" }

def prop_1_4.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_tensor_unconditional((i) 加法性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_tensor_unconditional") 3,
    .citation "[ABC3]" "prop_1_4_ii((ii) 下に一様に有界)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prop_1_4_ii") 3,
    .citation "[ABC3]" "htArith_sub_abs_le((iii) 計量だけ違えば差は一様に有界)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_sub_abs_le") 3,
    .citation "[ABC3]" "northcott_of_isAmple_coords_image((iv) 豊富なら Northcott、§9-961)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_isAmple_coords_image") 6,
    .implicitStep
      ("★★★★★逸脱 1(量化する対象): 原文を ∀ D : HeightTheoryData と読むと**偽である**" ++
       "——Check/GenEll/HeightAxiomGap.lean の prop_1_4_statement_false で機械検証済み。" ++
       "公理を足すのは check.mjs の B5 そのものなので、構成(ArcModel ＋ ArithCartier)に" ++
       "置き換えた") 6,
    .implicitStep
      ("★★逸脱 2((iii) の形): 原文は『生成ファイバーが同じなら』。" ++
       "因子表示では『因子が同じ(計量だけ違う)』であり、" ++
       "**垂直因子の差の分は含めていない**(Remark 1.4.1 の側で扱う)") 4,
    .implicitStep
      ("★★★★逸脱 3((iv) の形、2026-08-29 に変わった): " ++
       "以前は射影モデルをデータとして受けていたが、本 statement は" ++
       "**豊富性から射影モデルを作る**(§9-914〜920)ので、受けるのは原文どおり IsAmple である。" ++
       "★点の側は点とその同次座標だけ(同次座標は §9-941 で必ず取れる)。" ++
       "★★結論は正規化座標の**像**の有限性である" ++
       "——hinj は §9-960 で測ったとおり結論を添字集合へ移すためだけの条件であり、" ++
       "Northcott の内容ではない(点の族としての有限性が要るときは §9-952 の形に足す)") 8,
    .implicitStep
      ("★逸脱 4(Green 関数の連続性): 原文の |−|_L が連続なのはエルミート計量の定義から。" ++
       "本実装の GreenFn は任意の実数値関数なので連続性は利用者が与える" ++
       "(Prop16Proper.lean と同じ扱い)") 3 ]

end ABC3.Found.GenEll
