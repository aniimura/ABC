/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottCoord
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★`hinj` の要らない Northcott（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★これは何か —— `hinj` は最後に足す条件である

`§9-877`（`northcott_of_projModel`）は `hinj`（正規化座標が点を分ける）を
仮定として受けていた。★しかし**その仮定は結論を添字集合 `P` に移すためだけ**に使われている
——中身は「**座標の組の集合が有限**」であり、そこに単射性は要らない。

★★★★本ファイルはその**像の側**を取り出す:

    `{p | ht(p) ≤ C}` の**正規化座標の像**は有限

★★これが Northcott の内容そのものであり、`hinj`（`§9-936` の `hpt`）は
「点の族 `P` が重複を持たない」という**添字づけの条件**に過ぎないことが型で見える。

## ★★★機構 —— もともと証明の中にあった

`§9-877` の `finite_of_injOn_boundedAlg` の証明は

1. `himg : (f '' S).Finite`（次数・高さが有界な組は有限）
2. `Set.Finite.of_finite_image himg hinj`

の 2 段であり、★**単射性は 2 段目でしか使われていない**。
本ファイルは 1 段目を独立の定理として取り出し、鎖の上まで運ぶ。

## ★これで `Proposition 1.4, (iv)` は

★★★`hinj`／`hpt` なしで「**高さで有界な点の座標は有限個**」が言える。
点の族としての有限性が要るときだけ `hinj` を足せばよい。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★像の有限性（単射性なし） -/

/-- ★★★★**次数有界・高さ有界の組へ写る集合の像は有限**——単射性は要らない。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-877` の `finite_of_injOn_boundedAlg` の証明の**第 1 段**である。 -/
theorem finite_image_boundedAlg {α : Type*} {ι : Type*} [Finite ι] {d : ℕ} {B : ℝ}
    (f : α → (ι → ℂ)) (S : Set α)
    (hf : ∀ x ∈ S, ∀ i, f x i ∈ boundedAlg d B) : (f '' S).Finite := by
  refine Set.Finite.subset (finite_pi_boundedAlg (ι := ι) d B) ?_
  rintro y ⟨x, hx, rfl⟩
  exact fun i => hf x hx i

/-- ★★★★★**射影モデルでの Northcott（像の側）**——`hinj` は要らない。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-877` から単射性の仮定だけを落とし、結論を**像**にしたものである。 -/
theorem northcott_of_projModel_image {P : Type*} {ι : Type*} [Finite ι] (d : ℕ)
    (ht : P → ℝ)
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (crd : ∀ p, ι → (fld p)) (idx : ι) (const : ℝ)
    (hcmp : ∀ p, haveI := hnf p; Height.mulHeight (crd p) ≤ Real.exp (ht p + const))
    (C : ℝ) :
    ((fun (p : P) (i : ι) => ((crd p i / crd p idx : fld p) : ℂ)) '' {p | ht p ≤ C}).Finite := by
  refine finite_image_boundedAlg (d := d) (B := Real.exp (C + const))
    (fun (p : P) (i : ι) => ((crd p i / crd p idx : fld p) : ℂ)) _ ?_
  intro p hp i
  haveI := hnf p
  refine coord_mem_boundedAlg d (Real.exp (C + const)) (fld p) (hdeg p) (crd p) ?_ i idx
  refine le_trans (hcmp p) (Real.exp_le_exp.2 ?_)
  have h : ht p ≤ C := hp
  linarith

/-! ## ★★★★★★★★★★★★★★★★★★高さが素朴高さのときの像 -/

/-- ★★★★★★★★★★★★★★★★★★**高さが座標の素朴高さで書けるなら、
座標の像は有限**——`hinj` なし。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★`§9-934`（`northcott_of_log_mulHeight`）の像版である。
`H(x) = exp([F:ℚ]·ht)` と `[F:ℚ] ≤ d`、`ht ≥ 0` から `H(x) ≤ exp(d·ht)` を作る段は同じ。 -/
theorem northcott_of_log_mulHeight_image (d : ℕ) {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    {ι : Type} [Finite ι] (x : ∀ p, ι → NumberField.RingOfIntegers (fld p))
    (ht : P → ℝ)
    (hht : ∀ p, haveI := hnf p; ht p
      = Real.log (Height.mulHeight (fun k => ((x p k : fld p))))
        / (Module.finrank ℚ (fld p) : ℝ))
    (idx : ι) (C : ℝ) :
    ((fun (p : P) (k : ι) => ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ))
      '' {p | ht p ≤ C}).Finite := by
  have hcmp : ∀ p, haveI := hnf p;
      Height.mulHeight (fun k => ((x p k : fld p)))
        ≤ Real.exp ((d : ℝ) * ht p + 0) := by
    intro p
    haveI := hnf p
    have hpos : (0:ℝ) < Height.mulHeight (fun k => ((x p k : fld p))) :=
      lt_of_lt_of_le zero_lt_one (Height.one_le_mulHeight _)
    have hfr : (0:ℝ) < (Module.finrank ℚ (fld p) : ℝ) := by exact_mod_cast Module.finrank_pos
    have hlognn : 0 ≤ Real.log (Height.mulHeight (fun k => ((x p k : fld p)))) :=
      Real.log_nonneg (Height.one_le_mulHeight _)
    have htnn : 0 ≤ ht p := by rw [hht p]; exact div_nonneg hlognn hfr.le
    have hlog : (Module.finrank ℚ (fld p) : ℝ) * ht p
        = Real.log (Height.mulHeight (fun k => ((x p k : fld p)))) := by
      rw [hht p, mul_div_cancel₀ _ hfr.ne']
    have hle : (Module.finrank ℚ (fld p) : ℝ) * ht p ≤ (d : ℝ) * ht p :=
      mul_le_mul_of_nonneg_right (by exact_mod_cast hdeg p) htnn
    rw [add_zero]
    calc Height.mulHeight (fun k => ((x p k : fld p)))
        = Real.exp (Real.log (Height.mulHeight (fun k => ((x p k : fld p))))) :=
          (Real.exp_log hpos).symm
      _ = Real.exp ((Module.finrank ℚ (fld p) : ℝ) * ht p) := by rw [hlog]
      _ ≤ Real.exp ((d : ℝ) * ht p) := Real.exp_le_exp.2 hle
  refine Set.Finite.subset
    (northcott_of_projModel_image d (fun p => (d : ℝ) * ht p) fld hnf hdeg
      (fun p k => ((x p k : fld p))) idx 0 hcmp ((d : ℝ) * C))
    (Set.image_mono ?_)
  intro p hp
  have hp' : ht p ≤ C := hp
  exact mul_le_mul_of_nonneg_left hp' (Nat.cast_nonneg d)

/-! ## ★出典の紐付け(`.src`) -/

def finite_image_boundedAlg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(次数有界・高さ有界の組へ写る集合の像は有限)",
    sectionId := "genell-prop-1-4" }

def northcott_of_projModel_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(射影モデルでの Northcott——像の側、hinj なし)",
    sectionId := "genell-prop-1-4" }

def northcott_of_log_mulHeight_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(高さが素朴高さなら座標の像は有限——hinj なし)",
    sectionId := "genell-prop-1-4" }

def northcott_of_log_mulHeight_image.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "finite_pi_boundedAlg(次数・高さが有界な組は有限、§9-877)"
      (.inProject "ABC3" "ABC3.Found.GenEll.finite_pi_boundedAlg") 2,
    .citation "[ABC3]" "coord_mem_boundedAlg(射影高さから正規化座標へ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.coord_mem_boundedAlg") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-877 の hinj(正規化座標が点を分ける)は" ++
       "**結論を添字集合 P に移すためだけ**に使われていた" ++
       "——証明は (1) 像が有限 (2) 単射性で戻す、の 2 段で、" ++
       "単射性は 2 段目でしか要らない。" ++
       "★本ファイルは 1 段目を独立させ、鎖の上まで運んだ") 5,
    .implicitStep
      ("★★これで Proposition 1.4, (iv) は hinj／hpt なしで" ++
       "「高さで有界な点の**座標**は有限個」が言える。" ++
       "点の族としての有限性が要るときだけ hinj を足せばよい" ++
       "——hinj は添字づけの条件であって Northcott の内容ではない") 4 ]

end ABC3.Found.GenEll
