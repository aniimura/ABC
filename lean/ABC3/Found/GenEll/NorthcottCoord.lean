import ABC3.Found.GenEll.NorthcottTuple
import ABC3.Found.GenEll.HeightExtension

/-!
# [GenEll] Proposition 1.4, (iv) —— **射影座標へ降ろす 2 段**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★何が残っていたか

`NorthcottTuple.lean` は「**次数 `≤ d`・高さ `≤ B` の代数的数の組は有限個**」を取った。
★残るのは「**射影高さが `B` 以下なら、正規化座標の高さも `B` 以下**」という橋である。

★★2 段でできる:

| 段 | 使うもの |
|---|---|
| `H₁(x_i/x_j) ≤ H(x)` | mathlib の `mulHeight₁_div_eq_mulHeight` + `mulHeight_comp_le` |
| `H_L(z) ≤ B ⟹ H_E(z) ≤ B`(`E ⊆ L`) | `HeightExtension.lean` の `mulHeight₁_extension` |

★★★2 段目が要るのは、`NorthcottClassical` が高さを **`ℚ⟮α⟯` の上で**測っているからである
——大きい体 `L` で測った高さを、生成する部分体まで降ろさねばならない。
`H_L = H_E^{[L:E]}` と `H_E ≥ 1` から**直ちに出る**。
-/

namespace ABC3.Found.GenEll

open NumberField IntermediateField

/-! ## ★1 段目 —— 座標の比は射影高さで抑えられる -/

/-- ★★**`H₁(x_i / x_j) ≤ H(x)`**。

★mathlib の `mulHeight₁_div_eq_mulHeight`(比の高さは 2 成分の高さ)と
`mulHeight_comp_le`(部分組の高さは全体以下)を繋いだだけである。 -/
theorem mulHeight₁_div_le_mulHeight {K : Type*} [Field K] [NumberField K]
    {ι : Type*} [Finite ι] (x : ι → K) (i j : ι) :
    Height.mulHeight₁ (x i / x j) ≤ Height.mulHeight x := by
  rw [Height.mulHeight₁_div_eq_mulHeight]
  have h : ![x i, x j] = x ∘ ![i, j] := by
    funext k
    fin_cases k <;> rfl
  rw [h]
  exact Height.mulHeight_comp_le _ x

/-! ## ★★2 段目 —— 大きい体で測った高さを部分体へ降ろす -/

/-- ★★★**`H_L(z) ≤ B` なら `H_K(z) ≤ B`**(`K ⊆ L`、`z ∈ K`)。

★`H_L = H_K^{[L:K]}` と `H_K ≥ 1`・`[L:K] ≥ 1` から出る。 -/
theorem mulHeight₁_le_of_extension (K L : Type*) [Field K] [NumberField K] [Field L]
    [NumberField L] [Algebra K L] (z : K) (B : ℝ)
    (hB : Height.mulHeight₁ (algebraMap K L z) ≤ B) :
    Height.mulHeight₁ z ≤ B := by
  rcases eq_or_ne z 0 with rfl | hz
  · rw [Height.mulHeight₁_zero]
    simpa using hB
  · have hzu : IsUnit z := isUnit_iff_ne_zero.2 hz
    have hext := mulHeight₁_extension K L hzu.unit
    have hcoe : ((unitsExt K L hzu.unit : Lˣ) : L) = algebraMap K L z := by
      simp [unitsExt]
    rw [hcoe] at hext
    have h1 : (1:ℝ) ≤ Height.mulHeight₁ z := Height.one_le_mulHeight₁ z
    have hpos : 0 < Module.finrank K L := by
      haveI : FiniteDimensional K L := Module.Finite.of_restrictScalars_finite ℚ K L
      exact Module.finrank_pos
    have hle : Height.mulHeight₁ z ≤ Height.mulHeight₁ z ^ (Module.finrank K L) := by
      calc Height.mulHeight₁ z = Height.mulHeight₁ z ^ 1 := (pow_one _).symm
        _ ≤ Height.mulHeight₁ z ^ (Module.finrank K L) := by
            exact pow_le_pow_right₀ h1 hpos
    have hu : (hzu.unit : K) = z := rfl
    rw [hu] at hext
    linarith [hext ▸ hB, hle]

/-! ## ★★★★2 段を繋ぐ —— 正規化座標は `boundedAlg` に入る -/

/-- ★★★★★**射影高さが `B` 以下なら、正規化座標は `boundedAlg d B` に入る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★体は `ℂ` の中の中間体として取る——`X(ℚ̄)` の点は `ℂ` への埋め込みつきで来るからである。
★★これで `NorthcottTuple.lean` の `finite_of_injOn_boundedAlg` に渡せる形になる。 -/
theorem coord_mem_boundedAlg {ι : Type*} [Finite ι] (d : ℕ) (B : ℝ)
    (L : IntermediateField ℚ ℂ) [NumberField L]
    (hd : Module.finrank ℚ L ≤ d)
    (x : ι → L) (hB : Height.mulHeight x ≤ B) (i j : ι) :
    ((x i / x j : L) : ℂ) ∈ boundedAlg d B := by
  have hint : IsIntegral ℚ ((x i / x j : L) : ℂ) := by
    have h0 : IsIntegral ℚ (x i / x j) := Algebra.IsIntegral.isIntegral _
    exact h0.map (IsScalarTower.toAlgHom ℚ L ℂ)
  have hmemL : ((x i / x j : L) : ℂ) ∈ L := (x i / x j).2
  have hle : ℚ⟮((x i / x j : L) : ℂ)⟯ ≤ L :=
    (IntermediateField.adjoin_simple_le_iff).2 hmemL
  haveI hfd := IntermediateField.adjoin.finiteDimensional hint
  haveI : NumberField ℚ⟮((x i / x j : L) : ℂ)⟯ := ⟨⟩
  -- ★★`letI` でなければならない——`haveI` だと包む定義が消えて `rfl` が通らない
  letI : Algebra ℚ⟮((x i / x j : L) : ℂ)⟯ L :=
    (IntermediateField.inclusion hle).toRingHom.toAlgebra
  haveI : IsScalarTower ℚ ℚ⟮((x i / x j : L) : ℂ)⟯ L :=
    IsScalarTower.of_algebraMap_eq fun _ => rfl
  haveI : FiniteDimensional ℚ⟮((x i / x j : L) : ℂ)⟯ L :=
    Module.Finite.of_restrictScalars_finite ℚ _ L
  refine ⟨hint, ?_, ?_⟩
  · -- ★次数: mathlib の `finrank_le_of_le_right` そのもの
    exact le_trans (IntermediateField.finrank_le_of_le_right hle) hd
  · -- ★高さ: `L` で測ってから生成する部分体へ降ろす
    have hB1 : Height.mulHeight₁ (x i / x j) ≤ B :=
      le_trans (mulHeight₁_div_le_mulHeight x i j) hB
    refine mulHeight₁_le_of_extension ℚ⟮((x i / x j : L) : ℂ)⟯ L
      (⟨((x i / x j : L) : ℂ),
        IntermediateField.mem_adjoin_simple_self ℚ ((x i / x j : L) : ℂ)⟩) B ?_
    have hmap : (algebraMap ℚ⟮((x i / x j : L) : ℂ)⟯ L)
        (⟨((x i / x j : L) : ℂ),
          IntermediateField.mem_adjoin_simple_self ℚ ((x i / x j : L) : ℂ)⟩)
        = x i / x j := Subtype.ext rfl
    rw [hmap]
    exact hB1

/-! ## ★★★★★★閉じた形 -/

/-- ★★★★★★**`ℙⁿ(L)` の点で射影高さが `B` 以下のものは有限個**(`[L:ℚ] ≤ d`)。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★正規化座標(`j` 番目で割ったもの)を `ℂ` の中で見ているので、
**`L` が動いても同じ `boundedAlg d B` の中に入る**——それが次数を `d` で押さえた理由である。 -/
theorem finite_normalizedCoord {ι : Type*} [Finite ι] (d : ℕ) (B : ℝ)
    (L : IntermediateField ℚ ℂ) [NumberField L] (hd : Module.finrank ℚ L ≤ d) (jj : ι) :
    {y : ι → ℂ | ∃ x : ι → L, Height.mulHeight x ≤ B ∧
        ∀ i, y i = ((x i / x jj : L) : ℂ)}.Finite := by
  refine Set.Finite.subset (finite_pi_boundedAlg (ι := ι) d B) ?_
  rintro y ⟨x, hx, hy⟩ i
  rw [hy i]
  exact coord_mem_boundedAlg d B L hd x hx i jj

/-! ## ★出典の紐付け(`.src`) -/

def mulHeight₁_div_le_mulHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(座標の比の高さは射影高さ以下)",
    sectionId := "genell-prop-1-4" }

def mulHeight₁_le_of_extension.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(大きい体で測った高さを部分体へ降ろす)",
    sectionId := "genell-prop-1-4" }

def coord_mem_boundedAlg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(射影高さから正規化座標の Northcott 性へ)",
    sectionId := "genell-prop-1-4" }

def finite_normalizedCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ℙⁿ の次数有界 Northcott)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
