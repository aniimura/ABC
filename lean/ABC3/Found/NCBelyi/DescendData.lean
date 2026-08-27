import ABC3.Found.NCBelyi.NestedInduction

/-!
# [NCBelyi] Lemma 2.4 —— 1 段下げるデータ(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> applying an automorphism (with rational coeﬃcients!) as in Lemma 2.3 and then

## ★★★★★★★★★★何を取るか

`NestedInduction.lean` の `nested_induction_descend` は
**「1 段下げる操作」があれば全体が回る**ことを言う。その操作を実際に作る。

原文 (NCBelyi p.5):
> every α ∈f0(S0) [since f0(x),f

| 作るもの | 中身 |
|---|---|
| `x₀` | 最大の reduced degree を達成する点(`Finset.exists_mem_eq_sup`) |
| `f₀` | `x₀` の**最小多項式**(★原文が最小多項式に取る理由がここで効く) |
| `g` | `x ↦ f₀(x)` |
| `E` | `f₀(S₀)`——`S₀` は `f₀′` の根の集まり |

★★★**測度が下がる根拠 3 つ**(`measure_lt` の仮定)が、いずれもここで揃う:

- `redDeg (g x) ≤ redDeg x` —— `ℚ(g(x)) ⊆ ℚ(x)`(`redDeg_aeval_le`)
- `∀ y ∈ E, redDeg y < m(S)` —— `f₀′` の次数が `≤ d₀ − 1`(`redDeg_eval_root_derivative_lt`)
- `redDeg (g x₀) < m(S)` —— ★**`f₀(α₀) = 0`**、すなわち `redDeg 0 = 0`

★★★★**3 つめが原文の要である。** `f₀` を `α₀` の最小多項式に取るのは、
まさに `f₀(α₀) = 0` を得るためであり、原文はその理由を書いていない。

## ★★残っているもの

本ファイルは**集合と測度の段**である。`Lemma 2.4` の結論はさらに
**`ℙ¹` 上の有理関数としての組み立て**(不分岐性の帳簿)を要求する。
★出力は多項式ではなく `f(x) ∈ ℚ(x)` なので、`Separation.lean` の `P1C` が要る。
-/

namespace ABC3.Found.NCBelyi

open Polynomial IntermediateField

/-! ## ★最小多項式の次数と reduced degree -/

/-- ★**`deg (minpoly ℚ x) = redDeg x + 1`**。 -/
theorem natDegree_minpoly_eq_redDeg_succ {x : ℂ} (hx : IsIntegral ℚ x) :
    (minpoly ℚ x).natDegree = redDeg x + 1 := by
  have h1 : Module.finrank ℚ ℚ⟮x⟯ = (minpoly ℚ x).natDegree :=
    IntermediateField.adjoin.finrank hx
  have h2 : 1 ≤ Module.finrank ℚ ℚ⟮x⟯ := one_le_finrank_adjoin x hx
  rw [redDeg]
  omega

/-- ★**`redDeg 0 = 0`** —— `0` は有理点だから。 -/
theorem redDeg_zero : redDeg (0 : ℂ) = 0 := by
  have h : ((0 : ℚ) : ℂ) = (0 : ℂ) := by norm_num
  have := redDeg_ratCast (0 : ℚ)
  rwa [h] at this

/-! ## ★★★★★★★★★★1 段下げるデータ -/

/-- ★★★★★★★★★★**[NCBelyi] Lemma 2.4 の帰納段のデータ**。

原文 (NCBelyi p.5):
> every α ∈f0(S0) [since f0(x),f

`m(S) > 0` なら、最大層の点 `x₀` とその最小多項式 `f₀` を取ることで
`measure_lt` の仮定 3 つをすべて満たす `(E, g, x₀)` が作れる。

★★`E ≝ f₀(S₀)`(`S₀` は `f₀′` の `ℂ` での根)であり、原文の `S′ ≝ f₀(S) ∪ f₀(S₀)` は
`S.image g ∪ E` にあたる。 -/
theorem exists_descend_data (S : Finset ℂ) (hint : ∀ x ∈ S, IsIntegral ℚ x)
    (hS : 0 < maxRedDeg S) :
    ∃ (E : Finset ℂ) (x₀ : ℂ),
      x₀ ∈ S ∧ redDeg x₀ = maxRedDeg S
      ∧ (∀ x ∈ S, redDeg (Polynomial.aeval x (minpoly ℚ x₀)) ≤ redDeg x)
      ∧ (∀ y ∈ E, redDeg y < maxRedDeg S)
      ∧ redDeg (Polynomial.aeval x₀ (minpoly ℚ x₀)) < maxRedDeg S
      ∧ E = (((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)).roots.toFinset).image
              (fun w => Polynomial.aeval w (minpoly ℚ x₀)) := by
  classical
  have hne : S.Nonempty := by
    rw [← Finset.card_pos]
    by_contra hc
    have : S = ∅ := Finset.card_eq_zero.1 (by omega)
    rw [this] at hS
    simp [maxRedDeg] at hS
  obtain ⟨x₀, hx₀S, hx₀v⟩ := Finset.exists_mem_eq_sup S hne redDeg
  have hx₀top : redDeg x₀ = maxRedDeg S := hx₀v.symm
  have hx₀int : IsIntegral ℚ x₀ := hint x₀ hx₀S
  have hdeg : (minpoly ℚ x₀).natDegree = redDeg x₀ + 1 :=
    natDegree_minpoly_eq_redDeg_succ hx₀int
  have hdeg0 : 0 < (minpoly ℚ x₀).natDegree := by omega
  refine ⟨(((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)).roots.toFinset).image
      (fun w => Polynomial.aeval w (minpoly ℚ x₀)), x₀, hx₀S, hx₀top, ?_, ?_, ?_, rfl⟩
  · -- ★像は reduced degree を上げない
    intro x hx
    exact redDeg_aeval_le (hint x hx) _
  · -- ★`f₀(S₀)` の点は最大層より下
    intro y hy
    obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 hy
    rw [Multiset.mem_toFinset] at hw
    have hd0 : (derivative (minpoly ℚ x₀)) ≠ 0 := by
      intro hc
      have := (Polynomial.derivative_eq_zero (p := minpoly ℚ x₀)).1 hc
      omega
    have hmapne : ((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)) ≠ 0 :=
      (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
    have hwr : Polynomial.aeval w (derivative (minpoly ℚ x₀)) = 0 := by
      have := (Polynomial.mem_roots hmapne).1 hw
      simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using this
    have := redDeg_eval_root_derivative_lt hdeg0 hwr
    omega
  · -- ★★`f₀(α₀) = 0` —— ここが原文の要である
    rw [minpoly.aeval, redDeg_zero]
    omega

/-! ## ★★★★★★★★★★★不変量つきの帰納 -/

/-- ★★**不変量 `Q` を持ち回る形の入れ子帰納法**。

★`Lemma 2.4` の帰納では「`S` の点がすべて代数的である」を持ち回る必要がある
——`redDeg` の性質はどれも `IsIntegral ℚ x` を要求するからである。 -/
theorem nested_induction_descend' {P Q : Finset ℂ → Prop}
    (base : ∀ S : Finset ℂ, maxRedDeg S = 0 → Q S → P S)
    (descend : ∀ S : Finset ℂ, 0 < maxRedDeg S → Q S →
      ∃ (E : Finset ℂ) (g : ℂ → ℂ) (x₀ : ℂ),
        (∀ x ∈ S, redDeg (g x) ≤ redDeg x)
        ∧ (∀ y ∈ E, redDeg y < maxRedDeg S)
        ∧ x₀ ∈ S ∧ redDeg x₀ = maxRedDeg S ∧ redDeg (g x₀) < maxRedDeg S
        ∧ Q (S.image g ∪ E)
        ∧ (P (S.image g ∪ E) → P S)) :
    ∀ S : Finset ℂ, Q S → P S := by
  classical
  refine nested_induction (P := fun S => Q S → P S) (fun S h hQ => base S h hQ) ?_
  intro S hS ihm ihd hQ
  obtain ⟨E, g, x₀, hno, hE, hx₀S, hx₀top, hx₀drop, hQ', himp⟩ := descend S hS hQ
  refine himp ?_
  rcases measure_lt S E g hno hE hx₀S hx₀top hx₀drop hS with h | ⟨h1, h2⟩
  · exact ihm _ h hQ'
  · exact ihd _ h1 h2 hQ'

/-! ## ★★★★★★★★★★★★有理の場合への帰着(多項式の段) -/

/-- ★★★★★★★★★★★★**[NCBelyi] Lemma 2.4 の集合の段** ——
**代数的数の有限集合は `ℚ[x]` の非定数多項式 1 つで `ℚ` の中へ写せる**。

原文 (NCBelyi p.4):
> Lemma 2.4.

★原文の `Lemma 2.4` は `ℙ¹` 上の**有理関数** `f(x) ∈ ℚ(x)` を出すが、
**帰納の各段で使うのは多項式(最小多項式)である**——
有理関数になるのは正規化(`Lemma 2.3` の自己同型)のところだけである。
★★本定理はその**多項式の段を端から端まで**取ったものである。

★★★中身は `nested_induction_descend'` に `exists_descend_data` を渡すだけ:
- 基底段(`m(S) = 0`)では `S` はすべて有理点なので `h ≝ X` で足りる
- 帰納段では `h ≝ h′ ∘ f₀`(`f₀` は `x₀` の最小多項式)

★★★★**`f₀(α₀) = 0` が測度を下げ、合成が結論を運ぶ**——原文の
『replacing S by S′, β by f₀(β), applying the induction hypothesis,
and composing the resulting morphisms』そのものである。 -/
theorem exists_poly_image_rat (S : Finset ℂ) (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    ∃ h : ℚ[X], 0 < h.natDegree
      ∧ ∀ x ∈ S, ∃ q : ℚ, Polynomial.aeval x h = (q : ℂ) := by
  classical
  refine nested_induction_descend'
    (P := fun S => ∃ h : ℚ[X], 0 < h.natDegree
      ∧ ∀ x ∈ S, ∃ q : ℚ, Polynomial.aeval x h = (q : ℂ))
    (Q := fun S => ∀ x ∈ S, IsIntegral ℚ x) ?_ ?_ S hint
  · -- ★基底段: `m(S) = 0` ならすべて有理点
    intro T hT hQ
    refine ⟨X, by simp, fun x hx => ?_⟩
    obtain ⟨q, hq⟩ := (redDeg_eq_zero_iff (hQ x hx)).1 (maxRedDeg_eq_zero_iff.1 hT x hx)
    exact ⟨q, by simpa using hq⟩
  · -- ★★帰納段
    intro T hT hQ
    obtain ⟨E, x₀, hx₀S, hx₀top, hno, hE, hdrop, hEdef⟩ := exists_descend_data T hQ hT
    refine ⟨E, fun x => Polynomial.aeval x (minpoly ℚ x₀), x₀, hno, hE, hx₀S, hx₀top,
      hdrop, ?_, ?_⟩
    · -- ★`S′` の点も代数的
      intro y hy
      rcases Finset.mem_union.1 hy with hy | hy
      · obtain ⟨x, hx, rfl⟩ := Finset.mem_image.1 hy
        exact isIntegral_aeval (hQ x hx) _
      · rw [hEdef] at hy
        obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 hy
        rw [Multiset.mem_toFinset] at hw
        have hx₀int : IsIntegral ℚ x₀ := hQ x₀ hx₀S
        have hdeg : (minpoly ℚ x₀).natDegree = redDeg x₀ + 1 :=
          natDegree_minpoly_eq_redDeg_succ hx₀int
        have hd0 : (derivative (minpoly ℚ x₀)) ≠ 0 := by
          intro hc
          have := (Polynomial.derivative_eq_zero (p := minpoly ℚ x₀)).1 hc
          omega
        have hmapne : ((derivative (minpoly ℚ x₀)).map (algebraMap ℚ ℂ)) ≠ 0 :=
          (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
        have hwr : Polynomial.aeval w (derivative (minpoly ℚ x₀)) = 0 := by
          have := (Polynomial.mem_roots hmapne).1 hw
          simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using this
        exact isIntegral_aeval (isIntegral_of_aeval_eq_zero hd0 hwr) _
    · -- ★★★合成が結論を運ぶ
      rintro ⟨h', hdeg', hval'⟩
      refine ⟨h'.comp (minpoly ℚ x₀), ?_, fun x hx => ?_⟩
      · rw [natDegree_comp]
        refine Nat.mul_pos hdeg' ?_
        have hdeg : (minpoly ℚ x₀).natDegree = redDeg x₀ + 1 :=
          natDegree_minpoly_eq_redDeg_succ (hQ x₀ hx₀S)
        omega
      · rw [Polynomial.aeval_comp]
        exact hval' _ (Finset.mem_union_left _ (Finset.mem_image.2 ⟨x, hx, rfl⟩))

/-! ## ★出典の紐付け(`.src`) -/

def exists_descend_data.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(帰納段のデータ——α₀ の最小多項式 f₀ と S′ ≝ f₀(S) ∪ f₀(S₀))",
    sectionId := "ncbelyi-lemma-2-4" }

def exists_poly_image_rat.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.4(Reduction to the Rational Case——多項式の段)",
    sectionId := "ncbelyi-lemma-2-4" }

end ABC3.Found.NCBelyi
