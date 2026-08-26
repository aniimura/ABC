import ABC3.Found.NCBelyi.Separation

/-!
# [NCBelyi] Lemma 2.3 の有理版と、Lemma 2.4 の正規化(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.4–p.5。

原文 (NCBelyi p.4):
> if β ∈Q, then one may take λ ∈Q.

## ★★★本ファイルが取るもの

`Separation.lean` は `Lemma 2.3` の本体(`λ ∈ ℂ`)を取ったが、
★**「`β ∈ ℚ` なら `λ ∈ ℚ` に取れる」という原文の追記を落としていた**。
`Lemma 2.4` は『**with rational coefficients!**』と強調してこの追記を使うので、ここで取る。

原文 (NCBelyi p.5):
> multiplying by some positive rational number, we may assume that |α| ≤1, for all

`Lemma 2.4` の正規化はこの 2 段である:

| 原文 | ここでの形 |
|---|---|
| 「`Lemma 2.3` のような自己同型(有理係数!)」 | `lemma_2_3_rat` |
| 「正の有理数を掛けて `|α| ≤ 1`、`|β| ≥ C`」 | `exists_rat_normalization` |

## ★★★★有理数の余地を作るために `2C` で引く

★`λ` を決めたあとの倍率 `c` は、`|c·f(α)| ≤ 1` と `|c·f(β)| ≥ C` の**両方**を要求する。
`Lemma 2.3` を `C` で引くと `c = C/M` **ちょうど**しか許されず、有理数に取れない。

★★そこで **`2C` で引く**。すると `c ∈ (C/M, 2C/M)` の**開区間**が空でなくなり、
`ℚ` の稠密性(`exists_rat_btwn`)で有理数の `c` が取れる。
★★★原文が「for some sufficiently large C」と書くのは、この余地のことである。
-/

namespace ABC3.Found.NCBelyi

/-! ## ★`absInvShift` の基本性質 -/

theorem absInvShift_nonneg (lam : ℂ) (p : P1C) : 0 ≤ absInvShift lam p := by
  cases p with
  | none => simp
  | some x => simp

theorem absInvShift_pos_of_ne {lam b : ℂ} (h : lam ≠ b) :
    0 < absInvShift lam (some b) := by
  rw [absInvShift_some]
  refine norm_pos_iff.2 ?_
  refine div_ne_zero one_ne_zero ?_
  exact sub_ne_zero.2 (Ne.symm h)

/-! ## ★★`S` の点と `β` の距離の下界 -/

/-- ★**`S` の有限点はすべて `β` から `δ` 以上離れている**ような `δ > 0` が取れる。

★`Separation.lean` の `lemma_2_3` の証明の前半をそのまま切り出したものである
(`∞` は距離 `1` に潰して有限集合の `min'` を取る)。 -/
theorem exists_dist_lower_bound (S : Finset P1C) (b : ℂ) (hb : (some b) ∉ S) :
    ∃ delta : ℝ, 0 < delta ∧ ∀ p ∈ S, ∀ x : ℂ, p = some x → delta ≤ ‖x - b‖ := by
  classical
  set g : P1C → ℝ := fun p => match p with | none => 1 | some x => ‖x - b‖ with hg
  set D : Finset ℝ := S.image g with hD
  have hgpos : ∀ p ∈ S, 0 < g p := by
    intro p hp
    cases p with
    | none => simp [hg]
    | some x =>
      have hxb : x ≠ b := by
        intro hcon; rw [hcon] at hp; exact hb hp
      simp only [hg]
      exact norm_pos_iff.2 (sub_ne_zero.2 hxb)
  set delta : ℝ := if h : D.Nonempty then D.min' h else 1 with hdel
  have hdelta : 0 < delta := by
    simp only [hdel]
    split_ifs with h
    · obtain ⟨p, hp, hpe⟩ := Finset.mem_image.1 (D.min'_mem h)
      rw [← hpe]
      exact hgpos p hp
    · norm_num
  refine ⟨delta, hdelta, ?_⟩
  intro p hp x hpx
  have hmem : g p ∈ D := Finset.mem_image_of_mem g hp
  have hne : D.Nonempty := ⟨g p, hmem⟩
  have hle : delta ≤ g p := by simp only [hdel, dif_pos hne]; exact D.min'_le _ hmem
  rw [hpx] at hle
  simpa [hg] using hle

/-! ## ★★★★★`Lemma 2.3` の有理版 -/

/-- ★★★★★**[NCBelyi] Lemma 2.3 の追記** —— `β ∈ ℚ` なら `λ ∈ ℚ` に取れる。

原文 (NCBelyi p.4):
> if β ∈Q, then one may take λ ∈Q.

★`separation_core` は `λ ≝ β + ε` を返し、`ε` は
`0 < ε ≤ min(δ/4, δ/(4C))` を満たす**任意の実数**でよい。
★★だから **`ε` を有理数に取ればそのまま `λ ∈ ℚ`** である
——原文の「one may take」はこの自由度のことである。 -/
theorem lemma_2_3_rat (S : Finset P1C) (C : ℝ) (hC : 0 < C) (b : ℚ)
    (hb : (some (b : ℂ)) ∉ S) :
    ∃ lam : ℚ, (lam : ℂ) ≠ (b : ℂ) ∧ (some (lam : ℂ) ∉ S)
      ∧ ∀ p ∈ S, C * absInvShift (lam : ℂ) p
          ≤ absInvShift (lam : ℂ) (some (b : ℂ)) := by
  obtain ⟨delta, hdelta, hSle⟩ := exists_dist_lower_bound S (b : ℂ) hb
  have hmin : (0 : ℝ) < min (delta / 4) (delta / (4 * C)) :=
    lt_min (by positivity) (by positivity)
  obtain ⟨q, hq0, hqlt⟩ := exists_rat_btwn hmin
  refine ⟨b + q, ?_⟩
  have heps : (0 : ℝ) < (q : ℝ) := hq0
  have h1 : (q : ℝ) ≤ delta / 4 := le_trans hqlt.le (min_le_left _ _)
  have h2 : C * (q : ℝ) ≤ delta / 4 := by
    have hq2 : (q : ℝ) ≤ delta / (4 * C) := le_trans hqlt.le (min_le_right _ _)
    calc C * (q : ℝ) ≤ C * (delta / (4 * C)) := mul_le_mul_of_nonneg_left hq2 hC.le
      _ = delta / 4 := by field_simp
  have hcast : (((b + q : ℚ)) : ℂ) = (b : ℂ) + ((q : ℝ) : ℂ) := by push_cast; ring
  rw [hcast]
  exact separation_core S C hC (b : ℂ) delta (q : ℝ) hdelta heps hSle h1 h2

/-! ## ★★★★★★★`Lemma 2.4` の正規化 -/

/-- ★★★★★★★**正規化** —— `|f(α)| ≤ 1`(`∀ α ∈ S`)かつ `|f(β)| ≥ C`。

原文 (NCBelyi p.5):
> multiplying by some positive rational number, we may assume that |α| ≤1, for all

`f(x) = c/(x − λ)`(`λ, c ∈ ℚ`、`c > 0`)が取れる。

★★★★**`2C` で `Lemma 2.3` を引くのが要点**である。
`C` で引くと倍率は `C/M` ちょうどしか許されず有理数に取れないが、
`2C` で引けば `c ∈ (C/M, 2C/M)` の開区間が空でなくなり `ℚ` の稠密性が効く。 -/
theorem exists_rat_normalization (S : Finset P1C) (C : ℝ) (hC : 0 < C) (b : ℚ)
    (hb : (some (b : ℂ)) ∉ S) :
    ∃ lam c : ℚ, 0 < c ∧ (lam : ℂ) ≠ (b : ℂ) ∧ (some (lam : ℂ) ∉ S)
      ∧ (∀ p ∈ S, (c : ℝ) * absInvShift (lam : ℂ) p ≤ 1)
      ∧ C ≤ (c : ℝ) * absInvShift (lam : ℂ) (some (b : ℂ)) := by
  obtain ⟨lam, hne, hnotmem, hsep⟩ := lemma_2_3_rat S (2 * C) (by linarith) b hb
  have hM : 0 < absInvShift (lam : ℂ) (some (b : ℂ)) := absInvShift_pos_of_ne hne
  have hlt : C / absInvShift (lam : ℂ) (some (b : ℂ))
      < 2 * C / absInvShift (lam : ℂ) (some (b : ℂ)) := by
    rw [div_lt_div_iff_of_pos_right hM]
    linarith
  obtain ⟨c, hc1, hc2⟩ := exists_rat_btwn hlt
  have hcpos : (0 : ℝ) < (c : ℝ) := lt_trans (by positivity) hc1
  refine ⟨lam, c, by exact_mod_cast hcpos, hne, hnotmem, ?_, ?_⟩
  · intro p hp
    have hle : absInvShift (lam : ℂ) p
        ≤ absInvShift (lam : ℂ) (some (b : ℂ)) / (2 * C) := by
      rw [le_div_iff₀ (by linarith)]
      have := hsep p hp
      linarith
    have hcle : (c : ℝ) ≤ 2 * C / absInvShift (lam : ℂ) (some (b : ℂ)) := hc2.le
    calc (c : ℝ) * absInvShift (lam : ℂ) p
        ≤ (2 * C / absInvShift (lam : ℂ) (some (b : ℂ)))
            * (absInvShift (lam : ℂ) (some (b : ℂ)) / (2 * C)) :=
          mul_le_mul hcle hle (absInvShift_nonneg _ _) (by positivity)
      _ = 1 := by field_simp
  · have hlt2 := (div_lt_iff₀ hM).1 hc1
    linarith

/-! ## ★出典の紐付け(`.src`) -/

def lemma_2_3_rat.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.3(Moreover, if β ∈ ℚ, then one may take λ ∈ ℚ)",
    sectionId := "ncbelyi-lemma-2-3" }

def exists_rat_normalization.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(正規化——有理係数の自己同型と正の有理数倍)",
    sectionId := "ncbelyi-lemma-2-4" }

end ABC3.Found.NCBelyi
