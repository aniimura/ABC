/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ChartwiseUniform
import ABC3.Found.GenEll.HeightMetric

/-!
# [GenEll] Proposition 1.4, (iii) —— **生成ファイバーが同じなら**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★★★★★★★★逸脱が 1 つ消えた

`Skeleton/GenEll/Section1.lean` の `prop_1_4` は逸脱を 3 つ記録していた。
その 2 番目:

> 2. ★★**(iii) の形**: 原文は『生成ファイバーが同じなら』。因子表示では
>    『因子が同じ(計量だけ違う)』である。**垂直因子の差の分は含めていない**。

★★本ファイルはその差を含める:

> **因子が `ℤ[1/N]` の上で一致し**（各アフィンチャートで）、計量が連続なら、
> 高さの差は**すべての数体に一様に**有界である。

★★★これで (iii) は原文の『`L_ℚ` の同型類だけに依る』に一段近づいた
——「因子が完全に同じ」ではなく「生成ファイバーで同じ」で足りる。

## ★★★★★★機構は三角不等式 1 本

`E′ ≝ ⟨E.divisor, D.green⟩` を挟む:

| 段 | 何が違うか | 上界 | 出どころ |
|---|---|---|---|
| `ht_D` vs `ht_{E′}` | **因子だけ** | `log(N^m)` | `ChartwiseUniform.lean`（本セッション） |
| `ht_{E′}` vs `ht_E` | **計量だけ** | `C₂` | `HeightMetric.lean` の `htArith_sub_abs_le` |

★★どちらの上界も**数体に依らない**ので、和も依らない。

## ★残っているもの（明示）

★`Proposition 1.4` 全体には (i)(ii)(iv) が要る。
★★とくに (iv) は「`L_ℚ` が ample なら射影埋め込みが得られる」という段を
**データとして受けている**（逸脱 3）——`ample-and-projective-embedding` の在庫不足。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits NumberField

/-! ## ★★★★★★★数体に一様な上界 -/

/-- ★★★★★★★**層の水準で `n`-比較できるなら、すべての数体に一様な上界**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★`ChartwiseUniform.lean` の到達点は `BDeq`（数体ごと）だったが、
定数 `log n` は**数体に依らない**ので、一様な形で述べ直せる。
★★これが `htArith_sub_abs_le`（計量側）と足し合わせられる形である。 -/
theorem htArith_sub_abs_le_of_sheaf_comparable {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (D E : ArithCartier X) (n : ℕ) (hn : n ≠ 0)
    (h1 : (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(n : ℤ)})).comap f * D.divisor
        ≤ E.divisor)
    (h2 : (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(n : ℤ)})).comap f * E.divisor
        ≤ D.divisor)
    (hD0 : ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
      pullbackIdeal F D.divisor xF ≠ 0)
    (hE0 : ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
      pullbackIdeal F E.divisor xF ≠ 0)
    (harc : ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
      (archADiv F D.green xF).sum (fun _ r => r)
        = (archADiv F E.green xF).sum (fun _ r => r)) :
    ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
      |htArith F D xF - htArith F E xF| ≤ Real.log n := by
  intro F _ _ xF
  rw [htArith_sub_eq_degNormalized_sub F D E xF xF (harc F xF)]
  exact abs_degNormalized_idealADiv_sub_le F _ _ (hD0 F xF) (hE0 F xF) n hn
    (span_natCast_mul_pullbackIdeal_le F f n D.divisor E.divisor h1 xF)
    (span_natCast_mul_pullbackIdeal_le F f n E.divisor D.divisor h2 xF)

/-! ## ★★★★★★★★★★到達点 —— 逸脱 2 を含めた `(iii)` -/

/-- ★★★★★★★★★★**[GenEll] Proposition 1.4, (iii) —— 生成ファイバーが同じ場合**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

因子が各アフィンチャートで `N` を反転して一致し、計量がともに連続なら、
高さの差は**すべての数体に一様に**有界である。

★★`Skeleton/GenEll/Section1.lean` の `prop_1_4` の (iii) は
「因子が**完全に**同じ」を仮定していた（逸脱 2）。
★★★本定理はそれを「**生成ファイバーで**同じ」に弱める。

★★★★機構は `E′ ≝ ⟨E.divisor, D.green⟩` を挟む三角不等式 1 本である:
因子だけの差は `ChartwiseUniform.lean`、計量だけの差は `HeightMetric.lean`。 -/
theorem prop_1_4_iii_of_chartwise {X : Scheme.{0}} {V : Type}
    [NormedAddCommGroup V] [NormedSpace ℂ V] [FiniteDimensional ℂ V]
    (M : ArcModel X V) [Nonempty (complexPoints X)]
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (D E : ArithCartier X) (N : ℕ) (hN : N ≠ 0)
    {ι : Type} [Fintype ι] [DecidableEq ι] (U : ι → X.affineOpens)
    (hcov : ⨆ i, ((U i : X.affineOpens) : X.Opens) = ⊤)
    (hnoeth : ∀ i, IsNoetherianRing Γ((U i).1.toScheme, ⊤))
    (hDE : ∀ i,
      ((D.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤))))))
        ≤ ((E.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤)))))))
    (hED : ∀ i,
      ((E.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤))))))
        ≤ ((D.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤)))))))
    (hgD : @Continuous _ _ M.topology _ D.green)
    (hgE : @Continuous _ _ M.topology _ E.green)
    (hD0 : ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
      pullbackIdeal F D.divisor xF ≠ 0)
    (hE0 : ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
      pullbackIdeal F E.divisor xF ≠ 0) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
      (xF : specRingOfIntegers F ⟶ X), |htArith F D xF - htArith F E xF| ≤ C := by
  obtain ⟨m₁, h₁⟩ :=
    exists_uniform_sheaf_le_of_localization f D.divisor E.divisor (N : ℤ) U hcov hnoeth hDE
  obtain ⟨m₂, h₂⟩ :=
    exists_uniform_sheaf_le_of_localization f E.divisor D.divisor (N : ℤ) U hcov hnoeth hED
  set E' : ArithCartier X := ⟨E.divisor, D.green⟩ with hE'
  have hsub : ∀ k : ℕ, k ≤ max m₁ m₂ →
      (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(N : ℤ) ^ max m₁ m₂})).comap f
        ≤ (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(N : ℤ) ^ k})).comap f := by
    intro k hk
    refine Scheme.IdealSheafData.comap_mono f (specIdealSheaf_mono ?_)
    rw [Ideal.span_singleton_le_span_singleton]
    exact ⟨(N : ℤ) ^ (max m₁ m₂ - k), by rw [← pow_add]; congr 1; omega⟩
  have hcast : ((N ^ max m₁ m₂ : ℕ) : ℤ) = (N : ℤ) ^ max m₁ m₂ := by push_cast; ring
  have hA := htArith_sub_abs_le_of_sheaf_comparable f D E' (N ^ max m₁ m₂)
    (pow_ne_zero _ hN)
    (by rw [hcast]; exact le_trans (mul_le_mul' (hsub m₁ (le_max_left _ _)) le_rfl) h₁)
    (by rw [hcast]; exact le_trans (mul_le_mul' (hsub m₂ (le_max_right _ _)) le_rfl) h₂)
    hD0 hE0 (fun F _ _ xF => rfl)
  obtain ⟨C₂, hC₂, hB⟩ := htArith_sub_abs_le M E' E rfl hgD hgE
  have hlog : (0 : ℝ) ≤ Real.log ((N ^ max m₁ m₂ : ℕ) : ℝ) := by
    apply Real.log_nonneg
    exact_mod_cast Nat.one_le_iff_ne_zero.2 (pow_ne_zero (max m₁ m₂) hN)
  refine ⟨Real.log ((N ^ max m₁ m₂ : ℕ) : ℝ) + C₂, by linarith, fun F _ _ xF => ?_⟩
  calc |htArith F D xF - htArith F E xF|
      ≤ |htArith F D xF - htArith F E' xF| + |htArith F E' xF - htArith F E xF| :=
        abs_sub_le _ _ _
    _ ≤ Real.log ((N ^ max m₁ m₂ : ℕ) : ℝ) + C₂ := add_le_add (hA F xF) (hB F xF)

/-! ### ★出典の紐付け(`.src`)

★★**条つきである。** `Proposition 1.4` 全体には (i)(ii)(iv) が要り、
とくに (iv) は射影埋め込みをデータとして受けている（逸脱 3）。 -/

def htArith_sub_abs_le_of_sheaf_comparable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(層の水準で n-比較できる場合——数体に一様)",
    sectionId := "genell-prop-1-4" }

def prop_1_4_iii_of_chartwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(生成ファイバーが同じ場合。計量の差も込み)",
    sectionId := "genell-prop-1-4" }

def prop_1_4_iii_of_chartwise.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]"
      "htArith_bdeq_of_chartwise_localization(因子だけの差——垂直な差の一様有界性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_bdeq_of_chartwise_localization") 6,
    .citation "[ABC3]" "htArith_sub_abs_le(計量だけの差)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_sub_abs_le") 6,
    .implicitStep
      ("★★★★★**Skeleton の逸脱 2 が消えた**" ++
       "——(iii) は「因子が完全に同じ」ではなく「生成ファイバーで同じ」で足りる。" ++
       "★機構は E′ ≝ ⟨E.divisor, D.green⟩ を挟む三角不等式 1 本である") 6,
    .implicitStep
      ("★★残る逸脱 3: (iv) は「L_ℚ が ample なら射影埋め込みが得られる」段を" ++
       "データとして受けている(ample-and-projective-embedding の在庫不足)。" ++
       "★したがって Proposition 1.4 全体の .src はまだ置けない") 6 ]

end ABC3.Found.GenEll
