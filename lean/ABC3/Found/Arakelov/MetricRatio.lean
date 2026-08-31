/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.Prop14ii
import ABC3.Found.GenEll.Prop14Proper
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★計量を取り替えても高さの BD-類は変わらない（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> The BD-class of htL depends only on the isomorphism class of the line bundle LQ on XQ.

## ★★★★★★★★★★これは何か —— `Proposition 1.4, (iii)` の**計量の側**

原文 (iii) は「高さの BD-類は `L_ℚ` の同型類だけに依る」と言う。
★これは 2 つのことを含む:

| 条 | 内容 | 本ファイル |
|---|---|---|
| (a) | **計量**を取り替えても BD-類は同じ | ★★★**閉じた** |
| (b) | `L_ℚ ≅ M_ℚ` なら BD-類は同じ（2 つの `ℤ`-モデルの比較） | 残る |

★★(a) が本ファイルである。原文の証明は (i)(ii) を `L̄ ⊗ M̄⁻¹` に当てるが、
**計量だけを動かす場合はもっと直接に出る**——比が `X^arc` 上の正の連続関数だからである。

## ★★★★★★★★機構 —— 比 `h'/h` はチャートに依らない

計量は自明化ごとの基準ノルム `h_{V,e}` で書かれているので、
2 つの計量の比 `h'_{V,e}(p) / h_{V,e}(p)` を考えることができる。

★★**この比はチャート `(V,e)` に依らない**（`hRatio_chart_indep`）
——遷移単元 `u` は計量に依らないので、`compat` の欄で
分子・分母が**同じ `‖u‖` 倍**され、比では消える。

★★★したがって `hRatio m m' : X^arc → ℝ` が大域に定まり、

    `|s|_{m'}(p) = hRatio(p) · |s|_m(p)`   （`norm_ofMetric_eq`）

となる。★`X` が固有なら `X^arc` はコンパクトで、`hRatio` は正の連続関数なので
`|log hRatio|` は一様に有界である。

## ★★★★★★★引き戻しとの両立

    `hRatio (f^*m) (f^*m') (p) = hRatio m m' (p ≫ f)`   （`hRatio_pullback`）

★これがあるので、**定数は `F` にも `x_F` にも依らない**。
★★機構は `pullback_h_eq_of_chart`（`§9-802`）——引き戻した `h` は
`‖遷移単元‖⁻¹ · m.h(p ≫ f)` の形であり、`‖遷移単元‖` は計量に依らないので比で消える。

## ★★★★★★到達点

    exists_htMetricU_sub_abs_le : ∃C, ∀F ∀x_F, |ht'(x_F) − ht(x_F)| ≤ C
    heightBDClass_ofMetric_eq   : ★★BD-類が**等しい**

## ★残っている段（明示）

★(b)——`L_ℚ ≅ M_ℚ` から BD-類の一致を出す段。原文は (i)(ii) を `L̄ ⊗ M̄⁻¹` に
当てるが、そこでは「生成ファイバーが自明なら、ある捻りの後で大域切断を持つ」
という段が要る（因子表示では `Found/GenEll/VerticalTwist.lean` が持つ）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

variable {X : Scheme.{0}} {F : X.PresheafOfModules}

/-! ## ★★★(1) 比はチャートに依らない -/

/-- ★★**同じ開集合の上で自明化を取り替えても比は変わらない**。

★遷移単元 `u` は計量に依らないので、分子・分母が同じ `‖u‖` 倍される。 -/
theorem hRatio_triv_indep (m m' : LocalMetric X F) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) :
    m'.h V e' p / m.h V e' p = m'.h V e p / m.h V e p := by
  have h1 := m.compat V e e' p hp
  have h2 := m'.compat V e e' p hp
  have hm : m.h V e' p ≠ 0 := (m.pos V e' p).ne'
  have hm' : m.h V e p ≠ 0 := (m.pos V e p).ne'
  field_simp
  rw [← h1, ← h2]
  ring

theorem hRatio_restrict (m m' : LocalMetric X F) {V W : X.Opens} (hWV : W ≤ V)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hpW : p ⁻¹ᵁ W = ⊤) :
    m'.h W (trivialOfLe F hWV e) p / m.h W (trivialOfLe F hWV e) p
      = m'.h V e p / m.h V e p := by
  rw [m.restrict hWV e p hpW, m'.restrict hWV e p hpW]

/-- ★★★★★**比はチャートに全く依らない**（共通の細分へ降りる）。 -/
theorem hRatio_chart_indep (m m' : LocalMetric X F) {V V' : X.Opens}
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (e' : (restrictPresheafFunctor X V').obj F ≅ 𝟙_ (PresheafModulesOn X V'))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (hp' : p ⁻¹ᵁ V' = ⊤) :
    m'.h V e p / m.h V e p = m'.h V' e' p / m.h V' e' p := by
  have hpW : p ⁻¹ᵁ (V ⊓ V') = ⊤ := by
    show p ⁻¹ᵁ V ⊓ p ⁻¹ᵁ V' = ⊤
    rw [hp, hp', inf_idem]
  rw [← hRatio_restrict m m' (inf_le_left : V ⊓ V' ≤ V) e p hpW,
    ← hRatio_restrict m m' (inf_le_right : V ⊓ V' ≤ V') e' p hpW]
  exact hRatio_triv_indep m m' _ _ _ p hpW

/-! ## ★★★★★★(2) 大域の比 -/

open scoped Classical in
/-- ★★★★★★**2 つの計量の比** `h'/h : X^arc → ℝ`（チャートに依らない）。 -/
noncomputable def hRatio (m m' : LocalMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  if hc : Nonempty (NormChart X F p) then
    m'.h hc.some.V hc.some.e p / m.h hc.some.V hc.some.e p else 1

open scoped Classical in
@[simp] theorem hRatio_eq (m m' : LocalMetric X F) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) :
    hRatio m m' p = m'.h V e p / m.h V e p := by
  have hc : Nonempty (NormChart X F p) := ⟨⟨V, hp, e⟩⟩
  show (if hc : Nonempty (NormChart X F p) then
      m'.h hc.some.V hc.some.e p / m.h hc.some.V hc.some.e p else 1) = _
  rw [dif_pos hc]
  exact hRatio_chart_indep m m' hc.some.e e p hc.some.hp hp

theorem hRatio_pos (m m' : LocalMetric X F) (htriv : IsLocallyTrivial X F)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : 0 < hRatio m m' p := by
  obtain ⟨c⟩ := nonempty_normChart htriv p
  rw [hRatio_eq m m' c.V c.e p c.hp]
  exact div_pos (m'.pos _ _ _) (m.pos _ _ _)

/-- ★★★★★**比は連続である**（両方の計量が連続なら）。 -/
theorem continuous_hRatio (m m' : LocalMetric X F) (htriv : IsLocallyTrivial X F)
    (hm : m.IsContinuous) (hm' : m'.IsContinuous) :
    @Continuous _ ℝ (arcTopology X) _ (fun p => hRatio m m' p) := by
  letI := arcTopology X
  rw [continuous_iff_continuousAt]
  intro p₀
  obtain ⟨c⟩ := nonempty_normChart htriv p₀
  have hcont : @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => hRatio m m' p) (arcOpenSet c.V) := by
    refine ((hm' c.V c.e).div (hm c.V c.e) (fun p _ => (m.pos _ _ _).ne')).congr (fun p hp => ?_)
    show hRatio m m' p = m'.h c.V c.e p / m.h c.V c.e p
    exact hRatio_eq m m' c.V c.e p hp
  exact hcont.continuousAt ((isOpen_arcOpenSet c.V).mem_nhds c.hp)

/-- ★★★★★★★**`X` が固有なら `|log(h'/h)|` は一様に有界である**。 -/
theorem exists_abs_log_hRatio_le [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    (m m' : LocalMetric X F) (htriv : IsLocallyTrivial X F)
    (hm : m.IsContinuous) (hm' : m'.IsContinuous) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ p : Spec (CommRingCat.of ℂ) ⟶ X,
      |Real.log (hRatio m m' p)| ≤ C := by
  letI := arcTopology X
  haveI := compactSpace_arc hval
  obtain ⟨C, hC, hlo, hhi⟩ := exists_abs_bound_of_continuous
    (fun p => Real.log (hRatio m m' p))
    ((continuous_hRatio m m' htriv hm hm').log (fun p => (hRatio_pos m m' htriv p).ne'))
  exact ⟨C, hC, fun p => abs_le.2 ⟨hlo p, hhi p⟩⟩

/-! ## ★★★★★★(3) ノルムは比の分だけずれる -/

/-- ★同じ前層に 2 つの計量を載せたもの。 -/
def ofMetric (F : X.PresheafOfModules) (triv : IsLocallyTrivial X F)
    (m : LocalMetric X F) : AMetric X := ⟨F, triv, m⟩

/-- ★★★★★★**`|s|_{m'} = (h'/h) · |s|_m`**。 -/
theorem norm_ofMetric_eq (triv : IsLocallyTrivial X F) (m m' : LocalMetric X F)
    (s : (F.obj (op ⊤) : Type)) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (ofMetric F triv m').norm s p = hRatio m m' p * (ofMetric F triv m).norm s p := by
  obtain ⟨c⟩ := nonempty_normChart triv p
  rw [AMetric.norm_eq (ofMetric F triv m') s p c.V c.e c.hp,
    AMetric.norm_eq (ofMetric F triv m) s p c.V c.e c.hp,
    hRatio_eq m m' c.V c.e p c.hp]
  show trivSecNorm F c.V c.e p c.hp s * m'.h c.V c.e p
    = m'.h c.V c.e p / m.h c.V c.e p * (trivSecNorm F c.V c.e p c.hp s * m.h c.V c.e p)
  have hm : m.h c.V c.e p ≠ 0 := (m.pos _ _ _).ne'
  field_simp

/-! ## ★★★★★★★(4) 引き戻しとの両立 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**比は引き戻しと両立する** `hRatio (f^*m) (f^*m') p = hRatio m m' (p ≫ f)`。

★これがあるので定数は `F` にも `x_F` にも依らない。 -/
theorem hRatio_pullback {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (hL : IsLocallyTrivial Y L) (m m' : LocalMetric Y L)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    hRatio (LocalMetric.pullback f hL m) (LocalMetric.pullback f hL m') p
      = hRatio m m' (p ≫ f) := by
  obtain ⟨c⟩ := nonempty_normChart (isLocallyTrivial_pullbackPre f L hL) p
  rw [hRatio_eq _ _ c.V c.e p c.hp]
  obtain ⟨d⟩ := nonempty_pullChart f hL c.V p c.hp
  rw [pullback_h_eq_of_chart f hL m' d.hV'V d.hV'W d.eW c.e p d.hpV',
    pullback_h_eq_of_chart f hL m d.hV'V d.hV'W d.eW c.e p d.hpV',
    hRatio_eq m m' d.W d.eW (p ≫ f) d.comp_preimage_eq_top]
  have hu : ‖evalOn p d.V' d.hpV' (transUnit ((pullbackPre f).obj L) d.V'
      (pullTrivOfBase f L d.W d.eW d.hV'W)
      (trivialOfLe ((pullbackPre f).obj L) d.hV'V c.e))‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p d.V' d.hpV'
      (isUnit_transUnit _ d.V' _ _))
  have hm : m.h d.W d.eW (p ≫ f) ≠ 0 := (m.pos _ _ _).ne'
  field_simp

/-! ## ★★★★★★★★(5) 次数と高さの差 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**アルキメデス部分の差は `C` 以内**。 -/
theorem archDeg_ofMetric_sub_abs_le (F : Type) [Field F] [NumberField F]
    {P : (Spec (CommRingCat.of (𝓞 F))).PresheafOfModules}
    (triv : IsLocallyTrivial _ P) (m m' : LocalMetric _ P)
    (s : (P.obj (op ⊤) : Type))
    (hs : ∀ σ : F →+* ℂ, (ofMetric P triv m).norm s (embSpecPoint F σ) ≠ 0)
    (C : ℝ) (hb : ∀ p, |Real.log (hRatio m m' p)| ≤ C) :
    |archDeg F (ofMetric P triv m') s - archDeg F (ofMetric P triv m) s| ≤ C := by
  have hn : (0:ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have hlog : ∀ σ : F →+* ℂ,
      Real.log ((ofMetric P triv m').norm s (embSpecPoint F σ))
        = Real.log (hRatio m m' (embSpecPoint F σ))
          + Real.log ((ofMetric P triv m).norm s (embSpecPoint F σ)) := by
    intro σ
    rw [norm_ofMetric_eq triv m m' s, Real.log_mul (hRatio_pos m m' triv _).ne' (hs σ)]
  have hA : |∑ σ : F →+* ℂ, Real.log (hRatio m m' (embSpecPoint F σ))|
      ≤ (Module.finrank ℚ F : ℝ) * C := by
    refine le_trans (Finset.abs_sum_le_sum_abs _ _) ?_
    have h1 : (∑ σ : F →+* ℂ, |Real.log (hRatio m m' (embSpecPoint F σ))|)
        ≤ (∑ _σ : F →+* ℂ, C) := Finset.sum_le_sum (fun σ _ => hb _)
    have h2 : (∑ _σ : F →+* ℂ, C) = (Fintype.card (F →+* ℂ) : ℝ) * C := by
      rw [Finset.sum_const, nsmul_eq_mul, Finset.card_univ]
    rw [h2, NumberField.Embeddings.card F ℂ] at h1
    exact h1
  show |(-(∑ σ : F →+* ℂ, Real.log ((ofMetric P triv m').norm s (embSpecPoint F σ)))
      / (Module.finrank ℚ F : ℝ))
    - (-(∑ σ : F →+* ℂ, Real.log ((ofMetric P triv m).norm s (embSpecPoint F σ)))
      / (Module.finrank ℚ F : ℝ))| ≤ C
  simp only [hlog, Finset.sum_add_distrib]
  rw [show ∀ a b n : ℝ, -(a + b)/n - (-b/n) = -a/n from fun a b n => by ring]
  rw [abs_div, abs_neg, abs_of_pos hn, div_le_iff₀ hn]
  linarith [hA]

set_option maxHeartbeats 1000000 in
theorem degArithPre_ofMetric_sub_abs_le (F : Type) [Field F] [NumberField F]
    {P : (Spec (CommRingCat.of (𝓞 F))).PresheafOfModules}
    (triv : IsLocallyTrivial _ P) (m m' : LocalMetric _ P)
    (inv inv' : AMetric (Spec (CommRingCat.of (𝓞 F))))
    (hi : Isometric ((ofMetric P triv m) * inv) 1)
    (hi' : Isometric ((ofMetric P triv m') * inv') 1)
    (s : (P.obj (op ⊤) : Type))
    (hs : ∀ σ : F →+* ℂ, (ofMetric P triv m).norm s (embSpecPoint F σ) ≠ 0)
    (C : ℝ) (hb : ∀ p, |Real.log (hRatio m m' p)| ≤ C) :
    |degArithPre F ⟨ofMetric P triv m', inv', hi'⟩ s
      - degArithPre F ⟨ofMetric P triv m, inv, hi⟩ s| ≤ C := by
  have hfin : degFinPre (⟨ofMetric P triv m', inv', hi'⟩ : AInv _) s
      = degFinPre (⟨ofMetric P triv m, inv, hi⟩ : AInv _) s := rfl
  show |degFinPre (⟨ofMetric P triv m', inv', hi'⟩ : AInv _) s / (Module.finrank ℚ F : ℝ)
      + archDeg F (ofMetric P triv m') s
      - (degFinPre (⟨ofMetric P triv m, inv, hi⟩ : AInv _) s / (Module.finrank ℚ F : ℝ)
        + archDeg F (ofMetric P triv m) s)| ≤ C
  rw [hfin, show ∀ a b c : ℝ, a + b - (a + c) = b - c from fun a b c => by ring]
  exact archDeg_ofMetric_sub_abs_le F triv m m' s hs C hb

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**`deg_F` の差は `C` 以内**（有限部分は同じ、アルキメデス部分だけがずれる）。 -/
theorem degAInv_ofMetric_sub_abs_le (F : Type) [Field F] [NumberField F]
    {P : (Spec (CommRingCat.of (𝓞 F))).PresheafOfModules}
    (triv : IsLocallyTrivial _ P) (m m' : LocalMetric _ P)
    (inv inv' : AMetric (Spec (CommRingCat.of (𝓞 F))))
    (hi : Isometric ((ofMetric P triv m) * inv) 1)
    (hi' : Isometric ((ofMetric P triv m') * inv') 1)
    (C : ℝ) (hb : ∀ p, |Real.log (hRatio m m' p)| ≤ C) :
    |degAInv F ⟨ofMetric P triv m', inv', hi'⟩
      - degAInv F ⟨ofMetric P triv m, inv, hi⟩| ≤ C := by
  set L : AInv (Spec (CommRingCat.of (𝓞 F))) := ⟨ofMetric P triv m, inv, hi⟩ with hL
  set L' : AInv (Spec (CommRingCat.of (𝓞 F))) := ⟨ofMetric P triv m', inv', hi'⟩ with hL'
  obtain ⟨s, hs⟩ := exists_ne_zero_gammaPre L
  rw [degAInv_eq F L' s hs, degAInv_eq F L s hs]
  exact degArithPre_ofMetric_sub_abs_le F triv m m' inv inv' hi hi' s
    (fun σ => norm_ne_zero_of_ne_zero F L s hs σ) C hb

/-! ## ★★★★★★★★★★(6) 高さと BD-類 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**計量を取り替えても高さの差は一様に有界**。

原文 (GenEll p.6):
> The BD-class of htL depends only on the isomorphism class of the line bundle LQ on XQ.

★★定数は `F` にも `x_F` にも依らない——`hRatio_pullback` がそれを保証する。 -/
theorem exists_htMetricU_sub_abs_le {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    {P : X.PresheafOfModules} (triv : IsLocallyTrivial X P) (m m' : LocalMetric X P)
    (hm : m.IsContinuous) (hm' : m'.IsContinuous)
    (inv inv' : AMetric X)
    (hi : Isometric ((ofMetric P triv m) * inv) 1)
    (hi' : Isometric ((ofMetric P triv m') * inv') 1) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
      (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X),
      |htMetricU F ⟨ofMetric P triv m', inv', hi'⟩ xF
        - htMetricU F ⟨ofMetric P triv m, inv, hi⟩ xF| ≤ C := by
  obtain ⟨C, hC, hb⟩ := exists_abs_log_hRatio_le hval m m' triv hm hm'
  refine ⟨C, hC, fun F _ _ xF => ?_⟩
  have hb' : ∀ q : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (𝓞 F)),
      |Real.log (hRatio (LocalMetric.pullback xF triv m)
        (LocalMetric.pullback xF triv m') q)| ≤ C := by
    intro q
    rw [hRatio_pullback xF triv m m' q]
    exact hb _
  exact degAInv_ofMetric_sub_abs_le F (isLocallyTrivial_pullbackPre xF P triv)
    (LocalMetric.pullback xF triv m) (LocalMetric.pullback xF triv m')
    (AMetricPullback xF inv) (AMetricPullback xF inv') _ _ C hb'

set_option maxHeartbeats 1000000 in
theorem htMetricAlg_sub_abs_le {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    {P : X.PresheafOfModules} (triv : IsLocallyTrivial X P) (m m' : LocalMetric X P)
    (hm : m.IsContinuous) (hm' : m'.IsContinuous)
    (inv inv' : AMetric X)
    (hi : Isometric ((ofMetric P triv m) * inv) 1)
    (hi' : Isometric ((ofMetric P triv m') * inv') 1) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ x : AlgPointAnyClass X,
      |htMetricAlg (⟨ofMetric P triv m', inv', hi'⟩ : AInv X) x
        - htMetricAlg (⟨ofMetric P triv m, inv, hi⟩ : AInv X) x| ≤ C := by
  obtain ⟨C, hC, h⟩ := exists_htMetricU_sub_abs_le hval triv m m' hm hm' inv inv' hi hi'
  refine ⟨C, hC, fun x => ?_⟩
  induction x using AlgPointAnyClass.ind with
  | _ p => exact @h p.fld p.instField p.instNF p.map

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**計量を取り替えても BD-類は等しい**——
`Proposition 1.4, (iii)` の**計量の側**。

原文 (GenEll p.6):
> The BD-class of htL depends only on the isomorphism class of the line bundle LQ on XQ. -/
theorem heightBDClass_ofMetric_eq {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    {P : X.PresheafOfModules} (triv : IsLocallyTrivial X P) (m m' : LocalMetric X P)
    (hm : m.IsContinuous) (hm' : m'.IsContinuous)
    (inv inv' : AMetric X)
    (hi : Isometric ((ofMetric P triv m) * inv) 1)
    (hi' : Isometric ((ofMetric P triv m') * inv') 1)
    (S : Set (AlgPointAnyClass X)) :
    heightBDClass (⟨ofMetric P triv m', inv', hi'⟩ : AInv X) S
      = heightBDClass (⟨ofMetric P triv m, inv, hi⟩ : AInv X) S := by
  obtain ⟨C, hC, h⟩ := htMetricAlg_sub_abs_le hval triv m m' hm hm' inv inv' hi hi'
  refine (BDClass.mk_eq_mk _ _).2 ⟨C, fun x => ?_⟩
  exact h (x : AlgPointAnyClass X)

/-! ## ★出典の紐付け(`.src`) -/

def hRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(2 つの計量の比はチャートに依らない)",
    sectionId := "genell-prop-1-4" }

def hRatio_pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(比は引き戻しと両立する——定数が F に依らない理由)",
    sectionId := "genell-prop-1-4" }

def heightBDClass_ofMetric_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iii)(計量を取り替えても BD-類は等しい——計量の側のみ)",
    sectionId := "genell-prop-1-4" }

def heightBDClass_ofMetric_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "AMetric.continuous_norm / LocalMetric.IsContinuous(計量の連続性、§9-802)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.continuous_norm") 3,
    .citation "[ABC3]" "compactSpace_arc(固有性 ⟹ X^arc コンパクト)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.compactSpace_arc") 7,
    .citation "[ABC3]" "pullback_h_eq_of_chart(引き戻した h は固定チャートで書ける、§9-802)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullback_h_eq_of_chart") 3,
    .citation "[ABC3]" "BDClass / BDClass.mk_eq_mk(BD-類、Definition 1.2, (ii))"
      (.inProject "ABC3" "ABC3.Found.GenEll.BDClass") 5,
    .implicitStep
      ("★本ファイルは (iii) の**計量の側**だけである。残るのは " ++
       "「L_ℚ ≅ M_ℚ なら BD-類が同じ」——原文は (i)(ii) を L̄ ⊗ M̄⁻¹ に当てるが、" ++
       "そこでは『生成ファイバーが自明なら、ある捻りの後で大域切断を持つ』段が要る" ++
       "(因子表示では Found/GenEll/VerticalTwist.lean が持つ)") 6 ]

end ABC3.Found.Arakelov
