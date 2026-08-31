/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ConjMetricClosed
import ABC3.Found.Arakelov.ArcOpenSetOpen
import ABC3.Found.Arakelov.ArcEmbedding
import ABC3.Found.Arakelov.ArcEvalCont
import ABC3.Found.Arakelov.ArcMapCont
import ABC3.Found.Arakelov.AMetricNorm
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★計量は `X^arc` 上**連続**である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> pair consisting of a line bundle L on X and a hermitian metric

## ★★★★★★★★★★これは何か——「hermitian metric」の**もう一つの条項**

原文は「`X^arc` 上の直線束 `L^arc` の **hermitian 計量**」と書く。

原文 (GenEll p.3):
> Xarc for the compact normal complex analytic space determined by X

★`X^arc` は**複素解析空間**であり、そこでの計量は**連続**である。
★★ところが `LocalMetric` の欄は `h` / `pos` / `compat` / `restrict` の 4 つで、
**連続性が無かった**（2026-08-28 の実測）。★★★本ファイルはそれを入れる。

## ★★★★★★★★★到達点——切断のノルムは `X^arc` 上連続

    `AMetric.continuous_norm : Continuous (fun p => |s|_L̄ (p))`

★これが原文の「hermitian metric on `L^arc`」が**切断の水準で言っていること**である。

## ★★★★★★★★閉性

| 演算 | 定理 |
|---|---|
| 単位元 `Ō_X` | `structLocalMetric_isContinuous` |
| テンソル積 `⊗` | `LocalMetric.tensor_isContinuous` |
| 引き戻し `φ^*` | `LocalMetric.pullback_isContinuous` |
| 逆元 | `AInv.inv_isContinuous`（`isContinuous_of_tensor` から） |
| 等長同型 | `AMetric.isContinuous_of_isIsometry` |

★★★これで **`APicA X = APicC X ⊓ APicK X`**（共役両立かつ連続な類）が
`APicM X` の部分群になる——原文の `APic(X)` である。

## ★★★★★機構——「チャートの選択」を局所化で殺す

`LocalMetric.tensor` / `LocalMetric.pullback` の `h` は
`Nonempty (…Chart …)` の**選択**で書かれているので、式のまま連続性は言えない。

★そこで `p₀` の近くでは**固定したチャート**で書けること
（`tensor_h_eq_of_chart` / `pullback_h_eq_of_chart`、`chartH_indep` から）を示し、
`arcOpenSet W` が開である（`isOpen_arcOpenSet`）ことで
`ContinuousAt` に落とす。★★これが本ファイルの配管である。

## ★測定の記録

在庫がそのまま効いた（2026-08-28）:

    continuous_evalOn        : 開集合上の評価は連続（`ArcEvalCont.lean`）
    continuousOn_range_of_comp : `V.ι` の像の上での連続性（`ArcEmbedding.lean`）
    arcOpenSet_eq_range / isOpen_arcOpenSet : 像としての表示と開性
    continuous_comp_scheme   : `p ↦ p ≫ f` は連続（`ArcMapCont.lean`）
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

variable {X : Scheme.{0}} {F : X.PresheafOfModules}

/-! ## ★★★★★★★★(1) 条項そのもの -/

/-- ★★★★★★★★★★**[GenEll] Definition 1.1, (i)** ——
計量が `X^arc` 上**連続**であること。

原文 (GenEll p.3):
> pair consisting of a line bundle L on X and a hermitian metric

★「hermitian metric on the line bundle `L^arc` … on `X^arc`」の
**解析的な側**の中身である——`X^arc` は複素解析空間だから計量は連続である。
★★`p ⁻¹ᵁ V = ⊤` の側（`arcOpenSet V`）でのみ課す。そこを外した `h V e p` は
チャートの外の値であって意味を持たない。 -/
def LocalMetric.IsContinuous (m : LocalMetric X F) : Prop :=
  ∀ (V : X.Opens) (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)),
    @ContinuousOn _ ℝ (arcTopology X) _ (fun p => m.h V e p) (arcOpenSet V)

/-- ★★算術直線束の計量が連続であること。 -/
def AMetric.IsContinuous (L : AMetric X) : Prop := L.metric.IsContinuous

/-- ★★★**`V` の中へ持ち上げると連続性は「全体で連続」になる**。 -/
theorem continuous_h_comp {m : LocalMetric X F}
    (hm : m.IsContinuous) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    @Continuous _ ℝ (arcTopology V.toScheme) _ (fun q => m.h V e (q ≫ V.ι)) := by
  letI := arcTopology X
  letI := arcTopology V.toScheme
  exact (hm V e).comp_continuous (isOpenEmbedding_comp_ι V).continuous
    (fun q => comp_preimage_eq_top V q)

/-! ## ★★★★(2) 単位元 `Ō_X` -/

open scoped Classical in
/-- ★★★★★**構造層の標準計量は連続である**。

★`h` は遷移単元のノルムの逆数であり、`continuous_evalOn` がそのまま効く。 -/
theorem structLocalMetric_isContinuous (X : Scheme.{0}) :
    (structLocalMetric X).IsContinuous := by
  intro V e
  letI := arcTopology V.toScheme
  rw [arcOpenSet_eq_range]
  refine continuousOn_range_of_comp V _ ?_
  have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme =>
      (structLocalMetric X).h V e (q ≫ V.ι))
      = (fun q => ‖evalOn (q ≫ V.ι) V (comp_preimage_eq_top V q)
          (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e)‖⁻¹) := by
    funext q
    show (if hp : (q ≫ V.ι) ⁻¹ᵁ V = ⊤ then
        ‖evalOn (q ≫ V.ι) V hp (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e)‖⁻¹
      else 1) = _
    rw [dif_pos (comp_preimage_eq_top V q)]
  rw [he]
  refine ((continuous_evalOn V _).norm).inv₀ (fun q => ?_)
  exact norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit _ V (comp_preimage_eq_top V q)
    (isUnit_transUnit _ V (baseTriv X V) e))

/-! ## ★★★★★★★★★(3) 切断のノルムは連続 -/

/-- ★★★★★★**チャートの上では切断のノルムは連続である**。 -/
theorem AMetric.continuousOn_norm {L : AMetric X} (hL : L.IsContinuous)
    (s : L.sheaf.obj (op ⊤)) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V)) :
    @ContinuousOn _ ℝ (arcTopology X) _ (fun p => L.norm s p) (arcOpenSet V) := by
  letI := arcTopology X
  letI := arcTopology V.toScheme
  rw [arcOpenSet_eq_range]
  refine continuousOn_range_of_comp V _ ?_
  have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme => L.norm s (q ≫ V.ι))
      = (fun q => ‖evalOn (q ≫ V.ι) V (comp_preimage_eq_top V q) (trivValue L.sheaf V e s)‖
          * L.metric.h V e (q ≫ V.ι)) := by
    funext q
    rw [AMetric.norm_eq L s _ V e (comp_preimage_eq_top V q)]
    rfl
  rw [he]
  exact ((continuous_evalOn V _).norm).mul (continuous_h_comp hL V e)

/-- ★★★★★★★★★★**切断のノルム `|s|_L̄` は `X^arc` 上連続である**。

原文 (GenEll p.3):
> pair consisting of a line bundle L on X and a hermitian metric

★★これが原文の「hermitian metric on `L^arc`」が**切断の水準で言っていること**である。
★★★局所自明性からどの点にもチャートが取れ、`arcOpenSet` は開なので貼り合う。 -/
theorem AMetric.continuous_norm {L : AMetric X} (hL : L.IsContinuous)
    (s : L.sheaf.obj (op ⊤)) :
    @Continuous _ ℝ (arcTopology X) _ (fun p => L.norm s p) := by
  letI := arcTopology X
  rw [continuous_iff_continuousAt]
  intro p
  obtain ⟨c⟩ := nonempty_normChart L.triv p
  exact (AMetric.continuousOn_norm hL s c.V c.e).continuousAt
    ((isOpen_arcOpenSet c.V).mem_nhds c.hp)

/-! ## ★★★★★★(4) テンソル積 -/

open scoped Classical in
/-- ★★★★**固定したチャートで書ける**（`chartH_indep` から）。 -/
theorem tensor_h_eq_of_chart {A B : X.PresheafOfModules}
    (hlA : IsLocallyTrivial X A) (hlB : IsLocallyTrivial X B)
    (mA : LocalMetric X A) (mB : LocalMetric X B) {V W : X.Opens} (hWV : W ≤ V)
    (eA : (restrictPresheafFunctor X W).obj A ≅ 𝟙_ (PresheafModulesOn X W))
    (eB : (restrictPresheafFunctor X W).obj B ≅ 𝟙_ (PresheafModulesOn X W))
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hpW : p ⁻¹ᵁ W = ⊤) :
    (LocalMetric.tensor hlA hlB mA mB).h V f p
      = ‖evalOn p W hpW (transUnit (A ⊗ B) W (tensorTriv eA eB)
          (trivialOfLe (A ⊗ B) hWV f))‖⁻¹ * (mA.h W eA p * mB.h W eB p) := by
  have hc : Nonempty (TensorChart A B V p) := ⟨⟨W, hWV, hpW, eA, eB⟩⟩
  show (if h : Nonempty (TensorChart A B V p) then chartH mA mB h.some f else 1) = _
  rw [dif_pos hc, chartH_indep mA mB hc.some ⟨W, hWV, hpW, eA, eB⟩ f]
  rfl

/-- ★★★★★★★★**テンソル積は連続性を保つ**。 -/
theorem LocalMetric.tensor_isContinuous {A B : X.PresheafOfModules}
    (hlA : IsLocallyTrivial X A) (hlB : IsLocallyTrivial X B)
    {mA : LocalMetric X A} {mB : LocalMetric X B}
    (hA : mA.IsContinuous) (hB : mB.IsContinuous) :
    (LocalMetric.tensor hlA hlB mA mB).IsContinuous := by
  letI := arcTopology X
  intro V f p₀ hp₀
  obtain ⟨c⟩ := nonempty_tensorChart hlA hlB V p₀ hp₀
  have hcont : @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => (LocalMetric.tensor hlA hlB mA mB).h V f p) (arcOpenSet c.W) := by
    letI := arcTopology c.W.toScheme
    rw [arcOpenSet_eq_range]
    refine continuousOn_range_of_comp c.W _ ?_
    have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ c.W.toScheme =>
        (LocalMetric.tensor hlA hlB mA mB).h V f (q ≫ c.W.ι))
        = (fun q => ‖evalOn (q ≫ c.W.ι) c.W (comp_preimage_eq_top c.W q)
            (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB) (trivialOfLe (A ⊗ B) c.hWV f))‖⁻¹
          * (mA.h c.W c.eA (q ≫ c.W.ι) * mB.h c.W c.eB (q ≫ c.W.ι))) :=
      funext (fun q => tensor_h_eq_of_chart hlA hlB mA mB c.hWV c.eA c.eB f _
        (comp_preimage_eq_top c.W q))
    rw [he]
    refine Continuous.mul (((continuous_evalOn c.W _).norm).inv₀ (fun q => ?_))
      ((continuous_h_comp hA c.W c.eA).mul (continuous_h_comp hB c.W c.eB))
    exact norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit _ c.W (comp_preimage_eq_top c.W q)
      (isUnit_transUnit _ c.W _ _))
  exact (hcont.continuousAt ((isOpen_arcOpenSet c.W).mem_nhds c.hpW)).continuousWithinAt

/-- ★★★★★★★**テンソル積とその一方の因子が連続なら、他方も連続**。

★`h_{A⊗B} = h_A · h_B` を `h_A > 0` で割るだけである。 -/
theorem LocalMetric.isContinuous_of_tensor {A B : X.PresheafOfModules}
    (hlA : IsLocallyTrivial X A) (hlB : IsLocallyTrivial X B)
    {mA : LocalMetric X A} {mB : LocalMetric X B} {m : LocalMetric X (A ⊗ B)}
    (ht : IsTensorOf mA mB m) (hmc : m.IsContinuous) (hA : mA.IsContinuous) :
    mB.IsContinuous := by
  letI := arcTopology X
  intro V eB p₀ hp₀
  obtain ⟨c⟩ := nonempty_tensorChart hlA hlB V p₀ hp₀
  set eB' := trivialOfLe B c.hWV eB with heB'
  have hdiv : @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => m.h c.W (tensorTriv c.eA eB') p / mA.h c.W c.eA p) (arcOpenSet c.W) :=
    (hmc c.W (tensorTriv c.eA eB')).div (hA c.W c.eA)
      (fun p _ => (mA.pos c.W c.eA p).ne')
  have hcont : @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => mB.h V eB p) (arcOpenSet c.W) := by
    refine hdiv.congr (fun p hp => ?_)
    have hpW : p ⁻¹ᵁ c.W = ⊤ := hp
    have hne := (mA.pos c.W c.eA p).ne'
    show mB.h V eB p = m.h c.W (tensorTriv c.eA eB') p / mA.h c.W c.eA p
    rw [ht c.W c.eA eB' p hpW, ← mB.restrict c.hWV eB p hpW]
    field_simp
    rfl
  exact (hcont.continuousAt ((isOpen_arcOpenSet c.W).mem_nhds c.hpW)).continuousWithinAt

/-! ## ★★★★★★(5) 引き戻し -/

open scoped Classical in
/-- ★★★★**固定したチャートで書ける**（`pullChartH_indep` から）。 -/
theorem pullback_h_eq_of_chart {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    (hL : IsLocallyTrivial Y L) (m : LocalMetric Y L) {V V' : X.Opens} {W : Y.Opens}
    (hV'V : V' ≤ V) (hV'W : V' ≤ (Opens.map f.base).obj W)
    (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    (e : (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L) ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hpV' : p ⁻¹ᵁ V' = ⊤) :
    (LocalMetric.pullback f hL m).h V e p
      = ‖evalOn p V' hpV' (transUnit ((pullbackPre f).obj L) V' (pullTrivOfBase f L W eW hV'W)
          (trivialOfLe ((pullbackPre f).obj L) hV'V e))‖⁻¹ * m.h W eW (p ≫ f) := by
  have hc : Nonempty (PullChart f L V p) := ⟨⟨W, V', hV'V, hV'W, hpV', eW⟩⟩
  show (if h : Nonempty (PullChart f L V p) then pullChartH f m h.some e else 1) = _
  rw [dif_pos hc, pullChartH_indep f m hc.some ⟨W, V', hV'V, hV'W, hpV', eW⟩ e]
  rfl

/-- ★★★★★**底の側の計量を引き戻した関数は連続である**。

★`continuous_comp_scheme`（`p ↦ p ≫ f` は連続）との合成である。 -/
theorem continuous_h_comp_base {X Y : Scheme.{0}} (f : X ⟶ Y) {L : Y.PresheafOfModules}
    {m : LocalMetric Y L} (hm : m.IsContinuous) {V' : X.Opens} {W : Y.Opens}
    (hV'W : V' ≤ (Opens.map f.base).obj W)
    (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W)) :
    @Continuous _ ℝ (arcTopology V'.toScheme) _
      (fun q => m.h W eW ((q ≫ V'.ι) ≫ f)) := by
  letI := arcTopology X
  letI := arcTopology Y
  letI := arcTopology V'.toScheme
  refine (hm W eW).comp_continuous
    ((continuous_comp_scheme f).comp (isOpenEmbedding_comp_ι V').continuous) (fun q => ?_)
  show ((q ≫ V'.ι) ≫ f) ⁻¹ᵁ W = ⊤
  rw [Scheme.Hom.comp_preimage]
  exact preimage_eq_top_of_le hV'W (comp_preimage_eq_top V' q)

/-- ★★★★★★★★**引き戻しは連続性を保つ**。 -/
theorem LocalMetric.pullback_isContinuous {X Y : Scheme.{0}} (f : X ⟶ Y)
    {L : Y.PresheafOfModules} (hL : IsLocallyTrivial Y L) {m : LocalMetric Y L}
    (hm : m.IsContinuous) : (LocalMetric.pullback f hL m).IsContinuous := by
  letI := arcTopology X
  intro V e p₀ hp₀
  obtain ⟨c⟩ := nonempty_pullChart f hL V p₀ hp₀
  have hcont : @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => (LocalMetric.pullback f hL m).h V e p) (arcOpenSet c.V') := by
    letI := arcTopology c.V'.toScheme
    rw [arcOpenSet_eq_range]
    refine continuousOn_range_of_comp c.V' _ ?_
    have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ c.V'.toScheme =>
        (LocalMetric.pullback f hL m).h V e (q ≫ c.V'.ι))
        = (fun q => ‖evalOn (q ≫ c.V'.ι) c.V' (comp_preimage_eq_top c.V' q)
            (transUnit ((pullbackPre f).obj L) c.V' (pullTrivOfBase f L c.W c.eW c.hV'W)
              (trivialOfLe ((pullbackPre f).obj L) c.hV'V e))‖⁻¹
          * m.h c.W c.eW ((q ≫ c.V'.ι) ≫ f)) :=
      funext (fun q => pullback_h_eq_of_chart f hL m c.hV'V c.hV'W c.eW e _
        (comp_preimage_eq_top c.V' q))
    rw [he]
    refine Continuous.mul (((continuous_evalOn c.V' _).norm).inv₀ (fun q => ?_))
      (continuous_h_comp_base f hm c.hV'W c.eW)
    exact norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit _ c.V' (comp_preimage_eq_top c.V' q)
      (isUnit_transUnit _ c.V' _ _))
  exact (hcont.continuousAt ((isOpen_arcOpenSet c.V').mem_nhds c.hpV')).continuousWithinAt

/-! ## ★★★★★★★(6) 等長同型で移ること -/

/-- ★★★★★★**連続性は等長同型で移る**。 -/
theorem AMetric.isContinuous_of_isIsometry {L M : AMetric X} {φ : L.sheaf ≅ M.sheaf}
    (hφ : IsIsometry L M φ) (hL : L.IsContinuous) : M.IsContinuous := by
  letI := arcTopology X
  intro V e
  exact (hL V (pullTriv φ V e)).congr (fun p hp => (hφ V e p hp).symm)

theorem AMetric.isContinuous_of_isometric {L M : AMetric X} (h : Isometric L M)
    (hL : L.IsContinuous) : M.IsContinuous :=
  h.elim fun _ hφ => AMetric.isContinuous_of_isIsometry hφ hL

/-! ## ★★★★★★★★★(7) 算術直線束の側と `APic(X)` -/

theorem AMetric.one_isContinuous (X : Scheme.{0}) : (1 : AMetric X).IsContinuous :=
  structLocalMetric_isContinuous X

theorem AMetric.mul_isContinuous {L M : AMetric X} (hL : L.IsContinuous)
    (hM : M.IsContinuous) : (L * M).IsContinuous :=
  LocalMetric.tensor_isContinuous L.triv M.triv hL hM

theorem AMetricPullback_isContinuous {X Y : Scheme.{0}} (f : X ⟶ Y) {L : AMetric Y}
    (hL : L.IsContinuous) : (AMetricPullback f L).IsContinuous :=
  LocalMetric.pullback_isContinuous f L.triv hL

theorem AMetric.isContinuous_of_mul {L M : AMetric X}
    (h : (L * M).IsContinuous) (hL : L.IsContinuous) : M.IsContinuous :=
  LocalMetric.isContinuous_of_tensor L.triv M.triv
    (isTensorOf_tensor L.triv M.triv L.metric M.metric) h hL

def AInv.IsContinuous (L : AInv X) : Prop := L.carrier.IsContinuous

/-- ★★★★★★**逆元は自動で連続になる**。 -/
theorem AInv.inv_isContinuous {L : AInv X} (h : L.IsContinuous) : L.inv.IsContinuous :=
  AMetric.isContinuous_of_mul
    (AMetric.isContinuous_of_isometric (isometric_symm L.isInv)
      (AMetric.one_isContinuous X)) h

theorem AInv.isContinuous_congr {L M : AInv X} (h : (AInv.setoid X).r L M) :
    L.IsContinuous = M.IsContinuous :=
  propext ⟨fun hL => AMetric.isContinuous_of_isometric h hL,
    fun hM => AMetric.isContinuous_of_isometric (isometric_symm h) hM⟩

def APicM.IsContinuous : APicM X → Prop :=
  Quotient.lift AInv.IsContinuous (fun _ _ h => AInv.isContinuous_congr h)

/-- ★★★★★★★★連続な計量を持つ類がなす部分群。 -/
noncomputable def APicK (X : Scheme.{0}) : Subgroup (APicM X) where
  carrier := {a | APicM.IsContinuous a}
  mul_mem' := by
    rintro a b ha hb
    induction a using Quotient.ind with
    | _ L =>
      induction b using Quotient.ind with
      | _ M => exact AMetric.mul_isContinuous ha hb
  one_mem' := AMetric.one_isContinuous X
  inv_mem' := by
    rintro a ha
    induction a using Quotient.ind with
    | _ L => exact AInv.inv_isContinuous ha

/-- ★★★★★★★★★★**[GenEll] Definition 1.1, (i)** —— 原文の群 `APic(X)`。

原文 (GenEll p.3):
> The isomorphism classes of arithmetic line bundles on X, together with the operation
> of tensor product, thus determine a group APic(X).

★原文の「arithmetic line bundle」は**共役両立かつ連続な hermitian 計量**を持つものである。
★★その同型類が `APicM X` の**部分群**をなす——それが本定義である。 -/
noncomputable def APicA (X : Scheme.{0}) : Subgroup (APicM X) := APicC X ⊓ APicK X

theorem AInv.pullback_isContinuous {X Y : Scheme.{0}} (f : X ⟶ Y) {L : AInv Y}
    (h : L.IsContinuous) : (AInv.pullback f L).IsContinuous :=
  AMetricPullback_isContinuous f h

/-! ## ★出典の紐付け(`.src`) -/

def LocalMetric.IsContinuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(計量は X^arc 上連続である——hermitian metric の解析的な側)",
    sectionId := "genell-def-1-1-i" }

def AMetric.continuous_norm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(切断のノルム |s|_L̄ は X^arc 上連続)",
    sectionId := "genell-def-1-1-i" }

def LocalMetric.tensor_isContinuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(テンソル積は連続性を保つ)",
    sectionId := "genell-def-1-1-i" }

def LocalMetric.pullback_isContinuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しは連続性を保つ)",
    sectionId := "genell-def-1-1-ii" }

def APicA.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(共役両立かつ連続な算術直線束の同型類がなす群 APic(X))",
    sectionId := "genell-def-1-1-i" }

def APicA.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "continuous_evalOn(開集合上の評価は連続)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.continuous_evalOn") 3,
    .citation "[ABC3]" "continuousOn_range_of_comp / isOpen_arcOpenSet"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isOpen_arcOpenSet") 3,
    .citation "[ABC3]" "continuous_comp_scheme(p ↦ p ≫ f は連続)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.continuous_comp_scheme") 3,
    .citation "[ABC3]" "APicC(共役両立な類がなす部分群、§9-801)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicC") 3,
    .implicitStep
      ("★原文は「hermitian metric on the line bundle Larc … on Xarc」と書くだけだが、" ++
       "★★そこで畳まれていたのは「計量は複素解析空間 X^arc 上で連続である」" ++
       "という条項であり、それが 1・⊗・逆元・引き戻し・等長で閉じることである") 3 ]

end ABC3.Found.Arakelov
